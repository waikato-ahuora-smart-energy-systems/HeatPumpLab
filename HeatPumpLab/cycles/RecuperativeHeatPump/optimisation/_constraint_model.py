"""GP models of continuous constraint margins for the RHP searches.

Every solved evaluation, feasible or not, reports margins that are >= 0 when a
constraint holds (compressor discharge and suction superheat, and stage
pressure-ratio headroom). A GP on these margins gives a smooth probability of
feasibility, which multiplies the acquisition. This replaces the bounded
distance penalty for constraint violations; the distance penalty stays only for
numerical failures, which report no margins.
"""

from __future__ import annotations

import warnings

import numpy as np
import torch
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from gpytorch.mlls import ExactMarginalLogLikelihood

__all__ = ["MARGIN_GROUPS", "ConstraintModel", "aggregate_margins"]

# (group name, key suffix, scale, clip). Scaled margins are clipped so that
# designs far inside or outside a constraint do not dominate the fit.
MARGIN_GROUPS = (
    ("discharge_superheat", "_discharge_superheat", 10.0, 3.0),
    ("suction_superheat", "_suction_superheat", 10.0, 3.0),
    ("pressure_ratio", "_pressure_ratio", 1.0, 3.0),
)


def aggregate_margins(margins):
    """Return the scaled worst-stage margin per group, or None if any is unknown."""
    if not margins:
        return None
    output = []
    for _, suffix, scale, clip in MARGIN_GROUPS:
        values = [value for key, value in margins.items() if key.endswith(suffix)]
        if not values or any(value is None or not np.isfinite(value) for value in values):
            return None
        output.append(float(np.clip(min(values) / scale, -clip, clip)))
    return output


class ConstraintModel(torch.nn.Module):
    """Independent GPs for each margin group; returns log P(all margins >= 0)."""

    def __init__(self, x, margins, *, fit=True, state=None):
        """state: hyperparameters from a previous fit (see ``state()``). They are
        loaded when the active groups match; the model is then refitted only if
        ``fit`` is true. Without a usable state the model is always fitted."""
        super().__init__()
        X = torch.as_tensor(np.asarray(x), dtype=torch.double)
        Y = torch.as_tensor(np.asarray(margins), dtype=torch.double)
        # Groups that never approach violation carry no information: drop them.
        self.active = [j for j in range(Y.shape[-1]) if float(Y[:, j].min()) < 0.5]
        self.model = None
        if not self.active:
            return
        Y = Y[:, self.active]
        self.model = SingleTaskGP(X, Y, train_Yvar=torch.full_like(Y, 1e-6)).double()
        warm = state is not None and state.get("active") == self.active
        if warm:
            with torch.no_grad():
                for module, saved in ((self.model.covar_module, state["covar"]),
                                      (self.model.mean_module, state["mean"])):
                    for name, value in saved.items():
                        module.get_parameter(name).copy_(torch.tensor(value, dtype=torch.double))
        self.fitted = bool(fit or not warm)
        if self.fitted:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fit_gpytorch_mll(
                    ExactMarginalLogLikelihood(self.model.likelihood, self.model),
                    max_attempts=2,
                    optimizer_kwargs={"options": {"maxiter": 100}},
                )
        self.model.eval()

    def state(self):
        """Hyperparameters for warm starts (outcome scaling is not included)."""
        if self.model is None:
            return None
        return dict(active=list(self.active),
                    covar={k: v.detach().tolist() for k, v in self.model.covar_module.named_parameters()},
                    mean={k: v.detach().tolist() for k, v in self.model.mean_module.named_parameters()})

    def log_probability(self, X):
        """X: (..., 1, d) or (..., d). Returns (...) log P(feasible)."""
        if self.model is None:
            shape = X.shape[:-2] if X.dim() >= 3 else X.shape[:-1]
            return torch.zeros(shape, dtype=X.dtype)
        if X.dim() >= 3:
            X = X.squeeze(-2)
        posterior = self.model.posterior(X)
        mean = posterior.mean
        sd = posterior.variance.clamp_min(1e-12).sqrt()
        return torch.special.log_ndtr(mean / sd).sum(-1)

    def probability(self, X):
        with torch.no_grad():
            return self.log_probability(torch.as_tensor(X, dtype=torch.double)).exp().numpy()
