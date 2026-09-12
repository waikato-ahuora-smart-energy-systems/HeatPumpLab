"""Joint, dimension-independent BoTorch LogEI search over the RHP design domain."""

import warnings

import numpy as np
import torch
from botorch.acquisition import AcquisitionFunction, LogExpectedImprovement
from botorch.exceptions.errors import ModelFittingError
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.optim import optimize_acqf
from gpytorch.kernels import MaternKernel
from gpytorch.mlls import ExactMarginalLogLikelihood
from scipy.spatial.distance import cdist
from scipy.stats import qmc

from ._types import DEFAULT_BOUNDS


class DesignDomain:
    """Remove fixed coordinates and represent composition by n-1 fractions."""

    def __init__(self, fractions, max_lift, config):
        self.config = config
        self.fractions = np.asarray(fractions, dtype=float)
        self.max_lift = float(max_lift)
        self.base = np.r_[
            config.fixed_dT_rhx1_cold,
            config.fixed_dT_rhx2_cold,
            config.fixed_comp2_power_frac,
            self.fractions,
        ]
        self.axes = [
            i
            for i, enabled in enumerate(
                (
                    config.vary_dT_rhx1_cold,
                    config.vary_dT_rhx2_cold,
                    config.vary_comp2_power_frac,
                )
            )
            if enabled
        ]
        # A zero residual lift allowance makes the nominal lift axes fixed.
        fixed_lift = sum(self.base[i] for i in (0, 1) if i not in self.axes)
        if not np.isfinite(self.max_lift) or fixed_lift > self.max_lift:
            raise ValueError("Fixed recuperator lifts exceed the combined lift limit")
        if fixed_lift == self.max_lift:
            self.axes = [i for i in self.axes if i == 2]
            for i, enabled in enumerate(
                (config.vary_dT_rhx1_cold, config.vary_dT_rhx2_cold)
            ):
                if enabled:
                    self.base[i] = 0.0
        self.nc = len(fractions) - 1 if config.vary_composition else 0
        self.dimension = len(self.axes) + self.nc
        # Free lifts start at zero; fixed variables retain their prescribed value.
        for i in self.axes:
            self.base[i] = 0.0 if i < 2 else config.fixed_comp2_power_frac
        if not self.contains(self.encode(self.base)):
            raise ValueError("Initial design is outside the feasible domain")

    def encode(self, x):
        x = np.asarray(x, dtype=float)
        return np.r_[
            (x[self.axes] - DEFAULT_BOUNDS[self.axes, 0])
            / np.diff(DEFAULT_BOUNDS[self.axes], axis=1).ravel(),
            x[3 : 3 + self.nc],
        ]

    def decode(self, u):
        u = np.asarray(u, dtype=float)
        x = self.base.copy()
        k = len(self.axes)
        x[self.axes] = (
            DEFAULT_BOUNDS[self.axes, 0]
            + u[:k] * np.diff(DEFAULT_BOUNDS[self.axes], axis=1).ravel()
        )
        if self.nc:
            x[3:-1] = u[k:]
            x[-1] = 1.0 - sum(u[k:])
        return x

    def contains(self, u, tolerance=1e-12):
        u = np.asarray(u)
        if u.shape != (self.dimension,) or not np.isfinite(u).all():
            return False
        if np.any(u < -tolerance) or np.any(u > 1 + tolerance):
            return False
        x = self.decode(u)
        return bool(
            min(x[3:]) >= -tolerance and sum(x[:2]) <= self.max_lift + 100 * tolerance
        )

    def constraints(self):
        constraints = []
        indices = [j for j, i in enumerate(self.axes) if i < 2]
        fixed = sum(self.base[i] for i in (0, 1) if i not in self.axes)
        if indices:
            constraints.append(
                (
                    torch.tensor(indices),
                    -torch.ones(len(indices), dtype=torch.double),
                    -(self.max_lift - fixed) / 100,
                )
            )
        if self.nc:
            indices = torch.arange(len(self.axes), self.dimension)
            constraints.append(
                (indices, -torch.ones(self.nc, dtype=torch.double), -1.0)
            )
        return constraints or None


class FeasibilityWeightedLogEI(AcquisitionFunction):
    """Distance-based failure avoidance, not a calibrated feasibility probability."""

    def __init__(self, model, best_f, history, success):
        super().__init__(model=model)
        self.log_ei = LogExpectedImprovement(model, best_f=best_f)
        self.register_buffer("history", history)
        self.register_buffer("success", success)

    def forward(self, X):
        value = self.log_ei(X)
        if bool(self.success.all()):
            return value
        distance = (
            (X.squeeze(-2).unsqueeze(-2) - self.history).square().sum(-1) + 1e-20
        ).sqrt()
        good = distance[..., self.success].amin(-1)
        bad = distance[..., ~self.success].amin(-1)
        return value + torch.log((bad + 0.02) / (good + bad + 0.04))


class BayesianDesignSearch:
    """Origin, boundary-aware bootstrap, then joint GP LogEI and maximin exploration."""

    def __init__(self, domain, config):
        self.domain, self.config = domain, config
        self.x, self.y, self.success = [], [], []
        self.diagnostics = []
        torch.set_num_threads(1)
        self.initial = [domain.encode(domain.base)]
        if domain.dimension:
            u = self.initial[0]
            for axis in range(domain.dimension):
                for step in (0.05, -0.05):
                    p = u.copy()
                    p[axis] += step
                    if domain.contains(p, tolerance=0):
                        self.initial.append(p)
            if domain.nc:
                for fractions in np.eye(domain.nc + 1):
                    x = domain.base.copy()
                    x[3:] = fractions
                    self.initial.append(domain.encode(x))
            self.initial = list({tuple(p): p for p in self.initial}.values())
        self.initial = self.initial[: config.initial_evaluations]

    def observed(self, x):
        u = self.domain.encode(x)
        return bool(
            self.x and (self.domain.dimension == 0 or cdist([u], self.x).min() <= 1e-7)
        )

    def tell(self, x, objective, success):
        u = self.domain.encode(x)
        if not self.domain.contains(u) or not np.allclose(
            self.domain.decode(u), x, atol=1e-10, rtol=0
        ):
            raise ValueError(
                "Observation changed fixed variables or violates the domain"
            )
        if self.observed(x):
            raise ValueError("Duplicate BO observation")
        self.x.append(u.tolist())
        valid = bool(success and np.isfinite(objective))
        self.success.append(valid)
        self.y.append(float(objective) if valid else 0.0)

    def _candidates(self):
        d = self.domain.dimension
        seed = self.config.seed + 104729 * len(self.x)
        size = self.config.candidate_pool_size
        pool = qmc.Sobol(d, scramble=True, seed=seed).random_base2(
            int(np.ceil(np.log2(size)))
        )[:size]
        if self.domain.nc:
            # Uniform simplex samples avoid exponentially wasteful cube rejection.
            rng = np.random.default_rng(seed)
            pool[:, len(self.domain.axes) :] = rng.dirichlet(
                np.ones(self.domain.nc + 1), size
            )[:, :-1]
        parts = [pool]
        for axis in range(d):
            for bound in (0.0, 1.0):
                p = pool[: max(1, size // 8)].copy()
                p[:, axis] = bound
                parts.append(p)
        if self.domain.nc:
            for fractions in np.eye(self.domain.nc + 1):
                p = pool[: max(1, size // 8)].copy()
                p[:, len(self.domain.axes) :] = fractions[:-1]
                parts.append(p)
            p = pool.copy()
            comp = p[:, len(self.domain.axes) :]
            comp /= comp.sum(1, keepdims=True)
            parts.append(p)
        # Include the zero-lift face for joint composition/loading proposals.
        p = pool.copy()
        for j, axis in enumerate(self.domain.axes):
            if axis < 2:
                p[:, j] = 0.0
        parts.append(p)
        candidates = np.unique(np.vstack(parts), axis=0)
        candidates = np.asarray(
            [p for p in candidates if self.domain.contains(p, tolerance=0)]
        )
        if not len(candidates):
            raise RuntimeError("No feasible BO candidates")
        distances = cdist(candidates, self.x).min(1)
        return candidates[distances > 1e-7], distances[distances > 1e-7]

    def _fit(self):
        good = np.asarray(self.success)
        x = torch.tensor(np.asarray(self.x)[good], dtype=torch.double)
        values = np.asarray(self.y)[good]
        y = torch.tensor(
            -(values - values.mean()) / max(values.std(), 1e-8), dtype=torch.double
        )[:, None]
        model = SingleTaskGP(
            x,
            y,
            train_Yvar=torch.full_like(y, 1e-6),
            covar_module=MaternKernel(nu=2.5, ard_num_dims=self.domain.dimension),
            outcome_transform=None,
        ).double()
        fit_gpytorch_mll(
            ExactMarginalLogLikelihood(model.likelihood, model),
            max_attempts=1,
            optimizer_kwargs={"options": {"maxiter": 100}},
        )
        model.eval()
        return FeasibilityWeightedLogEI(
            model,
            y.max(),
            torch.tensor(self.x, dtype=torch.double),
            torch.tensor(good, dtype=torch.bool),
        )

    @staticmethod
    def _scores(acquisition, pool):
        with torch.no_grad():
            return torch.cat(
                [
                    acquisition(chunk[:, None])
                    for chunk in torch.tensor(pool, dtype=torch.double).split(256)
                ]
            ).numpy()

    def ask(self):
        for u in self.initial:
            x = self.domain.decode(u)
            if not self.observed(x):
                return x
        if not self.domain.dimension:
            return None
        pool, distances = self._candidates()
        if not len(pool):
            raise RuntimeError("BO exhausted the distinct candidate pool")
        fallback = pool[np.argmax(distances)]
        if len(self.x) < self.config.initial_evaluations or sum(self.success) < max(
            3, self.domain.dimension + 1
        ):
            return self.domain.decode(fallback)
        if (
            len(self.x) - self.config.initial_evaluations
        ) % self.config.exploration_interval == self.config.exploration_interval - 1:
            return self.domain.decode(fallback)
        chosen = fallback
        with (
            torch.random.fork_rng(devices=[]),
            warnings.catch_warnings(record=True) as caught,
        ):
            torch.manual_seed(self.config.seed + 104729 * len(self.x))
            try:
                acquisition = self._fit()
                scores = self._scores(acquisition, pool)
                if not np.isfinite(scores).any():
                    raise ValueError("No finite acquisition scores")
                index = int(np.nanargmax(scores))
                chosen, score = pool[index], scores[index]
                # Joint acquisition optimisation is cheap; it makes no cycle calls.
                try:
                    proposals, _ = optimize_acqf(
                        acquisition,
                        bounds=torch.tensor(
                            [
                                [0.0] * self.domain.dimension,
                                [1.0] * self.domain.dimension,
                            ],
                            dtype=torch.double,
                        ),
                        q=1,
                        num_restarts=self.config.acquisition_restarts,
                        raw_samples=128,
                        inequality_constraints=self.domain.constraints(),
                        options={"maxiter": 60, "batch_limit": 3},
                        return_best_only=False,
                        retry_on_optimization_warning=False,
                    )
                    proposals = (
                        proposals.detach().numpy().reshape(-1, self.domain.dimension)
                    )
                    proposals = np.asarray(
                        [
                            p
                            for p in proposals
                            if self.domain.contains(p, tolerance=0)
                            and cdist([p], self.x).min() > 1e-7
                        ]
                    )
                    if len(proposals):
                        scores = self._scores(acquisition, proposals)
                        if np.isfinite(scores).any() and np.nanmax(scores) > score:
                            chosen = proposals[np.nanargmax(scores)]
                except (ModelFittingError, RuntimeError, ValueError) as exc:
                    self.diagnostics.append(
                        dict(
                            evaluation=len(self.x) + 1,
                            stage="acquisition",
                            message=str(exc),
                        )
                    )
            except (ModelFittingError, RuntimeError, ValueError) as exc:
                self.diagnostics.append(
                    dict(
                        evaluation=len(self.x) + 1,
                        stage="GP fallback",
                        message=str(exc),
                    )
                )
            self.diagnostics.extend(
                dict(
                    evaluation=len(self.x) + 1, stage="warning", message=str(w.message)
                )
                for w in caught
            )
        return self.domain.decode(chosen)
