"""BoTorch single-objective LogEI with the existing RHP exploration policies."""

import warnings

import numpy as np
import torch
from botorch.acquisition import AcquisitionFunction, LogExpectedImprovement
from botorch.exceptions.errors import ModelFittingError
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.optim import optimize_acqf
from gpytorch.constraints import Interval
from gpytorch.kernels import MaternKernel
from gpytorch.means import ZeroMean
from gpytorch.mlls import ExactMarginalLogLikelihood
from scipy.spatial.distance import cdist

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.origin_search import OriginSurrogateSearch


class FeasibilityWeightedLogEI(AcquisitionFunction):
    """log(EI * existing distance-based feasibility heuristic), for q=1."""

    def __init__(self, model, best_f, history, success):
        super().__init__(model=model)
        self.log_ei = LogExpectedImprovement(model, best_f=best_f)
        self.register_buffer("history", history)
        self.register_buffer("success", success)

    def forward(self, X):
        value = self.log_ei(X)
        if bool(self.success.all()):
            return value
        delta = X.squeeze(-2).unsqueeze(-2) - self.history
        distances = torch.sqrt(torch.sum(delta.square(), dim=-1) + 1e-20)
        dg = distances[..., self.success].amin(dim=-1)
        db = distances[..., ~self.success].amin(dim=-1)
        return value + torch.log((db + 0.02)/(dg + db + 0.04))


class BotorchLogEISearch(OriginSurrogateSearch):
    """A controlled framework comparison, not BoTorch's default GP settings.

    Retains the custom BO bootstrap, Matérn GP assumptions, failure heuristic,
    explicit boundary candidate pool and every-eighth maximin step. BoTorch
    supplies GP fitting, LogEI and continuous acquisition optimisation.
    """

    def __init__(self, *, max_lift, seed, origin=(0., 0., 0.5)):
        super().__init__(method="gp_ei", max_lift=max_lift, seed=seed, origin=origin)
        torch.set_num_threads(1)
        self.seed = int(seed)
        self.diagnostics = []

    def _model(self, x, y, lengths):
        kernel = MaternKernel(nu=2.5, ard_num_dims=3,
                              lengthscale_constraint=Interval(0.02, 3.)).double()
        kernel.lengthscale = torch.tensor(lengths, dtype=torch.double)
        return SingleTaskGP(x, y, train_Yvar=torch.full_like(y, 1e-6),
                            covar_module=kernel, mean_module=ZeroMean(),
                            outcome_transform=None).double()

    def _fit_model(self):
        good = np.asarray(self.success)
        x = torch.tensor(np.asarray(self.x)[good], dtype=torch.double)
        values = np.asarray(self.y)[good]
        # BoTorch maximises: train on negative standardised SEC.
        y = torch.tensor(-(values-values.mean())/max(values.std(), 1e-8),
                         dtype=torch.double).unsqueeze(-1)
        fits = []
        if len(self.x) == len(self.initial) or len(self.x) % 5 == 0:
            for lengths in (self.lengthscales, np.full(3, 0.15)):
                model = self._model(x, y, lengths)
                mll = ExactMarginalLogLikelihood(model.likelihood, model)
                try:
                    with warnings.catch_warnings(record=True) as caught:
                        fit_gpytorch_mll(mll, max_attempts=1,
                                         optimizer_kwargs={"options": {"maxiter": 100}})
                    if caught:
                        self.diagnostics.append(dict(evaluation=len(self.x)+1, stage="GP warnings",
                                                     messages=[str(w.message) for w in caught]))
                    model.train()
                    with torch.no_grad():
                        loss = float(-mll(model(x), y.squeeze(-1)))
                    if np.isfinite(loss):
                        fits.append((loss, model))
                except (ModelFittingError, RuntimeError) as exc:
                    self.diagnostics.append(dict(evaluation=len(self.x)+1, stage="GP fit", message=str(exc)))
            if fits:
                model = min(fits, key=lambda pair: pair[0])[1]
                self.lengthscales = model.covar_module.lengthscale.detach().numpy().reshape(3).copy()
            else:
                model = self._model(x, y, self.lengthscales)
        else:
            model = self._model(x, y, self.lengthscales)
        model.eval()
        model.likelihood.eval()
        return model, y.max()

    @staticmethod
    def _scores(acquisition, points):
        with torch.no_grad():
            batches = [acquisition(chunk.unsqueeze(-2)) for chunk in
                       torch.tensor(points, dtype=torch.double).split(512)]
        return torch.cat(batches).numpy()

    def ask(self):
        if len(self.x) < len(self.initial):
            return self.decode(self.initial[len(self.x)])
        candidates, distances = self._candidates()
        good = np.asarray(self.success)
        if good.sum() < 5 or (len(self.x)-len(self.initial)) % 8 == 7:
            return self.decode(candidates[np.argmax(distances.min(axis=1))])
        # Scoped RNG makes restarts and checkpoint reconstruction deterministic.
        with torch.random.fork_rng(devices=[]):
            torch.manual_seed(self.seed + 104729 * len(self.x))
            model, best_f = self._fit_model()
            acquisition = FeasibilityWeightedLogEI(
                model, best_f, torch.tensor(np.asarray(self.x), dtype=torch.double),
                torch.tensor(good, dtype=torch.bool))
            pool_scores = self._scores(acquisition, candidates)
            if not np.isfinite(pool_scores).any():
                self.diagnostics.append(dict(evaluation=len(self.x)+1, stage="nonfinite acquisition"))
                return self.decode(candidates[np.argmax(distances.min(axis=1))])
            best_index = int(np.nanargmax(pool_scores))
            chosen, chosen_score = candidates[best_index], pool_scores[best_index]
            try:
                with warnings.catch_warnings(record=True) as caught:
                    proposed, _ = optimize_acqf(
                        acquisition, bounds=torch.tensor([[0., 0., 0.], [1., 1., 1.]], dtype=torch.double),
                        q=1, num_restarts=12, raw_samples=512,
                        inequality_constraints=[(torch.tensor([0, 1]),
                                                 torch.tensor([-1., -1.], dtype=torch.double),
                                                 -self.max_lift/100)],
                        options={"maxiter": 100, "batch_limit": 4}, return_best_only=False,
                        retry_on_optimization_warning=False)
                if caught:
                    self.diagnostics.append(dict(evaluation=len(self.x)+1, stage="acquisition warnings",
                                                 messages=[str(w.message) for w in caught]))
                pool = proposed.detach().numpy().reshape(-1, 3)
                # Optimiser tolerance can place results just outside a constraint.
                # Reject these cheaply; never silently repair an expensive point.
                pool = pool[np.all(np.isfinite(pool), axis=1)]
                pool = pool[np.all((pool >= 0) & (pool <= 1), axis=1)]
                pool = self._admissible(pool)
                if len(pool):
                    pool = pool[cdist(pool, self.x).min(axis=1) > 1e-7]
                if len(pool):
                    scores = self._scores(acquisition, pool)
                    index = int(np.nanargmax(scores))
                    if scores[index] > chosen_score:
                        chosen = pool[index]
            except (RuntimeError, ValueError) as exc:
                self.diagnostics.append(dict(evaluation=len(self.x)+1, stage="acquisition fallback", message=str(exc)))
            return self.decode(chosen)
