"""Four-dimensional COP/total-VHC Bayesian optimisation with BoTorch LogEHVI."""

import warnings

import numpy as np
import torch
from botorch.acquisition import AcquisitionFunction
from botorch.acquisition.multi_objective.logei import qLogExpectedHypervolumeImprovement
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.utils.multi_objective.box_decompositions.non_dominated import (
    NondominatedPartitioning,
)
from gpytorch.mlls import ExactMarginalLogLikelihood
from scipy.spatial.distance import cdist
from scipy.stats import qmc


class FeasibilityWeightedLogEHVI(AcquisitionFunction):
    """LogEHVI times the existing distance-based feasibility heuristic (q=1)."""

    def __init__(self, model, values, history, success, seed):
        super().__init__(model=model)
        ref = torch.zeros(2, dtype=torch.double)
        self.acquisition = qLogExpectedHypervolumeImprovement(
            model=model,
            ref_point=ref.tolist(),
            partitioning=NondominatedPartitioning(ref_point=ref, Y=values),
            sampler=SobolQMCNormalSampler(torch.Size([64]), seed=seed),
        )
        self.register_buffer("history", history)
        self.register_buffer("success", success)

    def forward(self, X):
        value = self.acquisition(X)
        if bool(self.success.all()):
            return value
        distances = (
            (X.squeeze(-2).unsqueeze(-2) - self.history).square().sum(-1) + 1e-20
        ).sqrt()
        good = distances[..., self.success].amin(-1)
        bad = distances[..., ~self.success].amin(-1)
        return value + torch.log((bad + 0.02) / (good + bad + 0.04))


class ParetoBOSearch:
    """Optimise RHX1, RHX2, compressor split and first-component mass fraction together."""

    def __init__(self, max_lift, seed=20260910):
        torch.set_num_threads(1)
        self.max_lift = float(max_lift)
        self.seed = int(seed)
        self.x, self.y, self.success = [], [], []
        self.kernel_state = None
        self.mean_state = None
        self.diagnostics = []
        anchors = [
            [0.0, 0.0, split, fraction]
            for fraction in (0.0, 0.25, 0.5, 0.75, 1.0)
            for split in (0.0, 0.5, 1.0)
        ]
        sobol = self._admissible(qmc.Sobol(4, scramble=True, seed=seed).random_base2(7))
        self.initial = np.asarray(anchors + sobol[:17].tolist())

    def _admissible(self, points, tolerance=0.0):
        p = np.asarray(points)
        return p[
            np.isfinite(p).all(1)
            & (p >= -tolerance).all(1)
            & (p <= 1 + tolerance).all(1)
            & (100 * p[:, :2].sum(1) <= self.max_lift + 100 * tolerance)
        ]

    @staticmethod
    def decode(u):
        a, b, split, fraction = u
        return np.asarray([100 * a, 100 * b, 0.1 + 0.8 * split, fraction, 1 - fraction])

    @staticmethod
    def encode(x):
        return np.asarray([x[0] / 100, x[1] / 100, (x[2] - 0.1) / 0.8, x[3]])

    def tell(self, x, cop, vhc, success):
        u = self.encode(x)
        # Encoding a physically valid boundary point can add a few ulps. This
        # check does not repair observations; proposal screening stays strict.
        if len(self._admissible(u[None], tolerance=1e-12)) != 1:
            raise ValueError("Observation is outside the design domain")
        if self.x and cdist(u[None], self.x).min() <= 1e-7:
            raise ValueError("Duplicate observation")
        valid = bool(
            success
            and cop is not None
            and vhc is not None
            and np.isfinite([cop, vhc]).all()
            and min(cop, vhc) > 0
        )
        self.x.append(u.tolist())
        self.y.append([float(cop), float(vhc)] if valid else [0.0, 0.0])
        self.success.append(valid)

    def _candidates(self):
        n = len(self.x)
        pool = qmc.Sobol(4, scramble=True, seed=self.seed + 104729 * n).random_base2(9)
        parts = [pool]
        for axis in range(4):
            for side in (0.0, 1.0):
                boundary = pool[:64].copy()
                boundary[:, axis] = side
                parts.append(boundary)
        origin = pool[:64].copy()
        origin[:, :2] = 0.0
        parts.append(origin)
        pool = np.unique(self._admissible(np.vstack(parts)), axis=0)
        distance = cdist(pool, self.x).min(1)
        return pool[distance > 1e-7], distance[distance > 1e-7]

    def _model(self):
        good = np.asarray(self.success)
        X = torch.tensor(np.asarray(self.x)[good], dtype=torch.double)
        Y = torch.tensor(np.asarray(self.y)[good], dtype=torch.double)
        model = SingleTaskGP(X, Y, train_Yvar=torch.full_like(Y, 1e-8)).double()
        if self.kernel_state is not None:
            with torch.no_grad():
                for module, saved in (
                    (model.covar_module, self.kernel_state),
                    (model.mean_module, self.mean_state),
                ):
                    for name, value in saved.items():
                        module.get_parameter(name).copy_(
                            torch.tensor(value, dtype=torch.double)
                        )
        if self.kernel_state is None or len(self.x) % 5 == 0:
            try:
                with warnings.catch_warnings(record=True) as caught:
                    fit_gpytorch_mll(
                        ExactMarginalLogLikelihood(model.likelihood, model),
                        max_attempts=1,
                        optimizer_kwargs={"options": {"maxiter": 100}},
                    )
                self.diagnostics.extend(
                    dict(
                        evaluation=len(self.x) + 1, stage="fit", message=str(w.message)
                    )
                    for w in caught
                )
            except RuntimeError as exc:
                self.diagnostics.append(
                    dict(
                        evaluation=len(self.x) + 1,
                        stage="fit fallback",
                        message=str(exc),
                    )
                )
            self.kernel_state = {
                k: v.detach().tolist() for k, v in model.covar_module.named_parameters()
            }
            self.mean_state = {
                k: v.detach().tolist() for k, v in model.mean_module.named_parameters()
            }
        model.eval()
        return model, Y

    @staticmethod
    def _scores(acquisition, pool):
        with torch.no_grad():
            values = [
                acquisition(chunk.unsqueeze(-2))
                for chunk in torch.tensor(pool, dtype=torch.double).split(8)
            ]
        return torch.cat(values).numpy()

    def ask(self):
        if len(self.x) < len(self.initial):
            return self.decode(self.initial[len(self.x)])
        pool, distances = self._candidates()
        if sum(self.success) < 6 or (len(self.x) - len(self.initial)) % 8 == 7:
            return self.decode(pool[np.argmax(distances)])
        with torch.random.fork_rng(devices=[]):
            seed = self.seed + 104729 * len(self.x)
            torch.manual_seed(seed)
            model, values = self._model()
            acquisition = FeasibilityWeightedLogEHVI(
                model,
                values,
                torch.tensor(self.x, dtype=torch.double),
                torch.tensor(self.success, dtype=torch.bool),
                seed,
            )
            scores = self._scores(acquisition, pool)
            if not np.isfinite(scores).any():
                self.diagnostics.append(
                    dict(evaluation=len(self.x) + 1, stage="nonfinite acquisition")
                )
                return self.decode(pool[np.argmax(distances)])
            index = int(np.nanargmax(scores))
            chosen, score = pool[index], scores[index]
            try:
                with warnings.catch_warnings(record=True) as caught:
                    proposals, _ = optimize_acqf(
                        acquisition,
                        bounds=torch.tensor([[0.0] * 4, [1.0] * 4], dtype=torch.double),
                        q=1,
                        num_restarts=6,
                        raw_samples=128,
                        inequality_constraints=[
                            (
                                torch.tensor([0, 1]),
                                torch.tensor([-1.0, -1.0], dtype=torch.double),
                                -self.max_lift / 100,
                            )
                        ],
                        options={
                            "maxiter": 60,
                            "batch_limit": 1,
                            "init_batch_limit": 8,
                        },
                        return_best_only=False,
                        retry_on_optimization_warning=False,
                    )
                self.diagnostics.extend(
                    dict(
                        evaluation=len(self.x) + 1,
                        stage="acquisition",
                        message=str(w.message),
                    )
                    for w in caught
                )
                proposals = self._admissible(proposals.detach().numpy().reshape(-1, 4))
                if len(proposals):
                    proposals = proposals[cdist(proposals, self.x).min(1) > 1e-7]
                if len(proposals):
                    new_scores = self._scores(acquisition, proposals)
                    if np.isfinite(new_scores).any() and np.nanmax(new_scores) > score:
                        chosen = proposals[np.nanargmax(new_scores)]
            except (RuntimeError, ValueError) as exc:
                self.diagnostics.append(
                    dict(
                        evaluation=len(self.x) + 1,
                        stage="acquisition fallback",
                        message=str(exc),
                    )
                )
            return self.decode(chosen)
