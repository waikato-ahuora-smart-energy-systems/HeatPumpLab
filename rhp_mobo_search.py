"""COP/total-VHC Bayesian optimisation with optional fixed compressor split."""

import warnings

import numpy as np
import torch
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._torch_env import ensure_ninja_on_path

ensure_ninja_on_path()  # before BoTorch first loads the fused kernel
from botorch.acquisition import AcquisitionFunction
from botorch.acquisition.multi_objective.logei import qLogExpectedHypervolumeImprovement
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.utils.multi_objective.box_decompositions.non_dominated import NondominatedPartitioning
from gpytorch.mlls import ExactMarginalLogLikelihood
from scipy.spatial.distance import cdist
from scipy.stats import qmc

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import COMP2_POWER_FRAC_BOUNDS
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._constraint_model import (
    ConstraintModel, aggregate_margins)


CONSTRAINT_REFIT_EVERY = 5


class FeasibilityWeightedLogEHVI(AcquisitionFunction):
    """LogEHVI + log P(constraint margins >= 0) + distance penalty (q=1).

    The constraint GP covers solved designs that violate a margin; the bounded
    distance penalty covers numerical failures only (no margins reported).
    """

    def __init__(self, model, values, history, success, seed, numerical_failure=None,
                 constraint=None, mc_samples=64):
        super().__init__(model=model)
        self.constraint = constraint
        if numerical_failure is None:
            numerical_failure = ~success
        self.register_buffer("numerical_failure", numerical_failure)
        ref = torch.zeros(2, dtype=torch.double)
        self.acquisition = qLogExpectedHypervolumeImprovement(
            model=model, ref_point=ref.tolist(),
            partitioning=NondominatedPartitioning(ref_point=ref, Y=values),
            sampler=SobolQMCNormalSampler(torch.Size([mc_samples]), seed=seed))
        self.register_buffer("history", history)
        self.register_buffer("success", success)

    def forward(self, X):
        value = self.acquisition(X)
        if self.constraint is not None:
            value = value+self.constraint.log_probability(X)
        if not bool(self.numerical_failure.any()) or not bool(self.success.any()):
            return value
        distances = ((X.squeeze(-2).unsqueeze(-2)-self.history).square().sum(-1)+1e-20).sqrt()
        good = distances[..., self.success].amin(-1)
        bad = distances[..., self.numerical_failure].amin(-1)
        return value+torch.log((bad+.02)/(good+bad+.04))


class ParetoBOSearch:
    """Optimise RHX1, RHX2, compressor split and propane mass fraction together."""

    def __init__(self, max_lift, seed=20260910, *, fixed_comp_frac=None):
        torch.set_num_threads(1)
        if fixed_comp_frac is not None and (not np.isfinite(fixed_comp_frac) or not .1 <= fixed_comp_frac <= .9):
            raise ValueError('fixed_comp_frac must be finite and between 0.1 and 0.9')
        self.fixed_comp_frac = None if fixed_comp_frac is None else float(fixed_comp_frac)
        self.dimension = 4 if self.fixed_comp_frac is None else 3
        self.max_lift = float(max_lift)
        self.seed = int(seed)
        self.x, self.y, self.success = [], [], []
        self.margins = []
        # Proposals handed out in the current batch but not yet told.
        self.pending = []
        # Leading records from ridge seeding; the initial design and step
        # schedule count from after them.
        self.offset = 0
        # Constraint-GP hyperparameters, keyed by training scope; refitted every
        # CONSTRAINT_REFIT_EVERY uses of a scope.
        self.constraint_states = {}
        self.kernel_state = None
        self.mean_state = None
        self.diagnostics = []
        anchors = [[0., 0., split, fraction] for fraction in (0., .25, .5, .75, 1.)
                   for split in (0., .5, 1.)]
        # Most of the front lies on a zero-lift face; seed both faces explicitly.
        faces = [[.3, 0.], [.6, 0.], [0., .3], [0., .6]]
        anchors += [[a, b, .5, fraction] for a, b in faces for fraction in (.25, .75)]
        if self.fixed_comp_frac is not None:
            anchors = [[0., 0., fraction] for fraction in (0., .25, .5, .75, 1.)]
            anchors += [[a, b, fraction] for a, b in faces for fraction in (.25, .75)]
        sobol = self._admissible(qmc.Sobol(self.dimension, scramble=True, seed=seed).random_base2(7))
        self.initial = np.asarray(anchors+sobol[:32-len(anchors)].tolist())

    def _admissible(self, points, tolerance=0.):
        p = np.asarray(points)
        return p[np.isfinite(p).all(1) & (p >= -tolerance).all(1) & (p <= 1+tolerance).all(1)
                 & (100*p[:, :2].sum(1) <= self.max_lift+100*tolerance)]

    def decode(self, u):
        if self.fixed_comp_frac is not None:
            a, b, fraction = u
            return np.asarray([100*a, 100*b, self.fixed_comp_frac, fraction, 1-fraction])
        a, b, split, fraction = u
        lo, hi = COMP2_POWER_FRAC_BOUNDS
        return np.asarray([100*a, 100*b, lo+(hi-lo)*split, fraction, 1-fraction])

    def encode(self, x):
        if self.fixed_comp_frac is not None:
            if not np.isclose(x[2], self.fixed_comp_frac, rtol=0, atol=1e-12):
                raise ValueError('Observation does not match fixed_comp_frac')
            return np.asarray([x[0]/100, x[1]/100, x[3]])
        lo, hi = COMP2_POWER_FRAC_BOUNDS
        return np.asarray([x[0]/100, x[1]/100, (x[2]-lo)/(hi-lo), x[3]])

    def tell(self, x, cop, vhc, success, margins=None):
        u = self.encode(x)
        self.pending = [p for p in self.pending if np.linalg.norm(np.subtract(p, u)) > 1e-9]
        # Encoding a physically valid boundary point can add a few ulps. This
        # check does not repair observations; proposal screening stays strict.
        if len(self._admissible(u[None], tolerance=1e-12)) != 1:
            raise ValueError("Observation is outside the design domain")
        if self.x and cdist(u[None], self.x).min() <= 1e-7:
            raise ValueError("Duplicate observation")
        valid = bool(success and cop is not None and vhc is not None
                     and np.isfinite([cop, vhc]).all() and min(cop, vhc) > 0)
        self.x.append(u.tolist())
        self.y.append([float(cop), float(vhc)] if valid else [0., 0.])
        self.success.append(valid)
        self.margins.append(aggregate_margins(margins))

    def numerical_failures(self):
        return torch.tensor([not s and m is None for s, m in zip(self.success, self.margins)],
                            dtype=torch.bool)

    def constraint_model(self, indices=None, key='global'):
        """GP on margins of solved designs; None until a violation is observed.

        Hyperparameters are warm-started per ``key`` and refitted every
        CONSTRAINT_REFIT_EVERY uses, as the objective GPs are.
        """
        known = [i for i, m in enumerate(self.margins) if m is not None
                 and (indices is None or i in indices)]
        if len(known) < self.dimension+2 or not any(min(self.margins[i]) < 0 for i in known):
            return None
        cached = self.constraint_states.get(key)
        uses = cached['uses'] if cached else 0
        try:
            model = ConstraintModel(np.asarray(self.x)[known], [self.margins[i] for i in known],
                                    fit=uses % CONSTRAINT_REFIT_EVERY == 0,
                                    state=cached['state'] if cached else None)
            self.constraint_states[key] = dict(uses=uses+1, state=model.state())
            return model
        except (RuntimeError, ValueError) as exc:
            self.diagnostics.append(dict(evaluation=len(self.x)+1, stage='constraint GP', message=str(exc)))
            return None

    def _candidates(self):
        n = len(self.x)
        pool = qmc.Sobol(self.dimension, scramble=True, seed=self.seed+104729*n).random_base2(9)
        parts = [pool]
        for axis in range(self.dimension):
            for side in (0., 1.):
                boundary = pool[:64].copy()
                boundary[:, axis] = side
                parts.append(boundary)
        origin = pool[:64].copy()
        origin[:, :2] = 0.
        parts.append(origin)
        pool = np.unique(self._admissible(np.vstack(parts)), axis=0)
        distance = cdist(pool, self.occupied()).min(1)
        return pool[distance > 1e-7], distance[distance > 1e-7]

    def occupied(self):
        """Observed and pending designs (encoded)."""
        return np.asarray(self.x + self.pending) if self.x or self.pending else np.empty((0, self.dimension))

    def step_index(self):
        return len(self.x)+len(self.pending)-self.offset

    def pending_penalty(self, points, radius=.05):
        """log-penalty keeping a batch spread out: 0 beyond radius of every pending point."""
        if not self.pending:
            return np.zeros(len(points))
        d = cdist(np.asarray(points), np.asarray(self.pending)).min(1)
        return np.log(np.clip(d/radius, 1e-3, 1.))

    def ask_batch(self, q):
        """q proposals for parallel evaluation; later ones avoid earlier pending ones."""
        batch = []
        for _ in range(q):
            x = self.ask()
            batch.append(x)
            self.pending.append(self.encode(x).tolist())
        return batch

    def _model(self):
        good = np.asarray(self.success)
        X = torch.tensor(np.asarray(self.x)[good], dtype=torch.double)
        Y = torch.tensor(np.asarray(self.y)[good], dtype=torch.double)
        model = SingleTaskGP(X, Y, train_Yvar=torch.full_like(Y, 1e-8)).double()
        if self.kernel_state is not None:
            with torch.no_grad():
                for module, saved in ((model.covar_module, self.kernel_state), (model.mean_module, self.mean_state)):
                    for name, value in saved.items():
                        module.get_parameter(name).copy_(torch.tensor(value, dtype=torch.double))
        if self.kernel_state is None or len(self.x) % 5 == 0:
            try:
                with warnings.catch_warnings(record=True) as caught:
                    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model),
                                     max_attempts=1, optimizer_kwargs={"options": {"maxiter": 100}})
                self.diagnostics.extend(dict(evaluation=len(self.x)+1, stage="fit", message=str(w.message))
                                        for w in caught)
            except RuntimeError as exc:
                self.diagnostics.append(dict(evaluation=len(self.x)+1, stage="fit fallback", message=str(exc)))
            self.kernel_state = {k: v.detach().tolist() for k, v in model.covar_module.named_parameters()}
            self.mean_state = {k: v.detach().tolist() for k, v in model.mean_module.named_parameters()}
        model.eval()
        return model, Y

    @staticmethod
    def _scores(acquisition, pool):
        with torch.no_grad():
            values = [acquisition(chunk.unsqueeze(-2)) for chunk in
                      torch.tensor(pool, dtype=torch.double).split(128)]
        return torch.cat(values).numpy()

    def face(self):
        """Every other global BO step is restricted to a zero-lift face (0: RHX1 = 0, 1: RHX2 = 0)."""
        return {1: 0, 3: 1}.get(self.step_index() % 4)

    def ask(self):
        n = self.step_index()
        if n < len(self.initial):
            u = np.asarray(self.initial[n])
            occupied = self.occupied()
            if not len(occupied) or cdist(u[None], occupied).min() > 1e-7:
                return self.decode(u)
            # Already observed (e.g. a ridge seed): use the most distant candidate instead.
            pool, distances = self._candidates()
            return self.decode(pool[np.argmax(distances)])
        pool, distances = self._candidates()
        if sum(self.success) < 6 or (n-len(self.initial)) % 8 == 7:
            return self.decode(pool[np.argmax(distances)])
        face = self.face()
        if face is not None:
            on_face = pool[:, face] == 0.
            if on_face.any():
                pool, distances = pool[on_face], distances[on_face]
            else:
                face = None
        with torch.random.fork_rng(devices=[]):
            seed = self.seed+104729*n
            torch.manual_seed(seed)
            model, values = self._model()
            acquisition = FeasibilityWeightedLogEHVI(model, values,
                torch.tensor(self.x, dtype=torch.double), torch.tensor(self.success, dtype=torch.bool), seed,
                numerical_failure=self.numerical_failures(), constraint=self.constraint_model())
            scores = self._scores(acquisition, pool)+self.pending_penalty(pool)
            if not np.isfinite(scores).any():
                self.diagnostics.append(dict(evaluation=n+1, stage="nonfinite acquisition"))
                return self.decode(pool[np.argmax(distances)])
            index = int(np.nanargmax(scores))
            chosen, score = pool[index], scores[index]
            try:
                with warnings.catch_warnings(record=True) as caught:
                    proposals, _ = optimize_acqf(acquisition,
                        bounds=torch.tensor([[0.]*self.dimension, [1.]*self.dimension], dtype=torch.double),
                        q=1, num_restarts=6, raw_samples=128,
                        inequality_constraints=[(torch.tensor([0, 1]), torch.tensor([-1., -1.], dtype=torch.double),
                                                 -self.max_lift/100)],
                        fixed_features=None if face is None else {face: 0.},
                        options={"maxiter": 60, "batch_limit": 3}, return_best_only=False,
                        retry_on_optimization_warning=False)
                self.diagnostics.extend(dict(evaluation=n+1, stage="acquisition", message=str(w.message))
                                        for w in caught)
                proposals = self._admissible(proposals.detach().numpy().reshape(-1, self.dimension))
                if face is not None and len(proposals):
                    proposals = proposals[proposals[:, face] == 0.]
                if len(proposals):
                    proposals = proposals[cdist(proposals, self.occupied()).min(1) > 1e-7]
                if len(proposals):
                    new_scores = self._scores(acquisition, proposals)+self.pending_penalty(proposals)
                    if np.isfinite(new_scores).any() and np.nanmax(new_scores) > score:
                        chosen = proposals[np.nanargmax(new_scores)]
            except (RuntimeError, ValueError) as exc:
                self.diagnostics.append(dict(evaluation=n+1, stage="acquisition fallback", message=str(exc)))
            return self.decode(chosen)
