"""Joint BO for COP, total VHC and an exactly known, discrete delivery temperature.

All temperatures share two GP outputs and a five-dimensional input space. The
128-point training subset is independent of the complete 3D hypervolume baseline.
"""
from __future__ import annotations

import time
import warnings

import numpy as np
import torch
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._torch_env import ensure_ninja_on_path

ensure_ninja_on_path()  # before BoTorch first loads the fused kernel
from botorch.acquisition import AcquisitionFunction
from botorch.acquisition.multi_objective.logei import qLogExpectedHypervolumeImprovement
from botorch.acquisition.multi_objective.objective import GenericMCMultiOutputObjective
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.optim import optimize_acqf
from botorch.sampling.normal import SobolQMCNormalSampler
from gpytorch.mlls import ExactMarginalLogLikelihood
from scipy.spatial import cKDTree
from scipy.stats import qmc

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._constraint_model import (
    ConstraintModel, aggregate_margins)

TEMPERATURES = np.arange(100., 221., 5.)
MILESTONES = (5000, 10000, 15000, 20000, 25000)
OUTCOME_SCALES = np.array([4., 10.])
CONSTRAINT_REFIT_EVERY = 5
PENDING_RADIUS = .05


def encode(x, temperature):
    if temperature not in TEMPERATURES:
        raise ValueError('Temperature must be on the 100:5:220 C grid')
    return np.array([x[0]/100, x[1]/100, (x[2]-.1)/.8, x[3], (temperature-100)/120])


def decode(u):
    u = np.asarray(u, dtype=float)
    temperature = 100.+5.*round(float(u[4])*24)
    if abs((temperature-100)/120-u[4]) > 1e-10:
        raise ValueError('Acquisition proposed an off-grid temperature')
    return np.array([100*u[0], 100*u[1], .1+.8*u[2], u[3], 1-u[3]]), temperature


def admissible(points, tolerance=1e-12):
    p = np.asarray(points, dtype=float).reshape(-1, 5)
    max_lift = (100+120*p[:, 4]-15)*2/3
    return (np.isfinite(p).all(1) & (p >= -tolerance).all(1)
            & (p <= 1+tolerance).all(1)
            & (abs(p[:, 4]*24-np.rint(p[:, 4]*24)) <= tolerance*24)
            & (100*p[:, :2].sum(1) <= max_lift+100*tolerance))


def objectives(record):
    return np.array([record['COP']/4, record['VHC_MJ_per_m3']/10,
                     (record['T_air_tar']-95)/125])


def usable(record):
    if not record['success'] or record.get('validation_rejected'):
        return False
    return bool(np.isfinite(objectives(record)).all() and (objectives(record) > 0).all())


def update_front(indices, values, new_index):
    """Incrementally retain the complete maximisation front, including ties."""
    value = values[new_index]
    old = np.asarray([values[i] for i in indices]).reshape(-1, 3)
    if np.any(np.all(old >= value, axis=1) & np.any(old > value, axis=1)):
        return indices
    dominated = np.all(value >= old, axis=1) & np.any(value > old, axis=1)
    return [i for i, remove in zip(indices, dominated) if not remove]+[new_index]


def front_indices(records):
    values, front = [], []
    for index, record in enumerate(records):
        values.append(objectives(record) if usable(record) else np.zeros(3))
        if usable(record):
            front = update_front(front, values, index)
    return front


def front_2d(points):
    points = np.unique(np.asarray(points, dtype=float).reshape(-1, 2), axis=0)
    if not len(points):
        return points
    order = np.lexsort((-points[:, 1], -points[:, 0]))
    selected, best = [], -np.inf
    for i in order:
        if points[i, 1] > best:
            selected.append(i)
            best = points[i, 1]
    return points[selected][::-1]  # Increasing COP, decreasing VHC.


def hypervolume_3d(points):
    """Exact volume against (COP=0, VHC=0, T=95 C), in fixed scaled units."""
    points = np.asarray(points, dtype=float).reshape(-1, 3)
    volume = 0.
    for level in range(1, 26):
        front = front_2d(points[points[:, 2] >= level/25-1e-12, :2])
        if len(front):
            volume += np.sum(np.diff(np.r_[0., front[:, 0]])*front[:, 1])/25
    return float(volume)


class TemperatureGridPartition:
    """Exact non-dominated boxes obtained from the 25 temperature slabs.

    The boxes cover the positive objective orthant. This avoids general 3D
    partitioning growth without omitting any observed Pareto points.
    """
    num_outcomes = 3

    def __init__(self, values):
        values = np.asarray(values, dtype=float).reshape(-1, 3)
        self.pareto_Y = torch.tensor(values, dtype=torch.double)
        lower, upper = [], []
        for level in range(1, 26):
            front = front_2d(values[values[:, 2] >= level/25-1e-12, :2])
            previous = 0.
            for cop, vhc in front:
                lower.append([previous, vhc, (level-1)/25])
                upper.append([cop, np.inf, level/25])
                previous = cop
            lower.append([previous, 0., (level-1)/25])
            upper.append([np.inf, np.inf, level/25])
        lower.append([0., 0., 1.])
        upper.append([np.inf, np.inf, np.inf])
        self.bounds = torch.tensor([lower, upper], dtype=torch.double)

    def get_hypercell_bounds(self):
        return self.bounds


def add_known_temperature(samples, X):
    """Attach the requested temperature to each MC draw with zero uncertainty."""
    known = (1+24*X[..., 4:5])/25
    return torch.cat((samples, known.expand(*samples.shape[:-1], 1)), dim=-1)


class JointAcquisition(AcquisitionFunction):
    """Joint qLogEHVI + log P(constraint margins >= 0) + distance penalty for numerical failures.

    Solved designs that violate a margin inform the constraint GP; the bounded distance
    penalty applies only to numerical failures, which report no margins.
    """

    def __init__(self, model, partition, history, success, seed, numerical_failure=None, constraint=None):
        super().__init__(model=model)
        self.constraint = constraint
        failure = ~np.asarray(success, dtype=bool) if numerical_failure is None else np.asarray(numerical_failure, dtype=bool)
        self.register_buffer('numerical_failure', torch.tensor(failure, dtype=torch.bool))
        self.ehvi = qLogExpectedHypervolumeImprovement(
            model=model, ref_point=[0., 0., 0.], partitioning=partition,
            objective=GenericMCMultiOutputObjective(add_known_temperature),
            sampler=SobolQMCNormalSampler(torch.Size([64]), seed=seed))
        self.register_buffer('history', torch.tensor(history, dtype=torch.double))
        self.register_buffer('success', torch.tensor(success, dtype=torch.bool))

    def forward(self, X):
        value = self.ehvi(X)
        if self.constraint is not None:
            value = value+self.constraint.log_probability(X)
        if not bool(self.numerical_failure.any()) or not bool(self.success.any()):
            return value
        # A sampling heuristic, not a claim that native failures are physical.
        distances = (X.squeeze(-2).unsqueeze(-2)-self.history).square().sum(-1).clamp_min(1e-20).sqrt()
        good = distances[..., self.success].amin(-1)
        bad = distances[..., self.numerical_failure].amin(-1)
        return value+torch.log((bad+.02)/(good+bad+.04))


class TemperatureParetoSearch:
    def __init__(self, records=(), *, seed=20260913, subset_size=128, state=None):
        torch.set_num_threads(1)
        self.seed, self.subset_size = int(seed), int(subset_size)
        self.records, self.x, self.y, self.success, self.front = [], [], [], [], []
        self.margins, self.pending = [], []
        for record in records:
            self.tell(record)
        # Ridge-seed records come first; adaptive step counts start after them.
        self.offset = sum(r.get('method') == 'ridge_seed' for r in self.records)
        self.state = state or dict(epoch=-1, centres=[], models={}, constraints={})
        self.diagnostic = {}
        # Explicit endpoint/composition anchors at every allowed temperature.
        anchors = [[0., 0., 0., fraction, level/24]
                   for level in range(25) for fraction in (0., .25, .5, .75, 1.)]
        pool = self.global_pool(1024, self.seed)
        self.initial = np.vstack((anchors, pool[:75]))

    def global_pool(self, count, seed):
        pool = qmc.Sobol(5, scramble=True, seed=seed).random_base2(int(np.ceil(np.log2(count))))
        pool[:, 4] = np.minimum(24, np.floor(pool[:, 4]*25))/24
        return pool[admissible(pool)]

    def tell(self, record):
        u = encode(record['x'], record['T_air_tar'])
        if not admissible(u[None])[0]:
            raise ValueError('Observation outside the temperature-dependent design domain')
        if self.x and cKDTree(self.x).query(u)[0] <= 1e-7:
            raise ValueError('Duplicate joint design/temperature observation')
        good = usable(record)
        self.pending = [p for p in self.pending if np.linalg.norm(np.subtract(p, u)) > 1e-9]
        self.records.append(record)
        self.x.append(u.tolist())
        self.y.append(objectives(record).tolist() if good else [0., 0., 0.])
        self.success.append(good)
        self.margins.append(aggregate_margins(record.get('constraint_margins')))
        if good:
            self.front = update_front(self.front, self.y, len(self.y)-1)

    def numerical_failures(self):
        return np.array([not s and m is None for s, m in zip(self.success, self.margins)], dtype=bool)

    def occupied(self):
        return np.asarray(self.x+self.pending) if self.x or self.pending else np.empty((0, 5))

    def step_index(self):
        """Adaptive step count: observed plus pending, excluding ridge seeds."""
        return len(self.x)+len(self.pending)-self.offset

    def pending_penalty(self, points):
        """log-penalty keeping a batch spread out: 0 beyond PENDING_RADIUS of every pending point."""
        if not self.pending:
            return np.zeros(len(points))
        d = cKDTree(self.pending).query(np.asarray(points))[0]
        return np.log(np.clip(d/PENDING_RADIUS, 1e-3, 1.))

    def ask_batch(self, q):
        """q proposals for parallel evaluation; later ones avoid earlier pending ones."""
        batch, diagnostics = [], []
        for _ in range(q):
            x, temperature = self.ask()
            batch.append((x, temperature))
            diagnostics.append(dict(self.diagnostic))
            self.pending.append(encode(x, temperature).tolist())
        self.batch_diagnostics = diagnostics
        return batch

    def constraint_model(self, centre, key):
        """GP on margins of the solved designs nearest the region centre; None until a violation.

        Hyperparameters are warm-started per region and refitted every CONSTRAINT_REFIT_EVERY uses.
        """
        known = np.array([i for i, m in enumerate(self.margins) if m is not None])
        if len(known) < 7:
            return None
        distance = np.linalg.norm(np.asarray(self.x)[known]-centre, axis=1)
        known = known[np.lexsort((known, distance))[:self.subset_size]]
        if not any(min(self.margins[i]) < 0 for i in known):
            return None
        states = self.state.setdefault('constraints', {})
        cached = states.get(key)
        uses = cached['uses'] if cached else 0
        try:
            model = ConstraintModel(np.asarray(self.x)[known], [self.margins[i] for i in known],
                                    fit=uses % CONSTRAINT_REFIT_EVERY == 0,
                                    state=cached['state'] if cached else None)
        except (RuntimeError, ValueError) as exc:
            warnings.warn(f'Constraint GP skipped: {exc}', stacklevel=2)
            return None
        states[key] = dict(uses=uses+1, state=model.state())
        return model

    def training_indices(self, centre):
        good = np.flatnonzero(self.success)
        distance = np.linalg.norm(np.asarray(self.x)[good]-centre, axis=1)
        return good[np.lexsort((good, distance))[:self.subset_size]]

    def _centres(self):
        epoch = max(0, (self.step_index()-len(self.initial))//48)
        if self.state['epoch'] != epoch:
            candidates = self.front or list(np.flatnonzero(self.success))
            y = np.asarray(self.y)[candidates]
            scaled = (y-y.min(0))/np.maximum(np.ptp(y, axis=0), 1e-12)
            selected = list(dict.fromkeys(int(np.argmax(scaled[:, i])) for i in range(3)))
            while len(selected) < min(8, len(candidates)):
                distances = np.linalg.norm(scaled[:, None]-scaled[selected], axis=-1).min(1)
                distances[selected] = -1
                selected.append(int(np.argmax(distances)))
            self.state = dict(epoch=epoch, centres=[self.x[candidates[i]] for i in selected], models={}, constraints={})
        return np.asarray(self.state['centres'])

    @staticmethod
    def scores(acquisition, pool):
        cells = len(acquisition.ehvi.cell_lower_bounds)
        batch_size = max(1, min(8, 4096//max(1, cells)))
        with torch.no_grad():
            return torch.cat([acquisition(chunk.unsqueeze(-2)) for chunk in
                              torch.tensor(pool, dtype=torch.double).split(batch_size)]).numpy()

    def ask(self):
        started, n = time.perf_counter(), self.step_index()
        occupied = self.occupied()
        tree = cKDTree(occupied) if len(occupied) else None
        seed = (self.seed+104729*n) % (2**31-1)
        if n < len(self.initial):
            # Initial points already observed (e.g. ridge-seed S = 0 designs) become maximin points.
            if tree is None or tree.query(self.initial[n])[0] > 1e-7:
                self.diagnostic = dict(method='initial', training_points=0)
                return decode(self.initial[n])
            pool = self.global_pool(2048, seed)
            self.diagnostic = dict(method='initial_replaced_maximin', training_points=0)
            return decode(pool[np.argmax(tree.query(pool)[0])])
        history = np.asarray(self.x)
        if sum(self.success) < 8 or (n-len(self.initial)) % 8 == 7:
            pool = self.global_pool(2048, seed)
            distance = tree.query(pool)[0]
            self.diagnostic = dict(method='global_maximin', training_points=0)
            return decode(pool[np.argmax(distance)])
        centres = self._centres()
        region = (n-len(self.initial)) % len(centres)
        centre = centres[region]
        lower, upper = np.maximum(centre[:4]-.08, 0), np.minimum(centre[:4]+.08, 1)
        local = qmc.Sobol(4, scramble=True, seed=seed).random_base2(4)
        local = lower+(upper-lower)*local
        # Enumerate temperatures: no continuous optimisation followed by rounding.
        pool = np.vstack([np.column_stack((local, np.full(len(local), level/24))) for level in range(25)])
        pool = np.vstack((pool, self.global_pool(256, seed+1)))
        pool = np.unique(pool[admissible(pool)], axis=0)
        pool = pool[tree.query(pool)[0] > 1e-7]
        indices = self.training_indices(centre)
        X = torch.tensor(history[indices], dtype=torch.double)
        Y = torch.tensor(np.asarray(self.y)[indices, :2], dtype=torch.double)
        key = str(region)
        old = self.state['models'].get(key)
        with torch.random.fork_rng(devices=[]), warnings.catch_warnings(record=True) as caught:
            torch.manual_seed(seed)
            model = SingleTaskGP(X, Y, train_Yvar=torch.full_like(Y, 1e-8)).double()
            if old:
                with torch.no_grad():
                    for module, saved in ((model.covar_module, old['kernel']), (model.mean_module, old['mean'])):
                        for name, value in saved.items():
                            module.get_parameter(name).copy_(torch.tensor(value, dtype=torch.double))
            visits = old['visits'] if old else 0
            fit_error = None
            if not old or visits % 5 == 0:
                try:
                    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model),
                                     max_attempts=1, optimizer_kwargs={'options': {'maxiter': 100}})
                except RuntimeError as exc:
                    fit_error = str(exc)
            self.state['models'][key] = dict(visits=visits+1,
                kernel={k: v.detach().tolist() for k, v in model.covar_module.named_parameters()},
                mean={k: v.detach().tolist() for k, v in model.mean_module.named_parameters()})
            model.eval()
            partition = TemperatureGridPartition(np.asarray(self.y)[self.front])
            constraint = self.constraint_model(centre, key)
            acquisition = JointAcquisition(model, partition, self.x, self.success, seed,
                                           numerical_failure=self.numerical_failures(), constraint=constraint)
            scores = self.scores(acquisition, pool)+self.pending_penalty(pool)
            if not np.isfinite(scores).any():
                raise RuntimeError('No finite joint acquisition scores')
            best = int(np.nanargmax(scores))
            chosen, score = pool[best], float(scores[best])
            # Refine the best two enumerated temperatures with temperature fixed.
            selected = []
            for index in np.argsort(-scores):
                temperature = float(pool[index, 4])
                if temperature in selected:
                    continue
                selected.append(temperature)
                try:
                    proposal, _ = optimize_acqf(acquisition,
                        bounds=torch.tensor([[0.]*5, [1.]*5], dtype=torch.double),
                        q=1, num_restarts=1, raw_samples=None,
                        batch_initial_conditions=torch.tensor(pool[index][None, None], dtype=torch.double),
                        fixed_features={4: temperature},
                        inequality_constraints=[(torch.tensor([0, 1]), torch.tensor([-1., -1.], dtype=torch.double),
                                                 -(100+120*temperature-15)*2/300)],
                        options={'maxiter': 35, 'batch_limit': 1, 'init_batch_limit': 8},
                        retry_on_optimization_warning=False)
                    proposed = proposal.detach().numpy().reshape(-1, 5)
                    proposed = proposed[admissible(proposed)]
                    if len(proposed) and tree.query(proposed)[0][0] > 1e-7:
                        value = float(self.scores(acquisition, proposed)[0]+self.pending_penalty(proposed)[0])
                        if np.isfinite(value) and value > score:
                            chosen, score = proposed[0], value
                except (ValueError, RuntimeError) as exc:
                    warnings.warn(f'Acquisition refinement skipped: {exc}', stacklevel=2)
                if len(selected) == 2:
                    break
            self.diagnostic = dict(method='joint_qLogEHVI', region=int(region), training_points=len(indices),
                training_indices=indices.tolist(), training_temperature_levels=sorted({self.records[i]['T_air_tar'] for i in indices}),
                frontier_points=len(self.front), hypercells=len(partition.bounds[0]), candidate_pool=len(pool),
                constraint_gp=constraint is not None and constraint.model is not None, pending=len(self.pending),
                selected_log_ehvi_weighted=score, fit_error=fit_error,
                warnings=list(dict.fromkeys(str(w.message) for w in caught)),
                proposal_seconds=time.perf_counter()-started)
        return decode(chosen)
