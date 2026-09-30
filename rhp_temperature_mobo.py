"""BO for the COP-VHC front at each delivery temperature on the 100:5:220 C grid.

Temperature is a condition, not an objective: the target is the separate 2D
(COP, VHC) Pareto front at every grid temperature, scored as the mean of the
per-temperature hypervolumes. All temperatures share one GP (five inputs,
temperature included, two outputs), so designs learnt at one temperature inform
the others. Each proposal targets one temperature on a fixed rotation and uses
qLogEHVI against that temperature's own front.

The 3D (COP, VHC, temperature) helpers below (front_indices, hypervolume_3d,
TemperatureGridPartition) are kept for reference and comparison with v19-v24.
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
from botorch.utils.multi_objective.box_decompositions.non_dominated import NondominatedPartitioning
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
LEVELS = len(TEMPERATURES)
SCHEDULE_STRIDE = 7        # adaptive step k targets level (7k) mod 25: every level once per 25 steps
CENTRES_PER_LEVEL = 4
LOCAL_POINTS, GLOBAL_POINTS = 64, 128


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


def level_of(u):
    return int(round(float(u[4])*24))


def per_temperature_front_indices(records):
    """Union of the 2D (COP, VHC) fronts at each temperature, ties kept."""
    by_level = {}
    for i, r in enumerate(records):
        if usable(r):
            by_level.setdefault(r['T_air_tar'], []).append(i)
    out = []
    for indices in by_level.values():
        values = {i: objectives(records[i])[:2] for i in indices}
        front = []
        for i in indices:
            front = update_front_2d(front, values, i)
        out.extend(front)
    return sorted(out)


def update_front_2d(indices, values, new_index):
    """Incrementally retain a 2D maximisation front, including ties."""
    value = np.asarray(values[new_index])
    old = np.asarray([values[i] for i in indices]).reshape(-1, 2)
    if np.any(np.all(old >= value, axis=1) & np.any(old > value, axis=1)):
        return indices
    dominated = np.all(value >= old, axis=1) & np.any(value > old, axis=1)
    return [i for i, remove in zip(indices, dominated) if not remove]+[new_index]


def hypervolume_2d(points):
    front = front_2d(points)
    if not len(front):
        return 0.
    return float(np.sum(np.diff(np.r_[0., front[:, 0]])*front[:, 1]))


def per_temperature_hypervolumes(records):
    """Scaled 2D hypervolume (COP/4, VHC/10 against the origin) at each grid temperature."""
    points = {float(t): [] for t in TEMPERATURES}
    for r in records:
        if usable(r):
            points[r['T_air_tar']].append(objectives(r)[:2])
    return {t: hypervolume_2d(p) for t, p in points.items()}


def mean_temperature_hypervolume(records):
    """Mean over the 25 grid temperatures of the per-temperature 2D hypervolume."""
    return float(np.mean(list(per_temperature_hypervolumes(records).values())))


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


class LevelAcquisition(AcquisitionFunction):
    """qLogEHVI against one temperature's 2D front + log P(margins >= 0) + numerical-failure penalty."""

    def __init__(self, model, front_values, history, success, seed, numerical_failure, constraint=None):
        super().__init__(model=model)
        self.constraint = constraint
        Y = torch.tensor(np.asarray(front_values, dtype=float).reshape(-1, 2), dtype=torch.double)
        if not len(Y):
            Y = torch.full((1, 2), 1e-9, dtype=torch.double)
        ref = torch.zeros(2, dtype=torch.double)
        self.ehvi = qLogExpectedHypervolumeImprovement(
            model=model, ref_point=ref.tolist(), partitioning=NondominatedPartitioning(ref_point=ref, Y=Y),
            sampler=SobolQMCNormalSampler(torch.Size([64]), seed=seed))
        self.register_buffer('history', torch.tensor(history, dtype=torch.double))
        self.register_buffer('success', torch.tensor(success, dtype=torch.bool))
        self.register_buffer('numerical_failure', torch.tensor(numerical_failure, dtype=torch.bool))

    def forward(self, X):
        value = self.ehvi(X)
        if self.constraint is not None:
            value = value+self.constraint.log_probability(X)
        if not bool(self.numerical_failure.any()) or not bool(self.success.any()):
            return value
        distances = (X.squeeze(-2).unsqueeze(-2)-self.history).square().sum(-1).clamp_min(1e-20).sqrt()
        good = distances[..., self.success].amin(-1)
        bad = distances[..., self.numerical_failure].amin(-1)
        return value+torch.log((bad+.02)/(good+bad+.04))


class TemperatureParetoSearch:
    def __init__(self, records=(), *, seed=20260913, subset_size=128, state=None):
        torch.set_num_threads(1)
        self.seed, self.subset_size = int(seed), int(subset_size)
        self.records, self.x, self.y, self.success = [], [], [], []
        self.fronts = {level: [] for level in range(LEVELS)}
        self.margins, self.pending = [], []
        for record in records:
            self.tell(record)
        # Ridge-seed records come first; adaptive step counts start after them.
        self.offset = sum(r.get('method') == 'ridge_seed' for r in self.records)
        self.state = state or dict(models={}, constraints={}, visits={})
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
            level = level_of(u)
            self.fronts[level] = update_front_2d(self.fronts[level], {i: self.y[i][:2] for i in
                                                  self.fronts[level]+[len(self.y)-1]}, len(self.y)-1)

    @property
    def front(self):
        """Union of the per-temperature fronts (observation indices)."""
        return sorted(i for f in self.fronts.values() for i in f)

    def hypervolume(self):
        """Mean scaled per-temperature 2D hypervolume."""
        return float(np.mean([hypervolume_2d(np.asarray(self.y)[f, :2]) if f else 0. for f in self.fronts.values()]))

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

    def level_centres(self, level):
        """Local-region centres for one temperature, spread along its 2D front.

        A sparse front borrows equipment designs from the nearest temperatures' fronts
        (lifts scaled into this temperature's lift limit), which is where the shared GP helps.
        """
        t = level/24
        candidates = [self.x[i] for i in self.fronts[level]]
        values = [self.y[i][:2] for i in self.fronts[level]]
        for other in sorted(range(LEVELS), key=lambda l: (abs(l-level), l)):
            if len(candidates) >= CENTRES_PER_LEVEL:
                break
            if other == level:
                continue
            for i in self.fronts[other]:
                u = np.array(self.x[i], dtype=float)
                u[4] = t
                limit = (100+120*t-15)*2/300
                if u[0]+u[1] > limit:
                    u[:2] *= limit/(u[0]+u[1])
                candidates.append(u.tolist())
                values.append(self.y[i][:2])
        if not candidates:
            return np.empty((0, 5))
        y = np.asarray(values)
        scaled = (y-y.min(0))/np.maximum(np.ptp(y, axis=0), 1e-12)
        selected = list(dict.fromkeys(int(np.argmax(scaled[:, i])) for i in range(2)))
        while len(selected) < min(CENTRES_PER_LEVEL, len(candidates)):
            distances = np.linalg.norm(scaled[:, None]-scaled[selected], axis=-1).min(1)
            distances[selected] = -1
            selected.append(int(np.argmax(distances)))
        return np.asarray([candidates[i] for i in selected])

    def level_pool(self, count, seed, level):
        pool = qmc.Sobol(4, scramble=True, seed=seed).random_base2(int(np.ceil(np.log2(count))))
        pool = np.column_stack((pool, np.full(len(pool), level/24)))
        return pool[admissible(pool)]

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
        k = n-len(self.initial)
        level = (SCHEDULE_STRIDE*k) % LEVELS
        t = level/24
        history = np.asarray(self.x)
        centres = self.level_centres(level)
        if sum(self.success) < 8 or k % 8 == 7 or not len(centres):
            pool = self.level_pool(1024, seed, level)
            self.diagnostic = dict(method='global_maximin', training_points=0, level=level,
                                   temperature_C=100.+5*level)
            return decode(pool[np.argmax(tree.query(pool)[0])])
        visits = self.state.setdefault('visits', {})
        visit = visits.get(str(level), 0)
        visits[str(level)] = visit+1
        region = visit % len(centres)
        centre = centres[region]
        lower, upper = np.maximum(centre[:4]-.08, 0), np.minimum(centre[:4]+.08, 1)
        local = qmc.Sobol(4, scramble=True, seed=seed).random_base2(int(np.log2(LOCAL_POINTS)))
        local = np.column_stack((lower+(upper-lower)*local, np.full(len(local), t)))
        pool = np.vstack((local, self.level_pool(GLOBAL_POINTS, seed+1, level)))
        pool = np.unique(pool[admissible(pool)], axis=0)
        pool = pool[tree.query(pool)[0] > 1e-7]
        if not len(pool):
            pool = self.level_pool(1024, seed, level)
            self.diagnostic = dict(method='global_maximin', training_points=0, level=level,
                                   temperature_C=100.+5*level)
            return decode(pool[np.argmax(tree.query(pool)[0])])
        indices = self.training_indices(centre)
        X = torch.tensor(history[indices], dtype=torch.double)
        Y = torch.tensor(np.asarray(self.y)[indices, :2], dtype=torch.double)
        key = f'{level}:{region}'
        old = self.state['models'].get(key)
        with torch.random.fork_rng(devices=[]), warnings.catch_warnings(record=True) as caught:
            torch.manual_seed(seed)
            model = SingleTaskGP(X, Y, train_Yvar=torch.full_like(Y, 1e-8)).double()
            if old:
                with torch.no_grad():
                    for module, saved in ((model.covar_module, old['kernel']), (model.mean_module, old['mean'])):
                        for name, value in saved.items():
                            module.get_parameter(name).copy_(torch.tensor(value, dtype=torch.double))
            fits = old['visits'] if old else 0
            fit_error = None
            if not old or fits % 5 == 0:
                try:
                    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model),
                                     max_attempts=1, optimizer_kwargs={'options': {'maxiter': 100}})
                except RuntimeError as exc:
                    fit_error = str(exc)
            self.state['models'][key] = dict(visits=fits+1,
                kernel={k: v.detach().tolist() for k, v in model.covar_module.named_parameters()},
                mean={k: v.detach().tolist() for k, v in model.mean_module.named_parameters()})
            model.eval()
            front_values = np.asarray(self.y)[self.fronts[level], :2] if self.fronts[level] else np.empty((0, 2))
            constraint = self.constraint_model(centre, key)
            acquisition = LevelAcquisition(model, front_values, self.x, self.success, seed,
                                           self.numerical_failures(), constraint=constraint)
            scores = self.scores(acquisition, pool)+self.pending_penalty(pool)
            if not np.isfinite(scores).any():
                raise RuntimeError('No finite acquisition scores')
            best = int(np.nanargmax(scores))
            chosen, score = pool[best], float(scores[best])
            try:
                proposal, _ = optimize_acqf(acquisition,
                    bounds=torch.tensor([[0.]*5, [1.]*5], dtype=torch.double),
                    q=1, num_restarts=1, raw_samples=None,
                    batch_initial_conditions=torch.tensor(pool[best][None, None], dtype=torch.double),
                    fixed_features={4: t},
                    inequality_constraints=[(torch.tensor([0, 1]), torch.tensor([-1., -1.], dtype=torch.double),
                                             -(100+120*t-15)*2/300)],
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
            self.diagnostic = dict(method='temperature_qLogEHVI', level=level, temperature_C=100.+5*level,
                region=int(region), training_points=len(indices),
                training_indices=indices.tolist(),
                training_temperature_levels=sorted({self.records[i]['T_air_tar'] for i in indices}),
                level_front_points=len(self.fronts[level]), candidate_pool=len(pool),
                selected_log_ehvi_weighted=score, fit_error=fit_error,
                constraint_gp=constraint is not None and constraint.model is not None, pending=len(self.pending),
                warnings=list(dict.fromkeys(str(w.message) for w in caught)),
                proposal_seconds=time.perf_counter()-started)
        return decode(chosen)
