"""Experimental local Pareto BO with full-data or nearest-neighbour GP models."""
import time
import warnings

import numpy as np
import torch
from botorch.fit import fit_gpytorch_mll
from botorch.models import SingleTaskGP
from botorch.optim import optimize_acqf
from gpytorch.mlls import ExactMarginalLogLikelihood
from scipy.spatial.distance import cdist
from scipy.stats import qmc

from rhp_mobo_search import ParetoBOSearch, FeasibilityWeightedLogEHVI


REFRESH_EVERY = 48
# Local proposal cost (the measured bottleneck is pool scoring and acquisition
# optimisation). A 0.05-radius box in 3-4 dimensions needs far fewer points
# than the global domain.
LOCAL_POOL_LOG2 = 8          # 256 Sobol points (was 512)
LOCAL_EDGE_POINTS = 32       # per box face (was 64)
LOCAL_MC_SAMPLES = 32        # qLogEHVI samples (was 64)
LOCAL_RESTARTS, LOCAL_RAW_SAMPLES = 4, 64   # acquisition optimisation (was 6, 128)


def max_lift_for(metadata):
    """Combined recuperator-lift limit for a run, from its base case."""
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _max_recuperator_cold_lift
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import DEFAULT_BASE_CASE
    return _max_recuperator_cold_lift(metadata.get('base_case') or dict(DEFAULT_BASE_CASE, T_air_tar=200.))


def restore_search(result):
    """Restore local-128 continuation state using the maintained search code.

    Anchor epochs and cached GP parameters retain the existing checkpoint
    format. Historical experiment runners keep their own frozen implementation.
    """
    phase = result['local128_phase']
    start = phase['start_evaluations']
    epoch = start + ((len(result['records']) - start) // REFRESH_EVERY) * REFRESH_EVERY
    baseline = dict(records=result['records'][:epoch], kernel_state=result.get('kernel_state'),
                    mean_state=result.get('mean_state'))
    search = LocalParetoSearch(baseline, 128, seed=phase['seed'] + 104729 * (epoch-start),
                              fixed_comp_frac=result.get('metadata', {}).get('fixed_comp_frac'),
                              max_lift=max_lift_for(result.get('metadata', {})),
                              components=len(result.get('metadata', {}).get('components') or [0, 0]),
                              min_fraction=result.get('metadata', {}).get('min_component_fraction', 0.))
    for record in result['records'][epoch:]:
        search.tell(record['x'], record['COP'], record['VHC_MJ_per_m3'],
                    record['success'] and not record.get('validation_rejected'),
                    record.get('constraint_margins'))
    saved = result.get('local128_state', {})
    if saved.get('epoch_start') == epoch:
        key = lambda k: int(k) if str(k).lstrip('-').isdigit() else k
        search.states = {key(k): value for k, value in saved['model_states'].items()}
        search.constraint_states = {key(k): v for k, v in saved.get('constraint_states', {}).items()}
    return search, epoch


def front_indices(values):
    """Maximise two finite positive objectives; keep exact duplicate optima."""
    y = np.asarray(values, dtype=float)
    if not len(y):
        return []
    valid = np.flatnonzero(np.isfinite(y).all(1) & (y > 0).all(1))
    order = sorted(valid, key=lambda i: (-y[i, 0], -y[i, 1]))
    best_v, best_c, front = -np.inf, -np.inf, []
    for i in order:
        c, v = y[i]
        if v > best_v:
            front.append(int(i))
            best_v, best_c = v, c
        elif v == best_v and c == best_c:
            front.append(int(i))
    return front


def hypervolume(values):
    y = np.asarray(values)
    indices = sorted(front_indices(y), key=lambda i: y[i, 1])
    previous, result = 0., 0.
    for i in indices:
        result += (y[i, 1]-previous)*y[i, 0]
        previous = y[i, 1]
    return float(result)


def _stratum(x):
    return 'rhx1_zero' if x[0] <= 1e-9 else ('rhx2_zero' if x[1] <= 1e-9 else 'interior')


def choose_anchors(records, count=6):
    """Pick local-region centres: one per zero-lift face, the rest spread along the front.

    The face anchors use each face's own nondominated set, so a branch of the
    front that the global search has only begun to explore still gets a local
    region. The remaining anchors keep the previous rule: the max-COP and
    max-VHC mixed-fluid front points, then maximin spacing in normalised
    objective space.

    Tested alternative (phase 8, Oct 2026): choosing the end and maximin anchors
    over the whole front, pure fluids included, refines the pure-fluid edge but
    lowers final hypervolume for every pair tested (Propane/Pentane, Propane/Butane,
    Isobutane/Pentane; 3 seeds with ridge seeding), because regions moved to that
    one-dimensional edge are taken from the mixed-fluid front.

    If fewer than `count` mixed designs are on the overall front, the next
    nondominated layers of the mixed designs fill the shortfall. A ValueError
    is raised only when there are fewer than `count` successful mixed designs.
    """
    values = np.asarray([[r['COP'], r['VHC_MJ_per_m3']] if r['success'] and
                         not r.get('validation_rejected') else [0., 0.] for r in records])
    mixed = np.array([max(r['x'][3:]) <= .99 for r in records])
    front = [i for i in front_indices(values) if mixed[i]]
    if len(front) < count:
        # Pure-fluid designs can dominate all but a few mixed designs (e.g. Neopentane
        # ternaries, where one evaluation removed five mixed front points at a refresh).
        # Top up from successive nondominated layers of the mixed designs alone.
        remaining = [i for i in range(len(records)) if mixed[i] and values[i, 0] > 0 and i not in front]
        while len(front) < count and remaining:
            layer = [remaining[k] for k in front_indices(values[remaining])]
            front += layer
            remaining = [i for i in remaining if i not in layer]
    if len(front) < count:
        raise ValueError('Too few mixed-fluid Pareto designs for the planned regions')
    scale = np.maximum(np.ptp(values[front], axis=0), 1e-12)
    low = values[front].min(0)
    norm = lambda idx: (values[idx]-low)/scale
    selected = []
    envelope = lambda i: max(values[j, 0] for j in front if values[j, 1] >= values[i, 1] - 1e-12) \
        if any(values[j, 1] >= values[i, 1] for j in front) else values[i, 0]
    for face in ('rhx1_zero', 'rhx2_zero'):
        members = [i for i in range(len(records)) if mixed[i] and values[i, 0] > 0
                   and _stratum(records[i]['x']) == face]
        if not members:
            continue
        sub = [members[k] for k in front_indices(values[members])]
        # The face-front design closest (relative COP) to the overall front.
        gap = [values[i, 0]/max(envelope(i), 1e-12) for i in sub]
        choice = sub[int(np.argmax(gap))]
        if choice not in selected:
            selected.append(choice)
    for i in (front[int(np.argmax(values[front, 0]))], front[int(np.argmax(values[front, 1]))]):
        if i not in selected and len(selected) < count:
            selected.append(i)
    while len(selected) < count:
        remaining = [i for i in front if i not in selected]
        if not remaining:
            break
        d = cdist(norm(remaining), norm(selected)).min(1)
        selected.append(remaining[int(np.argmax(d))])
    return sorted(selected, key=lambda i: values[i, 1])


class LocalParetoSearch(ParetoBOSearch):
    """Cycle through fixed small boxes, fitting on nearby successes or all data.

    The complete observed front always defines hypervolume improvement, and all
    failures remain available to the distance-based feasibility heuristic.
    """

    def __init__(self, baseline, subset_size, seed=20260911, radius=.05, *, fixed_comp_frac=None,
                 max_lift=None, components=2, min_fraction=0.):
        if max_lift is None:
            max_lift = max_lift_for({})
        super().__init__(max_lift, seed, fixed_comp_frac=fixed_comp_frac, components=components,
                         min_fraction=min_fraction)
        self.subset_size = subset_size
        self.radius = radius
        self.anchor_indices = choose_anchors(baseline['records'])
        self.centres = np.asarray([self.encode(baseline['records'][i]['x']) for i in self.anchor_indices])
        for r in baseline['records']:
            self.tell(r['x'], r['COP'], r['VHC_MJ_per_m3'], r['success'] and not r.get('validation_rejected'),
                      r.get('constraint_margins'))
        self.baseline_count = len(self.x)
        self.states = {}
        self.global_seed_state = (baseline.get('kernel_state'), baseline.get('mean_state'))
        self.last_diagnostic = {}

    def training_indices(self, centre):
        good = np.flatnonzero(self.success)
        if self.subset_size is None:
            return good
        distances = np.linalg.norm(np.asarray(self.x)[good]-centre, axis=1)
        order = np.lexsort((good, distances))
        return good[order[:self.subset_size]]

    def local_pool(self, lower, upper, seed):
        unit = qmc.Sobol(self.dimension, scramble=True, seed=seed).random_base2(LOCAL_POOL_LOG2)
        parts = [unit]
        for axis in range(self.dimension):
            for side in (0., 1.):
                edge = unit[:LOCAL_EDGE_POINTS].copy()
                edge[:, axis] = side
                parts.append(edge)
        points = lower+(upper-lower)*np.vstack(parts)
        if self.components > 2:
            # Composition edge with the last component absent, where it crosses this box.
            edge = self.simplex_edge(points[:LOCAL_EDGE_POINTS])
            points = np.vstack((points, edge[((edge >= lower-1e-12) & (edge <= upper+1e-12)).all(1)]))
        pool = np.unique(self._admissible(points), axis=0)
        distance = cdist(pool, self.occupied()).min(1)
        keep = distance > 1e-7
        return pool[keep], distance[keep]

    @staticmethod
    def scores(acquisition, pool):
        with torch.no_grad():
            values = [acquisition(t.unsqueeze(-2)) for t in
                      torch.tensor(pool, dtype=torch.double).split(64)]
        return torch.cat(values).numpy()

    def ask(self):
        started = time.perf_counter()
        step = len(self.x)+len(self.pending)-self.baseline_count
        region = step % len(self.centres)
        centre = self.centres[region]
        lower, upper = np.maximum(centre-self.radius, 0.), np.minimum(centre+self.radius, 1.)
        seed = self.seed+104729*step
        pool, distance = self.local_pool(lower, upper, seed)
        diag = dict(step=step, region=region, anchor_evaluation=self.anchor_indices[region]+1,
                    centre=centre.tolist(), bounds=[lower.tolist(), upper.tolist()],
                    candidate_pool=len(pool), observed_front_size=len(front_indices(self.y)))
        if step % 8 == 7:
            chosen = pool[np.argmax(distance)]
            diag.update(exploration=True, training_points=0, proposal_seconds=time.perf_counter()-started)
            self.last_diagnostic = diag
            return self.decode(chosen)
        indices = self.training_indices(centre)
        X = torch.tensor(np.asarray(self.x)[indices], dtype=torch.double)
        Y = torch.tensor(np.asarray(self.y)[indices], dtype=torch.double)
        cache_key = 'global' if self.subset_size is None else region
        old = self.states.get(cache_key)
        with torch.random.fork_rng(devices=[]), warnings.catch_warnings(record=True) as caught:
            torch.manual_seed(seed)
            model = SingleTaskGP(X, Y, train_Yvar=torch.full_like(Y, 1e-8)).double()
            initial = (old['kernel'], old['mean']) if old else self.global_seed_state
            if initial[0] is not None:
                with torch.no_grad():
                    for module, saved in ((model.covar_module, initial[0]), (model.mean_module, initial[1])):
                        for name, value in saved.items():
                            module.get_parameter(name).copy_(torch.tensor(value, dtype=torch.double))
            visits = old['visits'] if old else 0
            fit = old is None or visits % 5 == 0
            fit_start = time.perf_counter()
            if fit:
                try:
                    fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model), max_attempts=1,
                                     optimizer_kwargs={'options': {'maxiter': 100}})
                except RuntimeError as exc:
                    diag['fit_error'] = str(exc)
            diag['fit_seconds'] = time.perf_counter()-fit_start
            diag['hyperparameters_fitted'] = fit
            self.states[cache_key] = dict(visits=visits+1,
                kernel={k:v.detach().tolist() for k,v in model.covar_module.named_parameters()},
                mean={k:v.detach().tolist() for k,v in model.mean_module.named_parameters()})
            model.eval()
            # Crucial: use the full observed front, never just the training subset.
            observed_front = torch.tensor(np.asarray(self.y)[front_indices(self.y)], dtype=torch.double)
            # Constraint GP on the solved designs nearest this region's centre.
            solved = np.flatnonzero([m is not None for m in self.margins])
            if self.subset_size is not None and len(solved) > self.subset_size:
                near = np.linalg.norm(np.asarray(self.x)[solved]-centre, axis=1)
                solved = solved[np.lexsort((solved, near))[:self.subset_size]]
            constraint_start = time.perf_counter()
            constraint = self.constraint_model(set(solved.tolist()), key=cache_key)
            diag['constraint_seconds'] = time.perf_counter()-constraint_start
            diag['constraint_model'] = constraint is not None
            acquisition = FeasibilityWeightedLogEHVI(model, observed_front,
                torch.tensor(self.x, dtype=torch.double), torch.tensor(self.success, dtype=torch.bool), seed,
                numerical_failure=self.numerical_failures(), constraint=constraint,
                mc_samples=LOCAL_MC_SAMPLES)
            score_start = time.perf_counter()
            scores = self.scores(acquisition, pool)+self.pending_penalty(pool)
            diag['pool_score_seconds'] = time.perf_counter()-score_start
            if np.isfinite(scores).any():
                index = int(np.nanargmax(scores))
                chosen, score = pool[index], scores[index]
            else:
                chosen, score = pool[np.argmax(distance)], -np.inf
                diag['nonfinite_pool'] = True
            optimize_start = time.perf_counter()
            try:
                proposals, _ = optimize_acqf(acquisition, bounds=torch.tensor(np.vstack((lower, upper)), dtype=torch.double),
                    q=1, num_restarts=LOCAL_RESTARTS, raw_samples=LOCAL_RAW_SAMPLES,
                    inequality_constraints=self.inequality_constraints(),
                    options={'maxiter':60, 'batch_limit':1, 'init_batch_limit':8},
                    return_best_only=False, retry_on_optimization_warning=False)
                proposals = self._admissible(proposals.detach().numpy().reshape(-1,self.dimension))
                proposals = proposals[((proposals >= lower) & (proposals <= upper)).all(1)]
                if len(proposals):
                    proposals = proposals[cdist(proposals, self.occupied()).min(1) > 1e-7]
                if len(proposals):
                    more = self.scores(acquisition, proposals)+self.pending_penalty(proposals)
                    if np.isfinite(more).any() and np.nanmax(more) > score:
                        chosen, score = proposals[np.nanargmax(more)], float(np.nanmax(more))
            except (ValueError, RuntimeError) as exc:
                diag['acquisition_error'] = str(exc)
            diag['acquisition_opt_seconds'] = time.perf_counter()-optimize_start
            diag['warnings'] = list(dict.fromkeys(str(w.message) for w in caught))
            with torch.no_grad():
                posterior = model.posterior(torch.tensor(chosen[None], dtype=torch.double))
                diag['predicted_objectives'] = posterior.mean.squeeze(0).tolist()
                diag['predicted_sd'] = posterior.variance.sqrt().squeeze(0).tolist()
            diag['selected_log_ehvi_weighted'] = float(score) if np.isfinite(score) else None
        diag.update(training_points=len(indices), training_indices=indices.tolist(),
                    training_radius_max=float(np.linalg.norm(np.asarray(self.x)[indices]-centre, axis=1).max()),
                    proposal_seconds=time.perf_counter()-started)
        self.last_diagnostic = diag
        return self.decode(chosen)
