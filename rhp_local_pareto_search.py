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
                              fixed_comp_frac=result.get('metadata', {}).get('fixed_comp_frac'))
    for record in result['records'][epoch:]:
        search.tell(record['x'], record['COP'], record['VHC_MJ_per_m3'],
                    record['success'] and not record.get('validation_rejected'),
                    record.get('constraint_margins'))
    saved = result.get('local128_state', {})
    if saved.get('epoch_start') == epoch:
        search.states = {int(key): value for key, value in saved['model_states'].items()}
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


def choose_anchors(records, count=6):
    values = np.asarray([[r['COP'], r['VHC_MJ_per_m3']] if r['success'] and
                         not r.get('validation_rejected') else [0., 0.] for r in records])
    front = [i for i in front_indices(values) if .01 <= records[i]['x'][3] <= .99]
    if len(front) < count:
        raise ValueError('Too few mixed-fluid Pareto designs for the planned regions')
    points = values[front]
    scale = np.maximum(np.ptp(points, axis=0), 1e-12)
    points = (points-points.min(0))/scale
    selected = [int(np.argmax(points[:, 0])), int(np.argmax(points[:, 1]))]
    selected = list(dict.fromkeys(selected))
    while len(selected) < count:
        distances = cdist(points, points[selected]).min(1)
        distances[selected] = -1
        selected.append(int(np.argmax(distances)))
    return sorted([front[i] for i in selected], key=lambda i: values[i, 1])


class LocalParetoSearch(ParetoBOSearch):
    """Cycle through fixed small boxes, fitting on nearby successes or all data.

    The complete observed front always defines hypervolume improvement, and all
    failures remain available to the distance-based feasibility heuristic.
    """

    def __init__(self, baseline, subset_size, seed=20260911, radius=.05, *, fixed_comp_frac=None):
        super().__init__(123.33333333333333, seed, fixed_comp_frac=fixed_comp_frac)
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
        unit = qmc.Sobol(self.dimension, scramble=True, seed=seed).random_base2(9)
        parts = [unit]
        for axis in range(self.dimension):
            for side in (0., 1.):
                edge = unit[:64].copy()
                edge[:, axis] = side
                parts.append(edge)
        pool = np.unique(self._admissible(lower+(upper-lower)*np.vstack(parts)), axis=0)
        distance = cdist(pool, self.x).min(1)
        keep = distance > 1e-7
        return pool[keep], distance[keep]

    @staticmethod
    def scores(acquisition, pool):
        with torch.no_grad():
            values = [acquisition(t.unsqueeze(-2)) for t in
                      torch.tensor(pool, dtype=torch.double).split(8)]
        return torch.cat(values).numpy()

    def ask(self):
        started = time.perf_counter()
        step = len(self.x)-self.baseline_count
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
            constraint = self.constraint_model(set(solved.tolist()))
            diag['constraint_model'] = constraint is not None
            acquisition = FeasibilityWeightedLogEHVI(model, observed_front,
                torch.tensor(self.x, dtype=torch.double), torch.tensor(self.success, dtype=torch.bool), seed,
                numerical_failure=self.numerical_failures(), constraint=constraint)
            score_start = time.perf_counter()
            scores = self.scores(acquisition, pool)
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
                    q=1, num_restarts=6, raw_samples=128,
                    inequality_constraints=[(torch.tensor([0, 1]), torch.tensor([-1., -1.], dtype=torch.double), -self.max_lift/100)],
                    options={'maxiter':60, 'batch_limit':1, 'init_batch_limit':8},
                    return_best_only=False, retry_on_optimization_warning=False)
                proposals = self._admissible(proposals.detach().numpy().reshape(-1,self.dimension))
                proposals = proposals[((proposals >= lower) & (proposals <= upper)).all(1)]
                if len(proposals):
                    proposals = proposals[cdist(proposals, self.x).min(1) > 1e-7]
                if len(proposals):
                    more = self.scores(acquisition, proposals)
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
