"""Experimental COBYQA and TuRBO-1 searches; thermodynamic model is untouched.

TuRBO-1 follows Eriksson et al. (2019), with a SciPy exact Matérn GP and
joint Thompson sampling. Explicit domain rejection and failed-solve handling
are adaptations for HeatPumpLab. No upstream TuRBO code is vendored.
"""

import numpy as np
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import Bounds, LinearConstraint, minimize
from scipy.spatial.distance import cdist
from scipy.stats import qmc

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.surrogate_search import (
    JointSurrogateSearch, _matern52,
)


class TurboOne(JointSurrogateSearch):
    """Sequential single-region TuRBO with restart-local GP observations."""

    def __init__(self, *, max_lift, seed):
        # Initialise directly: the parent's bootstrap consumes a different RNG.
        self.max_lift = float(max_lift)
        self.rng = np.random.default_rng(seed)
        self.x, self.y, self.success = [], [], []
        self.restarts, self.events = [], []
        self._restart(first=True)

    def _restart(self, first=False):
        self.start = len(self.x)
        self.length = 0.8
        self.success_count = self.failure_count = 0
        self.hyper = np.log([0.5, 0.5, 0.5, 1.0, 0.005])
        self.initial = [np.array([0., 0., 0.5])] if first else []
        # Rejection preserves the box Latin-hypercube proposal distribution;
        # the accepted constrained subset is no longer a strict LHS.
        while len(self.initial) < 6:
            pool = qmc.LatinHypercube(3, seed=int(self.rng.integers(2**31))).random(12)
            for p in self._admissible(pool):
                if self.x and cdist([p], self.x).min() <= 1e-7:
                    continue
                self.initial.append(p)
                if len(self.initial) == 6:
                    break
        if not first:
            self.restarts.append(self.start)

    def tell(self, x, objective, success):
        old = np.asarray(self.y[self.start:])
        old_best = np.nanmin(old) if np.isfinite(old).any() else np.inf
        adaptive = len(old) >= 6
        super().tell(x, objective, success)
        if adaptive:
            improved = self.success[-1] and (
                not np.isfinite(old_best) or objective < old_best - 1e-3 * abs(old_best))
            self.success_count = self.success_count + 1 if improved else 0
            self.failure_count = 0 if improved else self.failure_count + 1
            before = self.length
            if self.success_count >= 3:
                self.length = min(2 * self.length, 1.6)
                self.success_count = 0
            elif self.failure_count >= 4:
                self.length /= 2
                self.failure_count = 0
            if before != self.length:
                self.events.append(dict(evaluation=len(self.x), length=self.length))
            if self.length < 0.5**7:
                self._restart()

    def _fit(self, x, y):
        def components(theta):
            lengths, amplitude, noise = np.exp(theta[:3]), np.exp(theta[3]), np.exp(theta[4])
            k = amplitude * _matern52(x, x, lengths)
            k.flat[::len(x) + 1] += noise + 1e-8
            chol = cho_factor(k, lower=True)
            ones = np.ones(len(x))
            ki_one = cho_solve(chol, ones)
            mean = float(ki_one @ y / (ones @ ki_one))
            alpha = cho_solve(chol, y - mean)
            return chol, mean, alpha

        def nll(theta):
            try:
                chol, mean, alpha = components(theta)
                return 0.5 * (y - mean) @ alpha + np.log(np.diag(chol[0])).sum()
            except (ValueError, np.linalg.LinAlgError):
                return 1e20

        bounds = [(np.log(0.005), np.log(2.))] * 3 + [
            (np.log(0.05), np.log(20.)), (np.log(5e-4), np.log(0.2))]
        fit = minimize(nll, self.hyper, method="L-BFGS-B", bounds=bounds,
                       options={"maxiter": 50})
        if np.isfinite(fit.fun) and fit.fun < 1e19:
            self.hyper = fit.x
        return components(self.hyper)

    @staticmethod
    def _posterior(x, candidates, y_mean, alpha, chol, hyper):
        lengths, amplitude, noise = np.exp(hyper[:3]), np.exp(hyper[3]), np.exp(hyper[4])
        cross = amplitude * _matern52(x, candidates, lengths)
        mean = y_mean + cross.T @ alpha
        v = solve_triangular(chol[0], cross, lower=True)
        covariance = amplitude * _matern52(candidates, candidates, lengths) - v.T @ v
        covariance = (covariance + covariance.T) / 2
        # Sample the noisy predictive distribution, as in the TuRBO reference.
        covariance.flat[::len(candidates) + 1] += noise + 1e-8
        return mean, covariance

    def ask(self):
        if len(self.x) - self.start < 6:
            return self.decode(self.initial[len(self.x) - self.start])
        local_x = np.asarray(self.x[self.start:])
        good = np.asarray(self.success[self.start:])
        if good.sum() < 3:
            pool = self._admissible(self.rng.random((1024, 3)))
            return self.decode(pool[np.argmax(cdist(pool, self.x).min(axis=1))])
        x = local_x[good]
        raw_y = np.asarray(self.y[self.start:])[good]
        spread = raw_y.std()
        y = (raw_y - np.median(raw_y)) / (spread if spread >= 1e-6 else 1.)
        center = x[np.argmin(raw_y)]
        chol, y_mean, alpha = self._fit(x, y)
        weights = np.exp(self.hyper[:3])
        weights /= np.exp(np.log(weights).mean())
        lower = np.maximum(0., center - weights * self.length / 2)
        upper = np.minimum(1., center + weights * self.length / 2)
        pool = qmc.Sobol(3, scramble=True, seed=int(self.rng.integers(2**31))).random_base2(10)
        candidates = self._admissible(lower + pool * (upper - lower))
        candidates = candidates[cdist(candidates, self.x).min(axis=1) > 1e-7][:300]
        if not len(candidates):
            raise RuntimeError("Empty TuRBO trust region")
        mean, covariance = self._posterior(x, candidates, y_mean, alpha, chol, self.hyper)
        sample = mean + np.linalg.cholesky(covariance) @ self.rng.standard_normal(len(candidates))
        # Failed solves are not SEC data. Use the same distance heuristic as BO,
        # applied additively to standardised Thompson scores (EI is multiplicative).
        all_good = np.asarray(self.success)
        if not all_good.all():
            distances = cdist(candidates, self.x)
            dg = distances[:, all_good].min(axis=1)
            db = distances[:, ~all_good].min(axis=1)
            feasibility = (db + 0.02) / (dg + db + 0.04)
            sample += 1 - feasibility
        return self.decode(candidates[np.argmin(sample)])


class _StartBudgetReached(Exception):
    pass


def multistart_cobyqa(evaluate, *, max_lift, seed, budget, per_start=80):
    """Minimise SEC with up to 80 unique evaluations per local start.

    evaluate(physical_x) returns a serialised outcome. Failures have a fixed
    scalar penalty of 1e6. Duplicate and cheap domain checks do not spend the
    expensive budget. Early convergence releases budget to additional starts.
    """
    rng = np.random.default_rng(seed)
    cache, starts = {}, []
    rejected = 0
    while len(cache) < budget:
        if not starts:
            u0 = np.array([0., 0., 0.5])
        else:
            pool = qmc.Sobol(3, scramble=True, seed=int(rng.integers(2**31))).random_base2(10)
            pool = pool[(pool[:, 0] + pool[:, 1]) * 100 <= max_lift]
            previous = JointSurrogateSearch.encode(np.asarray(list(cache)))
            u0 = pool[np.argmax(cdist(pool, previous).min(axis=1))]
        start_count = len(cache)
        limit = min(budget, start_count + per_start)

        def fun(u):
            nonlocal rejected
            if np.any(u < -1e-10) or np.any(u > 1 + 1e-10) or (u[0] + u[1]) * 100 > max_lift + 1e-8:
                rejected += 1
                return 1e6 + 1e6 * max(0., (u[0] + u[1]) * 100 - max_lift)
            physical = JointSurrogateSearch.decode(np.clip(u, 0, 1))
            key = tuple(np.round(physical, 8))
            if key in cache:
                return cache[key]
            if len(cache) >= limit:
                raise _StartBudgetReached
            item = evaluate(physical)
            cache[key] = float(item["SEC"]) if item["success"] else 1e6
            return cache[key]

        status, message = None, "per-start evaluation allowance reached"
        try:
            result = minimize(fun, u0, method="COBYQA", bounds=Bounds(np.zeros(3), np.ones(3)),
                              constraints=LinearConstraint([[1., 1., 0.]], -np.inf, max_lift / 100),
                              options=dict(maxfev=per_start * 5, maxiter=per_start * 10,
                                           initial_tr_radius=0.15, final_tr_radius=1e-5,
                                           feasibility_tol=1e-10, scale=False))
            status, message = int(result.status), str(result.message)
        except _StartBudgetReached:
            pass
        starts.append(dict(start_evaluation=start_count, end_evaluation=len(cache),
                           x0=JointSurrogateSearch.decode(u0).tolist(), status=status, message=message))
        if len(cache) == start_count:
            raise RuntimeError("COBYQA restart consumed no new evaluations")
    return dict(starts=starts, cheap_rejected_proposals=rejected)
