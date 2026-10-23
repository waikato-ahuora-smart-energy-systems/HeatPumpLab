"""Experimental fixed-composition, joint three-variable surrogate searches.

These pilots deliberately leave the production v11 optimiser unchanged. Only
successful SEC observations enter the response model; failed solves are handled
by a separate distance-based feasibility heuristic (not a calibrated classifier).
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RBFInterpolator
from scipy.linalg import cho_factor, cho_solve, solve_triangular
from scipy.optimize import minimize
from scipy.spatial.distance import cdist
from scipy.special import ndtr

from ._types import COMP2_POWER_FRAC_BOUNDS

_LOWER = np.array([0.0, 0.0, COMP2_POWER_FRAC_BOUNDS[0]])
_SPAN = np.array([100.0, 100.0, COMP2_POWER_FRAC_BOUNDS[1]-COMP2_POWER_FRAC_BOUNDS[0]])


class JointSurrogateSearch:
    """Ask/tell search in normalised physical coordinates, with explicit edges."""

    METHODS = ("rbf", "gp_ei", "space_filling")

    def __init__(self, *, method, max_lift, seed=20260908):
        if method not in self.METHODS:
            raise ValueError(f"Unknown method: {method}")
        if not np.isfinite(max_lift) or max_lift <= 0:
            raise ValueError("max_lift must be finite and positive")
        self.method = method
        self.max_lift = float(max_lift)
        self.rng = np.random.default_rng(seed)
        self.x = []
        self.y = []
        self.success = []
        self.lengthscales = np.full(3, 0.3)
        self.initial = self._initial_design()

    @staticmethod
    def decode(u):
        return np.asarray(u) * _SPAN + _LOWER

    @staticmethod
    def encode(x):
        return (np.asarray(x) - _LOWER) / _SPAN

    def _admissible(self, u):
        u = np.asarray(u)
        return u[(u[:, 0] + u[:, 1]) * 100 <= self.max_lift + 1e-10]

    def _initial_design(self):
        points = [[0, 0, s] for s in (0, 0.5, 1)]
        for lift in (0.5, 0.75):
            lift = min(lift, self.max_lift / 100)
            for s in (0, 0.5, 1):
                points.extend(([lift, 0, s], [0, lift, s]))
        points = np.unique(points, axis=0).tolist()
        # Identical initial design and random stream for every method.
        while len(points) < 24:
            p = self.rng.random(3)
            if (p[0] + p[1]) * 100 <= self.max_lift:
                points.append(p.tolist())
        return points

    def tell(self, x, objective, success):
        u = self.encode(x)
        if u.shape != (3,) or not np.all(np.isfinite(u)):
            raise ValueError("Expected three finite design coordinates")
        if np.any(u < -1e-10) or np.any(u > 1 + 1e-10):
            raise ValueError("Design outside bounds")
        if (u[0] + u[1]) * 100 > self.max_lift + 1e-8:
            raise ValueError("Design exceeds combined lift limit")
        if self.x and np.min(cdist([u], self.x)) < 1e-8:
            raise ValueError("Duplicate observation")
        self.x.append(u)
        self.success.append(bool(success and np.isfinite(objective)))
        self.y.append(float(objective) if self.success[-1] else np.nan)

    def _candidates(self):
        global_points = self.rng.random((2048, 3))
        # Boundary optima are common; a continuous random design misses them.
        for axis in range(3):
            start = axis * 256
            global_points[start:start + 128, axis] = 0
            global_points[start + 128:start + 256, axis] = 1
        pools = [global_points]
        good = np.flatnonzero(self.success)
        if len(good):
            order = good[np.argsort(np.asarray(self.y)[good])[:3]]
            for index in order:
                for scale in (0.02, 0.08, 0.25):
                    pools.append(np.clip(
                        self.x[index] + self.rng.normal(0, scale, (128, 3)), 0, 1,
                    ))
        candidates = np.unique(self._admissible(np.vstack(pools)), axis=0)
        distances = cdist(candidates, self.x)
        keep = distances.min(axis=1) > 1e-7
        return candidates[keep], distances[keep]

    def ask(self):
        if len(self.x) < len(self.initial):
            return self.decode(self.initial[len(self.x)])
        candidates, distances = self._candidates()
        nearest = distances.min(axis=1)
        good = np.asarray(self.success)
        # A maximin step every eighth adaptive call preserves global exploration.
        if (self.method == "space_filling" or good.sum() < 5
                or (len(self.x) - len(self.initial)) % 8 == 7):
            return self.decode(candidates[np.argmax(nearest)])
        train_x = np.asarray(self.x)[good]
        train_y = np.asarray(self.y)[good]
        y = (train_y - train_y.mean()) / max(train_y.std(), 1e-8)
        feasibility = np.ones(len(candidates))
        if not good.all():
            d_good = distances[:, good].min(axis=1)
            d_bad = distances[:, ~good].min(axis=1)
            feasibility = (d_bad + 0.02) / (d_good + d_bad + 0.04)
        try:
            if self.method == "rbf":
                model = RBFInterpolator(train_x, y, kernel="cubic", smoothing=1e-6)
                prediction = model(candidates)
                weight = (0.3, 0.5, 0.8, 0.95)[len(self.x) % 4]
                score = weight * _unit_range(prediction) - (1 - weight) * _unit_range(nearest)
                score += 1 - feasibility
                index = np.argmin(score)
            else:
                mean, std = self._gp_predict(train_x, y, candidates)
                improvement = y.min() - mean
                z = improvement / std
                ei = improvement * ndtr(z) + std * np.exp(-z*z/2) / np.sqrt(2*np.pi)
                index = np.argmax(ei * feasibility)
        except (ValueError, np.linalg.LinAlgError):
            index = np.argmax(nearest)
        return self.decode(candidates[index])

    def _gp_predict(self, x, y, candidates):
        def factor(log_lengths):
            kernel = _matern52(x, x, np.exp(log_lengths))
            kernel.flat[::len(x) + 1] += 1e-6
            return cho_factor(kernel, lower=True)

        def nll(log_lengths):
            try:
                chol = factor(log_lengths)
                return 0.5 * y @ cho_solve(chol, y) + np.log(np.diag(chol[0])).sum()
            except np.linalg.LinAlgError:
                return 1e20

        if len(self.x) == len(self.initial) or len(self.x) % 5 == 0:
            fits = [minimize(nll, np.log(start), method="L-BFGS-B",
                             bounds=[(np.log(0.02), np.log(3))] * 3)
                    for start in (self.lengthscales, np.full(3, 0.15))]
            fit = min(fits, key=lambda item: item.fun)
            if np.isfinite(fit.fun):
                self.lengthscales = np.exp(fit.x)
        chol = factor(np.log(self.lengthscales))
        cross = _matern52(x, candidates, self.lengthscales)
        mean = cross.T @ cho_solve(chol, y)
        v = solve_triangular(chol[0], cross, lower=True)
        std = np.sqrt(np.maximum(1 - np.sum(v*v, axis=0), 1e-12))
        return mean, std


def _unit_range(values):
    return (values - np.min(values)) / max(float(np.ptp(values)), 1e-12)


def _matern52(x, y, lengthscales):
    radius = np.sqrt(5) * cdist(x / lengthscales, y / lengthscales)
    return (1 + radius + radius**2 / 3) * np.exp(-radius)
