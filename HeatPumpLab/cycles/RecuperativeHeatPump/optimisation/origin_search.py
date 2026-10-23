"""Origin-start search variants, separate from the original cold-start pilot."""

import numpy as np

from .surrogate_search import JointSurrogateSearch


class OriginSurrogateSearch(JointSurrogateSearch):
    """Bootstrap a joint surrogate using only the origin and nearby probes.

    The five initial evaluations are part of the search budget. They give the
    response surface an affine-spanning design without a global prescreen.
    Subsequent adaptive proposals may explore the complete feasible domain.
    """

    def __init__(self, *, method, max_lift, origin, seed):
        self.origin = np.asarray(origin, dtype=float)
        super().__init__(method=method, max_lift=max_lift, seed=seed)

    def _initial_design(self):
        u = self.encode(self.origin)
        if u.shape != (3,) or np.any(~np.isfinite(u)):
            raise ValueError("Origin must contain three finite coordinates")
        if np.any(u < 0) or np.any(u > 1) or 100 * sum(u[:2]) > self.max_lift:
            raise ValueError("Origin is outside the feasible domain")
        points = [u.copy()]
        for axis in range(3):
            p = u.copy()
            p[axis] += 0.05 if u[axis] <= 0.5 else -0.05
            points.append(p)
        p = u.copy()
        p[2] += -0.05 if 0.05 <= u[2] <= 0.95 else (0.10 if u[2] < 0.05 else -0.10)
        points.append(p)
        points = np.asarray(points)
        if len(self._admissible(points)) != len(points):
            raise ValueError("Local initial design violates the combined lift limit")
        return points.tolist()


class OriginCMAES:
    """Official pycma with box transforms and rejection of the known lift constraint.

    Evaluate the origin first, then sample generations centred initially there.
    Native stopping triggers a fresh origin-centred restart until the external
    evaluation budget is spent. An incomplete final generation is still scored
    by the caller, but is not passed to CMA-ES's distribution update.
    """

    def __init__(self, *, max_lift, origin, seed):
        self.origin = np.asarray(origin, dtype=float)
        self.max_lift = max_lift
        self.seed = seed
        self.x = []
        self.pending = None
        self.generation_x = []
        self.generation_y = []
        self.restart_reasons = []
        self.rejected_candidates = 0
        self._start()

    def _start(self):
        import cma

        self.generation_best = []
        rng = np.random.RandomState(self.seed + 100003 * len(self.restart_reasons))
        self.es = cma.CMAEvolutionStrategy(
            JointSurrogateSearch.encode(self.origin), 0.25,
            dict(bounds=[0, 1], seed=self.seed, randn=rng.randn,
                 verbose=-9, verb_log=0, popsize=7),
        )

    def ask(self):
        if self.pending is not None:
            raise RuntimeError("tell must follow each ask")
        if not self.x:
            self.pending = self.origin.copy()
            return self.pending.copy()
        for _ in range(10000):
            u = self.es.ask(1)[0]
            if 100 * sum(u[:2]) > self.max_lift:
                self.rejected_candidates += 1
                continue
            x = JointSurrogateSearch.decode(u)
            if any(np.max(np.abs(x - previous)) < 1e-6 for previous in self.x):
                self.rejected_candidates += 1
                continue
            self.pending = x
            return x.copy()
        raise RuntimeError("CMA-ES could not propose a unique feasible candidate")

    def tell(self, x, objective, success):
        if self.pending is None or not np.allclose(x, self.pending, atol=1e-9, rtol=0):
            raise ValueError("Observation does not match the pending candidate")
        self.x.append(np.asarray(x).copy())
        self.pending = None
        if len(self.x) == 1:
            return
        self.generation_x.append(JointSurrogateSearch.encode(x))
        self.generation_y.append(float(objective) if success and np.isfinite(objective) else 1e6)
        if len(self.generation_x) == self.es.popsize:
            self.es.tell(self.generation_x, self.generation_y)
            self.generation_best.append(min(self.generation_y))
            self.generation_x = []
            self.generation_y = []
            stop = self.es.stop()
            # Native stopping can spend the remaining allowance polishing a
            # weak basin. Restart after ten generations without useful gain.
            if not stop and len(self.generation_best) >= 10:
                previous = min(self.generation_best[:-9])
                recent = min(self.generation_best[-9:])
                if previous - recent <= 1e-4 * max(1., abs(previous)):
                    stop = {"budget_stagnation": "10 generations; relative gain <= 1e-4"}
            if stop:
                self.restart_reasons.append({"evaluations": len(self.x), "reasons": str(stop)})
                self._start()
