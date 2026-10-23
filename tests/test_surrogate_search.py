"""Behavioural checks for the experimental joint surrogate search."""

import numpy as np
import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.surrogate_search import JointSurrogateSearch


@pytest.mark.parametrize("method", JointSurrogateSearch.METHODS)
def test_joint_search_respects_domain_and_does_not_repeat(method):
    search = JointSurrogateSearch(method=method, max_lift=123.333, seed=7)
    points = []
    for _ in range(32):
        x = search.ask()
        assert 0 <= x[0] <= 100 and 0 <= x[1] <= 100
        assert x[0] + x[1] <= 123.333 + 1e-9
        assert 0.1 <= x[2] <= 0.9
        assert not any(np.allclose(x, p, atol=1e-8, rtol=0) for p in points)
        points.append(x)
        # A failed region must not contaminate the smooth SEC response model.
        success = x[0] < 90
        search.tell(x, np.sum((search.encode(x) - [0.6, 0.2, 0.8])**2)
                    if success else 1e6, success)
    assert np.isnan(np.asarray(search.y)[~np.asarray(search.success)]).all()


@pytest.mark.parametrize("method", ("rbf", "gp_ei"))
def test_surrogate_improves_coupled_interior_solution(method):
    search = JointSurrogateSearch(method=method, max_lift=123.333, seed=42)
    scores = []
    for _ in range(48):
        x = search.ask()
        u = search.encode(x)
        # Optimum requires both lifts and a non-default compressor split.
        y = (u[0] + u[1] - 0.85)**2 + 2*(u[0] - u[1] - 0.35)**2 + (u[2] - 0.8)**2
        scores.append(y)
        search.tell(x, y, True)
    assert min(scores[24:]) < 0.25 * min(scores[:24])


@pytest.mark.parametrize("method", JointSurrogateSearch.METHODS)
def test_failure_only_history_still_explores(method):
    search = JointSurrogateSearch(method=method, max_lift=90, seed=3)
    for _ in range(27):
        x = search.ask()
        search.tell(x, 1e6, False)
    assert not any(search.success)


def test_shared_initial_design_and_reproducible_resume():
    searches = [JointSurrogateSearch(method=m, max_lift=123.333) for m in JointSurrogateSearch.METHODS]
    for _ in range(24):
        xs = [s.ask() for s in searches]
        assert all(np.array_equal(xs[0], x) for x in xs)
        for search, x in zip(searches, xs):
            search.tell(x, np.sum(search.encode(x)**2), True)
    original = searches[0]
    for _ in range(3):
        x = original.ask()
        original.tell(x, np.sum(original.encode(x)**2), True)
    restored = JointSurrogateSearch(method="rbf", max_lift=123.333)
    for u, y, success in zip(original.x, original.y, original.success):
        assert np.allclose(restored.ask(), restored.decode(u))
        restored.tell(restored.decode(u), y, success)
    assert np.array_equal(restored.ask(), original.ask())
