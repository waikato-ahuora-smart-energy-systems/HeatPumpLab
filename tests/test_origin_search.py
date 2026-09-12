"""Origin bootstrap, budget-length behaviour and reproducibility."""

import numpy as np
import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.origin_search import (
    OriginSurrogateSearch,
)
from main_RHPorigin_benchmark import make_search


def test_cma_stagnation_restart_is_replayable():
    pytest.importorskip("cma")
    search = make_search("cma_es", 123.3333, 29)
    restored = make_search("cma_es", 123.3333, 29)
    for i in range(100):
        x = search.ask()
        assert np.allclose(restored.ask(), x, atol=1e-9, rtol=0)
        for solver in (search, restored):
            solver.tell(x, 1.0 - i * 1e-8, True)
    assert any("budget_stagnation" in r["reasons"] for r in search.restart_reasons)
    assert search.restart_reasons == restored.restart_reasons


@pytest.mark.parametrize("method", ("gp_ei", "rbf", "cma_es"))
def test_origin_search_is_feasible_and_improves_coupled_objective(method):
    if method == "cma_es":
        pytest.importorskip("cma")
    search = make_search(method, 123.3333333, 19)
    history = []
    for _ in range(200):
        x = search.ask()
        assert np.all(x >= [0, 0, 0.1]) and np.all(x <= [100, 100, 0.9])
        assert x[0] + x[1] <= 123.3333333 + 1e-9
        assert not any(np.allclose(x, old, atol=1e-8, rtol=0) for old, _ in history)
        u = (x - [0, 0, 0.1]) / [100, 100, 0.8]
        y = np.sum((u - [0.6, 0.2, 0.8]) ** 2)
        history.append((x, y))
        search.tell(x, y, True)
    assert np.array_equal(history[0][0], [0, 0, 0.5])
    assert min(y for _, y in history) < 1e-3


@pytest.mark.parametrize("method", ("gp_ei", "rbf", "cma_es"))
def test_origin_search_restarts_deterministically_from_checkpoint_history(method):
    if method == "cma_es":
        pytest.importorskip("cma")
    original = make_search(method, 123.3333333, 29)
    observations = []
    for _ in range(23):
        x = original.ask()
        success = x[0] < 75
        y = float(sum(x * x)) if success else 1e6
        original.tell(x, y, success)
        observations.append((x, y, success))
    restored = make_search(method, 123.3333333, 29)
    for x, y, success in observations:
        assert np.allclose(restored.ask(), x, atol=1e-9, rtol=0)
        restored.tell(x, y, success)
    assert np.allclose(restored.ask(), original.ask(), atol=1e-9, rtol=0)


def test_surrogate_bootstrap_is_local_and_has_full_affine_rank():
    a = OriginSurrogateSearch(
        method="rbf", origin=[0, 0, 0.5], max_lift=123.333, seed=1
    )
    b = OriginSurrogateSearch(
        method="gp_ei", origin=[0, 0, 0.5], max_lift=123.333, seed=2
    )
    assert a.initial == b.initial
    assert len(a.initial) == 5
    offsets = np.asarray(a.initial) - a.encode([0, 0, 0.5])
    assert np.max(np.abs(offsets)) <= 0.05 + 1e-12
    assert np.linalg.matrix_rank(offsets) == 3
