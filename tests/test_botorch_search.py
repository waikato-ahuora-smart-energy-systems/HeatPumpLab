"""Check BoTorch's objective sign, gradients, failure handling and reproducibility."""

import numpy as np
import pytest
import torch

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.origin_search import OriginSurrogateSearch
from rhp_botorch_search import BotorchLogEISearch, FeasibilityWeightedLogEI


def score(x):
    u = BotorchLogEISearch.encode(x)
    return float((u[0] + u[1] - 0.85)**2 + 2*(u[0] - u[1] - 0.35)**2 + (u[2] - 0.8)**2)


def test_bootstrap_matches_custom_bo():
    a = BotorchLogEISearch(max_lift=123.333, seed=7)
    b = OriginSurrogateSearch(method="gp_ei", max_lift=123.333, seed=7, origin=(0., 0., 0.5))
    for _ in range(5):
        x = a.ask()
        assert np.array_equal(x, b.ask())
        a.tell(x, score(x), True)
        b.tell(x, score(x), True)


def test_botorch_minimises_sec_and_replays_with_failure():
    a = BotorchLogEISearch(max_lift=123.333, seed=42)
    history = []
    for i in range(24):
        x = a.ask()
        assert x[0] + x[1] <= 123.333 + 1e-8
        assert np.all(a.encode(x) >= 0) and np.all(a.encode(x) <= 1)
        # Force one failed observation to exercise the feasibility adapter.
        success = i != 6
        value = score(x) if success else 1e6
        history.append((x, value, success))
        a.tell(x, value, success)
    assert min(v for _, v, ok in history[5:] if ok) < 0.25 * min(v for _, v, _ in history[:5])
    assert np.isnan(a.y[6])
    b = BotorchLogEISearch(max_lift=123.333, seed=42)
    for x, value, success in history:
        assert np.allclose(b.ask(), x, rtol=0, atol=1e-9)
        b.tell(x, value, success)
    assert np.allclose(a.ask(), b.ask(), rtol=0, atol=1e-9)


def test_failure_only_history_remains_feasible_and_novel():
    a = BotorchLogEISearch(max_lift=90, seed=7)
    for _ in range(20):
        x = a.ask()
        a.tell(x, 1e6, False)
    assert not any(a.success)
    assert len({tuple(x) for x in a.x}) == 20


def test_log_feasibility_matches_existing_formula_and_has_finite_gradient():
    a = BotorchLogEISearch(max_lift=123.333, seed=7)
    for _ in range(6):
        x = a.ask()
        a.tell(x, score(x), True)
    x = a.ask()
    a.tell(x, 1e6, False)
    model, best = a._fit_model()
    acq = FeasibilityWeightedLogEI(model, best, torch.tensor(np.array(a.x), dtype=torch.double),
                                  torch.tensor(a.success))
    point = torch.tensor([[[0.2, 0.3, 0.4]]], dtype=torch.double, requires_grad=True)
    result = acq(point)
    delta = float((result-acq.log_ei(point)).detach())
    distance = np.linalg.norm(np.asarray(a.x)-[0.2, 0.3, 0.4], axis=1)
    expected = np.log((distance[-1]+0.02)/(distance[:-1].min()+distance[-1]+0.04))
    assert delta == pytest.approx(expected)
    result.sum().backward()
    assert torch.isfinite(point.grad).all()


def test_gp_posterior_uses_negative_sec_with_small_known_noise():
    a = BotorchLogEISearch(max_lift=123.333, seed=7)
    for _ in range(5):
        x = a.ask()
        a.tell(x, score(x), True)
    model, _ = a._fit_model()
    with torch.no_grad():
        predicted = model.posterior(torch.tensor(np.asarray(a.x), dtype=torch.double)).mean.squeeze(-1).numpy()
    expected = -(np.asarray(a.y)-np.mean(a.y))/np.std(a.y)
    assert np.allclose(predicted, expected, atol=1e-4, rtol=0)
    assert np.argmax(predicted) == np.argmin(a.y)
