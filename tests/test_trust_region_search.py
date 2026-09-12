"""Budget, restart, covariance and reproducibility checks without REFPROP."""

import numpy as np

from rhp_trust_region_search import TurboOne, multistart_cobyqa


def quadratic(x):
    u = TurboOne.encode(x)
    return float((u[0] + u[1] - 0.85)**2 + 2*(u[0] - u[1] - 0.35)**2 + (u[2] - 0.8)**2)


def test_cobyqa_uses_exact_budget_and_finds_coupled_minimum():
    points, values = [], []

    def evaluate(x):
        assert x[0] + x[1] <= 123.333 + 1e-8
        assert 0 <= x[0] <= 100 and 0 <= x[1] <= 100 and 0.1 <= x[2] <= 0.9
        points.append(tuple(np.round(x, 8)))
        values.append(quadratic(x))
        return dict(success=True, SEC=values[-1])

    info = multistart_cobyqa(evaluate, max_lift=123.333, seed=42, budget=180, per_start=60)
    assert len(points) == len(set(points)) == 180
    assert len(info["starts"]) >= 3
    assert min(values) < 1e-7


def test_cobyqa_failed_solves_count_and_trajectory_reproduces():
    histories = []
    for _ in range(2):
        history = []

        def evaluate(x):
            history.append(x.copy())
            return dict(success=False, SEC=None)

        multistart_cobyqa(evaluate, max_lift=90, seed=7, budget=35, per_start=15)
        histories.append(history)
    assert len(histories[0]) == 35
    assert np.array_equal(histories[0], histories[1])


def test_turbo_improves_and_replays_exactly_with_failed_region():
    search = TurboOne(max_lift=123.333, seed=42)
    history = []
    for _ in range(55):
        x = search.ask()
        assert x[0] + x[1] <= 123.333 + 1e-8
        assert np.all(search.encode(x) >= 0) and np.all(search.encode(x) <= 1)
        success = x[0] < 90
        y = quadratic(x) if success else 1e6
        history.append((x, y, success))
        search.tell(x, y, success)
    assert np.nanmin(search.y[6:]) < 0.1 * np.nanmin(search.y[:6])
    assert all(np.isnan(y) for y, success in zip(search.y, search.success) if not success)
    restored = TurboOne(max_lift=123.333, seed=42)
    for x, y, success in history:
        assert np.array_equal(restored.ask(), x)
        restored.tell(x, y, success)
    assert np.array_equal(restored.ask(), search.ask())


def test_turbo_shrinks_restarts_and_clears_local_history():
    search = TurboOne(max_lift=123.333, seed=3)
    for _ in range(34):
        x = search.ask()
        search.tell(x, 1e6, False)
    assert search.restarts == [34]
    assert search.start == 34 and search.length == 0.8
    assert not any(search.success)
    assert len(search.initial) == 6


def test_turbo_expands_after_three_meaningful_improvements():
    search = TurboOne(max_lift=123.333, seed=3)
    # Exercise controller independently of expensive GP fitting.
    for i in range(9):
        x = search.decode([0.01 * i, 0.1, 0.5])
        search.tell(x, 10. - i, True)
    assert search.length == 1.6


def test_turbo_draws_correlated_joint_posterior():
    search = TurboOne(max_lift=123.333, seed=3)
    x = np.array([[0., 0., 0.], [1., 1., 1.], [0.5, 0.5, 0.5]])
    y = np.array([1., 1., -1.])
    chol, mean, alpha = search._fit(x, y)
    candidates = np.array([[0.25, 0.25, 0.25], [0.25001, 0.25, 0.25]])
    _, covariance = search._posterior(x, candidates, mean, alpha, chol, search.hyper)
    assert covariance[0, 1] > 0
    assert covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1]) > 0.5
    assert np.linalg.eigvalsh(covariance).min() > 0
