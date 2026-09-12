"""Interrupted benchmark runs must resume without spending saved evaluations."""

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

import main_RHPbotorch_benchmark as benchmark


@pytest.mark.parametrize("method", benchmark.METHODS)
def test_checkpoint_resume_matches_uninterrupted_run(monkeypatch, tmp_path, method):
    state = dict(interrupt=True, search_calls=0, replay_calls=0)

    class FakeObjective:
        max_recuperator_lift = 123.333333333
        initial_mass_fractions = (0.3, 0.7)

        def __init__(self, max_evals, exact_evaluation_worker_max_tasks, **kwargs):
            self.max_evals = max_evals
            self.exact_evaluations = 0
            self.replay = exact_evaluation_worker_max_tasks == 1

        def evaluate(self, x):
            if state["interrupt"] and state["search_calls"] == 5:
                raise RuntimeError("simulated interruption")
            assert self.exact_evaluations < self.max_evals
            self.exact_evaluations += 1
            state["replay_calls" if self.replay else "search_calls"] += 1
            u = (np.asarray(x[:3]) - [0., 0., 0.1]) / [100., 100., 0.8]
            sec = float(1 + np.sum((u - [0.5, 0.3, 0.8])**2))
            return SimpleNamespace(x=tuple(x), success=True, sec=sec, cop=1/sec, pmax=1000.,
                                   report={}, elapsed_seconds=0., message="success", constraint_errors={})

        def close(self):
            pass

    monkeypatch.setattr(benchmark, "_OuterDesignObjective", FakeObjective)
    monkeypatch.setattr(benchmark, "load_case", lambda *args: ({}, {}))
    task = ("Propane;Pentane|0.3;0.7|mass", method, 7, 12, "unused", str(tmp_path / "resumed"))
    with pytest.raises(RuntimeError, match="simulated interruption"):
        benchmark.run_trial(task)
    state["interrupt"] = False
    resumed = json.loads(Path(benchmark.run_trial(task)).read_text())
    assert state["search_calls"] == 12 and state["replay_calls"] == 1
    assert resumed["complete"] and resumed["replay_matches"]
    baseline = json.loads(Path(benchmark.run_trial((*task[:-1], str(tmp_path / "fresh")))).read_text())
    assert np.array_equal([r["x"] for r in resumed["records"]], [r["x"] for r in baseline["records"]])
    assert [r["SEC"] for r in resumed["records"]] == [r["SEC"] for r in baseline["records"]]
    calls = state.copy()
    benchmark.run_trial(task)
    assert state == calls  # Completed checkpoint returns without oracle calls.
