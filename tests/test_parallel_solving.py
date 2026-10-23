from concurrent.futures import Future

from HeatPumpLab.utils import parallel_solving
from HeatPumpLab.utils.parallel_solving import solve_in_parallel


def test_solve_in_parallel_runs_serially_when_one_cpu_is_requested():
    calls = []

    def solve_fn(hp_cls, case, **kwargs):
        calls.append((hp_cls, case, kwargs))

    assert solve_in_parallel(
        "HP",
        [{"case": 1}, {"case": 2}],
        solve_fn,
        n_cpu=1,
        solve_kwargs={"results_dir": "custom"},
    )
    assert calls == [
        ("HP", {"case": 1}, {"results_dir": "custom"}),
        ("HP", {"case": 2}, {"results_dir": "custom"}),
    ]


def test_solve_in_parallel_submits_each_case_to_executor_factory():
    executors = []

    class FakeExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers
            self.submitted = []
            self.shutdown_called = False
            executors.append(self)

        def submit(self, fn, *args, **kwargs):
            self.submitted.append((fn, args, kwargs))
            future = Future()
            future.set_result(None)
            return future

        def shutdown(self):
            self.shutdown_called = True

    def solve_fn(hp_cls, case, **kwargs):
        return hp_cls, case

    cases = [{"case": 1}, {"case": 2}]

    assert solve_in_parallel(
        "HP",
        cases,
        solve_fn,
        n_cpu=2,
        executor_factory=FakeExecutor,
        solve_kwargs={"results_dir": "custom"},
    )

    assert len(executors) == 1
    assert executors[0].max_workers == 2
    assert executors[0].submitted == [
        (solve_fn, ("HP", {"case": 1}), {"results_dir": "custom"}),
        (solve_fn, ("HP", {"case": 2}), {"results_dir": "custom"}),
    ]
    assert executors[0].shutdown_called is True


def test_solve_in_parallel_uses_process_pool_executor_by_default(monkeypatch):
    executors = []

    class FakeProcessPoolExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers
            self.submitted = []
            self.shutdown_called = False
            executors.append(self)

        def submit(self, fn, *args, **kwargs):
            self.submitted.append((fn, args, kwargs))
            future = Future()
            future.set_result(None)
            return future

        def shutdown(self):
            self.shutdown_called = True

    def solve_fn(hp_cls, case, **kwargs):
        return hp_cls, case

    monkeypatch.setattr(
        parallel_solving,
        "ProcessPoolExecutor",
        FakeProcessPoolExecutor,
    )

    assert solve_in_parallel("HP", [{"case": 1}], solve_fn, n_cpu=2)
    assert len(executors) == 1
    assert executors[0].max_workers == 2
    assert executors[0].submitted == [(solve_fn, ("HP", {"case": 1}), {})]
    assert executors[0].shutdown_called is True


def test_solve_in_parallel_reports_worker_failures():
    class FailingExecutor:
        def __init__(self, max_workers):
            pass

        def submit(self, fn, *args, **kwargs):
            future = Future()
            try:
                future.set_result(fn(*args, **kwargs))
            except Exception as exc:
                future.set_exception(exc)
            return future

        def shutdown(self):
            pass

    def solve_fn(hp_cls, case):
        if case["case"] == 2:
            raise RuntimeError("REFPROP missing")

    assert not solve_in_parallel("HP", [{"case": 1}, {"case": 2}], solve_fn, n_cpu=2, executor_factory=FailingExecutor)


def test_default_n_cpu_is_at_least_one(monkeypatch):
    monkeypatch.setattr(parallel_solving.os, "cpu_count", lambda: 1)
    assert parallel_solving.default_n_cpu() == 1
    monkeypatch.setattr(parallel_solving.os, "cpu_count", lambda: None)
    assert parallel_solving.default_n_cpu() == 1
