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
