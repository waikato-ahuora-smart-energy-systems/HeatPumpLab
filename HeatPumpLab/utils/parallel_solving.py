"""Helper to submit case evaluations in parallel."""

from concurrent.futures import ProcessPoolExecutor
import os

__all__ = ["solve_in_parallel"]

#######################################################################################################
# Public API
#######################################################################################################


def solve_in_parallel(
    HP_fn,
    cases: list,
    f_solve,
    n_cpu: int | None = None,
    executor_factory=None,
    solve_kwargs: dict | None = None,
) -> bool:
    """Execute the provided solve function across cases using a process pool.

    Returns False (after reporting each failure) when any parallel case raised.
    """
    print(len(cases))
    if solve_kwargs is None:
        solve_kwargs = {}
    if n_cpu is None:
        n_cpu = default_n_cpu()

    if n_cpu == 1:
        for case in cases:
            f_solve(HP_fn, case, **solve_kwargs)
        print("All cases completed.")
        return True

    if executor_factory is None:
        executor_factory = ProcessPoolExecutor

    executor = executor_factory(max_workers=n_cpu)
    futures = [(case, executor.submit(f_solve, HP_fn, case, **solve_kwargs)) for case in cases]
    executor.shutdown()
    failures = [(case, exc) for case, future in futures if (exc := future.exception()) is not None]
    for case, exc in failures:
        print(f"Case failed: {case}: {type(exc).__name__}: {exc}")
    if failures:
        print(f"{len(failures)} of {len(cases)} cases failed.")
        return False
    print("All cases completed.")
    return True


def default_n_cpu() -> int:
    """One worker per CPU, leaving one free; at least one."""
    return max(1, (os.cpu_count() or 2) - 1)
