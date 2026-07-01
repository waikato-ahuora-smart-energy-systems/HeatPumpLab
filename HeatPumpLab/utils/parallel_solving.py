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
    n_cpu: int = os.cpu_count() - 1,
    executor_factory=None,
    solve_kwargs: dict | None = None,
) -> bool:
    """Execute the provided solve function across cases using a process pool."""
    print(len(cases))
    if solve_kwargs is None:
        solve_kwargs = {}

    if n_cpu == 1:
        for case in cases:
            f_solve(HP_fn, case, **solve_kwargs)
        print("All cases completed.")
        return True

    if executor_factory is None:
        executor_factory = ProcessPoolExecutor

    executor = executor_factory(max_workers=n_cpu)
    for case in cases:
        executor.submit(f_solve, HP_fn, case, **solve_kwargs)
    executor.shutdown()
    print("All cases completed.")
    return True
