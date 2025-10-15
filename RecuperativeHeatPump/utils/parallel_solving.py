"""Helper to submit case evaluations in parallel."""

from concurrent.futures import ProcessPoolExecutor as pool
import os
from .readWriteOutput import *

__all__ = ['solve_in_parallel']

#######################################################################################################
# Public API
#######################################################################################################

def solve_in_parallel(cases: list, f_solve, n_cpu: int = os.cpu_count() - 1) -> bool:
    """Execute the provided solve function across cases using a process pool."""
    # Use parallel processing distributed over all cores to speed up the analysis
    print(len(cases))
    executor: pool = pool(max_workers=n_cpu)
    [executor.submit(f_solve, case) for case in cases]
    executor.shutdown()
    print('All cases completed.')
    return True
 
