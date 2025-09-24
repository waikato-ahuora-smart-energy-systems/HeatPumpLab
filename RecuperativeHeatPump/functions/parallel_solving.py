from concurrent.futures import ProcessPoolExecutor as pool
import os
from .readWriteOutput import *

__all__ = ['solve_in_parallel']

#######################################################################################################
# Public API
#######################################################################################################

def solve_in_parallel(f_cases, f_solve, file_name: str = 'RHP_results', n_cpu: int = os.cpu_count() - 1) -> bool:
    # Use parallel processing distributed over all cores to speed up the analysis
    cases: list = f_cases('results/' + file_name)
    print(len(cases))
    executor: pool = pool(max_workers=n_cpu)
    [executor.submit(f_solve, case) for case in cases]
    executor.shutdown()
    print('All cases completed.')
    transfer_temp_files_to_results('results/' + file_name)
    return True
 