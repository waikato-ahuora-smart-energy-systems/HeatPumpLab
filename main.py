__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from functions import solve_single_case
from time import perf_counter as timer
# from functions import solve_single_case, solve_cases_in_parallel#, unsolved_case_analysis
# import type_enforced 

# @type_enforced.Enforcer
def main() -> None:
    tic: float = timer()
    solve_single_case()
    # unsolved_case_analysis()
    # solve_cases_in_parallel()
    toc: float = timer()
    print('Total time:', toc - tic)

if __name__=='__main__':
    main()