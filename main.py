__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from functions import analyse_single_case#, analyse_cases_in_parallel#, analyse_unsolved_cases
# import type_enforced 

# @type_enforced.Enforcer
def main() -> None:
    analyse_single_case()
    # analyse_unsolved_cases()
    # analyse_cases_in_parallel()

if __name__=='__main__':
    main()