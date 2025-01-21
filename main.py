__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2025'

from functions import analyse_cases, transfer_temp_files_to_results, optimise_cases, analyse_unsolved_cases, analyse_single_TEHP_case, analyse_single_case

def main():
    # transfer_temp_files_to_results()
    # analyse_single_TEHP_case()
    # analyse_single_case()
    # analyse_unsolved_cases(debug=False)
    analyse_cases()
    # optimise_cases()
    
 
if __name__=='__main__':
    main()