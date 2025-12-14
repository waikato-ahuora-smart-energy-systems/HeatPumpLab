__author__ = 'Tim Walmsley'
__credits__ = ['Keri Adamson', 'Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2025'

"""
Heat pump lab
By Tim Walmsley, 
University of Waikato, 
2021-2025 
tim.walmsley@waikato.ac.nz
"""

from HeatPumpLab import *
from aux import *
 
if __name__=='__main__':
    # analyse_single_case(RecuperativeHP)
    analyse_multiple_cases_in_parallel(RecuperativeHP, filename='RecuperatureHP_results', n_cpu=1)