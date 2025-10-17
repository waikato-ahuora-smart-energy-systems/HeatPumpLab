__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2025'

"""
High temperature recuperative heat pump model
By Tim Walmsley and Keri Adamson, 
University of Waikato, 
2021-2025 
tim.walmsley@waikato.ac.nz
"""

from HeatPumpLab import *
from aux import *
 
if __name__=='__main__':
    analyse_single_case()
    analyse_cases_in_parallel()
    # analyse_unsolved_cases(debug=False)
    # get_fluid_prop_overview_table()
    pass