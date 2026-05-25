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

import numpy as np
from HeatPumpLab import (
    RecuperativeHP, 
    analyse_multiple_cases_in_parallel,
)
from aux import *

_REFRIGERANT_LIST = [
        'Ammonia',
        'Propane', 
        'Butane', 
        'Isobutane',
        'Pentane',
        'Isopentane',
        'Neopentane',
        'Cyclopentane',        
        'Carbon dioxide', 
        'R1233zd(E)', 
        'R1336mzz(E)', 
        'R1336mzz(Z)', 
        'R1234ze(E)', 
        'R1234ze(Z)', 
        'R1234yf', 
        'R1243zf', 
        'R1224yd(Z)', 
        'PROPYLEN',
        'DME',
        'Methanol',
        # 'propane;butane|.5;.5|mol',
        # 'propane;Isobutane|.5;.5|mol',
        # 'propane;Pentane|.5;.5|mol',
        # 'propane;Isopentane|.5;.5|mol',
        # 'propane;Neopentane|.5;.5|mol',
        # 'Isobutane;butane|.5;.5|mol',
        # 'Isobutane;Pentane|.5;.5|mol',
        # 'Isobutane;Isopentane|.5;.5|mol',
        # 'Isobutane;Neopentane|.5;.5|mol',
        # 'butane;Pentane|.5;.5|mol',
        # 'butane;Isopentane|.5;.5|mol',
        # 'butane;Neopentane|.5;.5|mol',
    ]

def run_all_cases():  
    analyse_multiple_cases_in_parallel(
        HP_fn=RecuperativeHP, 
        filename='RecuperatureHP_results',
        refrigerant_ls=_REFRIGERANT_LIST,
        dT_rhx1_cold_ls=np.linspace(0, 100, num=41),
        dT_rhx2_cold_ls=np.linspace(0, 100, num=41),
        comp2_power_frac_ls=np.linspace(0.3, 0.7, num=9),
        # n_cpu=-1,
    )

if __name__=='__main__':
    run_all_cases()