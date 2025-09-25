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

from RecuperativeHeatPump import *
from aux import *

def analyse_single_case():
    args: list[float | str | bool]
    hp: RHP
    
    args = {'refrigerant': 'Propane',
            'dT_GC_min_tar': 20.0, 
            'T_air_tar': 200.0, 
            'T_air_inlet': 15.0,
            'T_source_in': 40.0,
            'dT_gc_out_min_offset': 0.0,
            'n_comp': 2,
            'ihx_on': 1,
            'dTmin_rhx': 5.0,
            'dT_rhx1_cold': 0.0,
            'dT_rhx2_cold': 70.0,
            'dT_rhx3_cold': 0.0,
            'comp2_power_frac': 0.6,
            'comp3_power_frac': 0.0,
            'comp_PR_uniform': False,
            'motor_eff': 1.00,
            'comp_model_index': 1,
            'debug': False,
            'sink': 'Air',
            'source': 'water',
            'Tref_evap_inlet' : 10,
            }

    hp = RHP(**args)
    hp.solve()
    
    if hp.converged:
        hp.get_report()

        if type(hp.output) is dict:
            print_result_to_terminal(output=hp.output)
            print_result_to_csv(output=hp.output, file_name='single-result')
        
        hp.get_driving_force_plot(hp.gas_coolers, hp.dT_GC_min_tar)
        hp.get_cycle_plot(x_axis='h', y_axis='T')
        hp.get_cycle_plot(x_axis='s', y_axis='T')
        hp.get_cycle_plot(x_axis='h', y_axis='P')
        hp.get_cycle_plot(x_axis='v', y_axis='P')
    elif not hp.feasible:
        print('Heat pump case is infeasible.')
    else:
        print('Heat pump model failed to solve.')

 
if __name__=='__main__':
    # analyse_single_case()
    analyse_cases_in_parallel()
    # analyse_unsolved_cases()
    pass