
from RecuperativeHeatPump import *

DEFAULT_ARGS: list[float | str | bool] = {
    'refrigerant': 'CO2',
    'dT_GC_min_tar': 20.0, 
    'T_air_tar': 200.0, 
    'T_air_inlet': 15.0,
    'T_source_in': 40.0,
    'dT_gc_out_min_offset': 8.43169831,
    'n_comp': 2,
    'ihx_on': 1,
    'dTmin_rhx': 5.0,
    'dT_rhx1_cold': 35.0,
    'dT_rhx2_cold': 0.0,
    'dT_rhx3_cold': 0.0,
    'comp2_power_frac': 0.5,
    'comp3_power_frac': 0.0,
    'comp_PR_uniform': False,
    'motor_eff': 0.96,
    'comp_model_index': 1,
    'debug': False,
    'sink': 'Air',
    'source': 'water',
    'Tref_evap_inlet' : 10,
}

def analyse_single_case(args: list[float | str | bool] = None):
    
    hp: RHP
    
    if args == None:
        args = DEFAULT_ARGS

    hp = RHP(**args)
    hp.solve()
    
    if hp.converged:
        hp.get_report()

        if type(hp.output) is dict:
            print_result_to_terminal(output=hp.output)
            print_result_to_csv(output=hp.output, file_name='single-result')
        
        refrigerant = args['refrigerant']

        hp.get_driving_force_plot(hp.gas_coolers, hp.dT_GC_min_tar, refrigerant)
        hp.get_cycle_plot(x_axis='h', y_axis='T', refrigerant=refrigerant)
        hp.get_cycle_plot(x_axis='s', y_axis='T', refrigerant=refrigerant)
        hp.get_cycle_plot(x_axis='h', y_axis='P', refrigerant=refrigerant, include_sink=False)
        hp.get_cycle_plot(x_axis='v', y_axis='P', refrigerant=refrigerant)
    elif not hp.feasible:
        print('Heat pump case is infeasible.')
    else:
        print('Heat pump model failed to solve.')

