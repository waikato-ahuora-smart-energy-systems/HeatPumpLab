__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from systems import TRHP, set_trhp_model
from functions.readWriteOutput import transfer_temp_files_to_results, print_result_to_terminal, read_csv_to_results
from functions.solver import solve_heat_pump_model
from functions.cases import get_cases_to_solve
from functions.optimiseTRHP import get_parallel_solve

def analyse_cases_in_parallel():
    file_name: str
    cases: list[list[float | str | bool]]

    file_name = 'results/TRHP_results'
    cases = get_cases_to_solve(file_name)
    get_parallel_solve(cases)
    transfer_temp_files_to_results(file_name=file_name)

def analyse_single_case():
    args: list[float | str | bool]
    hp: TRHP
    tic: float
    toc: float

    refrigerant_ls: list[str] = [
        'propane;butane|.5;.5|mol',
        'Butane', 
        'Carbon dioxide', 
        'Propane', 
        'Isobutane', 
        'R1233zd(E)', 
        'R1336mzz(Z)', 
        'R1234ze(E)', 
        'R1234ze(Z)', 
        'R1234yf', 
        'R1243zf', 
        'R1224yd(Z)', 
        'PROPYLEN',
        'DME',
        'Methanol'
        'Ammonia',
        ]

    args = {'refrigerant': 'propane;butane|.5;.5|mol',
            'dT_GC_min_tar': 20.0, 
            'T_air_tar': 250.0, 
            'T_air_inlet': 15.0, 
            'T_source_in': 40.0, 
            'T_gc_out_min': 35.0,
            'n_comp': 3, 
            'dTmin_phx': 5.0,
            'dT_phx1_cold': 20.0,
            'dT_phx2_cold': 30.0,
            'dT_phx3_cold': 20.0, 
            'comp2_power_frac': 0.5, 
            'comp3_power_frac': 0.0,
            'comp_PR_uniform': True,
            'eff_motor': 0.97,
            'debug': False,
            'sink': 'Air',
            'source': 'water',
            'Tref_evap_inlet' : 10,
            }
    
    hp = TRHP(**args)
    hp = set_trhp_model(hp)
    if hp.feasible:
        hp = solve_heat_pump_model(hp)
    
    if hp.converged:
        hp.get_summary()

        # if type(hp.output) is dict:
        #     print_result_to_terminal(output=hp.output)
        
        hp.get_driving_force_plot(hp.gas_coolers, hp.dT_GC_min_tar)
        hp.get_cycle_plot(x_axis='h', y_axis='T')
        # hp.get_cycle_plot(x_axis='s', y_axis='T')
        # hp.get_cycle_plot(x_axis='h', y_axis='P')
        # hp.get_cycle_plot(x_axis='v', y_axis='P')
    else:
        print('Heat model failed to solve.')

def analyse_unsolved_cases():
    # unsolved = read_csv_to_results(file_name='results/Cases-Unsolved')
    
    # for case in unsolved:
    #     args = {'refrigerant': case['Refrigerant:'],
    #             'dT_GC_min_tar': case['dT_min Gas Cooler:'], 
    #             'T_air_tar': case['T_air_tar:'],
    #             'T_air_inlet': case['T_air_inlet:'], 
    #             'T_source_in': case['T_source_in:'], 
    #             'n_comp': case['n_comp:'], 
    #             'dT_phx1_cold': case['dT phx max 1:'],
    #             'dT_phx2_cold': case['dT phx max 2:'],
    #             'dT_phx3_cold': case['dT phx max 3:'], 
    #             'comp2_power_frac': case['Comp 2 power frac:'], 
    #             'comp3_power_frac': case['Comp 3 power frac:'],
    #             'comp_PR_uniform': case['Comp_PR_uniform:'],
    #             'debug': False,
    #             }
    #     hp = TRHP(**args)
    #     hp.solve()
        
    #     # hp.get_cycle_plot(x_axis='s', y_axis='T')
    #     # hp.get_cycle_plot(x_axis='h', y_axis='P')
    #     # hp.get_cycle_plot(x_axis='v', y_axis='P')
        None
        