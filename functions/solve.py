__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from systems import HeatPump, set_trhp_model
from functions.readWriteOutput import transfer_temp_files_to_results
from functions.solver import solve_heat_pump_model
from functions.cases import get_cases_to_solve
from functions.optimiseTRHP import get_parallel_solve
from time import perf_counter as timer

def solve_cases_in_parallel():
    file_name: str
    cases: list[list[float | str | bool]]

    file_name = 'results/TRHP_results'
    cases = get_cases_to_solve(file_name)
    get_parallel_solve(cases)
    transfer_temp_files_to_results(file_name=file_name)

def solve_single_case():
    args: list[float | str | bool]
    hp: HeatPump
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

    args = {'refrigerant': 'R1336mzz(Z)',
            'dT_GC_min_tar': 20.0, 
            'T_air_tar': 150.0, 
            'T_air_inlet': 15.0, 
            'T_source_in': 40.0, 
            'n_comp': 2, 
            'dT_phx1_cold': 30.0,
            'dT_phx2_cold': 30.0,
            'dT_phx3_cold': 0.0, 
            'comp2_power_frac': 0.5, 
            'comp3_power_frac': 0.0,
            'comp_PR_uniform': True,
            'display': False,
            }

    tic = timer()
    hp = HeatPump(**args)
    hp = set_trhp_model(hp, 2)
    hp = solve_heat_pump_model(hp)
    toc = timer()
    hp.get_summary()
    # print_result_to_csv
    # print_result_to_terminal(hp.output)
    print('Solve time:', toc - tic)
    print('Discharge pressure:', hp.compressors[-1].flow_out.pressure)
    print('')
    # hp.get_driving_force_plot()

    tic = timer()
    hp = HeatPump(**args)
    hp = set_trhp_model(hp)
    hp.solve()
    toc = timer()
    hp.get_summary()
    # print_result_to_terminal(hp.output)
    print('Solve time:', toc - tic)
    print('Discharge pressure:', hp.compressors[-1].flow_out.pressure)
    # hp.get_driving_force_plot()
    None
    # if hp.converged:
    #     hp.get_summary()
    #     hp.output.update({
    #         'Obj function:' : '',
    #         'Iterations:' : 0,
    #         'Message:' : '',
    #         'Success:' : hp.converged
    #     })

    #     if type(hp.output) is dict:
    #         print_result_to_terminal(output=hp.output)
    #         print('Solve time:', toc - tic)
        
    #     hp.get_cycle_plot(x_axis='h', y_axis='T')
    #     hp.get_cycle_plot(x_axis='h', y_axis='P')
    #     hp.get_cycle_plot(x_axis='v', y_axis='P')

# def unsolved_case_analysis():
#     unsolved = read_csv_to_results(file_name='results/Cases-Unsolved')
    
#     for case in unsolved:
#         args = {'refrigerant': case['Refrigerant:'],
#                 'dT_GC_min_tar': case['dT_min Gas Cooler:'], 
#                 'T_air_tar': case['T_air_tar:'],
#                 'T_air_inlet': case['T_air_inlet:'], 
#                 'T_source_in': case['T_source_in:'], 
#                 'n_comp': case['n_comp:'], 
#                 'dT_phx1_cold': case['dT phx max 1:'],
#                 'dT_phx2_cold': case['dT phx max 2:'],
#                 'dT_phx3_cold': case['dT phx max 3:'], 
#                 'comp2_power_frac': case['Comp 2 power frac:'], 
#                 'comp3_power_frac': case['Comp 3 power frac:'],
#                 'comp_PR_uniform': case['Comp_PR_uniform:'],
#                 'display': False,
#                 }
#         hp = TRHP(**args)
#         hp.solve()
        
#         if hp.converged is False:
#             hp.print_failed_to_solve_report('results/Cases-Unsolved-rechecked')
        
#         # hp.get_cycle_plot(x_axis='s', y_axis='T')
#         # hp.get_cycle_plot(x_axis='h', y_axis='P')
#         # hp.get_cycle_plot(x_axis='v', y_axis='P')
#         None
        