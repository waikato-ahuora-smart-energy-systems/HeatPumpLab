__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from systems import TRHP, set_trhp_model
from streams.propertiesPackage import PropertyPackage
from functions.solver import solve_heat_pump_model
from functions.cases import get_cases_to_solve
from functions.readWriteOutput import print_result_to_terminal, print_result_to_csv, transfer_temp_files_to_results, read_csv_to_results
from functions.decorator_functions import time_func
from concurrent.futures import ProcessPoolExecutor as pool
from scipy.optimize import dual_annealing as glb_min, LinearConstraint, minimize as loc_min, approx_fprime
# from openbox import Optimizer, space as sp
# from bayes_opt import BayesianOptimization
import numpy as np, os

def analyse_single_case():
    args: list[float | str | bool]
    hp: TRHP

    refrigerant_ls = [
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
            'dT_gc_out_min_offset': 0.0,
            'n_comp': 3, 
            'dTmin_phx': 5.0,
            'dT_phx1_cold': 20.0,
            'dT_phx2_cold': 20.0,
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
    
    hp, dt = solve_single_case(args)
    print('Time to solve case:', dt)
    if hp.converged:
        hp.get_summary()

        # if type(hp.output) is dict:
        #     print_result_to_terminal(output=hp.output)
            # print_result_to_csv(output=hp.output)
        
        # hp.get_driving_force_plot(hp.gas_coolers, hp.dT_GC_min_tar)
        # hp.get_cycle_plot(x_axis='h', y_axis='T')
        # hp.get_cycle_plot(x_axis='s', y_axis='T')
        # hp.get_cycle_plot(x_axis='h', y_axis='P')
        # hp.get_cycle_plot(x_axis='v', y_axis='P')
    else:
        print('Heat model failed to solve.')

    return True

def analyse_unsolved_cases(debug: bool =True):
    transfer_temp_files_to_results('results/Cases-Resolved')
    unsolved_cases = read_csv_to_results(file_name='results/Cases-Unsolved')
    resolved_cases = read_csv_to_results(file_name='results/Cases-Resolved')

    if resolved_cases != None:
        for i in reversed(range(len(unsolved_cases))):
            for sol in resolved_cases:
                if (
                    sol['Refrigerant:'] == unsolved_cases[i]['Refrigerant:'] and 
                    int(sol['n comp:']) == int(unsolved_cases[i]['n comp:']) and
                    bool(sol['Comp PR uniform:']) == bool(unsolved_cases[i]['Comp PR uniform:']) and
                    abs(float(sol['dT GC min tar:']) - float(unsolved_cases[i]['dT GC min tar:'])) < 0.01 and 
                    abs(float(sol['T air inlet:']) - float(unsolved_cases[i]['T air inlet:'])) < 0.01 and 
                    abs(float(sol['T source in:']) - float(unsolved_cases[i]['T source in:'])) < 0.01 and 
                    abs(float(sol['T evap inlet:']) - float(unsolved_cases[i]['T evap inlet:'])) < 0.01
                    ):
                    unsolved_cases.pop(i)
                    break
                    
    if unsolved_cases != []:
        for case in unsolved_cases:
            args = {'refrigerant': str(case['Refrigerant:']),
                    'dT_GC_min_tar': float(case['dT GC min tar:']), 
                    'T_air_tar': float(case['T air tar:']),
                    'T_air_inlet': float(case['T air inlet:']), 
                    'T_source_in': float(case['T source in:']), 
                    'n_comp': int(case['n comp:']), 
                    'dT_phx1_cold': float(case['dT phx max 1:']),
                    'dT_phx2_cold': float(case['dT phx max 2:']),
                    'dT_phx3_cold': float(case['dT phx max 3:']), 
                    'comp2_power_frac': float(case['Comp 2 power frac:']), 
                    'comp3_power_frac': float(case['Comp 3 power frac:']),
                    'comp_PR_uniform': bool(case['Comp PR uniform:']),
                    'debug': bool(debug),
                    }
            hp, _ = solve_single_case(args)
            if hp.converged:
                hp.get_summary()
                fsuffix = _get_unique_f_suffix(args) + '-RESOLVED'
                _print_results(hp.output, fsuffix)
        
            # hp.get_cycle_plot(x_axis='s', y_axis='T')
            # hp.get_cycle_plot(x_axis='h', y_axis='P')
            # hp.get_cycle_plot(x_axis='v', y_axis='P')
            None
        transfer_temp_files_to_results('results/Cases-Resolved')
        
def analyse_cases_in_parallel(file_name: str = 'results/TRHP_results', n_cpu: int = os.cpu_count() - 1) -> None:
    # Use parallel processing distributed over all cores to speed up the analysis
    cases: list        = get_cases_to_solve(file_name)
    print(len(cases))
    # _optimise_case(cases[0])
    executor: pool     = pool(max_workers=n_cpu)
    results: list[any] = [executor.submit(_optimise_case, case) for case in cases]
    executor.shutdown()
    transfer_temp_files_to_results(file_name)

def _optimise_case(case: dict) -> None:
    tol: float = 1e-5
    dt: float = 0.0
    print_result_to_terminal(case)
    vars_input = _get_variable_input_data(case)
    case['additional_obj']['vars_input'] = vars_input
    hp, _ = solve_single_case( case | _get_var_kwargs(_get_initial_val(vars_input), case) )
    case['additional_obj']['COP_0'] = hp.COP
    case['additional_obj']['P_0'] = hp.comp_discharge.pressure
    result, dt = _get_global_optimum(case, vars_input, tol)
    result = _get_results_report(case, result, dt)

def _get_variable_input_data(param: dict, x: list =None) -> dict:
    #Define the bounds of the optimisation variables 
    vars_input = dict()
    dT_upper: float = 100.0
    i: int = 0

    if param['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    vars_input['dT_phx1_cold'] = {
        'xi' : 0.0 if x is None else x[i],
        'bnds' : (0.0, dT_upper),
        }
    i += 1
    if param['n_comp'] >= 2:
        vars_input['dT_phx2_cold'] = {
            'xi' : 0.0 if x is None else x[i],
            'bnds' : (0.0, dT_upper),
            }
        i += 1
    if param['n_comp'] == 3:
        vars_input['dT_phx3_cold'] = {
            'xi' : 0.0 if x is None else x[i],
            'bnds' : (0.0, dT_upper),
            }
        i += 1
    if param['comp_PR_uniform'] == False and param['n_comp'] > 1:
        if param['n_comp'] == 2:
            vars_input['comp2_power_frac'] = {
                'xi' : 0.5 if x is None else x[i],
                'bnds' : (0.2, 0.8),
                }
            i += 1
        if param['n_comp'] == 3:
            vars_input['comp2_power_frac'] = {
                'xi' : 1/3 if x is None else x[i],
                'bnds' : (0.1, 0.5),
                }
            i += 1
            vars_input['comp3_power_frac'] = {
                'xi' : 1/3 if x is None else x[i],
                'bnds' : (0.1, 0.5),
                }
            i += 1
    vars_input['dT_gc_out_min_offset'] = {
        'xi' : 0.0 if x is None else x[i],
        'bnds' : (0.0, dT_upper / 2),
        }

    return vars_input

def _get_initial_val(vars_input: dict) -> list:
    return [vars_input[key]['xi'] for key in vars_input]

def _update_initial_val(vars_input: dict, prev_res):
    vi = dict()
    for key in vars_input:
        vi[key] = {'xi'   : float(prev_res[key]),
                   'bnds' : vars_input[key]['bnds']}
    return vi

def _get_bounds(vars_input: dict) -> list:
    return [vars_input[key]['bnds'] for key in vars_input]

def _get_constraints(param: dict, vars_input: dict):
    m = np.array([1.0 if ('dT_phx' in key and '_cold' in key) else 0.0 for key in vars_input])
    return LinearConstraint(A=m, lb=0.0, ub=(param['T_air_tar'] - param['T_air_inlet']) * 2/3)   

@time_func
def _get_global_optimum(param: dict, vars_input: dict, tol: float) -> dict:
    return glb_min(func=lambda x: _get_obj(x, param), 
                   x0=_get_initial_val(vars_input),
                   bounds=_get_bounds(vars_input),
                   minimizer_kwargs={
                        'method' : 'L-BFGS-B',
                        'bounds' : _get_bounds(vars_input),
                        'jac'    : lambda x: _get_jac(x, param),
                   },
                   no_local_search=False,
                   maxiter=100,
                   )

@time_func
def _get_local_optimum(param: dict, vars_input: dict, tol: float) -> dict:
    return loc_min(fun=lambda x: _get_obj(x, param), 
                    x0=_get_initial_val(vars_input),
                    bounds=_get_bounds(vars_input),
                    # constraints=_get_constraints(param, vars_input),
                    method='SLSQP', #'SLSQP' 'L-BFGS-B'
                    # jac=lambda x: _get_jac(x, param),
                    tol=tol,
                    )

def _get_obj(vars: list[float], param: dict) -> float:
    hp: TRHP
    dt: float
    hp, dt = solve_single_case(param | _get_var_kwargs(vars, param))
    obj = _compute_obj(hp, param)
    print('Obj:', "{:.3f}".format(obj), '  Solve time:', "{:.3f}".format(dt), '  Vars:', ["{:.3f}".format(v) for v in vars], '  Ref:', param['refrigerant'])
    return -obj * 1000

def _compute_obj(hp: TRHP, param: dict) -> float:
    COP: float
    P: float
    COP_0: float = param['additional_obj']['COP_0']
    P_0: float   = param['additional_obj']['P_0']
    k: float     = param['additional_obj']['weight']
    COP, P = [hp.COP, hp.comp_discharge.pressure] if hp.feasible else [-np.inf, 0.0]
    return (COP / COP_0) * k - (P / P_0 - 2) * (1 - k) 

def _get_jac(vars: list[float], param: dict):
    hp: TRHP
    hp, _ = solve_single_case(param | _get_var_kwargs(vars, param))
    obj = _compute_obj(hp, param)
    if obj < 0:
        hp.compute_jac = True
        eps = np.sqrt(np.finfo(float).eps) * vars * 10
        return approx_fprime(vars, _compute_jac, eps, hp, param)
    else:
        return np.zeros(len(vars))

def _compute_jac(vars, hp, param):
    kwargs = _get_var_kwargs(vars, param)
    for key in kwargs:
        setattr(hp, key, kwargs[key])
    return _compute_obj(solve_heat_pump_model(hp), param)

def _get_var_kwargs(var: list[float], param: dict) -> dict[str, float]:
    i: int = 0
    d_var = dict()
    for key in param['additional_obj']['vars_input']:
        d_var[key] = var[i]
        i += 1
    return d_var

@time_func
def solve_single_case(args) -> TRHP:
    hp = TRHP(**args)
    hp = set_trhp_model(hp)
    if hp.feasible:
        hp = solve_heat_pump_model(hp)
    return hp

def _get_results_report(param: dict, result: dict, dt: float) -> list:
    if result.success:
        res = _get_results_summary(result, param, dt)
        fsuffix = _get_unique_f_suffix(param)
        _print_results(res, fsuffix)
        print('Time to optimise:', dt)

def _print_results(results: dict, fsuffix: str) -> None:
    print_result_to_terminal(results)
    print_result_to_csv(output=results, suffix=fsuffix)    

def _get_results_summary(result, param, dt: float =None) -> dict:
    hp, _ = solve_single_case(param | _get_var_kwargs(result.x, param))
    obj = _compute_obj(hp, param)
    hp.get_summary()
    hp.output.update({
        'Weight' : param['additional_obj']['weight'], 
        'COP_0' : param['additional_obj']['COP_0'],
        'P_0' : param['additional_obj']['P_0'],
        'Obj_function' : obj,
        'Iterations' : result.nfev,
        'Message' : result.message,
        'Success' : result.success,
        'Time_to_optimise' : dt,
    })
    return hp.output

def _get_unique_f_suffix(param) -> str:
    suffix: str = '-'
    
    r = PropertyPackage(param['refrigerant'])
    fld_comp = r.fld_name.replace('.FLD', '').split('*')
    if len(fld_comp) == 1:
        suffix = '-' + param['refrigerant']
    else:
        for i in range(len(fld_comp)):
            suffix += str(fld_comp[i]) + str(r.fld_composition[i]).replace('.', '_')
            if i + 1 < len(fld_comp): 
                suffix += '*'

    for ele in param:    
        if ele != 'additional_obj' and ele != 'refrigerant':
            suffix = suffix + '-' + "{:.0f}".format(param[ele])
    return suffix
