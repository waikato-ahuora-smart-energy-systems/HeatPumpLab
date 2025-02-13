__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from systems import TRHP, set_trhp_model, TEHP, set_tehp_model
from streams.propertiesPackage import PropertyPackage
from functions.solver import solve_heat_pump_model
from functions.cases import get_analysis_cases_to_solve, get_optimisation_cases_to_solve, get_optimisation_case_comparison
from functions.readWriteOutput import print_result_to_terminal, print_result_to_csv, transfer_temp_files_to_results, read_csv_to_results
from functions.decorator_functions import time_func
from concurrent.futures import ProcessPoolExecutor as pool
from scipy.optimize import dual_annealing as glb_min, LinearConstraint, minimize as loc_min, approx_fprime
import numpy as np, os, copy

def analyse_single_case():
    args: list[float | str | bool]
    hp: TRHP

    refrigerant_ls = [
        # 'propane;butane|.5;.5|mol',
        'Butane', 
        'Carbon dioxide', 
        'Propane', 
        'Isobutane', 
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
        'Methanol'
        'Ammonia',
        ]
    
    args = {'refrigerant': 'PROPYLEN',
            'dT_GC_min_tar': 20.0, 
            'T_air_tar': 200.0, 
            'T_air_inlet': 15.0,
            'T_source_in': 40.0,
            'dT_gc_out_min_offset': 0.0,
            'n_comp': 2,
            'ihx_on': 1,
            'dTmin_phx': 5.0,
            'dT_phx1_cold': 10.0,
            'dT_phx2_cold': 52.0,
            'dT_phx3_cold': 0.0,
            'comp2_power_frac': 0.9,
            'comp3_power_frac': 0.0,
            'comp_PR_uniform': False,
            # 'motor_eff': 1.00,
            'comp_model_index': 3,
            'debug': False,
            'sink': 'Air',
            'source': 'water',
            'Tref_evap_inlet' : 10,
            }
    
    hp, dt = solve_single_case(args)
    print('Time to solve case:', dt)
    if hp.converged:
        hp.get_summary()

        if type(hp.output) is dict:
            print_result_to_terminal(output=hp.output)
            print_result_to_csv(output=hp.output)
        
        hp.get_driving_force_plot(hp.gas_coolers, hp.dT_GC_min_tar)
        hp.get_cycle_plot(x_axis='h', y_axis='T')
        hp.get_cycle_plot(x_axis='s', y_axis='T')
        hp.get_cycle_plot(x_axis='h', y_axis='P')
        hp.get_cycle_plot(x_axis='v', y_axis='P')
    else:
        print('Heat model failed to solve.')

    return True

def analyse_single_TEHP_case():
    args: list[float | str | bool]
    hp: TEHP

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
            'dT_GC_min_tar': 5.0, 
            'T_air_tar': 200.0, 
            'T_air_inlet': 30.0, 
            'T_source_in': 30.0,
            'dT_gc_out_min_offset': 0.0,
            'n_comp': 2,
            'ihx_on': 1,
            'dTmin_ev_sm': 10.0,
            'dTmin_econ': 5.0,
            'y_econ': 0.5,
            'comp2_power_frac': 0.5,
            'comp_PR_uniform': False,
            'motor_eff': 1.0,
            'debug': False,
            'sink': 'Air',
            'source': 'water',
            'Tref_evap_inlet' : 20.0,
            }
    
    hp = TEHP(**args)
    hp = set_tehp_model(hp)
    if hp.feasible:
        hp = solve_heat_pump_model(hp)

    if hp.converged:
        hp.get_summary()

        if type(hp.output) is dict:
            print_result_to_terminal(output=hp.output)
            # print_result_to_csv(output=hp.output)
        
        hp.get_driving_force_plot(hp.gas_coolers, hp.dT_GC_min_tar)
        hp.get_cycle_plot(x_axis='h', y_axis='T')
        hp.get_cycle_plot(x_axis='s', y_axis='T')
        hp.get_cycle_plot(x_axis='h', y_axis='P')
        hp.get_cycle_plot(x_axis='v', y_axis='P')
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
                if get_optimisation_case_comparison(sol,unsolved_cases[i]):
                    unsolved_cases.pop(i)
                    break
                    
    if unsolved_cases != []:
        for case in unsolved_cases:
            args = {'refrigerant': str(case['refrigerant']),
                    'dT_GC_min_tar': float(case['dT_GC_min_tar']), 
                    'T_air_tar': float(case['T_air_tar']),
                    'T_air_inlet': float(case['T_air_inlet']), 
                    'T_source_in': float(case['T_source_in']), 
                    'n_comp': int(case['n_comp']), 
                    'dT_phx1_cold': float(case['dT_phx1_cold']),
                    'dT_phx2_cold': float(case['dT_phx2_cold']),
                    'dT_phx3_cold': float(case['dT_phx3_cold']), 
                    'comp2_power_frac': float(case['comp2_power_frac']), 
                    'comp3_power_frac': float(case['comp3_power_frac']),
                    'comp_PR_uniform': bool(case['comp_PR_uniform']),
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

def analyse_cases():
    _solve_in_parallel(
        f_cases=get_analysis_cases_to_solve,
        f_solve=_analyse_case,
        file_name='results/TRHP_results',
        # n_cpu=1,
    )

def optimise_cases():
    _solve_in_parallel(
        f_cases_to_solve=get_optimisation_cases_to_solve,
        f_solve=_optimise_case,
        # n_cpu=1,
    )

def _solve_in_parallel(f_cases, f_solve, file_name: str = 'results/TRHP_results', n_cpu: int = os.cpu_count() - 1) -> None:
    # Use parallel processing distributed over all cores to speed up the analysis
    cases: list        = f_cases(file_name)
    print(len(cases))
    executor: pool     = pool(max_workers=n_cpu)
    [executor.submit(f_solve, case) for case in cases]
    executor.shutdown()
    print('All cases completed.')
    transfer_temp_files_to_results(file_name)

def _analyse_case(case: dict) -> None:
    vars_input = { 'dT_gc_out_min_offset' : { 'xi' : 0.0, 'bnds' : (0.0, 40.0) } }
    case['additional_data'] = {'vars_input' : vars_input}

    result, dt = _get_local_optimum(
        param=case,
        vars_input=vars_input,
        tol=1e-4,
        f_obj=_compute_SEC_obj,
        )
    
    result = _get_results_report(
        param=case, 
        result=result, 
        f_obj=_compute_SEC_obj, 
        dt=dt, 
        print_to_terminal=False,
        )
    
def _optimise_case(case: dict) -> None:
    tol = 1e-5
    dt = 0.0
    # print_result_to_terminal(case)
    vars_input = _get_variable_input_data(case)
    case['additional_data']['vars_input'] = vars_input
    hp, _ = solve_single_case( case | _get_var_kwargs(_get_initial_val(vars_input), case) )
    case['additional_data']['SEC_0'] = hp.SEC
    case['additional_data']['P_0'] = hp.comp_discharge.pressure
    result, dt = _get_global_optimum(case, vars_input, tol)
    result = _get_results_report(case, result, dt)

def _get_variable_input_data(param: dict, x: list =None) -> dict:
    #Define the bounds of the optimisation variables 
    vars_input = dict()
    dT_upper: float = 100.0
    i: int = 0

    if param['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    # vars_input['dT_phx1_cold'] = {
    #     'xi' : 0.0 if x is None else x[i],
    #     'bnds' : (0.0, dT_upper),
    #     }
    # i += 1
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
    return glb_min(func=lambda x: _get_obj(x, param, _compute_multi_obj), 
                   x0=_get_initial_val(vars_input),
                   bounds=_get_bounds(vars_input),
                #    constraints=_get_constraints(param, vars_input),
                #    polish=True,
                   minimizer_kwargs={
                        'method' : 'L-BFGS-B',
                        'bounds' : _get_bounds(vars_input),
                        'jac'    : lambda x: _get_jac(x, param),
                   },
                   no_local_search=False,
                #    maxiter=30,
                   )

@time_func
def _get_local_optimum(param: dict, vars_input: dict, tol: float, f_obj) -> dict:
    return loc_min(fun=lambda x: _get_obj(x, param, f_obj), 
                    x0=_get_initial_val(vars_input),
                    bounds=_get_bounds(vars_input),
                    # constraints=_get_constraints(param, vars_input),
                    method='L-BFGS-B', #'SLSQP'
                    # jac=lambda y: _get_jac(y, param, f_obj),
                    tol=tol,
                    )

def _get_obj(vars: list[float], param: dict, f_obj) -> float:
    hp: TRHP
    hp, _ = solve_single_case(param | _get_var_kwargs(vars, param))
    obj = f_obj(hp, param)
    return obj * 1000

def _compute_SEC_obj(hp: TRHP, param: dict) -> float:
    return hp.SEC

def _compute_multi_obj(hp: TRHP, param: dict) -> float:
    SEC_0 = param['additional_data']['SEC_0']
    P_0   = param['additional_data']['P_0']
    k     = param['additional_data']['weight']
    SEC, P = [hp.SEC, hp.comp_discharge.pressure] if hp.feasible else [np.inf, np.inf]
    return (SEC / SEC_0) * k + (P / P_0) * (1 - k)

def _get_jac(vars: list[float], param: dict, f_obj):
    hp: TRHP
    hp, _ = solve_single_case(param | _get_var_kwargs(vars, param))
    _ = f_obj(hp, param)
    if hp.converged:
        hp.compute_jac = True
        eps = np.sqrt(np.finfo(float).eps) * vars * 10
        return approx_fprime(vars, _compute_jac, eps, hp, param) * 1000
    else:
        return np.inf(len(vars))

def _compute_jac(vars, hp, param, f_obj):
    kwargs = _get_var_kwargs(vars, param)
    for key in kwargs:
        setattr(hp, key, kwargs[key])
    return f_obj(solve_heat_pump_model(hp), param)

def _get_var_kwargs(var: list[float], param: dict) -> dict[str, float]:
    i: int = 0
    d_var = dict()
    for key in param['additional_data']['vars_input']:
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

def _get_results_report(param: dict, result: list[dict], f_obj: float, dt: float, print_to_terminal: bool=False) -> list:
    if result.success:
        res = _get_results_summary(result, param, f_obj, dt)
        fsuffix = _get_unique_f_suffix(param)
        _print_results(res, fsuffix, print_to_terminal)

def _print_results(results: dict, fsuffix: str, print_to_terminal: bool) -> None:
    if print_to_terminal:
        print_result_to_terminal(results)
    print_result_to_csv(output=results, suffix=fsuffix)    

def _get_results_summary(result, param, f_obj, dt: float =None) -> dict:
    hp: TRHP
    hp, _ = solve_single_case(param | _get_var_kwargs(result.x, param))
    obj = f_obj(hp, param)

    def _get_dQdW(hp: TRHP) -> float:
        hp1: TRHP
        dtemp = 0.001
        args = copy.deepcopy(hp.args)
        args['T_air_tar'] = hp.T_air_tar + dtemp
        args['P_H_est'] = hp.comp_discharge.pressure
        args['T_gc_out'] = hp.gas_cooler1.source_out.temperature
        hp1, _ = solve_single_case(args)
        if hp1.feasible:
            tcop = [(p.T_air_tar - p.T_air_inlet) / p.COP if p.COP > 0 else np.inf for p in [hp, hp1]]
            dqdw = dtemp / (tcop[1] - tcop[0]) if tcop[1] - tcop[0] != 0 and tcop[0] != np.inf and tcop[1] != np.inf else 0.0  
        else:
            dqdw = 0.0
        return dqdw   

    hp.dQdW = _get_dQdW(hp)

    hp.get_summary()
    if 'weight' in param['additional_data']:
        hp.output.update({
            'Weight' : param['additional_data']['weight'], 
            'SEC_0' : param['additional_data']['SEC_0'],
            'P_0' : param['additional_data']['P_0'],
        })
    hp.output.update({
        'Obj_function' : obj,
        'Function evaluations' : result.nfev,
        'Message' : result.message,
        'Success' : result.success,
        'Time_to_optimise' : dt,
    })
    return hp.output

def _get_unique_f_suffix(param, fld_only=False) -> str:
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
    if fld_only:
        for ele in param:    
            if ele != 'additional_data' and ele != 'refrigerant':
                suffix = suffix + '-' + "{:.0f}".format(param[ele])
    return suffix
