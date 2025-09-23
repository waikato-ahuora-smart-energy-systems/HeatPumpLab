import numpy as np
from scipy.optimize import dual_annealing as glb_min, LinearConstraint, minimize as loc_min, approx_fprime
from ..systems import RHP
from .solver import solve_heat_pump_model
from .decorator_functions import time_func
from .analysis import solve_single_case, get_results_report, solve_in_parallel
from .optimisation_cases import get_optimisation_cases_to_solve

def optimise_cases():
    solve_in_parallel(
        f_cases_to_solve=get_optimisation_cases_to_solve,
        f_solve=_optimise_case,
        # n_cpu=1,
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
    result = get_results_report(case, result, dt)

def _get_variable_input_data(param: dict, x: list =None) -> dict:
    #Define the bounds of the optimisation variables 
    vars_input = dict()
    dT_upper: float = 100.0
    i: int = 0

    if param['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    # vars_input['dT_rhx1_cold'] = {
    #     'xi' : 0.0 if x is None else x[i],
    #     'bnds' : (0.0, dT_upper),
    #     }
    # i += 1
    if param['n_comp'] >= 2:
        vars_input['dT_rhx2_cold'] = {
            'xi' : 0.0 if x is None else x[i],
            'bnds' : (0.0, dT_upper),
            }
        i += 1
    if param['n_comp'] == 3:
        vars_input['dT_rhx3_cold'] = {
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
    m = np.array([1.0 if ('dT_rhx' in key and '_cold' in key) else 0.0 for key in vars_input])
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
    hp: RHP
    hp, _ = solve_single_case(param | _get_var_kwargs(vars, param))
    obj = f_obj(hp, param)
    return obj * 1000

def _compute_SEC_obj(hp: RHP, param: dict) -> float:
    return hp.SEC

def _compute_multi_obj(hp: RHP, param: dict) -> float:
    SEC_0 = param['additional_data']['SEC_0']
    P_0   = param['additional_data']['P_0']
    k     = param['additional_data']['weight']
    SEC, P = [hp.SEC, hp.comp_discharge.pressure] if hp.feasible else [np.inf, np.inf]
    return (SEC / SEC_0) * k + (P / P_0) * (1 - k)

def _get_jac(vars: list[float], param: dict, f_obj):
    hp: RHP
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
