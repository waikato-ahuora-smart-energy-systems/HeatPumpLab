"""Analysis utilities for running and post-processing heat pump cases."""

from typing import Tuple
from scipy.optimize import minimize as loc_min, approx_fprime
import numpy as np, copy, os
from ..cycles import *
from ..streams import *
from .parallel_solving import *
from .cases import *
from .readWriteOutput import *
from .decorator_functions import *

__all__ = ['analyse_cases_in_parallel', 'analyse_unsolved_cases']

#######################################################################################################
# Public API
#######################################################################################################

def analyse_cases_in_parallel(
        filename: str = 'RecuperativeHP_results',
        refrigerant_ls: list = None, 
        dT_rhx1_cold_ls: list = None,
        dT_rhx2_cold_ls: list = None, 
        comp2_power_frac_ls: list = None,
        n_cpu: int = os.cpu_count() - 1,
        ):
    """Solve a batch of cases concurrently and persist intermediate results."""
    file_name = filename
    transfer_temp_files_to_results(file_name)
    cases = get_case_list_to_solve(
        file_name, 
        refrigerant_ls, 
        dT_rhx1_cold_ls, 
        dT_rhx2_cold_ls, 
        comp2_power_frac_ls
    )
    solve_in_parallel(
        cases=cases,
        f_solve=_analyse_case,
        n_cpu=n_cpu,
    )
    transfer_temp_files_to_results(file_name)

def analyse_unsolved_cases(debug: bool =True):
    """Retry previously unsolved cases and optionally emit debug output."""
    hp: RecuperativeHP
    unsolved_cases = get_unsolved_case_list_to_resolve()
                    
    if unsolved_cases is not None and unsolved_cases != []:
        for case in unsolved_cases:
            args = {'refrigerant': str(case['refrigerant']),
                    'dT_GC_min_tar': float(case['dT_GC_min_tar']), 
                    'T_air_tar': float(case['T_air_tar']),
                    'T_air_inlet': float(case['T_air_inlet']), 
                    'T_source_in': float(case['T_source_in']), 
                    'n_comp': int(case['n_comp']), 
                    'ihx_on': int(case['ihx_on']),
                    'dTmin_rhx': float(case['dTmin_rhx']),
                    'dT_rhx1_cold': float(case['dT_rhx1_cold']),
                    'dT_rhx2_cold': float(case['dT_rhx2_cold']),
                    'dT_rhx3_cold': float(case['dT_rhx3_cold']), 
                    'comp2_power_frac': float(case['comp2_power_frac']), 
                    'comp3_power_frac': float(case['comp3_power_frac']),
                    'comp_PR_uniform': bool(case['comp_PR_uniform']),
                    'comp_model_index': int(case['comp_model_index']),
                    'motor_eff': float(case['motor_eff']),
                    'sink': str(case['sink']),
                    'source': str(case['source']),
                    'Tref_evap_inlet' : float(case['Tref_evap_inlet']),                
                    'debug': bool(debug),
                    }
            
            hp, _ = _analyse_case(args)
            if isinstance(hp, RecuperativeHP):
                if hp.converged and debug:
                    hp.get_report()
                    fsuffix = _get_unique_f_suffix(args) + '-RESOLVED'
                    _print_results(hp.output, fsuffix)
            
                    # hp.get_cycle_plot(x_axis='h', y_axis='T')
                    # hp.get_cycle_plot(x_axis='s', y_axis='T')
                    # hp.get_cycle_plot(x_axis='h', y_axis='P')
                    # hp.get_cycle_plot(x_axis='v', y_axis='P')
        transfer_temp_files_to_results('Cases-Resolved')

@timer_func
def solve_single_case_with_timer(args) -> RecuperativeHP:
    """Instantiate and solve a single RecuperativeHP case while timing execution."""
    hp = RecuperativeHP(**args)
    hp.solve()
    return hp

def get_results_report(hp: RecuperativeHP, param: dict, result: list[dict], obj: float, dt: float, print_to_terminal: bool=False) -> list:
    """Format, persist, and optionally print a solved case summary."""
    res = _get_results_summary(hp, result, param, obj, dt)
    fsuffix = _get_unique_f_suffix(param)
    _print_results(res, fsuffix, print_to_terminal)

#######################################################################################################
# Helper functions
#######################################################################################################

def _analyse_case(case: dict) -> Tuple[RecuperativeHP, float]:
    """Optimise the gas-cooler offset for a single case and solve it."""
    hp: RecuperativeHP
    dt_ls = []
    
    dT_gc_out_min_offset = 0.0 if 'dT_gc_out_min_offset' not in case else case['dT_gc_out_min_offset']

    vars_input = { 
        'dT_gc_out_min_offset' : 
            { 
                'xi' : dT_gc_out_min_offset, 
                'bnds' : (0.0, 40.0) 
            } 
        }
    case['additional_data'] = {
        'vars_input' : vars_input
        }
    
    try:
        # Attempt to calculate the gradient at the boundary, if it fails, 
        # assume that further optimisation could be possible
        fprime, dt = _get_fprime(
            param=case, 
            f_obj=_compute_SEC_obj,
        )
        dt_ls.append(dt)
    except:
        fprime = 1

    if fprime > 0.0:
        result, dt = _get_local_optimum(
            param=case,
            vars_input=vars_input,
            tol=1e-4,
            f_obj=_compute_SEC_obj,
            method='COBYQA',
        )
        dt_ls.append(dt)
        dT_gc_out_min_offset = result.x[0]
        is_solved = result.success
    elif fprime <= 0.0:
        dT_gc_out_min_offset = 0.0
        result = {
            'nfev': 2,
            'message': 'Success',
            'success': True,
        }
        is_solved = True
    else:
        is_solved = False
    
    if is_solved:
        hp, dt = solve_single_case_with_timer(case | _get_var_kwargs(dT_gc_out_min_offset, case))
        dt_ls.append(dt)        
        obj = _compute_SEC_obj(hp, case)

        hp.dQdW = _get_dQdW(hp)
        get_results_report(
            hp=hp,
            param=case, 
            result=result, 
            obj=obj, 
            dt=sum(dt_ls), 
            print_to_terminal=False,
        )
        return hp, dt
    else:
        return None, dt
    
def _get_initial_val(vars_input: dict) -> list:
    """Return the initial optimisation variables in argument order."""
    return [vars_input[key]['xi'] for key in vars_input]

def _get_bounds(vars_input: dict) -> list:
    """Return bounds corresponding to the optimisation variables."""
    return [vars_input[key]['bnds'] for key in vars_input]

@timer_func
def _get_local_optimum(param: dict, vars_input: dict, tol: float, f_obj, method: str) -> dict:
    """Run a local optimisation over the specified variable set."""
    return loc_min(
        fun=lambda x: _get_obj(x, param, f_obj), 
        x0=_get_initial_val(vars_input),
        bounds=_get_bounds(vars_input),
        method=method,
        tol=tol,
    )

@timer_func
def _get_fprime(param: dict, f_obj) -> dict:
    """Estimate the derivative of the objective at the lower bound."""
    return approx_fprime(
        xk=0.0,
        f=lambda x: _get_obj(x, param, f_obj),
        epsilon=0.1
    ) * -1

def _get_obj(vars: list[float], param: dict, f_obj) -> float:
    """Evaluate the objective function for a candidate variable vector."""
    hp: RecuperativeHP
    if not isinstance(vars, list):
        vars = [vars]
    hp, _ = solve_single_case_with_timer(param | _get_var_kwargs(vars, param))
    obj = f_obj(hp, param)
    return obj * 1000

def _compute_SEC_obj(hp: RecuperativeHP, param: dict) -> float:
    """Objective function: specific electricity consumption."""
    return hp.SEC

def _get_var_kwargs(vars: list[float], param: dict) -> dict[str, float]:
    """Map optimisation vector back to keyword arguments."""
    d_var = dict()
    if not isinstance(vars, list):
        vars = [vars]    
    for i, key in enumerate(param['additional_data']['vars_input']):
        d_var[key] = vars[i]
    return d_var

def _print_results(results: dict, fsuffix: str, print_to_terminal: bool =True) -> None:
    """Persist results to CSV and optionally print to stdout."""
    if print_to_terminal:
        print_result_to_terminal(results)
    print_result_to_csv(output=results, suffix=fsuffix, file_name='temp/R')    

def _get_dQdW(hp0: RecuperativeHP) -> float:
    """Finite-difference estimate of heating capacity to work sensitivity."""
    hp1: RecuperativeHP
    dtemp = 0.001
    args = copy.deepcopy(hp0.args)
    args['T_air_tar'] = hp0.T_air_tar + dtemp
    args['P_H_est'] = hp0.comp_discharge.pressure
    args['T_gc_out'] = hp0.gas_cooler1.source_out.temperature
    hp1, _ = solve_single_case_with_timer(args)
    if hp1.feasible:
        tcop = [(p.T_air_tar - p.T_air_inlet) / p.COP if p.COP > 0 else np.inf for p in [hp0, hp1]]
        dqdw = dtemp / (tcop[1] - tcop[0]) if tcop[1] - tcop[0] != 0 and tcop[0] != np.inf and tcop[1] != np.inf else 0.0  
    else:
        dqdw = 0.0
    return dqdw

def _get_results_summary(hp: RecuperativeHP, result, param, obj, dt: float =None) -> dict:
    """Combine optimisation metadata with the heat pump report."""
    hp.get_report()

    if 'weight' in param['additional_data']:
        hp.output.update({
            'Weight' : param['additional_data']['weight'], 
            'SEC_0' : param['additional_data']['SEC_0'],
            'P_0' : param['additional_data']['P_0'],
        })

    hp.output.update({
        'Obj_function' : obj,
        'Function evaluations' : result["nfev"],
        'Message' : result["message"],
        'Success' : result["success"],
        'Time_to_optimise' : dt,
    })
    return hp.output

def _get_unique_f_suffix(param, fld_only=False) -> str:
    """Build a unique filename suffix encoding refrigerant and key inputs."""
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
                suffix = suffix + '-' + '{:.0f}'.format(param[ele])
    return suffix
