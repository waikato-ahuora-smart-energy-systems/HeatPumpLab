__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from systems import TRHP, set_trhp_model
from streams.propertiesPackage import PropertyPackage
from functions.solver import solve_heat_pump_model
from functions.cases import get_cases_to_solve
from functions.readWriteOutput import print_result_to_terminal, print_result_to_csv, transfer_temp_files_to_results, read_csv_to_results
from functions.decorator_functions import time_func
from concurrent.futures import ProcessPoolExecutor as pool
from scipy.optimize import shgo as glb_opt_mthd, LinearConstraint, minimize
import numpy as np, os

def analyse_single_case():
    args: list[float | str | bool]
    hp: TRHP

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
    
    hp, dt = solve_single_case(args)
    print('Time to solve case:', dt)
    if hp.converged:
        hp.get_summary()

        if type(hp.output) is dict:
            print_result_to_terminal(output=hp.output)
            # print_result_to_csv(output=hp.output)
        
        # hp.get_driving_force_plot(hp.gas_coolers, hp.dT_GC_min_tar)
        # hp.get_cycle_plot(x_axis='h', y_axis='T')
        # hp.get_cycle_plot(x_axis='s', y_axis='T')
        # hp.get_cycle_plot(x_axis='h', y_axis='P')
        # hp.get_cycle_plot(x_axis='v', y_axis='P')
    else:
        print('Heat model failed to solve.')

def analyse_unsolved_cases(debug: bool =True):
    unsolved = read_csv_to_results(file_name='results/Cases-Unsolved')
    for case in unsolved:
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
        
def analyse_cases_in_parallel(file_name: str = 'results/TRHP_results', n_cpu: int = os.cpu_count()) -> None:
    # Use parallel processing distributed over all cores to speed up the analysis
    cases: list        = get_cases_to_solve(file_name)
    print(len(cases))
    executor: pool     = pool(max_workers=n_cpu)
    results: list[any] = [executor.submit(_optimise_case_set, case) for case in cases]
    executor.shutdown()
    transfer_temp_files_to_results(file_name)

def _optimise_case_set(case: dict) -> None:
    tol: float = 1e-3
    T_air_tar_ls: list = case['T_air_tar']
    results_ls: list = case['results']
    case.pop('results')
    start: int = len(results_ls)
    end: int = len(T_air_tar_ls)
    param: dict = case
    for i in range(start, end):
        param['T_air_tar'] = T_air_tar_ls[i]
        bnds, var0 = _get_bounds(param)
        hp, _ = solve_single_case(
            _get_hp_args(var0, param)
            )
        param['result_parameters']['COP_0'] = hp.COP
        param['result_parameters']['P_0'] = hp.comp_discharge.pressure

        if i == 0:    
            result = _get_global_optimum(param, bnds, tol) 
        else:
            method = 'trust-constr' #'SLSQP'
            xi = _get_initial(results_ls[i-1], param)
            result = _get_local_optima(method, param, bnds, tol, xi) 
            if result.success and -result.fun < 1.0:
                result = _get_local_optima(method, param, bnds, tol, var0)

        if result.success:
            res = _get_results_summary(result, param)
            fsuffix = _get_unique_f_suffix(param)
            _print_results(res, fsuffix)
            results_ls.append(res)

def _get_initial(res: dict, param: dict) -> list:
    if param['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")
    var = [res['dT phx max 1:']]
    if param['n_comp'] > 1:
        var.append(res['dT phx max 2:'])
    if param['n_comp'] > 2:
        var.append(res['dT phx max 3:'])
    if param['comp_PR_uniform'] == False and param['n_comp'] > 1:
        if param['n_comp'] == 2:
            var.append(res['Comp 2 power frac:'])
        elif param['n_comp'] == 3:
            var.append(res['Comp 2 power frac:'])
            var.append(res['Comp 3 power frac:'])
    return var  

def _get_bounds(param: dict, x: list =None) -> list:
    #Define the bounds of the optimisation variables 
    dT_upper = 100.0
    if param['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    bnds = [(0.0, dT_upper)]
    xi = [0.0] if x is None else [x[0]]
    if param['n_comp'] >= 2:
        bnds += [(0.0, dT_upper)]
        xi += [0.0] if x is None else [x[1]]
    if param['n_comp'] == 3:
        bnds += [(0.0, dT_upper)]
        xi += [0.0] if x is None else [x[2]] 
    if param['comp_PR_uniform'] == False and param['n_comp'] > 1:
        if param['n_comp'] >= 2:
            bnds += [(0.25, 0.75)]
            xi += [0.5]
        if param['n_comp'] == 3:
            bnds += [(0.25, 0.75)]
            xi += [0.5] 

    return bnds, xi

def _get_global_optimum(param: dict, bnds: tuple, tol: float) -> dict:
    return glb_opt_mthd(func=_get_obj_val, 
                        args=[param], 
                        bounds=bnds,
                        constraints=LinearConstraint(A=np.array([1.0] * len(bnds)), 
                                                    lb=0.0, 
                                                    ub=(param['T_air_tar'] - param['T_air_inlet']) * 2/3),
                        workers=1, 
                        iters=3, 
                        n=40, 
                        sampling_method='sobol', 
                        options={'tol': tol}, 
                        minimizer_kwargs={'method' : 'SLSQP',
                                        'options': {'tol': tol}}
                        )    

def _get_local_optima(method: str, param: dict, bnds: tuple, tol: float, xi: float) -> dict:
    return minimize(fun=_get_obj_val, 
                    x0=xi,
                    args=param, 
                    bounds=bnds,
                    constraints=LinearConstraint(A=np.array([1.0] * len(bnds)), 
                                                lb=0.0, 
                                                ub=(param['T_air_tar'] - param['T_air_inlet']) * 2/3),
                    method=method,
                    tol=tol / 10,
                    )

def _get_obj_val(var: dict, param: dict) -> float:
    _, obj, dt = _get_obj(var, param)
    print('Obj:', "{:.4f}".format(-obj), '  Solve time:', "{:.2f}".format(dt), '  Vars:', ["{:.6f}".format(v) for v in var], '  Ref:', param['refrigerant'])    
    return obj 

def _get_obj(var: dict, param: dict) -> any:
    hp: TRHP
    COP: float
    P: float
    COP_0: float = param['result_parameters']['COP_0']
    P_0: float   = param['result_parameters']['P_0']
    k: float     = param['result_parameters']['weight']
    
    hp, dt = solve_single_case(
        _get_hp_args(var, param)
        )
    COP, P = hp.COP, hp.comp_discharge.pressure
    obj = -(COP / COP_0) * k + (P / P_0 - 2) * (1 - k) 
    return hp, obj, dt

def _get_hp_args(var: list[float], param) -> dict[str, float]:
    if param['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    d_var = {'dT_phx1_cold': var[0]}
    if param['n_comp'] > 1:
        d_var['dT_phx2_cold'] = var[1]
    if param['n_comp'] > 2:
        d_var['dT_phx3_cold'] = var[2]
    if param['comp_PR_uniform'] == False and param['n_comp'] > 1:
        if param['n_comp'] == 2:
            d_var['comp2_power_frac'] = var[-1]
        elif param['n_comp'] == 3:
            d_var['comp2_power_frac'] = var[-2]
            d_var['comp3_power_frac'] = var[-1]

    return param | d_var

@time_func
def solve_single_case(args) -> TRHP:
    hp = TRHP(**args)
    hp = set_trhp_model(hp)
    if hp.feasible:
        hp = solve_heat_pump_model(hp)
    return hp

def _get_results_summary(result, param) -> dict:
    hp: TRHP
    hp, obj, _ = _get_obj(result.x, param)
    hp.get_summary()
    hp.output.update({
        'Weight:' : param['result_parameters']['weight'], 
        'COP_0:' : param['result_parameters']['COP_0'],
        'P_0:' : param['result_parameters']['P_0'],
        'Obj function:' : obj
    })
    hp.output.update({
        'Iterations:' : result.nfev,
        'Message:' : result.message,
        'Success:' : result.success
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

    # self.fld_name, self.fld_composition, self.n_components, self.iMass

    for ele in param:    
        if ele != 'result_parameters' and ele != 'refrigerant':
            suffix = suffix + '-' + "{:.0f}".format(param[ele])
    return suffix

def _print_results(results: dict, fsuffix: str) -> None:
    print_result_to_terminal(results)
    print_result_to_csv(output=results, suffix=fsuffix)    

