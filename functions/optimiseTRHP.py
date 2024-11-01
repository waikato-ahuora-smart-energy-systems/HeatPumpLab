__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from scipy.optimize import shgo, LinearConstraint, minimize
from systems import TRHP_model_functions
from functions.readWriteOutput import print_result_to_terminal, print_result_to_csv
from time import perf_counter as timer
from concurrent.futures import ProcessPoolExecutor as pool
import numpy as np, os

def dict_var(var: list[float], parameters) -> dict[str, float]:
    if parameters['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    d_var = {'dT_phx1_cold': var[0]}
    if parameters['n_comp'] > 1:
        d_var['dT_phx2_cold'] = var[1]
    if parameters['n_comp'] > 2:
        d_var['dT_phx3_cold'] = var[2]
    if parameters['comp_PR_uniform'] == False and parameters['n_comp'] > 1:
        if parameters['n_comp'] == 2:
            d_var['comp2_power_frac'] = var[-1]
        elif parameters['n_comp'] == 3:
            d_var['comp2_power_frac'] = var[-2]
            d_var['comp3_power_frac'] = var[-1]

    return d_var

def get_initial(res, parameters):
    if parameters['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    var = [res['dT phx max 1:']]
    if parameters['n_comp'] > 1:
        var.append(res['dT phx max 2:'])
    if parameters['n_comp'] > 2:
        var.append(res['dT phx max 3:'])
    if parameters['comp_PR_uniform'] == False and parameters['n_comp'] > 1:
        if parameters['n_comp'] == 2:
            var.append(res['Comp 2 power frac:'])
        elif parameters['n_comp'] == 3:
            var.append(res['Comp 2 power frac:'])
            var.append(res['Comp 3 power frac:'])

    return var  

def get_bounds(parameters, x=None):
    #Define the bounds of the optimisation variables 
    dT_upper = 100.0
    if parameters['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    bnds = [(0.0, dT_upper)]
    xi = [0.0] if x is None else [x[0]]
    if parameters['n_comp'] >= 2:
        bnds += [(0.0, dT_upper)]
        xi += [0.0] if x is None else [x[1]]
    if parameters['n_comp'] == 3:
        bnds += [(0.0, dT_upper)]
        xi += [0.0] if x is None else [x[2]] 
    if parameters['comp_PR_uniform'] == False and parameters['n_comp'] > 1:
        if parameters['n_comp'] >= 2:
            bnds += [(0.25, 0.75)]
            xi += [0.5]
        if parameters['n_comp'] == 3:
            bnds += [(0.25, 0.75)]
            xi += [0.5] 

    return bnds, xi

def get_max_multiObj(var, parameters):
    k = parameters['result_parameters']['weight']
    COP_0 = parameters['result_parameters']['COP_0']
    P_0 = parameters['result_parameters']['P_0']
    
    tic = timer()
    d_var = dict_var(var, parameters)
    m = TRHP_model_functions(**parameters, **d_var)
    COP, P = m.solve()
    obj = -(COP / COP_0) * k + (P / P_0 - 2) * (1 - k)
    toc = timer()
    
    print('Obj:', "{:.4f}".format(-obj), '  Solve time:', "{:.2f}".format(toc - tic), '  Vars:', ["{:.6f}".format(v) for v in var])    
    return obj 

def get_results_summary(result, parameters):
    d_var = dict_var(result.x, parameters)
    m = TRHP_model_functions(**parameters, **d_var)    
    COP, P = m.solve()
    m.get_summary()

    k = parameters['result_parameters']['weight']
    COP_0 = parameters['result_parameters']['COP_0']
    P_0 = parameters['result_parameters']['P_0']
    COP, P = m.solve()
    obj = -(COP / COP_0) * k + (P / P_0 - 2) * (1 - k)

    m.output.update({
        'Weight:' : parameters['result_parameters']['weight'], 
        'COP_0:' : parameters['result_parameters']['COP_0'],
        'P_0:' : parameters['result_parameters']['P_0'],
        'Obj function:' : obj
    })

    if result != None:
        m.output.update({
            'Iterations:' : result.nfev,
            'Message:' : result.message,
            'Success:' : result.success
        })

    for ele in parameters:
        if ele == 'refrigerant':
            suffix = '-' + str(parameters[ele])
        else:
            if ele != 'result_parameters':
                suffix = suffix + '-' + "{:.0f}".format(parameters[ele])

    print_result_to_terminal(m.output)
    print_result_to_csv(output=m.output, suffix=suffix)

    return m.output

def get_global_optimum(parameters, bnds, tol):
    return shgo(func=get_max_multiObj, 
                args=[parameters], 
                bounds=bnds,
                constraints=LinearConstraint(A=np.array([1.0] * len(bnds)), 
                                             lb=0.0, 
                                             ub=(parameters['T_air_tar'] - parameters['T_air_inlet']) * 2/3),
                workers=1, 
                iters=3, 
                n=40, 
                sampling_method='sobol', 
                options={'tol': tol}, 
                minimizer_kwargs={'method' : 'SLSQP',
                                  'options': {'tol': tol}}
                )    

def get_local_optima(method, parameters, bnds, tol, xi):
    return minimize(fun=get_max_multiObj, 
                x0=xi,
                args=parameters, 
                bounds=bnds,
                constraints=LinearConstraint(A=np.array([1.0] * len(bnds)), 
                                             lb=0.0, 
                                             ub=(parameters['T_air_tar'] - parameters['T_air_inlet']) * 2/3),
                method=method,
                tol=tol / 10,
                )

def get_TRHP_optimal(parameters):
    tol = 1e-3
    T_air_tar_ls = parameters['T_air_tar']
    results_ls = parameters['results']
    parameters.pop('results')
    for i in range(len(results_ls), len(T_air_tar_ls)):
        parameters['T_air_tar'] = T_air_tar_ls[i]
        bnds, x0 = get_bounds(parameters)
        d_x0 = dict_var(x0, parameters)
        m = TRHP_model_functions(**parameters, **d_x0)
        COP_0, P_0 = m.solve()
        parameters['result_parameters']['COP_0'] = COP_0
        parameters['result_parameters']['P_0'] = P_0

        if i == 0 or 1:    
            result = get_global_optimum(parameters, bnds, tol) 
        else:
            xi = get_initial(results_ls[i-1], parameters)
            result = get_local_optima('trust-constr', parameters, bnds, tol, xi) #'SLSQP',
            if result.success and -result.fun < 1.0:
                result = get_local_optima('trust-constr', parameters, bnds, tol, x0)

        if result.success:
            res = get_results_summary(result, parameters)
            results_ls.append(res)

def get_parallel_solve(cases):
    # Use parallel processing distributed over all cores to speed up the analysis
    n_cpu: int
    executor: pool
    results: list[any]

    n_cpu = os.cpu_count()
    executor = pool(max_workers=n_cpu)
    results = [executor.submit(get_TRHP_optimal, case) for case in cases]
    executor.shutdown()

