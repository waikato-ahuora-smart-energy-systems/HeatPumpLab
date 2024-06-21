__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from scipy.optimize import shgo, LinearConstraint
from heat_pumps import TRHP
from outputs import print_result_to_terminal, print_result_to_csv
from timeit import default_timer as timer
import numpy as np

def dict_var(var, parameters):
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
        xi += [0.0] if x is None else [x[1]] 
    if parameters['comp_PR_uniform'] == False and parameters['n_comp'] > 1:
        if parameters['n_comp'] >= 2:
            bnds += [(0.4, 0.6)]
            xi += [0.5]
        if parameters['n_comp'] == 3:
            bnds += [(0.4, 0.6)]
            xi += [0.5] 

    return bnds, xi

def get_max_multiObj(var, parameters):
    w = parameters['result_parameters']['weight']
    COP_0 = parameters['result_parameters']['COP_0']
    P_0 = parameters['result_parameters']['P_0']
    
    TIC = timer()
    d_var = dict_var(var, parameters)
    m = TRHP(**parameters, **d_var)
    COP, P = m.solve()
    obj = (COP / COP_0) * w - (P / P_0) * (1 - w)
    TOC = timer()
    
    print('Obj:', "{:.4f}".format(obj), '  Solve time:', "{:.2f}".format(TOC - TIC), '  Vars:', ["{:.6f}".format(v) for v in var])    
    return -obj 

def get_results_summary(result, x, parameters):
    d_var = dict_var(x, parameters)
    m = TRHP(**parameters, **d_var)    
    COP, P = m.solve()
    m.get_summary()

    w = parameters['result_parameters']['weight']
    COP_0 = parameters['result_parameters']['COP_0']
    P_0 = parameters['result_parameters']['P_0']
    COP, P = m.solve()
    obj = (COP / COP_0) * (w) + (1 - P / P_0) * (1 - w)

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

def get_TRHP_optimal(parameters):
    tol = 1e-4
    parameters.update({'comp_PR_uniform': False})

    bnds, xi = get_bounds(parameters)
    d_xi = dict_var(xi, parameters)
    m = TRHP(**parameters, **d_xi)
    COP_0, P_0 = m.solve()
    parameters['result_parameters']['COP_0'] = COP_0
    parameters['result_parameters']['P_0'] = P_0
    
    if COP_0 > 0.0:        
        result = shgo(func=get_max_multiObj, 
                      args=[parameters], 
                      bounds=bnds,
                      constraints=LinearConstraint(A=np.array([1.0] * len(xi)), 
                                                   lb=0.0, 
                                                   ub=(parameters['T_air_tar'] - parameters['T_air_inlet']) * 2/3),
                      workers=1, 
                      iters=4, 
                      #n=50, 
                      sampling_method='sobol', 
                      options={'tol': tol}, 
                      minimizer_kwargs={'method' : 'SLSQP',
                                        'options': {'tol': tol}})

        print('Obj:', -result.fun, 
              '  Iter:', result.nfev, 
              '  Vars:' ,result.x)

        if result.success:
            get_results_summary(result, result.x, parameters)

    # m.get_cycle_plot(x_axis='h', y_axis='P')
    # m.get_cycle_plot(x_axis='s', y_axis='T')
    # m.get_cycle_plot(x_axis='h', y_axis='T')
    # m.get_cycle_plot(x_axis='v', y_axis='P')