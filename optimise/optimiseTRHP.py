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
    
    TIC = timer()
    d_var = dict_var(var, parameters)
    m = TRHP(**parameters, **d_var)
    COP, P = m.solve()
    obj = -(COP / COP_0) * k + (P / P_0 - 2) * (1 - k)
    TOC = timer()
    
    print('Obj:', "{:.4f}".format(-obj), '  Solve time:', "{:.2f}".format(TOC - TIC), '  Vars:', ["{:.6f}".format(v) for v in var])    
    return obj 

def get_results_summary(result, x, parameters):
    d_var = dict_var(x, parameters)
    m = TRHP(**parameters, **d_var)    
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
                      n=50, 
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

    # Obj: 1.183411077993845   Iter: 456   Vars: [15.83295351 73.8423853   0.6       ]
    # COP: 2.2564099485166604
    # Tol: 6.320988177321851e-13
    # Converged: True
    # refrigerant: Butane
    # n_comp: 2
    # T_source_in: 40.0
    # dT_min Gas Cooler: 20.0
    # T_air_inlet: 15.0
    # T_air_tar: 250.0
    # Actual air T: 249.99999999984198
    # Evap T: 10.0
    # Evap P: 148.45189997734684
    # Comp 1 P out: 638.4188066475795
    # Comp 2 P out: 5693.5814731072205
    # Comp 3 P out: 5693.5814731072205
    # Comp 1 power frac: 0.4
    # Comp 2 power frac: 0.6
    # Comp 3 power frac: 0.0
    # dT phx max 1: 15.832953513302762
    # dT phx max 2: 73.84238529969029
    # dT phx max 3: 0.0
    # PR 1: 4.300509503381227
    # PR 2: 8.918254621922811
    # PR 3: 1.0
    # Overall P ratio: 38.353038755152596
    # Gas cooler(s) inlet T: 269.99999999984243
    # Gas cooler(s) outlet T: 35.0
    # Weight: 1.0
    # COP_0: 1.906700038960085
    # P_0: 41017.075411147205
    # Obj function: 1.1834110779932776
    # Iterations: 456
    # Message: Optimization terminated successfully.
    # Success: True