__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from scipy.optimize import minimize, differential_evolution as DE_opt, dual_annealing as DA_opt, shgo, direct as DI_opt, LinearConstraint
from heat_pumps import TRHP
from outputs import print_result_to_terminal, print_result_to_csv
from timeit import default_timer as timer
from streams.properties import get_property, Properties
from math import log
import numpy as np

def dict_var(var, parameters):
    if parameters['n_comp'] > 3:
        raise ValueError("Number of compressors must be 1, 2 or 3.")     
    
    d_var = {'dT_phx1_cold': var[0]}
    if parameters['n_comp'] == 2:
        d_var['dT_phx2_cold'] = var[1]
    if parameters['n_comp'] == 3:
        d_var['dT_phx3_cold'] = var[2]
    if parameters['comp_PR_uniform'] == False and parameters['n_comp'] > 1:
        if parameters['n_comp'] == 2:
            d_var['comp2_power_frac'] = var[-1]
        elif parameters['n_comp'] == 3:
            d_var['comp2_power_frac'] = var[-2]
            d_var['comp3_power_frac'] = var[-1]

    return d_var

def get_bounds_constraints(parameters, x=None, include_P=False):
    #Define the bounds of the optimisation variables 
    dT_upper = 100.0 
    dT_int = (parameters['T_air_tar'] - parameters['T_air_inlet']) * 2/3  
    if parameters['comp_PR_uniform'] == True:
        if parameters['n_comp'] == 1:
            bnds = [(0.0, dT_upper)]
            xi = [0.0]
            lc = LinearConstraint([[1]], 0.0, dT_int)

        elif parameters['n_comp'] == 2:
            bnds = [(0.0, dT_upper), (0.0, dT_upper)]
            xi = [0.0, 0.0] if x is None else [x[0], x[1]]
            lc = LinearConstraint([[1, 1]], 0.0, dT_int)

        elif parameters['n_comp'] == 3:
            bnds = [(0.0, dT_upper), (0.0, dT_upper), (0.0, dT_upper)]
            xi = [0.0, 0.0, 0.0] if x is None else [x[0], x[1], x[2]]
            lc = LinearConstraint([[1, 1, 1]], 0.0, dT_int)
        
        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.") 

    else:
        if parameters['n_comp'] == 1:
            bnds = [(0.0, dT_upper)]
            xi = [0.0]
            lc = LinearConstraint([[1]], 0.0, dT_int)

        elif parameters['n_comp'] == 2:
            bnds = [(0.0, dT_upper), (0.0, dT_upper), (0.4, 0.6)]
            xi = [0.0, 0.0, 50.0] if x is None else [x[0], x[1], 0.5]
            lc = LinearConstraint([[1, 1, 0]], 0.0, dT_int)

        elif parameters['n_comp'] == 3:
            bnds = [(0.0, dT_upper), (0.0, dT_upper), (0.0, dT_upper), (0.2333, 0.4333), (0.2333, 0.4333)]
            xi = [0.0, 0.0, 0.0, 1/3, 1/3] if x is None else [x[0], x[1], x[2], 1/3, 1/3]
            lc = LinearConstraint([[1, 1, 1, 0, 0]], 0.0, dT_int)
        
        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.") 

    if include_P:
        RP = Properties().init_properties()
        P_crit = get_property('PCRIT', parameters['refrigerant'], 'TQ', 100, 0, RP)
        P_low = get_property('P', parameters['refrigerant'], 'TQ', parameters['Tsat_evap'], 1.0, RP)
        Log_P_h_min = log(P_low * 2)
        Log_P_h_max = log(P_crit * 2.5)
        bnds = bnds + [(Log_P_h_min, Log_P_h_max)]

    return bnds, lc, xi

def get_max_multiObj(var, parameters):
    w = parameters['result_parameters']['weight']
    COP_0 = parameters['result_parameters']['COP_0']
    P_0 = parameters['result_parameters']['P_0']
    
    TIC = timer()
    d_var = dict_var(var, parameters)
    m = TRHP(**parameters, **d_var)
    COP, P = m.solve()
    obj = (COP / COP_0) * (w) + (1 - P / P_0) * (1 - w)
    TOC = timer()    
    
    print('Obj:', "{:.4f}".format(obj), '  Solve time:', "{:.3f}".format(TOC - TIC), '  Vars:', ["{:.8f}".format(v) for v in var])   
    return obj * -1000

def get_results_summary(result, x, parameters):
    if parameters['comp_PR_uniform'] == True:
        if parameters['n_comp'] == 1:
            m = TRHP(
                **parameters,
                dT_phx1_cold=x[0]
                )
        elif parameters['n_comp'] == 2:   
            m = TRHP(
                **parameters,
                dT_phx1_cold=x[0], 
                dT_phx2_cold=x[1]
                )        
        elif parameters['n_comp'] == 3:
            m = TRHP(
                **parameters,
                dT_phx1_cold=x[0], 
                dT_phx2_cold=x[1], 
                dT_phx3_cold=x[2]
                )
            
    else:  
        if parameters['n_comp'] == 1:
            m = TRHP(
                **parameters,
                dT_phx1_cold=x[0]
                )
        elif parameters['n_comp'] == 2:   
            m = TRHP(
                **parameters,
                dT_phx1_cold=x[0], 
                dT_phx2_cold=x[1],  
                comp2_power_frac=x[2]
                )        
        elif parameters['n_comp'] == 3:
            m = TRHP(
                **parameters,
                dT_phx1_cold=x[0], 
                dT_phx2_cold=x[1], 
                dT_phx3_cold=x[2], 
                comp2_power_frac=x[3], 
                comp3_power_frac=x[4]
                )
    
    m.solve()
    m.get_summary()
    if result != None:
        m.output.update({
            'weight' : parameters['result_parameters']['weight'], 
            'COP_max' : parameters['result_parameters']['COP_0'],
            'P_max' : parameters['result_parameters']['P_0'],
            'Obj function:' : -result.fun / 1000,
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
    parameters.update({'comp_PR_uniform': True})
    bnds, lc, xi = get_bounds_constraints(parameters)
    d_xi = dict_var(xi, parameters)

    m = TRHP(**parameters, **d_xi)
    COP_0, P_0 = m.solve()
    # m.get_cycle_plot(x_axis='h', y_axis='P')
    # m.get_cycle_plot(x_axis='s', y_axis='T')
    # m.get_cycle_plot(x_axis='h', y_axis='T')
    # m.get_cycle_plot(x_axis='v', y_axis='P')
    
    # if COP_0 > 0.0:
    #     get_results_summary(None, xi, parameters)    

    parameters['result_parameters'] = {'weight': 1.0, 
                                       'COP_0':  COP_0,
                                       'P_0':    P_0}
    for j in range(10):
        TIC = timer()
        result = DE_opt(func=get_max_multiObj, args=[parameters], bounds=bnds, constraints=lc, workers=-1, polish=False, tol=tol) 
        TOC = timer()
        dt = TOC - TIC

        print('SH - Obj function:', -result.fun / 1000, '  time taken:', dt,  '  iter:', result.nfev, result.x)

        if result.success:
            get_results_summary(result, result.x, parameters)

        parameters['result_parameters']['weight'] -= 0.05

    # TIC = timer()
    # result = shgo(func=get_max_multiObj, args=[parameters], bounds=bnds, workers=1, iters=1, constraints=lc, options={'tol': tol}, minimizer_kwargs={'options': {'tol': tol}}) 
    # TOC = timer()
    # dt = TOC - TIC 

    # TIC = timer()
    # result1 = DE_opt(func=get_max_multiObj, args=[parameters], bounds=bnds, constraints=lc, workers=-1, popsize=30, init='sobol', polish=False, tol=tol)
    # TOC = timer()
    # dt1 = TOC - TIC

    # TIC = timer()
    # result2 = DE_opt(func=get_max_multiObj, args=[parameters], bounds=bnds, constraints=lc, workers=-1, init='sobol', polish=False, tol=tol) 
    # # result2 = minimize(method='COBYLA', fun=get_max_multiObj, args=parameters, bounds=bnds, constraints=lc, x0=result1.x, tol=1e-4) 
    # TOC = timer()
    # dt2 = TOC - TIC 

    # TIC = timer()
    # result3 = DE_opt(func=get_max_multiObj, args=[parameters], bounds=bnds, constraints=lc, workers=-1, popsize=30, polish=False, tol=tol) 
    # # result3 = minimize(method='SLSQP', fun=get_max_multiObj, args=parameters, bounds=bnds, x0=result1.x, constraints=lc, tol=1e-6) 
    # TOC = timer()
    # dt3 = TOC - TIC     

    # TIC = timer()
    # result4 = DE_opt(func=get_max_multiObj, args=[parameters], bounds=bnds, constraints=lc, workers=-1, polish=False, tol=tol) 
    # # result4 = minimize(method='l-bfgs-b', fun=get_max_multiObj, args=parameters, bounds=bnds, x0=result1.x, tol=1e-6) 
    # TOC = timer()
    # dt4 = TOC - TIC    

    # TIC = timer()
    # result5 = minimize(method='nelder-mead', fun=get_max_multiObj, args=parameters, bounds=bnds, x0=xi, tol=1e-4) 
    # TOC = timer()
    # dt5 = TOC - TIC 

    # TIC = timer()
    # result6 = minimize(method='powell', fun=get_max_multiObj, args=parameters, bounds=bnds, x0=xi, tol=1e-4) 
    # TOC = timer()
    # dt6 = TOC - TIC 

    # print('Case:', parameters)
    # print('SH - Obj function:', -result.fun / 1000, '  time taken:', dt,  '  iter:', result.nfev, result.x)
    # print('1 - Obj function:', -result1.fun / 1000, '  time taken:', dt1,  '  iter:', result1.nfev, result1.x)
    # print('2 - Obj function:', -result2.fun / 1000, '  time taken:', dt2,  '  iter:', result2.nfev, result2.x)
    # print('3 - Obj function:', -result3.fun / 1000, '  time taken:', dt3,  '  iter:', result3.nfev, result3.x)
    # print('4 - Obj function:', -result4.fun / 1000, '  time taken:', dt4,  '  iter:', result4.nfev, result4.x)
    # print('NM - Obj function:', -result5.fun / 1000, '  time taken:', dt5,  '  iter:', result5.nfev, result5.x)
    # print('PO - Obj function:', -result6.fun / 1000, '  time taken:', dt6,  '  iter:', result6.nfev, result6.x)
    
    # if result.fun < 0.0:
    #     return get_results_summary(result, parameters)

    # SH - Obj function: 1.0920403900035898   time taken: 291.5808960830036   iter: 337 [84.64207665  2.0065959   2.02504565]
    # DE - Obj function: 1.0967210448531715   time taken: 205.5725983750017   iter: 282 [ 6.75369068 57.06274948  9.5667201 ]
    # CO - Obj function: 1.081442095804895   time taken: 69.20527524998761   iter: 85 [ 8.11770086 14.74049575 27.44669489]
    # SQ - Obj function: 1.087246330258263   time taken: 70.97270229199785   iter: 98 [43.84438408 35.99317723  0.        ]
    # BF - Obj function: 1.0595740417601036   time taken: 19.200484708009753   iter: 32 [ 0.         56.97558247 64.93886961]
    # NM - Obj function: 1.0817586292773815   time taken: 133.22085754200816   iter: 169 [ 6.29425085 14.57161934 27.94951824]
    # PO - Obj function: 1.0892481746100446   time taken: 203.31527549997554   iter: 198 [8.86124212e+01 5.36344546e-05 1.41985683e-03]
    # DE - Obj function: 1.0944394837116338   time taken: 206.553638124984   iter: 281 [10.48118176 56.37081546 13.74189146]
    # SQ - Obj function: 1.096722148624305   time taken: 68.31465733301593   iter: 90 [ 9.57320536 57.88912903 11.5953677 ]
    # DE - Obj function: 1.0966082201741973   time taken: 196.75819420898915   iter: 265 [ 9.50080573 57.88589885 11.75163457]
    # BF - Obj function: 1.0966375611071375   time taken: 47.731355333002284   iter: 60 [ 9.50852004 57.90703301 11.73876879] 

    # 1 - Obj function: 1.0915244080082995   time taken: 57.852384124998935   iter: 643 [87.05316307 25.72662237  3.06535077]
    # 2 - Obj function: 1.091004199524171   time taken: 22.572867000009865   iter: 255 [86.84299188 47.91791973  3.41171639]
    # 3 - Obj function: 1.0867617963993943   time taken: 28.0631263749965   iter: 288 [86.62229785 18.16989943 16.15693834]
    # 4 - Obj function: 1.0917191333980216   time taken: 26.13374516699696   iter: 254 [87.04834626  2.03127733  2.40630719]

    # 1 - Obj function: 1.0909664643762378   time taken: 55.44563954201294   iter: 599 [87.0855029  16.16188493  5.09221192]
    # 2 - Obj function: 1.0914616019811294   time taken: 35.29208095901413   iter: 359 [87.26783893 38.60218709  3.89127883]
    # 3 - Obj function: 1.0912989273159475   time taken: 35.56015595901408   iter: 400 [87.60361832 29.5333796   1.40794312]
    # 4 - Obj function: 1.0861068135564136   time taken: 15.084947208000813   iter: 157 [87.09915641 33.63796361 22.0068382 ]  
    # 4 - Obj function: 1.0920981096228377   time taken: 93.3139755000011   iter: 845 [87.19684774 20.05521968  1.94944369] 