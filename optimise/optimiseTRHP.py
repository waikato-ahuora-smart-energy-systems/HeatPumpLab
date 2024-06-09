__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from scipy.optimize import minimize, differential_evolution as DE_opt, dual_annealing as DA_opt, shgo, direct as DI_opt, basinhopping as BH_opt
from heat_pumps import TRHP
from outputs import print_result_to_terminal, print_result_to_csv
from timeit import default_timer as timer
from streams.properties import get_property
from math import log, exp

def create_m(var, parameters):
    if parameters['comp_PR_uniform'] == True:
        if parameters['n_comp'] == 1:
            m = TRHP(
                **parameters,
                dT_phx1_cold=var[0]
                )
            
        elif parameters['n_comp'] == 2:
            m = TRHP(
                **parameters,
                dT_phx1_cold=var[0], 
                dT_phx2_cold=var[1]
                )

        elif parameters['n_comp'] == 3:
            m = TRHP(
                **parameters,
                dT_phx1_cold=var[0], 
                dT_phx2_cold=var[1], 
                dT_phx3_cold=var[2]
            )

        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.") 
        
    else:
        if parameters['n_comp'] == 1:
            m = TRHP(
                **parameters,
                dT_phx1_cold=var[0]
                )
            
        elif parameters['n_comp'] == 2:
            m = TRHP(
                **parameters,
                dT_phx1_cold=var[0], 
                dT_phx2_cold=var[1], 
                comp2_power_frac=var[2]/100
                )

        elif parameters['n_comp'] == 3:
            m = TRHP(
                **parameters,
                dT_phx1_cold=var[0], 
                dT_phx2_cold=var[1], 
                dT_phx3_cold=var[2], 
                comp2_power_frac=var[3]/100, 
                comp3_power_frac=var[4]/100
            )
    
        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.")    
    
    return m

def get_bounds_constraints(parameters, x=None, include_P=False):
    #Define the bounds of the optimisation variables    
    if parameters['comp_PR_uniform'] == True:
        if parameters['n_comp'] == 1:
            bnds = [(0.0, 60.0)]
            xi = [10.0] 

        elif parameters['n_comp'] == 2:
            bnds = [(0.0, 60.0), (0.0, 60.0)]
            xi = [10.0, 10.0] if x is None else [x[0], x[1]]

        elif parameters['n_comp'] == 3:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 60.0)]
            xi = [10.0, 10.0, 10.0] if x is None else [x[0], x[1], x[2]]
        
        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.") 

    else:
        if parameters['n_comp'] == 1:
            bnds = [(0.0, 60.0)]
            xi = [10.0]

        elif parameters['n_comp'] == 2:
            bnds = [(0.0, 60.0), (0.0, 60.0), (25.0, 75.0)]
            xi = [10.0, 10.0, 50.0] if x is None else [x[0], x[1], 100/2]

        elif parameters['n_comp'] == 3:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 60.0), (15.0, 45.0), (15.0, 45.0)]
            xi = [10.0, 10.0, 10.0, 100/3, 100/3] if x is None else [x[0], x[1], x[2], 100/3, 100/3]
        
        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.") 

    if include_P:
        P_crit = get_property('PCRIT', parameters['refrigerant'], 'TQ', 100, 0)
        P_low = get_property('P', parameters['refrigerant'], 'TQ', parameters['Tsat_evap'], 1.0)
        Log_P_h_min = log(P_low * 2)
        Log_P_h_max = log(P_crit * 2.5)
        bnds = bnds + [(Log_P_h_min, Log_P_h_max)]

    cons = ({'type': 'ineq', 'fun': g, 'args': [parameters]})

    return bnds, cons, xi

def get_TRHP_obj_fun(var, parameters):
    TIC = timer()
    m = create_m(var, parameters)
    COP = m.solve()
    TOC = timer()
    print('COP:', "{:.4f}".format(COP), '  dt:', "{:.1f}".format(TOC - TIC), '  ', ["{:.3f}".format(v) for v in var])
    return COP * -1000

def get_results_summary(result, parameters):  
    if parameters['n_comp'] == 1:
        m = TRHP(
            **parameters,
            dT_phx1_cold=result.x[0]
            )
    elif parameters['n_comp'] == 2:   
        m = TRHP(
            **parameters,
            dT_phx1_cold=result.x[0], 
            dT_phx2_cold=result.x[1],  
            comp2_power_frac=result.x[2]
            )        
    elif parameters['n_comp'] == 3:
        m = TRHP(
            **parameters,
            dT_phx1_cold=result.x[0], 
            dT_phx2_cold=result.x[1], 
            dT_phx3_cold=result.x[2], 
            comp2_power_frac=result.x[3], 
            comp3_power_frac=result.x[4]
            )
    
    m.solve()
    m.get_summary()
    m.output.update({
        'Obj function:' : -result.fun / 1000,
        'Iterations:' : result.nfev,
        'Message:' : result.message,
        'Success:' : result.success
    })

    for ele in parameters:
        if ele == 'refrigerant':
            suffix = '-' + str(parameters[ele])
        else:
            suffix = suffix + '-' + "{:.0f}".format(parameters[ele])

    print_result_to_terminal(m.output)
    print_result_to_csv(output=m.output, suffix=suffix)

    return m.output

def get_TRHP_optimal(parameters):
    
    parameters.update({'comp_PR_uniform': True})
    bnds, null, xi = get_bounds_constraints(parameters)

    # Stage 1 results
    # TNC:         Obj function: 2.4460048742934912   time taken: 83.31311983399792   iter: 76  [10.33852704 11.10224061 12.06421683]
    # COBYLA:      Obj function: 2.4472800163803337   time taken: 67.44320199999493   iter: 61  [10.59877176 12.06284308 11.6946918 ]
    # L-BFGS-B:    Obj function: 2.450975749198941    time taken: 93.72743958298815   iter: 80  [ 8.89610068 13.70109655 11.76477823]
    # DA           Obj function: 2.483999500956222    time taken: 285.50347345799673  iter: 249 [39.02741184 20.5775166   5.68172743]
    # SHGO:        Obj function: 2.4928093772501727   time taken: 189.85609366599238  iter: 187 [60.          6.28234605  0.        ]
    # SLSQP (-3):  Obj function: 2.5078354680816624   time taken: 98.8724535000074    iter: 69  [ 5.77521535 41.68543104  3.52243245]
    # Powell(-2):  Obj function: 2.508021193791541    time taken: 134.78243087499868  iter: 103 [17.4705338  36.69965176  1.53318787]
    # Powell(-3):  Obj function: 2.508032855605729    time taken: 169.6446572910063   iter: 128 [17.46920884 36.70298654  1.53294837]
    # Powell:      Obj function: 2.508034717676362    time taken: 368.26358054199954  iter: 271 [17.46971595 36.70311604  1.53268726]
    # DF:          Obj function: 2.5084170892152353   time taken: 976.7336959999957   iter: 695 [ 6.63906037 41.6696337   3.14823521]
    # nelder (-3): Obj function: 2.509937861912262    time taken: 748.6417926249997   iter: 464 [12.21400893 39.73371679  2.14883141]
    # nelder-mead: Obj function: 2.5099380782021825   time taken: 786.1336665419949   iter: 489 [12.21459005 39.73342848  2.14857862]

    TIC = timer()
    result = minimize(method='SLSQP', fun=get_TRHP_obj_fun, args=parameters, bounds=bnds, x0=xi, tol=1e-3)
    TOC = timer()
    dt5 = TOC - TIC

    # Powell (-3): Obj function: 2.5080337810370046   time taken: 138.55814145899785  iter: 105 [17.4695599  36.70293412  1.53282015 33.33329619 33.33317214]
    # nelder (-3): Obj function: 2.5366470558544436   time taken: 1283.897543625004   iter: 917 [30.60922086 41.4015971   0.         45.         44.99983587]
    # SLSQP (-3):  Obj function: 2.5078355843904285   time taken: 21.571080500012613  iter: 13  [ 5.77522075 41.6854673   3.52240467 33.33329395 33.33331956]
    # SLSQP (-4):  Obj function: 2.5078355843904285   time taken: 22.44235362499603   iter: 13  [ 5.77522075 41.6854673   3.52240467 33.33329395 33.33331956]
    # SLSQP (-5):  Obj function: 2.5246633642281338   time taken: 1125.0843905830116  iter: 598 [ 8.74022909 53.9815596   1.73066445 37.32273063 45.        ]
    
    parameters.update({'comp_PR_uniform': False})
    bnds, null, xi = get_bounds_constraints(parameters, result.x)

    TIC = timer()
    result = minimize(method='SLSQP', fun=get_TRHP_obj_fun, args=parameters, bounds=bnds, x0=xi, tol=1e-3)
    TOC = timer()
    dt4 = TOC - TIC

    # 3.337977405464159
    # print('Case:', parameters)
    # print('1a: Obj function:',  -result3.fun / 1000, '  time taken:', dt3 , '  iter:', result3.nfev, result3.x)
    # print('2a: Obj function:',  -result1.fun / 1000,  ' time taken:', dt1,  '  iter:', result1.nfev, result1.x)
    # print('2b: Obj function:', -result0.fun / 1000, '  time taken:', dt0,  '  iter:', result0.nfev, result0.x)

    if result.fun < 0.0:
        return get_results_summary(result, parameters)

def f(var, parameters):
    m = create_m(var, parameters)
    dT = m.get_solve_state(P_comp_out=exp(var[len(var) - 1]))
    print('COP:', "{:.4f}".format(m.COP), '  T_sink_out:', "{:.4f}".format(m.gas_cooler4.sink_out.temperature), '  dT:', "{:.4f}".format(dT))
    return m.COP * -1000

def g(var, parameters):
    m = create_m(var, parameters)
    dT = m.get_solve_state(P_comp_out=exp(var[len(var) - 1]))
    return dT * 1000

def get_TRHP_optimal_v2(parameters):
    parameters.update({'comp_PR_uniform': True})
    bnds, cons, null = get_bounds_constraints(parameters)
    
    TIC = timer()
    result = shgo(func=f, args=[parameters], bounds=bnds, constraints=cons, workers=1, #sampling_method='sobol', iters=5,
                  minimizer_kwargs={'method': 'SLSQP', 'bounds': bnds, 'options': {'ftol': 1e-3}, 'constraints': cons})
    TOC = timer()
    dt = TOC - TIC

    print('Obj function:', -result.fun / 1000, '  time taken:', dt, '  iter:', result.nfev, result.success)
    print(-f(result.x, parameters)/1000)
    print(g(result.x, parameters)/100)
    
    # 3.337977405464159
