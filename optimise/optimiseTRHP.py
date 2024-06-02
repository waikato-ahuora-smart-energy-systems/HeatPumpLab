__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from scipy.optimize import minimize, differential_evolution as DE_opt, dual_annealing as DA_opt, shgo, direct as DI_opt, basinhopping as BH_opt
from heat_pumps import TRHP
from outputs import print_result_to_terminal, print_result_to_csv
from timeit import default_timer as timer
from streams.properties import get_property

def get_TRHP_obj_fun(var, parameters):
    TIC = timer()
    COP = 0.0
    
    if parameters['n_comp'] == 1:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0]
            )
        valid = True
        
    elif parameters['n_comp'] == 2:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0], 
            dT_phx2_cold=var[1], 
            comp2_power_frac=var[2]
            )
        valid = True if var[2] < 1.0 and var[2] > 0.0 else False

    elif parameters['n_comp'] == 3:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0], 
            dT_phx2_cold=var[1], 
            dT_phx3_cold=var[2], 
            comp2_power_frac=var[3], 
            comp3_power_frac=var[4]
            )
        valid = True if var[3] + var[4] < 1.0 and var[3] > 0.0 and var[4] > 0.0 else False
    
    else:
        raise ValueError("Number of compressors must be 1, 2 or 3.")

    if valid:
        # print(var, parameters)
        m.solve()
        COP = m.COP
        TOC = timer()
        print('COP:', "{:.4f}".format(COP), '  dt:', "{:.1f}".format(TOC - TIC), '  ', var)
        # print('COP:', "{:.4f}".format(COP), '  dt:', "{:.1f}".format(TOC - TIC), '  ', parameters, var)
    else:
        COP = 0
    return COP * -1000

def get_TRHP_optimal_old(parameters):
    try:
        #Define the bounds of the optimisation variables
        if parameters['n_comp'] == 1:
            bnds = [(0.0, 60.0)]
            xi = [0.0]

        elif parameters['n_comp'] == 2:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.25, 0.75)]
            xi = [0.0, 0.0, 1/2]    

        elif parameters['n_comp'] == 3:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 60.0), (0.15, 0.45), (0.15, 0.45)]
            xi = [0.0, 0.0, 0.0, 1/3, 1/3]

        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.")
        
        # Perform the optimisation
        # TIC = timer()
        # result4 = DI_opt(func=get_TRHP_obj_fun, args=[parameters], bounds=bnds)
        # TOC = timer()
        # dt4 = TOC - TIC

        # TIC = timer()
        # result3 = DA_opt(func=get_TRHP_obj_fun, args=[parameters], bounds=bnds, x0=xi)
        # TOC = timer()
        # dt3 = TOC - TIC

        TIC = timer()
        result2 = shgo(func=get_TRHP_obj_fun, args=[parameters], bounds=bnds, workers=-1, iters=3 )
        TOC = timer()
        dt2 = TOC - TIC

        TIC = timer()
        result = DE_opt(func=get_TRHP_obj_fun, args=[parameters], bounds=bnds, x0=xi, polish=False, workers=12, popsize=2, tol=1e-3)
        result = minimize(method='SLSQP', fun=get_TRHP_obj_fun, args=parameters, bounds=bnds, x0=result.x, tol=1e-4)
        TOC = timer()
        dt = TOC - TIC
        # 3.337977405464159
        print('DE: Obj function:', -result.fun / 1000,  '  time taken:', dt , '  iter:', result.nfev)
        print('SH: Obj function:', -result2.fun / 1000, '  time taken:', dt2, '  iter:', result2.nfev)
        # print('DA: Obj function:', -result3.fun / 1000, '  time taken:', dt3, '  iter:', result3.nfev)
        # print('DI: Obj function:', -result4.fun / 1000, '  time taken:', dt4, '  iter:', result4.nfev)

        if result.fun < 0.0:  
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

    except:
        return 0 

def f(var, parameters):
    COP = 0.0
 
    if parameters['n_comp'] == 1:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0]
            )
        valid = True
        
    elif parameters['n_comp'] == 2:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0], 
            dT_phx2_cold=var[1], 
            comp2_power_frac=var[2]
            )
        valid = True if var[2] < 1.0 and var[2] > 0.0 else False

    elif parameters['n_comp'] == 3:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0], 
            dT_phx2_cold=var[1], 
            dT_phx3_cold=var[2], 
            comp2_power_frac=var[3], 
            comp3_power_frac=var[4]
            )
        valid = True if var[3] + var[4] < 1.0 and var[3] > 0.0 and var[4] > 0.0 else False
    
    else:
        raise ValueError("Number of compressors must be 1, 2 or 3.")

    if valid:
        COP = m.f_cop(P_comp_out=var[len(var) - 1])
        print('COP:', "{:.4f}".format(COP), '  T_sink_out:', m.gas_cooler4.sink_out.temperature)

    return COP * -1000

def g(var, parameters):
    dT = 0.0

    if parameters['n_comp'] == 1:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0]
            )
        valid = True
        
    elif parameters['n_comp'] == 2:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0], 
            dT_phx2_cold=var[1], 
            comp2_power_frac=var[2]
            )
        valid = True if var[2] < 1.0 and var[2] > 0.0 else False

    elif parameters['n_comp'] == 3:
        m = TRHP(
            **parameters,
            dT_phx1_cold=var[0], 
            dT_phx2_cold=var[1], 
            dT_phx3_cold=var[2], 
            comp2_power_frac=var[3], 
            comp3_power_frac=var[4]
            )
        valid = True if var[3] + var[4] < 1.0 and var[3] > 0.0 and var[4] > 0.0 else False
    
    else:
        raise ValueError("Number of compressors must be 1, 2 or 3.")

    if valid:
        dT = m.g_dt(P_comp_out=var[len(var) - 1])

    return dT

def get_TRHP_optimal(parameters):
    try:
        #Define the bounds of the optimisation variables
        P_crit = get_property('PCRIT', parameters['refrigerant'], 'TQ', 100, 0)
        P_low = get_property('P', parameters['refrigerant'], 'TQ', parameters['Tsat_evap'], 1.0)
        PR_max = P_crit / P_low * 2.0

        if parameters['n_comp'] == 1:
            bnds = [(0.0, 60.0), (1.1, PR_max)]

        elif parameters['n_comp'] == 2:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.25, 0.75), (1.1, PR_max)]

        elif parameters['n_comp'] == 3:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 60.0), (0.15, 0.45), (0.15, 0.45), (1.1, PR_max)]

        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.")

        p = [parameters]
        cons = ({'type': 'ineq', 'fun': g, 'args': p})
        
        TIC = timer()
        result = shgo(func=f, args=parameters, bounds=bnds, constraints=cons, workers=1, iters=3 )
        TOC = timer()
        dt = TOC - TIC

        # 3.337977405464159
        print('Obj function:', -result.fun / 1000, '  time taken:', dt, '  iter:', result.nfev)

    except:
        return 0 


