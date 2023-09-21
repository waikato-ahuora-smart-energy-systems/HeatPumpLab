__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from scipy.optimize import minimize, LinearConstraint, differential_evolution as global_minimize
from heat_pumps import TRHP
import gc
from outputs import print_result_to_terminal

def get_TRHP_obj_fun(var, parameters):
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
        m.solve()
        COP = m.COP
        print(COP, var)
    
    return COP * -1000

def get_TRHP_optimal(parameters):
    try:
        #Define the bounds of the optimisation variables
        if parameters['n_comp'] == 1:
            bnds = [(0.0, 60.0)]
            lc = LinearConstraint([[0]], [0], [1])
            xi = [0]

        elif parameters['n_comp'] == 2:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 1.0)]
            lc = LinearConstraint([[0, 0, 1]], [0], [1])
            xi = [0, 0, 1/parameters['n_comp']]    

        elif parameters['n_comp'] == 3:
            bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 60.0), (0.0, 1.0), (0.0, 1.0)]
            lc = LinearConstraint([[0, 0, 0, 1, 1]], [0], [1])
            xi = [0, 0, 0, 1/parameters['n_comp'], 1/parameters['n_comp']]

        else:
            raise ValueError("Number of compressors must be 1, 2 or 3.")
        
        # Perform the optimisation
        # pool = multiprocessing.Pool(processes=1)
        # result = global_minimize(func=get_TRHP_obj_fun, args=[parameters], bounds=bnds, x0=xi, constraints=lc, polish=False, updating='deferred', workers=pool.map, popsize=4, tol=0.001)
        result = global_minimize(func=get_TRHP_obj_fun, args=[parameters], bounds=bnds, x0=xi, constraints=lc, polish=False, updating='immediate', workers=1, popsize=4, tol=0.001)
        # pool.close()
        result = minimize(method='SLSQP', fun=get_TRHP_obj_fun, args=parameters, bounds=bnds, x0=result.x, constraints=lc, tol=0.001)
        

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
                    comp2_power_frac=result.x[3]
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

            print_result_to_terminal(m.output)

            res = m.output
            del result
            del m
            gc.collect()
            return res

    except:
        return None 
