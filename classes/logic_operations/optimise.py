__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import gc
from scipy.optimize import minimize, LinearConstraint, differential_evolution, dual_annealing
from classes.heat_pumps.TRHP import TRHP


def get_optimal(parameters):

    def get_obj_fun(var, parameters):
        COP = 0
        
        if var[3]+ var[4] < 1:
            if parameters['n_comp'] == 1:
                m = TRHP(
                    **parameters,
                    dT_phx1_max = var[0]
                    )

            elif parameters['n_comp'] == 2:
                m = TRHP(
                    **parameters,
                    dT_phx1_max = var[0], 
                    dT_phx2_max = var[1], 
                    comp2_power_frac = var[2]
                    )

            elif parameters['n_comp'] == 3:
                m = TRHP(
                    **parameters,
                    dT_phx1_max = var[0], 
                    dT_phx2_max = var[1], 
                    dT_phx3_max = var[2], 
                    comp2_power_frac = var[3], 
                    comp3_power_frac = var[4]
                    )
            
            else:
                raise ValueError("Number of compressors must be 1, 2 or 3.")

            try:
                m.solve()
                COP = m.COP
                print(COP, m.converged, var) 
                del m
                collected = gc.collect()
            except:
                pass

        
        return -COP


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
    result = differential_evolution(func=get_obj_fun, args=[parameters], bounds=bnds, x0=xi, constraints=lc, polish=False, workers=1, popsize=3, tol=0.001)
    result = minimize(method='SLSQP', fun=get_obj_fun, args=parameters, bounds=bnds, x0=result.x, constraints=lc, tol=0.001)

    if result.fun < 0.0:  
        if parameters['n_comp'] == 1:
            m = TRHP(
                **parameters,
                dT_phx1_max = result.x[0]
                )
        elif parameters['n_comp'] == 2:   
            m = TRHP(
                **parameters,
                dT_phx1_max = result.x[0], 
                dT_phx2_max = result.x[1],  
                comp2_power_frac = result.x[3]
                )        
        elif parameters['n_comp'] == 3:
            m = TRHP(
                **parameters,
                dT_phx1_max = result.x[0], 
                dT_phx2_max = result.x[1], 
                dT_phx3_max = result.x[2], 
                comp2_power_frac = result.x[3], 
                comp3_power_frac = result.x[4]
                )
        
        m.solve()

        m.output.update({
            'Obj function:' : result.fun,
            'Iterations:' : result.nfev,
            'Message:' : result.message,
            'Success:' : result.success
        })

        return m.output       

