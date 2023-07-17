__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import numpy as np, gc
from scipy.optimize import minimize, LinearConstraint, differential_evolution, dual_annealing, shgo
from classes.outputs.writeOutput import Write_Output

from classes.heat_pumps.TRHP import TRHP


class Optimise:
    def __init__(self) -> None:
        pass


    def get_obj_fun(self, var, parameters):

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

        m.solve()
        COP = m.COP
        #del m
        #collected = gc.collect()
        print(COP, var) 
        return -COP

        
    def get_optimal(self, parameters):
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
        #result = differential_evolution(func=self.get_obj_fun, bounds=bnds, x0=xi, args=[parameters], constraints=lc, updating='deferred', workers=-1, polish=True, popsize=10, atol=0.001) 
        result = dual_annealing(func=self.get_obj_fun, bounds=bnds, x0=xi, args=[parameters])
        #result = shgo(func=self.get_obj_fun, bounds=bnds, args=[parameters])

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

            o = Write_Output(m.output)
            o.terminal_print()
            print('Obj diff: ', abs(m.COP / result.fun + 1))
            print('Iterations: ', result.nfev)
            o.append_output_dict(result=result)
            o.append_csv()
    