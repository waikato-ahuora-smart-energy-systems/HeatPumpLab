__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import numpy as np, gc
from scipy.optimize import minimize, LinearConstraint, differential_evolution, dual_annealing, shgo
from functions.write import append_result_to_csv, append_result_to_output, terminal_print_result
from globalVariables import filename_output as fout

from classes import TRHP

def obj(var, parameters):

    m = TRHP(
        Refrigerant = parameters[0],
        dT_GC_min_tar = parameters[1],
        T_air_tar = parameters[2],
        T_air_inlet = parameters[3],
        T_air_exh = parameters[4],
        n_comp = parameters[5],
        dT_phx1_max = var[0], 
        dT_phx2_max = var[1], 
        dT_phx3_max = var[2], 
        comp2_power_frac = var[3], 
        comp3_power_frac = var[4]
        )
    
    try:
        m.solve()
        COP = m.COP
    except:
        COP = 0

    del m

    try:
        collected = gc.collect()
    except:
        pass

    print(COP, var) 
    return -COP

    
def optimiser_function(parameters):
    #Define the bounds of the optimisation variables
    if parameters[5] == 1:
        bnds = [(0.0, 60.0)]
        xi = [0]
    if parameters[5] == 2:
        bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 1.0)]
        lc = LinearConstraint([[0, 0, 1]], [0], [1])
        xi = [0, 0, 1/parameters[5]]    
    if parameters[5] == 3:
        bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 60.0), (0.0, 1.0), (0.0, 1.0)]
        lc = LinearConstraint([[0, 0, 0, 1, 1]], [0], [1])
        xi = [0, 0, 0, 1/parameters[5], 1/parameters[5]]
    
    # Perform the optimisation
    result = differential_evolution(func=obj, bounds=bnds, x0=xi, args=[parameters], updating='deferred', workers=-1, polish=True, constraints=lc, popsize=4, atol=0.01) 
    #result = dual_annealing(func=obj, bounds=bnds, x0=xi, args=[parameters])
    #result = direct(func=obj, bounds=bnds, args=[parameters])
    #result = shgo(func=obj, bounds=bnds, args=[parameters])

    if result.fun < 0.0:  
        if parameters[5] == 1:
            m = TRHP(
                Refrigerant = parameters[0],
                dT_GC_min_tar = parameters[1],
                T_air_tar = parameters[2],
                T_air_inlet = parameters[3],
                T_air_exh = parameters[4],
                n_comp = parameters[5],
                dT_phx1_max = result.x[0]
                )
        elif parameters[5] == 2:   
            m = TRHP(
                Refrigerant = parameters[0],
                dT_GC_min_tar = parameters[1],
                T_air_tar = parameters[2],
                T_air_inlet = parameters[3],
                T_air_exh = parameters[4],
                n_comp = parameters[5],
                dT_phx1_max = result.x[0], 
                dT_phx2_max = result.x[1],  
                comp2_power_frac = result.x[3]
                )        
        elif parameters[5] == 3:
            m = TRHP(
                Refrigerant = parameters[0],
                dT_GC_min_tar = parameters[1],
                T_air_tar = parameters[2],
                T_air_inlet = parameters[3],
                T_air_exh = parameters[4],
                n_comp = parameters[5],
                dT_phx1_max = result.x[0], 
                dT_phx2_max = result.x[1], 
                dT_phx3_max = result.x[2], 
                comp2_power_frac = result.x[3], 
                comp3_power_frac = result.x[4]
                )
        
        m.solve()
        terminal_print_result(m.output)
        print('Obj diff: ', abs(m.COP / result.fun + 1))
        print('Iterations: ', result.nfev)
        #output = append_result_to_output(output=m.output, result=result)
        #append_result_to_csv(output=output)
       