__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import numpy as np
from scipy.optimize import minimize, differential_evolution, Bounds
from functions.refprop import REFPROP
from functions.write import append_result_to_csv, append_result_to_output
from globalVariables import temp, filename_output as fout

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
        m.solve
        COP = - m.COP
        temp = m.output
    except:
        COP = 0
        temp = {'' : ''}

    del m 
    return COP

    
def optimiser_function(parameters):
    #Define the bounds of the optimisation variables
    pcrit = REFPROP('Pcrit', parameters[0], 'PH', 10, 200) * 1000
    p_RT = REFPROP('P', parameters[0], 'QT', 1, 5) * 1000
    Lower_bound = 1.0
    upperbound_PR2 = max(pcrit / p_RT / 2, 1.1)
    upperbound_PR3 = max(upperbound_PR2 / 2, 1.1)

    bnds = [(0.0, 60.0), (0.0, 60.0), (0.0, 60.0), (0.0, 1.0), (0.0, 1.0)]
    
    # Perform the optimisation
    result = differential_evolution(func=obj, bounds=bnds, args=[parameters], maxiter=200, updating='immediate', workers=1, polish=False)
    #DE_nfev = result.nfev
    #result = minimize(fun=obj, x0=result.x, args=parameters, method='SLSQP', bounds=bnds)
    
    #Covert the min obj result to a positive COP in W/W
    result.fun = - result.fun
    #result.nfev += DE_nfev

    output = append_result_to_output(output=temp, result=result)
    append_result_to_csv(output=output)