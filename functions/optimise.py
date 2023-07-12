__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import numpy as np
from scipy.optimize import minimize, differential_evolution, Bounds
from functions.refprop import REFPROP
from functions.output import write_output
from globalVariables import filename_output as fout
from classes import TRHP

def obj(var, indep_var):

    m = TRHP(
                Refrigerant = indep_var[0],
                dT_GC_min_tar = indep_var[1],
                T_air_tar = indep_var[2],
                T_air_inlet = indep_var[3],
                T_air_exh = indep_var[4],
                n_comp = indep_var[5], 
                dT_phx1_max = var[0], 
                dT_phx2_max = var[1], 
                PR_2 = var[2], 
                dT_phx3_max = var[3], 
                PR_3 = var[4]
            )
    
    try:
        m.solve
        COP = m.COP
    except:
        COP = 0

    del m 
    return COP * -1000

    
def optimiser_function(indep_var):
    #Define the bounds of the optimisation variables
    pcrit = REFPROP('Pcrit', indep_var[0], 'PH', 10, 200) * 1000
    p_RT = REFPROP('P', indep_var[0], 'QT', 1, 5) * 1000
    Lower_bound = 1.0
    upperbound_PR2 = max(pcrit / p_RT / 2, 1.1)
    upperbound_PR3 = max(upperbound_PR2 / 2, 1.1)

    bnds = [(0.0, 60.0), (0.0, 60.0), (Lower_bound, upperbound_PR2), (0.0, 60.0), (Lower_bound, upperbound_PR3)]
    
    # Perform the optimisation
    result = differential_evolution(func=obj, bounds=bnds, args=[indep_var], maxiter=200, updating='deferred', workers=-1, polish=False)
    DE_nfev = result.nfev
    #result = minimize(fun=obj, x0=result.x, args=indep_var, method='SLSQP', bounds=bnds)
    
    #Covert the min obj result to a positive COP in W/W
    result.fun = result.fun / -1000
    result.nfev += DE_nfev
    write_output(filename=fout, result=result, n_comp=indep_var[5])

    return result