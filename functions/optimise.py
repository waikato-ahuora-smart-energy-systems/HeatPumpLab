__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from scipy.optimize import minimize, differential_evolution
from functions.refprop import REFPROP
from functions.output import write_output
from globalVariables import filename_output as fout
from classes import TRHP

def obj(var, indep_var):
    if indep_var[5] == 1:
        var.extend([0.0, 1.0, 0.0, 1.0])
    
    if indep_var[5] == 2:
        var.extend([0.0, 1.0])

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

    '''
    #Constraint to IHX bounds
    if indep_var[0] == 'R1233zd(E)':
        PHX1_x0 = 39.9
        #PHX1_UB = 5.0
        #PHX2_x0 = 6.0
        #PHX2_UB =25.0
        #PHX3_x0 = 17.0
    if indep_var[0] == 'R1234yf':
        PHX1_x0 = 40.0
        #PHX1_UB = 5.0
        #PHX2_x0 = 26.0
        #PHX2_UB =25.0
        #PHX3_x0 = 15.0
    if indep_var[0] == 'R1234ze(Z)':
        PHX1_x0 = 40.0
        #PHX1_UB = 5.0
        #PHX2_x0 = 26.0
        #PHX2_UB =25.0
        #PHX3_x0 = 3.0
    if indep_var[0] == 'R1243zf':
        PHX1_x0 = 34.3
        #PHX1_UB = 5.0
        #PHX2_x0 = 22.0
        #PHX2_UB =25.0
        #PHX3_x0 = 2.0
    if indep_var[0] == 'R1336mzz(Z)':
        PHX1_x0 = 4.3
        #PHX1_UB = 5.0
        #PHX2_x0 = 6.0
        #PHX2_UB =25.0
        #PHX3_x0 = 39.0
    if indep_var[0] == 'R1224yd(Z)':
        PHX1_x0 = 40.0
        #PHX1_UB = 5.0
        #PHX2_x0 = 14.0
        #PHX2_UB =25.0
        #PHX3_x0 = 19.0
    if indep_var[0] == 'R1234ze(E)':
        PHX1_x0 = 30.0
        #PHX1_UB = 5.0
        #PHX2_x0 = 26.0
        #PHX2_UB =25.0
        #PHX3_x0 = 3.0
    if indep_var[0] == 'Carbon dioxide':
        PHX1_x0 = 4.1
        #PHX1_UB = 5.0
        #PHX2_x0 = 0.0
        #PHX2_UB =25.0
        #PHX3_x0 = 1.0
    if indep_var[0] == 'Ammonia':
        PHX1_x0 = 4.99
        #PHX1_UB = 5.0
        #PHX2_x0 = 4.0
        #PHX2_UB =25.0
        #PHX3_x0 = 1.0
    if indep_var[0] == 'Isobutane':
        PHX1_x0 = 40.0
        #PHX1_UB = 5.0
        #PHX2_x0 = 38.0
        #PHX2_UB =25.0
        #PHX3_x0 = 0.0
    if indep_var[0] == 'Butane':
        PHX1_x0 = 40.0
        #PHX1_UB = 5.0
        #PHX2_x0 = 11.0
        #PHX2_UB =25.0
        #PHX3_x0 = 23.0
    '''
    PHX1_x0 = 40.0
    
    n_comp = indep_var[5]
    x0 = [PHX1_x0]
    bnds = [[0.0, 60.0]]
    if n_comp > 1:
        bnds.append([0.0, 60.0]) #dT_phx2
        x0.append(0.0) #PHX2_x0
        bnds.append([Lower_bound, upperbound_PR2]) #PR_2
        x0.append(1.0)

    if n_comp > 2:
        bnds.append([0.0, 60.0]) #dT_phx3
        x0.append(0.0) #PHX2_x0
        bnds.append([Lower_bound, upperbound_PR3]) #PR_3 
        x0.append(1.0)
    
    # Perform the optimisation
    result = differential_evolution(func=obj, bounds=bnds, args=[indep_var], maxiter=200, updating='deferred', workers=-1, polish=False) #allows for multiprocessing
    DE_nfev = result.nfev
    #result = minimize(fun=obj, x0=result.x, args=indep_var, method='SLSQP', bounds=bnds)

    #Covert the min obj result to a positive COP in W/W
    result.fun = result.fun / -1000
    result.nfev += DE_nfev
    write_output(filename=fout, result=result, n_comp=n_comp)

    return result