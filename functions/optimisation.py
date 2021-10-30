__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

#import os, numpy as np, scipy as sp, math, glob, pandas, timeit, matplotlib.pyplot as plt, xlwings, copy
from classes import HeatPump
from scipy.optimize import differential_evolution as opt_method, slsqp #dual_annealing  differential_evolution shgo
#from scipy.optimize import dual_annealing as opt_method
from scipy.optimize import minimize
from globalFunctions import REFPROP
from multiprocessing import Lock

import threading 
lock = threading.Lock()
def obj_TRHP(var, indep_var):

    #Independent variables 
    Refrigerant = indep_var[0]
    dT_GC_min_tar = indep_var[1]
    T_air_tar = indep_var[2]
    T_air_inlet = indep_var[3]
    T_air_exh = indep_var[4]
    #RH_air_exh = indep_var[5]
    n_comp = indep_var[5]
    #T_ihx_in_max = indep_var[6]
    #T_evap_init = indep_var[7]
    
    T_gc_out_min = T_air_inlet + dT_GC_min_tar 
    T_ihx_in_max = T_air_exh - dT_GC_min_tar
    
    #Optimisation variables
    dT_phx1_max = var[0]
    dT_phx2_max = 0.0#var[1] if n_comp > 1 else 0
    dT_phx3_max = 0.0#var[2] if n_comp > 2 else 0
    PR_2 = 0.0#var[3] 
    PR_3 = 0.0#var[4] 

    #Fixed variables
    eff_phx = 0.90  
    eff_motor = 0.95
    try:
        heat_pump = HeatPump(Refrigerant=Refrigerant, T_air_inlet=T_air_inlet, T_air_exh=T_air_exh, T_air_tar=T_air_tar, dT_GC_min_tar=dT_GC_min_tar, 
                            T_gc_out_min=T_gc_out_min, T_ihx_in_max=T_ihx_in_max, eff_phx=eff_phx, dT_phx1_max=dT_phx1_max, dT_phx2_max=dT_phx2_max, dT_phx3_max=dT_phx3_max, 
                            PR_2_init=PR_2, PR_3_init=PR_3, eff_motor=eff_motor, n_comp=n_comp) #RH_air_exh=RH_air_exh

        obj_COP = heat_pump.solve()

        
        

        print(var, indep_var, (obj_COP))
        #print("dT_phx1_max",dT_phx1_max, "dT_phx2_max",dT_phx2_max, 'T_ihx_in_max',T_ihx_in_max, 'PR_offset1', PR_offset)
        #print('T_air_tar=',T_air_tar, 'dT_GC_min_tar=',dT_GC_min_tar, 'T_air_inlet = ',T_air_inlet, 'T_air_exh =',T_air_exh, 'T_gc_out_min = ',T_gc_out_min)
    except:
        obj_COP = 1
        pass
    return obj_COP * -1

    
def optimiser_function(indep_var):
    lock.acquire()
    Lower_bound=1.0
    #Define the bounds of the optimisation variables
    pcrit = REFPROP('Pcrit', indep_var[0], 'PH', 10000 / 1000, 200)*1000
    p_RT = REFPROP('P', indep_var[0], 'QT', 1,15)*1000
    upperbound_PR2=(pcrit / p_RT)
    
    #constraint to PR bounds
    '''
    if (pcrit / p_RT) > 25.0:
        upperbound_PR2 = 25.0 
    else:
        upperbound_PR2=(pcrit / p_RT)
    ''' 
    upperbound_PR3 = upperbound_PR2 /2
    
    
    if upperbound_PR3 < 1.0:
        Lower_bound = 0.1
    else:
        Lower_bound = 1.0

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
        
    bound_dT_phx1_max = [0.0, 40.0]
    bound_dT_phx2_max = [0.0, 40.0]
    bound_dT_phx3_max = [0.0, 40.0] 
    bound_PR_2 = [Lower_bound, upperbound_PR2]
    bound_PR_3 = [Lower_bound, upperbound_PR3]
    bnds = [bound_dT_phx1_max]#, bound_dT_phx2_max,  bound_dT_phx3_max,  bound_PR_2, bound_PR_3] #bound_dT_phx1_max, bound_dT_phx2_max,  bound_dT_phx3_max
    
    
    start =  [PHX1_x0]#, PHX2_x0, PHX3_x0,1.0, 1.0]
    
    # perform the dual annealing search
   #result = opt_method(func=obj_TRHP, bounds=bnds, args=[indep_var], maxiter=1000) #, maxiter=100
    
    result = minimize(fun=obj_TRHP, x0=start, args=indep_var, method= 'SLSQP', bounds= bnds)
    result.fun *= -1
    print(result.fun) 
    lock.release()

    #for base case
    #var = 0
    #r#esult = obj_TRHP(var, indep_var)
    #result *= -1
    #print(result) 

    output = []
    output.append(['Obj function:', result.fun])
    output.append(['iterations:', result.nfev])
    #output.append(['Obj function:', result])
    output.append(['Refrigerant:', indep_var[0]])
    output.append(['dT Gas Cooler:', indep_var[1]])
    output.append(['T_air_tar:', indep_var[2]])
    output.append(['T_air_inlet', indep_var[3]])
    #output.append(['T_ihx_in_max', indep_var[6]])
    output.append(['dT phx max 1:', result.x[0]])
    #output.append(['dT phx max 2:', result.x[1]])
    #output.append(['dT phx max 3:', result.x[2]])
    #output.append(['PR 2:', result.x[3]])
    #output.append(['PR 3:', result.x[4]])
    output.append(['Message:', result.message])
    output.append(['Success:', result.success])

    return output
   

    