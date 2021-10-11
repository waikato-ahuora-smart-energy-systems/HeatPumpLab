__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

#import os, numpy as np, scipy as sp, math, glob, pandas, timeit, matplotlib.pyplot as plt, xlwings, copy
from classes import HeatPump
from scipy.optimize import differential_evolution as opt_method #dual_annealing  differential_evolution shgo

Streams = []
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
    dT_phx2_max = var[1] if n_comp > 1 else 0
    dT_phx3_max = var[2] if n_comp > 2 else 0
    PR_2 = 0.0 #var[3]
    PR_3 = 0.0 #var[4]

    #Fixed variables
    eff_phx = 0.90  
    eff_motor = 0.95
    try:
        heat_pump = HeatPump(Refrigerant=Refrigerant, T_air_inlet=T_air_inlet, T_air_exh=T_air_exh, T_air_tar=T_air_tar, dT_GC_min_tar=dT_GC_min_tar, 
                            T_gc_out_min=T_gc_out_min, T_ihx_in_max=T_ihx_in_max, eff_phx=eff_phx, dT_phx1_max=dT_phx1_max, dT_phx2_max=dT_phx2_max, dT_phx3_max=dT_phx3_max, 
                            PR_2_init=PR_2, PR_3_init=PR_3, eff_motor=eff_motor, n_comp=n_comp) #RH_air_exh=RH_air_exh

        obj_COP = heat_pump.solve()

        global Streams 
        Streams= heat_pump.Stream_outputs()    

        print(var, indep_var, obj_COP)
        #print("dT_phx1_max",dT_phx1_max, "dT_phx2_max",dT_phx2_max, 'T_ihx_in_max',T_ihx_in_max, 'PR_offset1', PR_offset)
        #print('T_air_tar=',T_air_tar, 'dT_GC_min_tar=',dT_GC_min_tar, 'T_air_inlet = ',T_air_inlet, 'T_air_exh =',T_air_exh, 'T_gc_out_min = ',T_gc_out_min)
    except:
        obj_COP = 1
        pass
    return obj_COP * -1

    
def optimiser_function(lock, indep_var):
    #Define the bounds of the optimisation variables
    bound_dT_phx1_max = [0, 40.0]
    bound_dT_phx2_max = [0, 40.0]
    bound_dT_phx3_max = [0, 40.0] 
    bound_PR_2 = [1.0, 50.0]
    bound_PR_3 = [1.0, 50.0]
    bounds = [bound_dT_phx1_max, bound_dT_phx2_max] #, bound_dT_phx3_max]#, bound_PR_2]#, bound_PR_3] 

    # perform the dual annealing search
    result = opt_method(func=obj_TRHP, bounds=bounds, args=[indep_var], atol = 0.001, maxiter=100) #, maxiter=100
    result.fun *= -1
    print(result.fun) 

    #for base case
    #var = 0
    #result = obj_TRHP(var, indep_var)
    #result *= -1/1000
    #print(result) 

    output = []
    output.append(['Obj function:', result.fun])
    #output.append(['Obj function:', result])
    output.append(['Refrigerant:', indep_var[0]])
    #output.append(['dT Gas Cooler:', indep_var[1]])
    #output.append(['T_air_tar:', indep_var[2]])
    #output.append(['T_air_inlet', indep_var[3]])
    #output.append(['T_ihx_in_max', indep_var[6]])
    output.append(['dT phx max 1:', result.x[0]])
    output.append(['dT phx max 2:', result.x[1]])
    #output.append(['dT phx max 3:', result.x[2]])
    output.append(['P1:', Streams[0]])
    output.append(['P2', Streams[1]])
    output.append(['P3', Streams[2]])
    #output.append(['P4', Streams[3]])
    output.append(['Target air', Streams[4]])
    output.append(['Heat Capacity', Streams[5]])
    output.append(['Message:', result.message])
    output.append(['Success:', result.success])

    return output
   

    