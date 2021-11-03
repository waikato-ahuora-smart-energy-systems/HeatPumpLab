__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes import HeatPump
from scipy import optimize
from globalFunctions import REFPROP

def obj_TRHP(var, indep_var):
    try:
        #Independent variables 
        Refrigerant = indep_var[0]
        dT_GC_min_tar = indep_var[1]
        T_air_tar = indep_var[2]
        T_air_inlet = indep_var[3]
        T_air_exh = indep_var[4]
        n_comp = indep_var[5]
        
        T_gc_out_min = T_air_inlet + dT_GC_min_tar 
        T_ihx_in_max = T_air_exh - dT_GC_min_tar
        
        #Optimisation variables
        dT_phx1_max = var[0]
        dT_phx2_max = 0.0
        PR_2 = 1.0 
        dT_phx3_max = 0.0
        PR_3 = 1.0

        if n_comp > 1:
            dT_phx2_max = var[1]
            PR_2 = var[2]
            if n_comp > 2:
                dT_phx3_max = var[3]
                PR_3 = var[4]

        #Fixed variables
        eff_phx = 0.90
        eff_motor = 0.95
    
        heat_pump = HeatPump(Refrigerant=Refrigerant, T_air_inlet=T_air_inlet, T_air_exh=T_air_exh, T_air_tar=T_air_tar, dT_GC_min_tar=dT_GC_min_tar, 
                            T_gc_out_min=T_gc_out_min, T_ihx_in_max=T_ihx_in_max, eff_phx=eff_phx, dT_phx1_max=dT_phx1_max, dT_phx2_max=dT_phx2_max, dT_phx3_max=dT_phx3_max, 
                            PR_2_init=PR_2, PR_3_init=PR_3, eff_motor=eff_motor, n_comp=n_comp)

        obj_COP = heat_pump.solve()
       
        print(var, indep_var, (obj_COP))

    except: 
        obj_COP = 1 #if the model fails, always return a COP of 1
        pass

    heat_pump = 0 #Wipe the memory of the heat pump data
    return 1000 / obj_COP #Convert to W/Q in W/kW for minimisation

    
def optimiser_function(indep_var):

    #Define the bounds of the optimisation variables
    pcrit = REFPROP('Pcrit', indep_var[0], 'PH', 10, 200) * 1000
    p_RT = REFPROP('P', indep_var[0], 'QT', 1, 5) * 1000
    Lower_bound = 1.0
    upperbound_PR2 = max(pcrit / p_RT / 2, 1.1)
    upperbound_PR3 = max(upperbound_PR2 / 2, 1.1)

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
    if 1:
        result = optimize.differential_evolution(func=obj_TRHP, bounds=bnds, args=[indep_var], maxiter=200, init='sobol', updating='deferred', workers=-1, polish=False) #allows for multiprocessing
        result = optimize.minimize(fun=obj_TRHP, x0=result.x, args=indep_var, method='SLSQP', bounds=bnds)
    elif 0:
        result = optimize.dual_annealing(func=obj_TRHP, bounds=bnds, args=[indep_var], maxiter=200) #does not produce the global minimum value
    else:
        result = optimize.minimize(fun=obj_TRHP, x0=x0, args=indep_var, method='SLSQP', bounds=bnds)

    #Covert the min obj result to the conventional COP = Q/W 
    result.fun = 1000 / result.fun
    
    output = []
    output.append(['Obj function:', result.fun])
    output.append(['iterations:', result.nfev])
    output.append(['Refrigerant:', indep_var[0]])
    output.append(['dT Gas Cooler:', indep_var[1]])
    output.append(['T_air_tar:', indep_var[2]])
    output.append(['T_air_inlet', indep_var[3]])
    output.append(['T_air_exh', indep_var[4]])
    output.append(['n_comp', indep_var[5]])
    output.append(['dT phx max 1:', result.x[0]])
    if n_comp > 1:
        output.append(['dT phx max 2:', result.x[1]])
        output.append(['PR 2:', result.x[2]])
    if n_comp > 2:
        output.append(['dT phx max 3:', result.x[3]])
        output.append(['PR 3:', result.x[4]])
    output.append(['Message:', result.message])
    output.append(['Success:', result.success])

    result = 0

    return output