__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley, Keri Adamson and Alex Geary, University of Waikato, 2021 (tim.walmsley@waikato.ac.nz)
"""
import timeit, xlwings
from classes import HeatPump
from functions import optimiser_function, graph_heat_exchangers
from globalFunctions import write_output
#import multiprocessing


@xlwings.func
def Trans_HTHP(Refrigerant='Carbon dioxide', T_air_tar=200.0, dT_GC_min_tar=20.0, n_comp=2, dT_phx1_max=0.0, dT_phx2_max=0.0, dT_phx3_max=0.0, 
               PR_offset1=0.0, PR_offset2=0.0, T_air_inlet=15.0, T_air_exh=70.0, T_gc_out_min=0.0, T_ihx_in_max=46.0, eff_motor=0.95, eff_phx=0.9):

    TIC = timeit.default_timer()
    heat_pump = HeatPump(Refrigerant=Refrigerant, T_air_tar=T_air_tar, dT_GC_min_tar=dT_GC_min_tar, n_comp=n_comp,
               dT_phx1_max=dT_phx1_max, dT_phx2_max=dT_phx2_max, dT_phx3_max=dT_phx3_max, PR_offset1=PR_offset1, PR_offset2=PR_offset2,
               T_air_inlet=T_air_inlet, T_air_exh=T_air_exh, T_gc_out_min=T_gc_out_min, T_ihx_in_max=T_ihx_in_max, eff_motor=eff_motor, eff_phx=eff_phx)

    heat_pump.solve()
    TOC = timeit.default_timer()
    
    return heat_pump.Write_Outputs(TOC - TIC) 
    
def multi_optimal_main():

    output = ['Obj function:','Iterations:','Refrigerant:','dT Gas Cooler:','T_air_tar:','T_air_inlet:','T_air_exh:','n_comp:','dT phx max 1:','dT phx max 2:','PR 2:','dT phx max 3:','PR 3:','Message:','Success:']
    write_output('TRHP_refs_v5.csv', output, 'w')

    Refrigerant = ['R1233zd(E)', 'R1234yf','R1234ze(Z)', 'R1243zf','R1336mzz(Z)', 'R1224yd(Z)', 'R1234ze(E)', 'Carbon dioxide', 'Ammonia', 'Isobutane', 'Butane']
    dT_GC_min_tar = [i for i in range(15,30,5)]
    T_air_tar = [i for i in range(180,220,10)]
    T_air_inlet = [i for i in range(10,30,5)]
    T_air_exh = 65.0
    n_comp = [i + 1 for i in range(3)]
    
    [optimiser_function([Refrigerant[j], dT_GC_min_tar[1], T_air_tar[2], T_air_inlet[1], T_air_exh, n_comp[1]]) for j in range(len(Refrigerant))]

    #indep_var =  [[Refrigerant[0], dT_GC_min_tar[1], T_air_tar[i], T_air_inlet[1], T_air_exh, n_comp[1]] for i in range(len(T_air_tar))]
    #with multiprocessing.Pool() as pool:
    #    results = pool.map(optimiser_function,(indep_var))
    

def single_state_main():
    heat_pump =  HeatPump(Refrigerant='R1336mzz(Z)', T_air_inlet=15, T_air_tar=200, T_air_exh=60, dT_phx1_max=0.0, dT_phx2_max=0.0, dT_phx3_max=0.0, n_comp=3, dT_GC_min_tar=20.0)        
    obj_COP = heat_pump.solve()

    if 0:
        print(heat_pump.Refrigerant)
        print(obj_COP)
        print(heat_pump.evaporator.ColdIn.Pressure)
        print(heat_pump.comp1.FlowOut.Pressure)
        print(heat_pump.comp2.FlowOut.Pressure)
        print(heat_pump.gas_cooler4.HotIn.Pressure)
        print(heat_pump.gas_cooler4.ColdOut.Temperature)
        print(heat_pump.gas_cooler1.ColdIn.Temperature)
        print(heat_pump.dT_GC_min_tar)
        print(heat_pump.evaporator.HotIn.Temperature)
        print(heat_pump.T_evap_init)      
        print(heat_pump.dT_phx1_max)
        print(heat_pump.dT_phx2_max)
        print(heat_pump.dT_phx3_max)
        print(heat_pump.PR_2_init)
        print(heat_pump.PR_3_init)
        print(heat_pump.T_air_tar)
        print(heat_pump.gas_cooler4.HotIn.Pressure / heat_pump.evaporator.ColdIn.Pressure)
    else:
        print('COP:', obj_COP)
        print('Refrigerant:', heat_pump.Refrigerant)
        print('Comp 3 pressure out:', heat_pump.gas_cooler4.HotIn.Pressure)
        print('Comp 2 pressure out:', heat_pump.comp2.FlowOut.Pressure)
        print('Comp 1 pressure out:', heat_pump.comp1.FlowOut.Pressure)
        print('Evap pressure:', heat_pump.evaporator.ColdIn.Pressure)
        print('PHX1:', heat_pump.dT_phx1_max)
        print('PHX2:', heat_pump.dT_phx2_max)
        print('PHX3:', heat_pump.dT_phx3_max)
        print('PR2 init:', heat_pump.PR_2_init)
        print('PR3 init:', heat_pump.PR_3_init)
        print('Gas cooler(s) inlet temp:', heat_pump.gas_cooler4.ColdOut.Temperature)
        print('Gas cooler(s) outlet temp:', heat_pump.gas_cooler1.ColdIn.Temperature)
        print('Evap temp:', heat_pump.T_evap_init)
        print('Target air:', heat_pump.T_air_tar)
        print('Dtmin:', heat_pump.dT_GC_min_tar)
        print('Overall pressure ratio:', heat_pump.gas_cooler4.HotIn.Pressure / heat_pump.evaporator.ColdIn.Pressure)
        print("Dryer exh air temp:", heat_pump.evaporator.HotIn.Temperature)
        print("Air temp out:", heat_pump.gas_cooler4.ColdOut.Temperature)
    
    if 0:
        graph_heat_exchangers([heat_pump.gas_cooler1, heat_pump.gas_cooler2, heat_pump.gas_cooler3, heat_pump.gas_cooler4, heat_pump.phx1, heat_pump.phx2, heat_pump.phx3])
    if 0:
        graph_heat_exchangers([heat_pump.gas_cooler1, heat_pump.phx1, heat_pump.gas_cooler2, heat_pump.phx2, heat_pump.gas_cooler3, heat_pump.phx3, heat_pump.gas_cooler4])


if __name__=='__main__':
    TIC = timeit.default_timer()

    multi_optimal_main()
    #single_state_main()

    TOC = timeit.default_timer()
    print("Total time:",TOC - TIC) 