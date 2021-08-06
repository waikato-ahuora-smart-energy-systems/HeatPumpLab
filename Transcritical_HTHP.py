__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley, Keri Adamson and Alex Geary, University of Waikato, 2021 (tim.walmsley@waikato.ac.nz)
"""
import timeit, xlwings, multiprocessing
from classes import HeatPump
from functions import optimiser_function, graph_heat_exchangers
from globalFunctions import REFPROP
from multiprocessing import pool 



@xlwings.func
def Trans_HTHP(Refrigerant='Carbon dioxide', T_air_tar=200.0, dT_GC_min_tar=20.0, n_comp=2, dT_phx1_max=0.0, dT_phx2_max=0.0, dT_phx3_max=0.0, 
               PR_offset1=0.0, PR_offset2=0.0, T_air_inlet=15.0, T_air_exh=70.0, T_gc_out_min=0.0, T_ihx_in_max=46.0, eff_motor=0.95, eff_phx=0.9):

    TIC = timeit.default_timer()
    heat_pump_solution = HeatPump(Refrigerant=Refrigerant, T_air_tar=T_air_tar, dT_GC_min_tar=dT_GC_min_tar, n_comp=n_comp,
               dT_phx1_max=dT_phx1_max, dT_phx2_max=dT_phx2_max, dT_phx3_max=dT_phx3_max, PR_offset1=PR_offset1, PR_offset2=PR_offset2,
               T_air_inlet=T_air_inlet, T_air_exh=T_air_exh, T_gc_out_min=T_gc_out_min, T_ihx_in_max=T_ihx_in_max, eff_motor=eff_motor, eff_phx=eff_phx)

    heat_pump_solution.solve()
    TOC = timeit.default_timer()
    
    return heat_pump.Write_Outputs(TOC - TIC) 
    

#Run function for cmd line run
if __name__=='__main__':
    TIC = timeit.default_timer()
    if 1:
        T_air_inlet = 15.0
        T_air_exh = 70.0
        RH_air_exh = 0.15
        n_comp = 2
        T_air_tar = 200.0 #[i for i in range(180,205,5)]
        dT_GC_min_tar = 20 #[j for j in range(15,26,1)]
        
        Refrigerant = [ 'R1233zd(E)', 'R1234yf', 'R1234ze(Z)', 'R1243zf', 'R1336mzz(Z)', 'R1224yd(Z)', 'R1234ze(E)'] #'Carbon dioxide', ,'Ammonia', 'Isobutane', 'Butane'
        indep_var =  [[Refrigerant[j],dT_GC_min_tar, T_air_tar, T_air_inlet, T_air_exh, RH_air_exh, n_comp] for j in range(len(Refrigerant))]
        #            [[Refrigerant, float(dT_GC_min_tar[j]), T_air_tar, T_air_inlet, T_air_exh, RH_air_exh] for j in range(len(dT_GC_min_tar))]
        #            [dT_GC_min_tar, T_air_tar, T_air_inlet, T_air_exh, Refrigerant] 
        #            [[float(dT_GC_min_tar[j]), T_air_tar, T_air_inlet, T_air_exh, Refrigerant] for j in range(len(dT_GC_min_tar))]         

        with multiprocessing.Pool() as pool:
            results = pool.map(optimiser_function, indep_var)
            pool.start()
            pool.join()
        
        #results = optimiser_function(indep_var)
        
        print(results)
        f = open("TRHP_refs_v2.csv", "w")
        for result in results:
            for element in result:
                print(element, file=f)
            print('\n', file=f)
         
        f.close()        

    else:  
        heat_pump = HeatPump( Refrigerant = 'R1336mzz(Z)', T_air_tar=200.0,  dT_phx1_max=20, dT_phx2_max=35,  dT_phx3_max=0, n_comp = 3)     
        obj_COP = heat_pump.solve() 
        
        print(obj_COP)
        print(heat_pump.gas_cooler4.HotIn.Pressure)
        print(heat_pump.comp2.FlowOut.Pressure)
        print(heat_pump.comp1.FlowOut.Pressure)
        print(heat_pump.evaporator.ColdIn.Pressure)
        print(heat_pump.gas_cooler4.HotIn.Pressure / heat_pump.evaporator.ColdIn.Pressure)
        print("air temp out",heat_pump.gas_cooler4.ColdOut.Temperature)


        if 1:
            graph_heat_exchangers([heat_pump.gas_cooler1, heat_pump.gas_cooler2, heat_pump.gas_cooler3, heat_pump.gas_cooler4, heat_pump.phx1, heat_pump.phx2, heat_pump.phx3])
        else:
            graph_heat_exchangers([heat_pump.gas_cooler1, heat_pump.phx1, heat_pump.gas_cooler2, heat_pump.phx2, heat_pump.gas_cooler3, heat_pump.phx3, heat_pump.gas_cooler4])


    TOC = timeit.default_timer()
    print("time",TOC - TIC) 