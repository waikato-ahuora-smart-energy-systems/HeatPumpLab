__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley, Keri Adamson and Alex Geary, University of Waikato, 2021 (tim.walmsley@waikato.ac.nz)
"""

from timeit import default_timer as timer
from functions.optimise import optimiser_function, TRHP
from functions.write import terminal_print_result
from functions.graph import graph_heat_exchangers

if __name__=='__main__':
    TIC = timer()

    Refrigerant = ['R1233zd(E)', 'R1234yf','R1234ze(Z)', 'R1243zf','R1336mzz(Z)', 'R1224yd(Z)', 'R1234ze(E)', 'Carbon dioxide', 'Isobutane', 'Butane', 'Propane']
    dT_GC_min_tar = [j for j in range(20,40,5)]
    T_air_tar = [k for k in range(100,300,10)]
    T_air_inlet = [15.0]
    T_air_exh = [65.0]
    n_comp = [n + 1 for n in range(3)]

    if 0:
        [optimiser_function(
            parameters=[
                Refrigerant[0], 
                dT_GC_min_tar[2], 
                T_air_tar[12], 
                T_air_inlet[0],
                T_air_exh[0], 
                n_comp[1]
                ]
            )] # for j in range(len(Refrigerant))
        
    else:
        m = TRHP(
            Refrigerant = Refrigerant[0],
            dT_GC_min_tar = dT_GC_min_tar[2],
            T_air_tar = T_air_tar[12],
            T_air_inlet = T_air_inlet[0],
            T_air_exh = T_air_exh[0],
            n_comp = n_comp[2],
            dT_phx1_max = 60,
            dT_phx2_max = 60,
            dT_phx3_max = 20, 
            comp2_power_frac = 0.05, 
            comp3_power_frac = 0.8
        )
        m.solve()
        terminal_print_result(m.output)
        if 1:
            graph_heat_exchangers([m.gas_cooler1, m.gas_cooler2, m.gas_cooler3, m.gas_cooler4, m.phx1, m.phx2, m.phx3])
            #graph_heat_exchangers([m.gas_cooler1, m.phx1, m.gas_cooler2, m.phx2, m.gas_cooler3, m.phx3, m.gas_cooler4])  


    TOC = timer()
    print('Total time:',TOC - TIC) 

    