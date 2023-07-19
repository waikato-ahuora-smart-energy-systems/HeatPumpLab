__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley, Keri Adamson and Alex Geary, University of Waikato, 2021 (tim.walmsley@waikato.ac.nz)
"""

import concurrent.futures
from timeit import default_timer as timer
from classes.logic_operations.optimise import get_optimal
from classes.heat_pumps.TRHP import TRHP
from classes.outputs.writeOutput import terminal_print, append_csv
from classes.outputs.graph import Graph


if __name__=='__main__':
    TIC = timer()

    Refrigerant = ['R1233zd(E)', 'R1234yf','R1234ze(Z)', 'R1243zf','R1336mzz(Z)', 'R1224yd(Z)', 'R1234ze(E)', 'Carbon dioxide', 'Isobutane', 'Butane', 'Propane']
    dT_GC_min_tar = [j for j in range(20,40,5)]
    T_air_tar = [k for k in range(100,300,10)]
    T_air_inlet = [15.0]
    T_air_exh = [65.0]
    n_comp = [n + 1 for n in range(3)]

    if 0:
        args = [{
                'Refrigerant' : Refrigerant[j],
                'dT_GC_min_tar' : dT_GC_min_tar[k],
                'T_air_tar' : T_air_tar[10],
                'T_air_inlet' : T_air_inlet[0],
                'T_air_exh' : T_air_exh[0],
                'n_comp' : n_comp[n]
            } for j in range(len(Refrigerant)) for k in range(len(dT_GC_min_tar)) for n in range(len(n_comp))] 
        
        i = 0 
        with concurrent.futures.ProcessPoolExecutor() as executor:
            for result in executor.map(get_optimal, args):
                if type(result) is dict:
                    append_csv(output=result)
                    i += result
        
        print(i)
        print(len(args))


    else:
        m = TRHP(
            Refrigerant = Refrigerant[0],
            dT_GC_min_tar = dT_GC_min_tar[2],
            T_air_tar = T_air_tar[12],
            T_air_inlet = T_air_inlet[0],
            T_air_exh = T_air_exh[0],
            n_comp = n_comp[0],
            dT_phx1_max = 4.17231090e+01,
            dT_phx2_max = 3.89716524e+01,
            dT_phx3_max = 7.09109499e+00, 
            comp2_power_frac = 1.08420217e-19, 
            comp3_power_frac = 9.45968903e-01
            )
        
        m.solve()

        m.output.update({
            'Obj function:' : 0,
            'Iterations:' : 0,
            'Message:' : 0,
            'Success:' : 0
        })

        if type(m.output) is dict:
            terminal_print(output=m.output)

        if 0:
            g = Graph()
            g.heat_exchangers([m.gas_cooler1, m.gas_cooler2, m.gas_cooler3, m.gas_cooler4, m.phx1, m.phx2, m.phx3])
            g.heat_exchangers([m.gas_cooler1, m.phx1, m.gas_cooler2, m.phx2, m.gas_cooler3, m.phx3, m.gas_cooler4])  

    TOC = timer()
    print('Total time:',TOC - TIC) 

    