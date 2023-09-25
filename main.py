__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2023'

import concurrent.futures, os
from timeit import default_timer as timer
from heat_pumps import TRHP
from optimise import get_TRHP_optimal
from outputs import print_result_to_terminal, print_result_to_csv, read_csv_to_results
from outputs import Graph

if __name__=='__main__':
    TIC = timer()
    refrigerant = ['Isobutane', 'Butane', 'R1233zd(E)', 'R1234ze(E)', 'R1234ze(Z)', 'R1336mzz(Z)', 'Carbon dioxide', 'Propane']
    dT_GC_min_tar = [float(j) for j in range(20,35,5)]
    T_air_tar = [float(k) for k in range(100,255,10)]
    T_air_inlet = [float(m) for m in range(15,60,10)]
    T_source_in = [45.0]
    n_comp = [n + 1 for n in range(3)]

    if 1:
        file_name='TRHP_results_v1'
        prev_results = read_csv_to_results(file_name=file_name)
        n_cpu = os.cpu_count() - 2
        for j in range(len(dT_GC_min_tar)):
            for i in range(len(refrigerant)):
                for m in range(len(T_air_inlet)):
                    for n in reversed(range(len(n_comp))):
                        cases = [{
                                'refrigerant' : refrigerant[i],
                                'dT_GC_min_tar' : dT_GC_min_tar[j],
                                'T_air_tar' : T_air_tar[k],
                                'T_air_inlet' : T_air_inlet[m],
                                'T_source_in' : T_source_in[0],
                                'n_comp' : n_comp[n]
                                } for k in range(len(T_air_tar)) ]
                                                                        
                        if prev_results != None:
                            r = []
                            i = 0
                            for case in cases:
                                for res in prev_results:
                                    if (case['refrigerant'] == res['refrigerant:'] and 
                                        abs(case['dT_GC_min_tar'] - float(res['dT_min Gas Cooler:'])) < 0.01 and 
                                        abs(case['T_air_tar'] - float(res['T_air_tar:'])) < 0.01 and 
                                        abs(case['T_air_inlet'] - float(res['T_air_inlet:'])) < 0.01 and 
                                        abs(case['T_source_in'] - float(res['Dryer exh air T:'])) < 0.01 and 
                                        case['n_comp'] == int(res['n_comp:'])):
                                        r.append(i)
                                i = i + 1
                                    
                            for i in reversed(r):
                                cases.pop(i)

                        for c in cases:
                            print(c)
                        
                        if cases != None and cases != []:
                            results = [] 

                            executor = concurrent.futures.ProcessPoolExecutor(max_workers=n_cpu)
                            results = [executor.submit(get_TRHP_optimal, case) for case in cases]
                            # results = executor.map(get_TRHP_optimal, cases)
                            executor.shutdown()

                            for result in results:
                                if result != 0:
                                    print_result_to_csv(output=result, file_name=file_name)

    else:
        m = TRHP(
            refrigerant=refrigerant[0],
            dT_GC_min_tar=dT_GC_min_tar[0],
            T_air_tar=T_air_tar[5],
            T_air_inlet=T_air_inlet[0],
            T_source_in=T_source_in[0],
            n_comp=1,
            dT_phx1_cold=10.0,
            dT_phx2_cold=1.0,
            dT_phx3_cold=10.0, 
            comp2_power_frac=0.5, 
            comp3_power_frac=0.2,
            dTmin_phx=5.0
            )
        
        m.solve()
        
        if m.converged:
            m.get_summary()
            m.output.update({
                'Obj function:' : '',
                'Iterations:' : 0,
                'Message:' : '',
                'Success:' : m.converged
            })

            if type(m.output) is dict:
                print_result_to_terminal(output=m.output)

            if 0:
                g = Graph()
                g.heat_exchangers([m.gas_cooler1, m.gas_cooler2, m.gas_cooler3, m.gas_cooler4, m.phx1, m.phx2, m.phx3])
                g.heat_exchangers([m.gas_cooler1, m.phx1, m.gas_cooler2, m.phx2, m.gas_cooler3, m.phx3, m.gas_cooler4])  

    TOC = timer()
    print('Total time:',TOC - TIC) 

    