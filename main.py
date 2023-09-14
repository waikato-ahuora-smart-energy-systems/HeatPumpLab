__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2023'

import concurrent.futures
from timeit import default_timer as timer
from heat_pumps import TRHP
from optimise import get_TRHP_optimal
from outputs import print_result_to_terminal, print_result_to_csv, read_csv_to_results
from outputs import Graph

if __name__=='__main__':
    TIC = timer()
    refrigerant = ['Isobutane', 'Butane', 'Propane', 'Carbon dioxide', 'R1233zd(E)', 'R1234ze(E)', 'R1234ze(Z)', 'R1336mzz(Z)', 'R1224yd(Z)', 'R245fa', 'R365mfc']
    dT_GC_min_tar = [j for j in range(20,35,5)]
    T_air_tar = [k for k in range(100,255,5)]
    T_air_inlet = [15.0]
    T_source_in = [45.0]
    n_comp = [n + 1 for n in range(3)]

    if 1:
        prev_results = read_csv_to_results(file_name='TRHP_results_v4')

        for n in range(len(n_comp)):
            for i in range(len(refrigerant)):
                for j in range(len(dT_GC_min_tar)):
                    cases = [{
                            'refrigerant' : refrigerant[i],
                            'dT_GC_min_tar' : dT_GC_min_tar[j],
                            'T_air_tar' : T_air_tar[k],
                            'T_air_inlet' : T_air_inlet[0],
                            'T_source_in' : T_source_in[0],
                            'n_comp' : n_comp[n]
                        } for k in range(len(T_air_tar)) ]
                    
                    for c in cases:
                        print(c)
                    
                    completed = False
                    case = cases[0]
                    for res in prev_results:
                        if (case['refrigerant'] == res['refrigerant:'] and 
                            case['dT_GC_min_tar'] == int(res['dT_min Gas Cooler:']) and 
                            case['T_air_tar'] == int(res['T_air_tar:']) and 
                            case['T_air_inlet'] == float(res['T_air_inlet:']) and 
                            case['T_source_in'] == float(res['Dryer exh air T:']) and 
                            case['n_comp'] == int(res['n_comp:'])):
                            completed = True

                    if completed == False:
                        results = [] 
                        pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
                        futures = [pool.submit(get_TRHP_optimal, case) for case in cases]
                        for future in concurrent.futures.as_completed(futures):
                            result = future.result()
                            results.append(result)
                                
                        for result in results:
                            print_result_to_csv(output=result, file_name='TRHP_results_v4')
                        pool.shutdown()

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

    