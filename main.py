__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

import concurrent.futures, os
from timeit import default_timer as timer
from heat_pumps import TRHP
from optimise import get_TRHP_optimal
from outputs import print_result_to_terminal, read_csv_to_results, transfer_results, Graph

if __name__=='__main__':
    TIC = timer()

    if 1:
        refrigerant = ['Butane', 'Carbon dioxide', 'Propane', 'Isobutane', 'R1233zd(E)', 'R1336mzz(Z)', 'R1234ze(E)', 'R1234ze(Z)', 'R1234yf', 'R1243zf', 'R1224yd(Z)', 'Ammonia']

        dT_GC_min_tar = [float(j) for j in range(30,10,-10)]
        T_air_tar = [float(k) for k in range(250,110,-10)]
        T_air_inlet = [float(m) for m in range(35,10,-20)]
        T_source_in = [40.0]
        n_comp = [3 - n for n in range(3)]
        w = [float(o/100) for o in range(100, 50, -10)]
        Tsat_evap = [10.0]

        file_name='results/TRHP_results'
        transfer_results(file_name=file_name)

        prev_results = read_csv_to_results(file_name=file_name)
        n_cpu = os.cpu_count() 
        
        cases = [{
                    'refrigerant' : refrigerant[i],
                    'dT_GC_min_tar' : dT_GC_min_tar[j],
                    'T_air_tar' : T_air_tar[k],
                    'T_air_inlet' : T_air_inlet[l],
                    'T_source_in' : T_source_in[m],
                    'n_comp' : n_comp[n],
                    'Tsat_evap' : Tsat_evap[0],
                    'result_parameters' : {'weight': w(o)}
                    } 
                    for i in range(len(refrigerant))
                    for j in range(len(dT_GC_min_tar))
                    for k in range(len(T_air_tar)) 
                    for l in range(len(T_air_inlet))
                    for m in range(len(T_source_in))
                    for n in range(len(n_comp))
                    for o in range(len(w))
                ]

        if prev_results != None:
            c = []
            c_index = 0
            for case in cases:
                for res in prev_results:
                    if (
                        case['refrigerant'] == res['refrigerant:'] and 
                        abs(float(case['dT_GC_min_tar']) - float(res['dT_min Gas Cooler:'])) < 0.01 and 
                        abs(float(case['T_air_tar']) - float(res['T_air_tar:'])) < 0.01 and 
                        abs(float(case['T_air_inlet']) - float(res['T_air_inlet:'])) < 0.01 and 
                        abs(float(case['T_source_in']) - float(res['T_source_in:'])) < 0.01 and 
                        int(case['n_comp']) == int(res['n_comp:'])
                        ):
                        c.append(c_index)
                c_index = c_index + 1
                    
            for c_index in reversed(c):
                cases.pop(c_index)

        print('Cases remaining: ', len(cases))
        
        if cases != None and cases != []: 
            get_TRHP_optimal(parameters=cases[0])
            # executor = concurrent.futures.ProcessPoolExecutor(max_workers=n_cpu)
            # results = [executor.submit(get_TRHP_optimal, case) for case in cases]
            # executor.shutdown()

            # transfer_results(file_name=file_name)

    else:
        m = TRHP(
            refrigerant='Butane',
            dT_GC_min_tar=20.0,
            T_air_tar=250.0,
            T_air_inlet=15.0,
            T_source_in=40.0,
            n_comp=3,
            dT_phx1_cold=9.84375000,
            dT_phx2_cold=30.46875000,
            dT_phx3_cold=1.40625000, 
            comp_PR_uniform=True 
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

            # if 1:
            #     g = Graph()
            #     g.heat_exchangers([m.gas_cooler1, m.gas_cooler2, m.gas_cooler3, m.gas_cooler4, m.phx1, m.phx2, m.phx3])
            #     g.heat_exchangers([m.gas_cooler1, m.phx1, m.gas_cooler2, m.phx2, m.gas_cooler3, m.phx3, m.gas_cooler4])  

    TOC = timer()
    print('Total time:', TOC - TIC)

    