__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

import concurrent.futures, os
from timeit import default_timer as timer
from heat_pumps import TRHP
from optimise import get_TRHP_optimal
from outputs import print_result_to_terminal, read_csv_to_results, transfer_results

def optimisation_analysis():
    # Refrigerants analysed
    refrigerant = ['Butane', 
                   'Carbon dioxide', 
                   'Propane', 
                   'Isobutane', 
                   'R1233zd(E)', 
                   'R1336mzz(Z)', 
                   'R1234ze(E)', 
                   'R1234ze(Z)', 
                   'R1234yf', 
                   'R1243zf', 
                   'R1224yd(Z)', 
                   'Ammonia'
                   ]

    # Cycle parameters and their ranges tested
    dT_GC_min_tar = [20.0]
    T_air_tar = [float(k) for k in range(200, 90, -10)] 
    T_air_inlet = [15.0] 
    T_source_in = [40.0]
    n_comp = [2, 1]
    w = [float(o/100) for o in range(100, 40, -10)]
    Tsat_evap = [10.0]
    # comp_PR_uniform = [False]

    # Location of saved results
    file_name='results/TRHP_results'

    # Check if any past results in the temp folder need to be transferred and added to the complete results summary
    transfer_results(file_name=file_name)

    # Enumerate all the cases to be analysed
    cases = [{
                'refrigerant' : refrigerant[i],
                'dT_GC_min_tar' : dT_GC_min_tar[j],
                'T_air_tar' : T_air_tar[k],
                'T_air_inlet' : T_air_inlet[l],
                'T_source_in' : T_source_in[m],
                'n_comp' : n_comp[n],
                'Tsat_evap' : Tsat_evap[0],
                'result_parameters' : {'weight': w[o]},
                } 
                for m in range(len(T_source_in))
                for j in range(len(dT_GC_min_tar))
                for l in range(len(T_air_inlet))
                for i in range(len(refrigerant))
                for n in range(len(n_comp))
                for k in range(len(T_air_tar)) 
                for o in range(len(w))
            ]

    # Remove any cases that have already been analysed
    prev_results = read_csv_to_results(file_name=file_name)
    if prev_results != None:
        for res in prev_results:
            for i in range(len(cases), 0, -1):
                if (
                    cases[i-1]['refrigerant'] == res['Refrigerant:'] and 
                    abs(float(cases[i-1]['dT_GC_min_tar']) - float(res['dT_min Gas Cooler:'])) < 0.01 and 
                    abs(float(cases[i-1]['T_air_tar']) - float(res['T_air_tar:'])) < 0.01 and 
                    abs(float(cases[i-1]['T_air_inlet']) - float(res['T_air_inlet:'])) < 0.01 and 
                    abs(float(cases[i-1]['T_source_in']) - float(res['T_source_in:'])) < 0.01 and 
                    abs(float(cases[i-1]['result_parameters']['weight']) - float(res['Weight:'])) < 0.01 and 
                    int(cases[i-1]['n_comp']) == int(res['n_comp:'])
                    ):
                    cases.pop(i-1)
                    break

    # Disply the number of remaining cases that have no optimal result recorded
    print('Cases remaining: ', len(cases))

    # Use parallel processing distributed over all cores to speed up the analysis
    n_cpu = os.cpu_count() 
    if cases != None and cases != []: 
        # get_TRHP_optimal(parameters=cases[0])
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=n_cpu)
        results = [executor.submit(get_TRHP_optimal, case) for case in cases]
        executor.shutdown()
        transfer_results(file_name=file_name)

def single_case_analysis():
        args = {'refrigerant': 'Carbon dioxide', 
                'dT_GC_min_tar': 20.0, 
                'T_air_tar': 170.0, 
                'T_air_inlet': 15.0, 
                'T_source_in': 40.0, 
                'n_comp': 2, 
                'dT_phx1_cold': 0.0, 
                'dT_phx2_cold': 0.0, 
                'dT_phx3_cold': 0.0, 
                'comp2_power_frac': 0.5, 
                'comp3_power_frac': 0.0, 
                'comp_PR_uniform': True}
        m = TRHP(**args)
        TIC = timer()
        m.solve()
        TOC = timer()
        
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
                print('Solve time:', TOC - TIC)
            
            # m.get_cycle_plot(x_axis='h', y_axis='T')
            # m.get_cycle_plot(x_axis='h', y_axis='P')
            # m.get_cycle_plot(x_axis='v', y_axis='P')
            None

if __name__=='__main__':
    single_case_analysis()
    # optimisation_analysis()
    None
    