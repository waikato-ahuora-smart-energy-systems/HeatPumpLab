__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from heat_pumps import TRHP
from functions import print_result_to_terminal, print_result_to_csv
from optimise.optimiseTRHP import get_TRHP_optimal
from timeit import default_timer as timer
import concurrent.futures, os, csv, glob

def get_case_list(refrigerant,
                  dT_GC_min_tar,
                  T_air_tar,
                  T_air_inlet,
                  T_source_in,
                  n_comp,
                  Tsat_evap,
                  comp_PR_uniform,
                  w,
                  ):
    # Enumerate all the cases to be analysed
    cases = [{
                'refrigerant' : refrigerant[i],
                'dT_GC_min_tar' : dT_GC_min_tar[j],
                'T_air_tar' : [T_air_tar[k] for k in range(len(T_air_tar))],
                'T_air_inlet' : T_air_inlet[l],
                'T_source_in' : T_source_in[m],
                'n_comp' : n_comp[n],
                'Tsat_evap' : Tsat_evap[r],
                'comp_PR_uniform' : comp_PR_uniform[q],
                'result_parameters' : {'weight': w[o]},
                'results' : [],
             } 
                for r in range(len(Tsat_evap))
                for m in range(len(T_source_in))
                for j in range(len(dT_GC_min_tar))
                for l in range(len(T_air_inlet))
                for o in range(len(w))
                for q in range(len(comp_PR_uniform))
                for n in range(len(n_comp))
                for i in range(len(refrigerant))
            ]
    print('Cases remaining: ', len(cases))
    return cases

def transfer_results(file_name='results/untitled'):
        csv_file_list = glob.glob(os.getcwd() + '/results/temp/*.csv')
        for f in csv_file_list:
            output = read_csv_to_results(file_name=f, file_ext='')
            for o in output:
                print_result_to_csv(output=o, file_name=file_name)
            os.remove(f)

def read_csv_to_results(file_name='results/untitled', file_ext='.csv'):
    try:
        with open(file_name + file_ext, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            return list(reader)

    except:
        print('Failed to open CSV file.')

def get_previous_results(cases, file_name): 
    # Check if any past results in the temp folder need to be transferred and added to the complete results summary
    transfer_results(file_name=file_name)
    prev_results = read_csv_to_results(file_name=file_name)
    if prev_results != None:
        for res in prev_results:
            for case in cases:
                if (
                    case['refrigerant'] == res['Refrigerant:'] and 
                    bool(case['comp_PR_uniform']) == bool(res['Comp_PR_uniform:']) and
                    abs(float(case['dT_GC_min_tar']) - float(res['dT_min Gas Cooler:'])) < 0.01 and 
                    abs(float(case['T_air_inlet']) - float(res['T_air_inlet:'])) < 0.01 and 
                    abs(float(case['T_source_in']) - float(res['T_source_in:'])) < 0.01 and 
                    abs(float(case['Tsat_evap']) - float(res['Evap T:'])) < 0.01 and
                    abs(float(case['result_parameters']['weight']) - float(res['Weight:'])) < 0.01 and 
                    int(case['n_comp']) == int(res['n_comp:'])
                    ):
                    case['results'].append(res)
                    break
    return cases

def get_ordered_previous_results(cases):
    for case in cases:
        if len(case['results']) > 1:   
            case['results'].sort(reverse=True, key=lambda x: x['T_air_tar:'])         
    return cases

def get_parallel_solve(cases):
    # Use parallel processing distributed over all cores to speed up the analysis
    n_cpu = os.cpu_count()
    executor = concurrent.futures.ProcessPoolExecutor(max_workers=n_cpu)
    results = [executor.submit(get_TRHP_optimal, case) for case in cases]
    executor.shutdown()

def massive_parallel_case_analysis(parameters):
    cases = get_case_list(**parameters)
    # Location of saved results
    file_name='results/TRHP_results'
    cases = get_previous_results(cases, file_name)
    cases = get_ordered_previous_results(cases)
    get_parallel_solve(cases)
    transfer_results(file_name=file_name)

def single_case_analysis():
        args = {'refrigerant': 'Propane', #;butane|.5;.5|mass
                'dT_GC_min_tar': 20.0, 
                'T_air_tar': 200.0, 
                'T_air_inlet': 15.0, 
                'T_source_in': 40.0, 
                'n_comp': 1, 
                'dT_phx1_cold': 67.96875,
                'dT_phx2_cold': 0.0, 
                'dT_phx3_cold': 0.0, 
                'comp2_power_frac': 0.5, 
                'comp3_power_frac': 0.0,
                'comp_PR_uniform': True,
                'display': True,
                }
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
            
            m.get_cycle_plot(x_axis='h', y_axis='T')
            # m.get_cycle_plot(x_axis='h', y_axis='P')
            # m.get_cycle_plot(x_axis='v', y_axis='P')
            None
            