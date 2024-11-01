__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from functions.readWriteOutput import transfer_temp_files_to_results, read_csv_to_results

def get_cases_to_solve(file_name: str) -> list[list[float | str | bool]]:
    all_cases: list[list[float | str | bool]]
    unsolved_cases: list[list[float | str | bool]]
    
    all_cases = _get_all_case_list()
    unsolved_cases = _get_previous_results(all_cases, file_name)
    return _order_previous_results(unsolved_cases)

def _get_all_case_list() -> list[list[float | str | bool]]:
    refrigerant_ls: list[str] = [
        'propane;butane|.5;.5|mol',
        'Butane', 
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
        'PROPYLEN',
        'DME',
        'Methanol'
        'Ammonia',
        ]
    
    # Enumerate all the cases to be analysed
    return  [{
                'refrigerant' : [r],
                'dT_GC_min_tar' : [20.0],
                'T_air_tar' : [float(x) for x in range(250, 99, -5)],
                'T_air_inlet' : [15.0],
                'T_source_in' : [40.0],
                'n_comp' : [j + 1],
                'Tref_evap_inlet' : [10.0],
                'comp_PR_uniform' : [True],
                'result_parameters' : {'weight': [k / 100]},
                'results' : [],
             }
                for k in [k0 for k0 in range(100, 40, -10)]
                for j in range(2)
                for r in refrigerant_ls
            ]

def _get_previous_results(cases, file_name): 
    # Check if any past results in the temp folder need to be transferred and added to the complete results summary
    transfer_temp_files_to_results(file_name=file_name)
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
                    abs(float(case['Tref_evap_inlet']) - float(res['Evap T:'])) < 0.01 and
                    abs(float(case['result_parameters']['weight']) - float(res['Weight:'])) < 0.01 and 
                    int(case['n_comp']) == int(res['n_comp:'])
                    ):
                    case['results'].append(res)
                    break
    return cases

def _order_previous_results(cases):
    for case in cases:
        if len(case['results']) > 1:   
            case['results'].sort(reverse=True, key=lambda x: x['T_air_tar:'])         
    return cases
