__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from functions.readWriteOutput import transfer_temp_files_to_results, read_csv_to_results

def get_cases_to_solve(file_name: str) -> list[list[float | str | bool]]:
    return  _remove_completed_cases(
                _order_all_cases(
                    _order_previous_results(
                        _get_previous_results(
                            _get_all_case_list(), file_name
                            )
                        )
                    )
                )

def _get_all_case_list() -> list[list[float | str | bool]]:
    refrigerant_ls: list[str] = [
        'Butane', 
        # 'Carbon dioxide', 
        # 'Propane', 
        # 'Isobutane', 
        # 'R1233zd(E)', 
        # 'R1336mzz(Z)', 
        # 'R1234ze(E)', 
        # 'R1234ze(Z)', 
        # 'R1234yf', 
        # 'R1243zf', 
        # 'R1224yd(Z)', 
        # 'PROPYLEN',
        # 'DME',
        # 'Methanol',
        # 'Ammonia',
        # 'propane;butane|.5;.5|mol',
        # 'Methanol;butane|.5;.5|mol',
        # 'R1234ze(E);R1336mzz(Z)|.5;.5|mol',
        # 'propane;Isobutane|.5;.5|mol',
        # 'Propane;DME|.5;.5|mol',
        ]
    
    # Enumerate all the cases to be analysed
    return  [{
                'refrigerant' : r,
                'dT_GC_min_tar' : 20.0,
                'T_air_tar' : [float(x) for x in range(250, 99, -5)],
                'T_air_inlet' : 15.0,
                'T_source_in' : 40.0,
                'n_comp' : j + 1,
                'Tref_evap_inlet' : 10.0,
                'comp_PR_uniform' : False,
                'additional_obj' : {'weight': k / 100},
                'results' : [],
             }
                for k in [k0 for k0 in range(100, 40, -10)]
                for j in range(3)
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
                    int(case['n_comp']) == int(res['n comp:']) and
                    bool(case['comp_PR_uniform']) == bool(res['Comp PR uniform:']) and
                    abs(float(case['dT_GC_min_tar']) - float(res['dT GC min tar:'])) < 0.01 and 
                    abs(float(case['T_air_inlet']) - float(res['T air inlet:'])) < 0.01 and 
                    abs(float(case['T_source_in']) - float(res['T source in:'])) < 0.01 and 
                    abs(float(case['Tref_evap_inlet']) - float(res['T evap inlet:'])) < 0.01 and
                    abs(float(case['additional_obj']['weight']) - float(res['Weight:'])) < 0.01
                    ):
                    case['results'].append(res)
                    break
    return cases

def _order_previous_results(cases):
    for case in cases:
        if len(case['results']) > 1:   
            case['results'].sort(reverse=True, key=lambda x: x['T air tar:'])       
    return cases

def _order_all_cases(cases):
    cases.sort(reverse=True, key=lambda x: x['additional_obj']['weight'])  
    cases.sort(reverse=True, key=lambda x: x['n_comp']) 
    cases.sort(reverse=True, key=lambda x: x['refrigerant']) 
    return cases

def _remove_completed_cases(cases_with_res):
    for i in reversed(range(len(cases_with_res))):
        if len(cases_with_res[i]['results']) == len(cases_with_res[i]['T_air_tar']):
            cases_with_res.pop(i)
    return cases_with_res
