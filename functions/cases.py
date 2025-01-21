__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from functions.readWriteOutput import transfer_temp_files_to_results, read_csv_to_results

def get_analysis_cases_to_solve(file_name: str ='results/untitled') -> list[list[float | str | bool]]:
    return  _remove_completed_cases(
            cases = _get_all_analysis_case_list(),
            file_name = file_name,
            f_comparison = get_analysis_case_comparison,
    )

def get_optimisation_cases_to_solve(file_name: str ='results/untitled') -> list[list[float | str | bool]]:
    return  _order_cases(
                _remove_completed_cases(
                    _get_all_optimisation_case_list(), 
                    file_name,
                    get_optimisation_case_comparison,
                    )
                )

def get_analysis_case_comparison(res, case):
    return (
        case['refrigerant'] == res['refrigerant'] and 
        int(case['n_comp']) == int(res['n_comp']) and
        abs(float(case['dT_GC_min_tar']) - float(res['dT_GC_min_tar'])) < 0.01 and 
        abs(float(case['T_air_inlet']) - float(res['T_air_inlet'])) < 0.01 and
        abs(float(case['T_air_tar']) - float(res['T_air_tar'])) < 0.01 and 
        abs(float(case['T_source_in']) - float(res['T_source_in'])) < 0.01 and 
        abs(float(case['Tref_evap_inlet']) - float(res['Tref_evap_inlet'])) < 0.01 and
        abs(float(case['dTmin_phx']) - float(res['dTmin_phx'])) < 0.01 and
        abs(float(case['dT_phx1_cold']) - float(res['dT_phx1_cold'])) < 0.01 and
        abs(float(case['dT_phx2_cold']) - float(res['dT_phx2_cold'])) < 0.01 and
        abs(float(case['dT_phx3_cold']) - float(res['dT_phx3_cold'])) < 0.01 and
        abs(float(case['comp2_power_frac']) - float(res['comp2_power_frac'])) < 0.01 and
        abs(float(case['comp3_power_frac']) - float(res['comp3_power_frac'])) < 0.01
        )

def get_optimisation_case_comparison(res, case):
    return (
        case['refrigerant'] == res['refrigerant'] and 
        int(case['n_comp']) == int(res['n_comp']) and
        str(case['comp_PR_uniform']) == str(res['comp_PR_uniform']) and
        abs(float(case['dT_GC_min_tar']) - float(res['dT_GC_min_tar'])) < 0.01 and 
        abs(float(case['T_air_inlet']) - float(res['T_air_inlet'])) < 0.01 and
        abs(float(case['T_air_tar']) - float(res['T_air_tar'])) < 0.01 and 
        abs(float(case['T_source_in']) - float(res['T_source_in'])) < 0.01 and 
        abs(float(case['Tref_evap_inlet']) - float(res['Tref_evap_inlet'])) < 0.01 and
        abs(float(case['additional_data']['weight']) - float(res['Weight'])) < 0.01
        )

def _get_all_analysis_case_list() -> list[list[float | str | bool]]:
    refrigerant_ls: list[str] = [
        'Propane', 
        'Butane', 
        'Isobutane',
        'Pentane',
        'Isopentane',
        'Neopentane',
        # 'Carbon dioxide', 
        'R1233zd(E)', 
        'R1336mzz(E)', 
        'R1336mzz(Z)', 
        'R1234ze(E)', 
        'R1234ze(Z)', 
        'R1234yf', 
        'R1243zf', 
        'R1224yd(Z)', 
        'PROPYLEN',
        'DME',
        'Methanol',
        'Ammonia',
        # 'propane;butane|.5;.5|mol',
        # 'Methanol;butane|.5;.5|mol',
        # 'R1234ze(E);R1336mzz(Z)|.5;.5|mol',
        # 'propane;Isobutane|.5;.5|mol',
        # 'Propane;DME|.5;.5|mol',
        ]
    
    # Enumerate all the cases to be analysed
    return  [{
                'refrigerant': r,
                'dT_GC_min_tar': 20.0, 
                'T_air_tar': 200.0, 
                'T_air_inlet': 15.0, 
                'T_source_in': 40.0,
                'n_comp': 2,
                'ihx_on': 1,
                'dTmin_phx': 5.0,
                'dT_phx1_cold': float(dt1),
                'dT_phx2_cold': float(dt2),
                'dT_phx3_cold': 0.0, 
                'comp2_power_frac': float(cpf2) / 100, 
                'comp3_power_frac': 0.0,
                'comp_PR_uniform': False,
                'sink': 'Air',
                'source': 'water',
                'Tref_evap_inlet' : 10.0,
             }
                for dt1 in range(0, 21, 2)
                for dt2 in range(0, 81, 1)
                for cpf2 in range(10, 91, 5)
                for r in refrigerant_ls
            ]

def _get_all_optimisation_case_list() -> list[list[float | str | bool]]:
    refrigerant_ls: list[str] = [
        'Butane', 
        # 'Isobutane', 
        # 'Carbon dioxide', 
        # 'Propane', 
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
                'T_air_tar' : float(x),
                'T_air_inlet' : 15.0,
                'T_source_in' : 40.0,
                'n_comp' : j + 1 + 1,
                'Tref_evap_inlet' : 10.0,
                'comp_PR_uniform' : False,
                'additional_data' : {'weight': 0.98}, #k / 100
             }
                # for k in [k0 for k0 in range(100, 40, -10)]
                for j in range(1)
                for r in refrigerant_ls
                for x in range(250, 119, -10)
            ]

def _remove_completed_cases(cases: list, file_name, f_comparison): 
    # Check if any past results in the temp folder need to be transferred and added to the complete results summary
    transfer_temp_files_to_results(file_name=file_name)
    prev_results = read_csv_to_results(file_name=file_name)
    if prev_results != None:
        for res in prev_results:
            for i in range(len(cases)):
                case = cases[i]
                if f_comparison(res, case):
                    cases.pop(i)
                    break
    return cases

def _order_cases(cases):
    cases.sort(reverse=True, key=lambda x: x['additional_data']['weight'])  
    cases.sort(reverse=False, key=lambda x: x['refrigerant']) 
    cases.sort(reverse=True, key=lambda x: x['n_comp'])
    cases.sort(reverse=True, key=lambda x: x['T_air_tar']) 
    return cases
