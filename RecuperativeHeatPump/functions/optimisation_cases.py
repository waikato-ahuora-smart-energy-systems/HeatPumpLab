from .cases import remove_completed_cases

def get_optimisation_cases_to_solve(file_name: str ='results/untitled') -> list[list[float | str | bool]]:
    return  _order_cases(
                remove_completed_cases(
                    _get_all_optimisation_case_list(), 
                    file_name,
                    _get_optimisation_case_comparison,
                    )
                )

def _get_optimisation_case_comparison(res, case):
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

def _order_cases(cases):
    cases.sort(reverse=True, key=lambda x: x['additional_data']['weight'])  
    cases.sort(reverse=False, key=lambda x: x['refrigerant']) 
    cases.sort(reverse=True, key=lambda x: x['n_comp'])
    cases.sort(reverse=True, key=lambda x: x['T_air_tar']) 
    return cases
