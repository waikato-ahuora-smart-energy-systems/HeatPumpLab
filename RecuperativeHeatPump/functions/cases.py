import numpy as np
from .readWriteOutput import *

__all__ = ['get_case_list_to_solve', 'get_unsolved_case_list_to_resolve', 'get_optimisation_cases_to_solve', 'preselected_refrigerant_ls']

#######################################################################################################
# Public API
#######################################################################################################

def get_case_list_to_solve(
        file_name: str ='results/untitled', 
        refrigerant_ls: list = None, 
        dT_rhx1_cold_ls: list = None,
        dT_rhx2_cold_ls: list = None, 
        comp2_power_frac_ls: list = None,
        ) -> list[list[float | str | bool]]:
    return  (
        _remove_completed_cases_from_list(
            cases = _get_all_analysis_case_list(
                refrigerant_ls, 
                dT_rhx1_cold_ls, 
                dT_rhx2_cold_ls, 
                comp2_power_frac_ls
            ),
            file_name = file_name,
            f_comparison = _check_if_analysis_cases_are_identical,
        )
    )

def get_unsolved_case_list_to_resolve(filename_unsolved: str = 'Cases-Unsolved', filename_resolved: str = 'Cases-Resolved'):
    transfer_temp_files_to_results(f'results/{filename_resolved}')
    unsolved_cases = read_csv_to_results(file_name=f'results/{filename_unsolved}')
    resolved_cases = read_csv_to_results(file_name=f'results/{filename_resolved}')
    if resolved_cases != None:
        for i in reversed(range(len(unsolved_cases))):
            for sol in resolved_cases:
                if _check_if_analysis_cases_are_identical(sol,unsolved_cases[i]):
                    unsolved_cases.pop(i)
                    break
    return unsolved_cases

def get_optimisation_cases_to_solve(file_name: str ='results/untitled') -> list[list[float | str | bool]]:
    return  _order_cases(
                _remove_completed_cases_from_list(
                    _get_all_optimisation_case_list(), 
                    file_name,
                    _check_if_optimisation_cases_are_identical,
                    )
                )

def preselected_refrigerant_ls():
    return [
        'Ammonia',
        'Propane', 
        'Butane', 
        'Isobutane',
        'Pentane',
        'Isopentane',
        'Neopentane',
        'Cyclopentane',        
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
        # 'propane;butane|.5;.5|mol',
        # 'propane;Isobutane|.5;.5|mol',
        # 'propane;Pentane|.5;.5|mol',
        # 'propane;Isopentane|.5;.5|mol',
        # 'propane;Neopentane|.5;.5|mol',
        # 'Isobutane;butane|.5;.5|mol',
        # 'Isobutane;Pentane|.5;.5|mol',
        # 'Isobutane;Isopentane|.5;.5|mol',
        # 'Isobutane;Neopentane|.5;.5|mol',
        # 'butane;Pentane|.5;.5|mol',
        # 'butane;Isopentane|.5;.5|mol',
        # 'butane;Neopentane|.5;.5|mol',
    ]

#######################################################################################################
# Helper functions
#######################################################################################################

def _get_all_analysis_case_list(
        refrigerant_ls: list = ['Butane'], 
        T_source_in_ls: list = [40.0],
        Tref_evap_inlet_ls: list = [10.0],
        T_air_inlet_ls: list = [15.0],
        T_air_tar_ls: list = [200.0],
        dT_GC_min_tar_ls: list = [20.0],
        dTmin_rhx_ls: list = [5.0],
        dT_rhx1_cold_ls: list = [0.0],
        dT_rhx2_cold_ls: list = [0.0],
        dT_rhx3_cold_ls: list = [0.0], 
        comp2_power_frac_ls: list = [0.5],
        comp3_power_frac_ls: list = [0.0],
        n_comp_ls: list = [2],
        motor_eff_ls: list = [0.96],
        ihx_on_ls: list = [True],
        comp_PR_uniform_ls: list = [False],
        comp_model_index_ls: list = [1]
    ) -> list[list[float | str | bool]]:
    
    refrigerant_ls: list[str] = preselected_refrigerant_ls()
    dT_rhx1_cold_ls = np.linspace(0, 100, num=21)
    dT_rhx2_cold_ls = np.linspace(0, 100, num=21)
    comp2_power_frac_ls = np.linspace(0.3, 0.7, num=9)
    
    # Enumerate all the cases to be analysed
    return  [
        {
            'refrigerant': refrigerant,
            'T_source_in': T_source_in,
            'Tref_evap_inlet' : Tref_evap_inlet,
            'T_air_inlet': T_air_inlet, 
            'T_air_tar': T_air_tar, 
            'dT_GC_min_tar': dT_GC_min_tar, 
            'dTmin_rhx': dTmin_rhx,
            'dT_rhx1_cold': dT_rhx1_cold,
            'dT_rhx2_cold': dT_rhx2_cold,
            'dT_rhx3_cold': dT_rhx3_cold, 
            'comp2_power_frac': comp2_power_frac, 
            'comp3_power_frac': comp3_power_frac,
            'comp_PR_uniform': comp_PR_uniform,
            'comp_model_index': comp_model_index,
            'n_comp': n_comp,
            'ihx_on': ihx_on,
            'motor_eff': motor_eff,
            'sink': 'Air',
            'source': 'water',
        }
        for dT_rhx1_cold in dT_rhx1_cold_ls
        for dT_rhx2_cold in dT_rhx2_cold_ls
        for dT_rhx3_cold in dT_rhx3_cold_ls
        for comp2_power_frac in comp2_power_frac_ls
        for comp3_power_frac in comp3_power_frac_ls
        for refrigerant in refrigerant_ls
        for n_comp in n_comp_ls
        for motor_eff in motor_eff_ls
        for Tref_evap_inlet in Tref_evap_inlet_ls
        for T_air_tar in T_air_tar_ls
        for T_air_inlet in T_air_inlet_ls
        for dT_GC_min_tar in dT_GC_min_tar_ls
        for T_source_in in T_source_in_ls
        for dTmin_rhx in dTmin_rhx_ls
        for comp_PR_uniform in comp_PR_uniform_ls
        for ihx_on in ihx_on_ls
        for comp_model_index in comp_model_index_ls
    ]

def _get_all_optimisation_case_list() -> list[list[float | str | bool]]:
    refrigerant_ls: list[str] = [
        'Butane', 
        'Isobutane', 
        'Carbon dioxide', 
        'Propane', 
        'R1233zd(E)', 
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

def _check_if_analysis_cases_are_identical(res: dict, case: dict) -> bool:
    return (
        case['refrigerant'] == res['refrigerant'] and 
        int(case['n_comp']) == int(res['n_comp']) and
        abs(float(case['dT_GC_min_tar']) - float(res['dT_GC_min_tar'])) < 0.01 and 
        abs(float(case['T_air_inlet']) - float(res['T_air_inlet'])) < 0.01 and
        abs(float(case['T_air_tar']) - float(res['T_air_tar'])) < 0.01 and 
        abs(float(case['T_source_in']) - float(res['T_source_in'])) < 0.01 and 
        abs(float(case['Tref_evap_inlet']) - float(res['Tref_evap_inlet'])) < 0.01 and
        abs(float(case['dTmin_rhx']) - float(res['dTmin_rhx'])) < 0.01 and
        abs(float(case['dT_rhx1_cold']) - float(res['dT_rhx1_cold'])) < 0.01 and
        abs(float(case['dT_rhx2_cold']) - float(res['dT_rhx2_cold'])) < 0.01 and
        abs(float(case['dT_rhx3_cold']) - float(res['dT_rhx3_cold'])) < 0.01 and
        abs(float(case['comp2_power_frac']) - float(res['comp2_power_frac'])) < 0.01 and
        abs(float(case['comp3_power_frac']) - float(res['comp3_power_frac'])) < 0.01
        )

def _check_if_optimisation_cases_are_identical(res: dict, case: dict) -> bool:
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

def _remove_completed_cases_from_list(cases: list, file_name: str, f_comparison) -> list: 
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
