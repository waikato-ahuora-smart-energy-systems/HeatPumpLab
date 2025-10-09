import numpy as np
import pandas as pd

from .interpolation import estimate_operating_targets_from_results
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
                refrigerant_ls=refrigerant_ls, 
                dT_rhx1_cold_ls=dT_rhx1_cold_ls, 
                dT_rhx2_cold_ls=dT_rhx2_cold_ls, 
                comp2_power_frac_ls=comp2_power_frac_ls,
            ),
            file_name = file_name,
            f_comparison = _check_if_analysis_cases_are_identical,
            f_signature = _analysis_case_signature,
        )
    )

def get_unsolved_case_list_to_resolve(filename_unsolved: str = 'Cases-Unsolved', filename_resolved: str = 'Cases-Resolved'):
    # transfer_temp_files_to_results(f'{filename_resolved}')
    unsolved_cases = read_csv_to_results(file_name=f'{filename_unsolved}')
    resolved_cases = read_csv_to_results(file_name=f'{filename_resolved}')
    if resolved_cases != None:
        for i in reversed(range(len(unsolved_cases))):
            for sol in resolved_cases:
                if _check_if_analysis_cases_are_identical(sol,unsolved_cases[i]):
                    unsolved_cases.pop(i)
                    break
    return unsolved_cases

def get_optimisation_cases_to_solve(file_name: str ='untitled') -> list[list[float | str | bool]]:
    return  _order_cases(
                _remove_completed_cases_from_list(
                    _get_all_optimisation_case_list(), 
                    file_name,
                    _check_if_optimisation_cases_are_identical,
                    _optimisation_case_signature,
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
        'Carbon dioxide', 
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
    dT_rhx1_cold_ls = np.linspace(0, 100, num=41)
    dT_rhx2_cold_ls = np.linspace(0, 100, num=41)
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
    return  [
        {
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
    return _analysis_case_signature(res) == _analysis_case_signature(case)

def _check_if_optimisation_cases_are_identical(res: dict, case: dict) -> bool:
    return _optimisation_case_signature(res) == _optimisation_case_signature(case)

def _remove_completed_cases_from_list(
    cases: list,
    file_name: str,
    f_comparison,
    f_signature=None,
) -> list: 
    # Check if any past results in the temp folder need to be transferred and added to the complete results summary
    transfer_temp_files_to_results(file_name=file_name)
    prev_results = read_csv_to_results(file_name=file_name)
    results_df = pd.DataFrame(prev_results) if prev_results not in (None, []) else None

    signature_lookup: set = set()
    if f_signature is not None and prev_results not in (None, []):
        for res in prev_results:
            try:
                signature = f_signature(res)
            except Exception:
                signature = None
            if signature is not None:
                signature_lookup.add(signature)

    remaining_cases: list[dict] = []
    for case in cases:
        case_signature = None
        if f_signature is not None:
            try:
                case_signature = f_signature(case)
            except Exception:
                case_signature = None
            if case_signature in signature_lookup:
                continue
        elif prev_results not in (None, []):
            if any(f_comparison(res, case) for res in prev_results):
                continue
        _apply_interpolation_if_possible(case, results_df)
        if f_signature is not None and case_signature is not None:
            signature_lookup.add(case_signature)
        remaining_cases.append(case)
    return remaining_cases

def _order_cases(cases):
    cases.sort(reverse=True, key=lambda x: x['additional_data']['weight'])  
    cases.sort(reverse=False, key=lambda x: x['refrigerant']) 
    cases.sort(reverse=True, key=lambda x: x['n_comp'])
    cases.sort(reverse=True, key=lambda x: x['T_air_tar']) 
    return cases


def _apply_interpolation_if_possible(case: dict, results_df: pd.DataFrame | None) -> None:
    """Populate initial estimates for a case using historical interpolation where available."""
    if results_df is None or results_df.empty:
        case.setdefault('dT_gc_out_min_offset', 0.0)
        return

    inputs = dict(case)
    comp1_frac = inputs.get('comp1_power_frac')
    if comp1_frac is None:
        comp2 = float(inputs.get('comp2_power_frac', 0.0))
        comp3 = float(inputs.get('comp3_power_frac', 0.0))
        comp1_frac = max(0.0, 1.0 - comp2 - comp3)
        inputs['comp1_power_frac'] = comp1_frac
    try:
        estimates = estimate_operating_targets_from_results(results_df, inputs)
    except (ValueError, TypeError):
        case.setdefault('dT_gc_out_min_offset', 0.0)
        return

    case['dT_gc_out_min_offset'] = float(estimates.get('dT_gc_out_min_offset', case.get('dT_gc_out_min_offset', 0.0)))
    if 'Pmax' in estimates and np.isfinite(estimates['Pmax']):
        case['P_H_est'] = float(estimates['Pmax']) / 1000.0
    if estimates.get('source', None):
        estimates.pop('source', None)
    case['interpolated_targets'] = estimates


def _normalize_str(value: object) -> str:
    return str(value).strip().lower()


def _round_float(value: object, ndigits: int = 4) -> float | None:
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return None


_ANALYSIS_FLOAT_FIELDS = [
    'dT_GC_min_tar',
    'T_air_inlet',
    'T_air_tar',
    'T_source_in',
    'Tref_evap_inlet',
    'dTmin_rhx',
    'dT_rhx1_cold',
    'dT_rhx2_cold',
    'dT_rhx3_cold',
    'comp2_power_frac',
    'comp3_power_frac',
]


def _analysis_case_signature(data: dict) -> tuple | None:
    try:
        refrigerant = _normalize_str(data['refrigerant'])
        n_comp = int(float(data['n_comp']))
    except (KeyError, TypeError, ValueError):
        return None

    float_parts = tuple(_round_float(data.get(field)) for field in _ANALYSIS_FLOAT_FIELDS)
    return (refrigerant, n_comp) + float_parts


_OPTIMISATION_FLOAT_FIELDS = [
    'dT_GC_min_tar',
    'T_air_inlet',
    'T_air_tar',
    'T_source_in',
    'Tref_evap_inlet',
]


def _optimisation_case_signature(data: dict) -> tuple | None:
    try:
        refrigerant = _normalize_str(data['refrigerant'])
        n_comp = int(float(data['n_comp']))
    except (KeyError, TypeError, ValueError):
        return None

    comp_pr_uniform = data.get('comp_PR_uniform', False)
    if isinstance(comp_pr_uniform, str):
        comp_pr_uniform_normalized = comp_pr_uniform.strip().lower()
        if comp_pr_uniform_normalized in {'false', '0'}:
            comp_pr_uniform_bool = False
        elif comp_pr_uniform_normalized in {'true', '1'}:
            comp_pr_uniform_bool = True
        else:
            comp_pr_uniform_bool = bool(comp_pr_uniform)
    else:
        comp_pr_uniform_bool = bool(comp_pr_uniform)

    float_parts = tuple(_round_float(data.get(field)) for field in _OPTIMISATION_FLOAT_FIELDS)

    weight = None
    if isinstance(data, dict):
        if 'additional_data' in data and isinstance(data['additional_data'], dict):
            weight = _round_float(data['additional_data'].get('weight'))
        elif 'Weight' in data:
            weight = _round_float(data.get('Weight'))

    return (refrigerant, n_comp, comp_pr_uniform_bool) + float_parts + (weight,)
