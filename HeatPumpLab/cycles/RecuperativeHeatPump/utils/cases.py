"""Helpers to enumerate, deduplicate, and prioritise recuperative heat pump cases."""

import pandas as pd

from ....utils.readWriteOutput import *

__all__ = ["get_RHP_case_list_to_solve", "get_unsolved_RHP_case_list_to_resolve"]

#######################################################################################################
# Public API
#######################################################################################################


def get_RHP_case_list_to_solve(
    file_name: str = "results/untitled",
    refrigerant_ls: list = None,
    dT_rhx1_cold_ls: list = None,
    dT_rhx2_cold_ls: list = None,
    comp2_power_frac_ls: list = None,
) -> list[list[float | str | bool]]:
    """Return pending analysis cases after removing already-solved entries."""
    return _remove_completed_cases_from_list(
        cases=_get_all_analysis_case_list(
            refrigerant_ls=refrigerant_ls,
            dT_rhx1_cold_ls=dT_rhx1_cold_ls,
            dT_rhx2_cold_ls=dT_rhx2_cold_ls,
            comp2_power_frac_ls=comp2_power_frac_ls,
        ),
        file_name=file_name,
        f_comparison=_check_if_analysis_cases_are_identical,
        f_signature=_analysis_case_signature,
    )


def get_unsolved_RHP_case_list_to_resolve(
    filename_unsolved: str = "Cases-Unsolved", filename_resolved: str = "Cases-Resolved"
):
    """Return unsolved cases after subtracting any newly resolved ones."""
    # transfer_temp_files_to_results(f'{filename_resolved}')
    unsolved_cases = read_csv_to_results(file_name=f"{filename_unsolved}")
    resolved_cases = read_csv_to_results(file_name=f"{filename_resolved}")
    if resolved_cases != None:
        for i in reversed(range(len(unsolved_cases))):
            for sol in resolved_cases:
                if _check_if_analysis_cases_are_identical(sol, unsolved_cases[i]):
                    unsolved_cases.pop(i)
                    break
    return unsolved_cases


#######################################################################################################
# Helper functions
#######################################################################################################


def _get_all_analysis_case_list(
    refrigerant_ls: list = ["Butane"],
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
    comp_model_index_ls: list = [1],
) -> list[list[float | str | bool]]:
    """Build the full factorial analysis case list before filtering."""
    # Enumerate all the cases to be analysed
    return [
        {
            "refrigerant": refrigerant,
            "T_source_in": T_source_in,
            "Tref_evap_inlet": Tref_evap_inlet,
            "T_air_inlet": T_air_inlet,
            "T_air_tar": T_air_tar,
            "dT_GC_min_tar": dT_GC_min_tar,
            "dTmin_rhx": dTmin_rhx,
            "dT_rhx1_cold": dT_rhx1_cold,
            "dT_rhx2_cold": dT_rhx2_cold,
            "dT_rhx3_cold": dT_rhx3_cold,
            "comp2_power_frac": comp2_power_frac,
            "comp3_power_frac": comp3_power_frac,
            "comp_PR_uniform": comp_PR_uniform,
            "comp_model_index": comp_model_index,
            "n_comp": n_comp,
            "ihx_on": ihx_on,
            "motor_eff": motor_eff,
            "sink": "Air",
            "source": "water",
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


def _check_if_analysis_cases_are_identical(res: dict, case: dict) -> bool:
    """Compare two analysis cases using a signature."""
    return _analysis_case_signature(res) == _analysis_case_signature(case)


def _remove_completed_cases_from_list(
    cases: list,
    file_name: str,
    f_comparison,
    f_signature=None,
) -> list:
    """Filter out cases already present in the results store."""
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
        # _apply_interpolation_if_possible(case, results_df)
        if f_signature is not None and case_signature is not None:
            signature_lookup.add(case_signature)
        remaining_cases.append(case)
    return remaining_cases


_ANALYSIS_FLOAT_FIELDS = [
    "dT_GC_min_tar",
    "T_air_inlet",
    "T_air_tar",
    "T_source_in",
    "Tref_evap_inlet",
    "dTmin_rhx",
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "dT_rhx3_cold",
    "comp2_power_frac",
    "comp3_power_frac",
]


def _analysis_case_signature(data: dict) -> tuple | None:
    """Signature tuple uniquely identifying an analysis case."""
    try:
        refrigerant = _normalize_str(data["refrigerant"])
        n_comp = int(float(data["n_comp"]))
    except (KeyError, TypeError, ValueError):
        return None

    float_parts = tuple(
        _round_float(data.get(field)) for field in _ANALYSIS_FLOAT_FIELDS
    )
    return (refrigerant, n_comp) + float_parts


#######################################################################################################
# Helper functions - common
#######################################################################################################


def _normalize_str(value: object) -> str:
    """Normalise strings for consistent signature generation."""
    return str(value).strip().lower()


def _round_float(value: object, ndigits: int = 4) -> float | None:
    """Safely round a value, returning None on failure."""
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return None
