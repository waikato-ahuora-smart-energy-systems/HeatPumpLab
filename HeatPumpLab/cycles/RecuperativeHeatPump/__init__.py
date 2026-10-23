"""Expose the RecuperativeHP model."""

from .rhp_base import RecuperativeHP
from .utils import (
    analyse_RHP_case,
    analyse_unsolved_RHP_cases,
    get_best_cases_summary,
    get_comp2_power_frac_comparison_analysis,
    get_fluid_prop_overview_table,
    get_heatmap_results_comparison,
    get_RHP_case_list_to_solve,
    get_unsolved_RHP_case_list_to_resolve,
)

__all__ = [
    "RecuperativeHP",
    "analyse_RHP_case",
    "analyse_unsolved_RHP_cases",
    "get_best_cases_summary",
    "get_comp2_power_frac_comparison_analysis",
    "get_fluid_prop_overview_table",
    "get_heatmap_results_comparison",
    "get_RHP_case_list_to_solve",
    "get_unsolved_RHP_case_list_to_resolve",
]
