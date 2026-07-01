"""Run the default recuperative heat-pump binary-mixture sweep and reports."""

import numpy as np
from HeatPumpLab.utils import get_parameter_sweep
from HeatPumpLab.cycles.RecuperativeHeatPump import (
    RecuperativeHP,
    get_fluid_prop_overview_table,
    get_heatmap_results_comparison,
    get_best_cases_summary,
    get_comp2_power_frac_comparison_analysis,
)

REFRIGERANT_LIST = [
    "Propane;Butane|0.5;0.5|mass",
]


def run():
    # Note: the sweep typically takes 30 hours to run on a standard computer,
    # so adjust the parameter lists as needed to reduce runtime
    param_kwargs = {
        "refrigerant_ls": REFRIGERANT_LIST,
        "dT_rhx1_cold_ls": np.linspace(0, 100, num=41),
        "dT_rhx2_cold_ls": np.linspace(0, 100, num=41),
        "comp2_power_frac_ls": np.linspace(0.3, 0.7, num=9),
    }
    get_parameter_sweep(
        HP_fn=RecuperativeHP,
        filename="RecuperatureHP_results",
        param_kwargs=param_kwargs,
    )


def results_analysis():
    get_fluid_prop_overview_table(
        refrigerant_ls=REFRIGERANT_LIST,
        dest_filename="results/Refrigerant props.csv",
    )
    get_heatmap_results_comparison(
        filepath="results/RecuperatureHP_results.csv",
        output_dir="results/heatmaps/",
    )
    get_best_cases_summary()
    get_comp2_power_frac_comparison_analysis()


if __name__ == "__main__":
    run()
    results_analysis()
