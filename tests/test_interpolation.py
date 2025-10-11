from pathlib import Path
import sys

import pandas as pd
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from RecuperativeHeatPump.utils.interpolation import (  # noqa: E402
    estimate_operating_targets_from_results,
)


def _base_case_row() -> dict:
    """Return default column values shared by the synthetic datasets."""
    return {
        "refrigerant": "R1234ze(E)",
        "n_comp": 2,
        "comp_PR_uniform": False,
        "ihx_on": 0,
        "comp_model_index": 1,
        "sink": "Air",
        "source": "Water",
        "T_air_inlet": 15.0,
        "T_air_tar": 30.0,
        "Tref_evap_inlet": 10.0,
        "dTmin_rhx": 5.0,
        "motor_eff": 0.9,
        "dT_rhx1_cold": 0.0,
        "dT_rhx2_cold": 0.0,
        "dT_rhx3_cold": 0.0,
        "comp1_power_frac": 1.0,
        "comp2_power_frac": 0.0,
        "comp3_power_frac": 0.0,
    }


def test_estimate_operating_targets_linear_branch():
    base = _base_case_row()
    rows = [
        dict(base, T_source_in=0.0, dT_GC_min_tar=0.0, Pmax=0.0, dT_gc_out_min_offset=10.0),
        dict(base, T_source_in=1.0, dT_GC_min_tar=0.0, Pmax=1.0, dT_gc_out_min_offset=13.0),
        dict(base, T_source_in=0.0, dT_GC_min_tar=1.0, Pmax=2.0, dT_gc_out_min_offset=11.0),
    ]
    df = pd.DataFrame(rows)

    inputs = dict(base, T_source_in=0.5, dT_GC_min_tar=0.25)

    estimates = estimate_operating_targets_from_results(df, inputs)

    assert estimates["Pmax"] == pytest.approx(1.0)
    assert estimates["dT_gc_out_min_offset"] == pytest.approx(11.75)


def test_estimate_operating_targets_nearest_branch():
    base = _base_case_row()
    rows = [
        dict(base, T_source_in=0.0, dT_GC_min_tar=0.0, Pmax=5.0, dT_gc_out_min_offset=7.0),
        dict(base, T_source_in=2.0, dT_GC_min_tar=2.0, Pmax=9.0, dT_gc_out_min_offset=11.0),
    ]
    df = pd.DataFrame(rows)

    inputs = dict(base, T_source_in=1.6, dT_GC_min_tar=1.7)

    estimates = estimate_operating_targets_from_results(df, inputs)

    assert estimates["Pmax"] == pytest.approx(9.0)
    assert estimates["dT_gc_out_min_offset"] == pytest.approx(11.0)
