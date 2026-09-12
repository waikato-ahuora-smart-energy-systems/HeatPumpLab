import math
import json
from concurrent.futures import ThreadPoolExecutor

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.utils import composition_pareto
from HeatPumpLab.cycles.RecuperativeHeatPump.utils.composition_pareto import (
    add_pareto_metrics,
    nondominated_mask,
)


def test_pareto_maximises_both_metrics_and_handles_ties_and_invalid_points():
    assert nondominated_mask([
        (2, 4), (3, 3), (2, 2), (3, 3), (math.nan, 10), (0, 1), (1, math.inf),
    ]) == [True, True, False, True, False, False, False]
    assert nondominated_mask([(2, 4), (3, 4), (3, 5)]) == [False, False, True]


def test_pareto_uses_paired_replay_cop_and_vhc_and_reuses_endpoints():
    rows = [
        {"case_signature": "a", "input_refrigerant": "Butane", "COP": 5, "success": True},
        {"case_signature": "b", "input_refrigerant": "Propane", "COP": 2, "success": True},
        {"case_signature": "a", "input_refrigerant": "Butane", "COP": 5, "success": True},
        {"case_signature": "c", "input_refrigerant": "Pentane", "COP": None, "success": False},
    ]
    add_pareto_metrics(rows, {
        "a": {"COP": 2.5, "VHC_MJ_per_m3": 1.0, "COP_change_vs_optimisation_percent": -50},
        "b": {"COP": 3.0, "VHC_MJ_per_m3": 2.0, "COP_change_vs_optimisation_percent": 50},
    })
    assert [row["non_dominated"] for row in rows] == [False, True, False, False]
    assert rows[0]["replay_COP"] == 2.5
    assert rows[2]["VHC_MJ_per_m3"] == 1.0
    assert rows[3]["VHC_MJ_per_m3"] is None


def test_successful_design_requires_vhc_replay():
    with pytest.raises(ValueError, match="Missing COP/VHC replay"):
        add_pareto_metrics([{"case_signature": "a", "success": True, "input_refrigerant": "Butane"}], {})


def test_replay_cache_refreshes_old_vhc_definition_and_reuses_current_version(tmp_path, monkeypatch):
    cache_path = tmp_path / "design_replays.json"
    cache_path.write_text(json.dumps({"a": {"VHC_MJ_per_m3": 6.}}))
    calls = []
    fresh = {"VHC_MJ_per_m3": 4., "optimiser_version": composition_pareto.OPTIMISER_VERSION}

    def replay(row, results_dir):
        calls.append(row["case_signature"])
        return fresh

    monkeypatch.setattr(composition_pareto, "ProcessPoolExecutor", ThreadPoolExecutor)
    monkeypatch.setattr(composition_pareto, "_replay_design", replay)
    rows = [{"case_signature": "a", "batch_success": True, "input_refrigerant": "Butane"}]
    assert composition_pareto.replay_designs(rows, tmp_path) == {"a": fresh}
    assert composition_pareto.replay_designs(rows, tmp_path) == {"a": fresh}
    assert json.loads(cache_path.read_text()) == {"a": fresh}
    assert calls == ["a"]
