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


# --- _replay_design -------------------------------------------------------

import matplotlib  # noqa: E402

matplotlib.use("Agg")

from types import SimpleNamespace  # noqa: E402

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _objective  # noqa: E402
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import DEFAULT_BASE_CASE  # noqa: E402


def _design_row(**overrides):
    row = {key: str(value) for key, value in DEFAULT_BASE_CASE.items()}
    row.update({
        "case_signature": "sig-1",
        "input_refrigerant": "Propane;Pentane|0.3;0.7|mass",
        "best_dT_rhx1_cold": "12.5",
        "best_dT_rhx2_cold": "7.25",
        "best_comp2_power_frac": "0.4",
        "best_mass_fractions": "0.3;0.7",
        "best_COP": "2.5",
    })
    row.update(overrides)
    return row


def _install_fake_objective(monkeypatch, record=None, error=None):
    calls = {"closed": 0}

    class FakeObjective:
        def __init__(self, **kwargs):
            calls["kwargs"] = kwargs

        def evaluate(self, x):
            calls["x"] = x
            if error is not None:
                raise error
            return record

        def close(self):
            calls["closed"] += 1

    monkeypatch.setattr(_objective, "_OuterDesignObjective", FakeObjective)
    return calls


def _record(**overrides):
    values = dict(
        success=True, cop=3.0, sec=1 / 3, pmax=2500.0, message="ok",
        report={"VHC": 4500.0, "T_comp_discharge": 151.0},
    )
    values.update(overrides)
    return SimpleNamespace(**values)


def test_replay_design_converts_units_and_reports_cop_change(monkeypatch, tmp_path):
    calls = _install_fake_objective(monkeypatch, _record())
    result = composition_pareto._replay_design(_design_row(), str(tmp_path))
    assert result == {
        "case_signature": "sig-1",
        "optimiser_version": composition_pareto.OPTIMISER_VERSION,
        "COP": 3.0,
        "VHC_MJ_per_m3": pytest.approx(4.5),  # kJ/m3 -> MJ/m3
        "SEC": pytest.approx(1 / 3),
        "Pmax_MPa": pytest.approx(2.5),  # kPa -> MPa
        "COP_change_vs_optimisation_percent": pytest.approx(20.0),
        "T_comp_discharge_C": 151.0,
    }
    assert calls["x"] == (12.5, 7.25, 0.4, 0.3, 0.7)
    assert calls["kwargs"]["refrigerant"] == "Propane;Pentane|0.3;0.7|mass"
    assert calls["kwargs"]["results_dir"] == str(tmp_path)
    assert calls["kwargs"]["max_evals"] == 1
    assert calls["closed"] == 1


def test_replay_design_parses_string_booleans_and_casts_base_case_types(monkeypatch, tmp_path):
    calls = _install_fake_objective(monkeypatch, _record())
    row = _design_row(comp_PR_uniform="True", ihx_on="False", debug="false",
                      n_comp="3", T_air_tar="180", sink="Air")
    composition_pareto._replay_design(row, str(tmp_path))
    base_case = calls["kwargs"]["base_case"]
    assert set(base_case) == set(DEFAULT_BASE_CASE)
    assert base_case["comp_PR_uniform"] is True
    assert base_case["ihx_on"] is False
    assert base_case["debug"] is False
    assert base_case["n_comp"] == 3 and type(base_case["n_comp"]) is int
    assert base_case["T_air_tar"] == 180.0 and type(base_case["T_air_tar"]) is float
    assert base_case["sink"] == "Air"


def test_replay_design_accepts_native_bool_values(monkeypatch, tmp_path):
    calls = _install_fake_objective(monkeypatch, _record())
    composition_pareto._replay_design(_design_row(ihx_on=False, comp_PR_uniform=True), str(tmp_path))
    assert calls["kwargs"]["base_case"]["ihx_on"] is False
    assert calls["kwargs"]["base_case"]["comp_PR_uniform"] is True


@pytest.mark.parametrize("record", [
    _record(success=False, message="solver failed"),
    _record(report={}),
    _record(report={"VHC": 0.0}),
    _record(report={"VHC": -10.0}),
    _record(report={"VHC": math.inf}),
], ids=["unsuccessful", "missing-vhc", "zero-vhc", "negative-vhc", "infinite-vhc"])
def test_replay_design_raises_on_failed_replay(monkeypatch, tmp_path, record):
    calls = _install_fake_objective(monkeypatch, record)
    with pytest.raises(RuntimeError, match="VHC replay failed for Propane;Pentane"):
        composition_pareto._replay_design(_design_row(), str(tmp_path))
    assert calls["closed"] == 1


def test_replay_design_closes_objective_when_evaluation_raises(monkeypatch, tmp_path):
    calls = _install_fake_objective(monkeypatch, error=OSError("worker died"))
    with pytest.raises(OSError, match="worker died"):
        composition_pareto._replay_design(_design_row(), str(tmp_path))
    assert calls["closed"] == 1


# --- plot_cop_vhc ---------------------------------------------------------

_FAMILIES = {
    "R290 / R600": "Butane",
    "R290 / R601A": "Isopentane",
    "R290 / R601": "Pentane",
    "R290 / R600A": "Isobutane",
}


def _pareto_rows():
    """Each family: pure propane, one 50% mixture, pure partner. Propane is shared."""
    propane = {"case_signature": "propane", "input_refrigerant": "Propane", "success": True,
               "R290_mass_percent": 100}
    mixture_points = {  # (COP, VHC) of the 50% blend
        "R290 / R600": (3.1, 3.0),  # dominated by R600A blend
        "R290 / R601A": (3.6, 2.0),  # front
        "R290 / R601": (3.0, 1.5),
        "R290 / R600A": (3.2, 4.0),  # front
    }
    pure_points = {"Butane": (3.1, 2.5), "Isopentane": (3.7, 1.0), "Pentane": (3.8, 0.8), "Isobutane": (2.9, 3.5)}
    replays = {"propane": {"COP": 2.5, "VHC_MJ_per_m3": 6.0, "COP_change_vs_optimisation_percent": 0}}
    rows = []
    for family, partner in _FAMILIES.items():
        rows.append(dict(propane, family=family))
        rows.append({"case_signature": f"{partner}-50", "input_refrigerant": f"Propane;{partner}",
                     "success": True, "R290_mass_percent": 50, "family": family})
        rows.append({"case_signature": partner, "input_refrigerant": partner, "success": True,
                     "R290_mass_percent": 0, "family": family})
        cop, vhc = mixture_points[family]
        replays[f"{partner}-50"] = {"COP": cop, "VHC_MJ_per_m3": vhc, "COP_change_vs_optimisation_percent": 0}
        cop, vhc = pure_points[partner]
        replays[partner] = {"COP": cop, "VHC_MJ_per_m3": vhc, "COP_change_vs_optimisation_percent": 0}
    return rows, replays


def _capture_figures(monkeypatch):
    from matplotlib.figure import Figure

    figures = []
    original = Figure.savefig

    def savefig(self, *args, **kwargs):
        figures.append(self)
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", savefig)
    return figures


def _front_line(figure):
    (ax,) = figure.axes
    (line,) = [line for line in ax.lines if line.get_label() == "Non-dominated sampled designs"]
    return list(zip(line.get_xdata(), line.get_ydata()))


def test_plot_cop_vhc_writes_files_and_highlights_only_non_dominated_designs(tmp_path, monkeypatch):
    rows, replays = _pareto_rows()
    add_pareto_metrics(rows, replays)
    figures = _capture_figures(monkeypatch)

    composition_pareto.plot_cop_vhc(rows, tmp_path, 200.0, status_note="draft")

    assert (tmp_path / "cop_vhc_pareto.png").stat().st_size > 0
    assert (tmp_path / "cop_vhc_pareto.svg").stat().st_size > 0
    # Max COP & max VHC: Pentane (3.8, 0.8), Isopentane (3.7, 1.0), R601A blend (3.6, 2.0),
    # R600A blend (3.2, 4.0), Propane (2.5, 6.0); sorted by VHC.
    assert _front_line(figures[0]) == [(0.8, 3.8), (1.0, 3.7), (2.0, 3.6), (4.0, 3.2), (6.0, 2.5)]
    assert "draft" in figures[0].axes[0].get_title()


def test_plot_cop_vhc_tolerates_failed_designs_without_replay(tmp_path, monkeypatch):
    rows, replays = _pareto_rows()
    rows.insert(1, {"case_signature": "Butane-70", "input_refrigerant": "Propane;Butane", "success": False,
                    "R290_mass_percent": 70, "family": "R290 / R600"})
    add_pareto_metrics(rows, replays)
    rows[1]["replay_COP"] = math.nan  # mix None and NaN as missing values
    figures = _capture_figures(monkeypatch)

    composition_pareto.plot_cop_vhc(rows, tmp_path, 200.0)

    assert (tmp_path / "cop_vhc_pareto.png").exists()
    assert rows[1]["non_dominated"] is False
    assert (math.nan, math.nan) not in _front_line(figures[0])
    assert len(_front_line(figures[0])) == 5
