import csv

import pytest

import main_RHPcomposition as sweep
from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import normalise_refrigerant_input
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._cases import _build_optimisation_cases


def test_sweep_covers_requested_mass_grid_and_reuses_pure_endpoints():
    points = sweep.composition_points()
    assert len(points) == 45
    assert len({point.refrigerant for point in points}) == 41
    assert points[0].refrigerant == "Butane"
    for family, partner in (
        ("R290 / R600", "Butane"),
        ("R290 / R601A", "Isopentane"),
        ("R290 / R601", "Pentane"),
        ("R290 / R600A", "Isobutane"),
    ):
        family_points = [p for p in points if p.family == family]
        assert [p.propane_mass_percent for p in family_points] == list(range(0, 101, 10))
        assert family_points[0].refrigerant == partner
        assert family_points[-1].refrigerant == "Propane"
        for point in family_points[1:-1]:
            spec = normalise_refrigerant_input(point.refrigerant)
            assert spec.components == ("Propane", partner)
            assert spec.basis == "mass"
            assert spec.fractions == pytest.approx(
                (point.propane_mass_percent / 100, 1 - point.propane_mass_percent / 100)
            )


def test_batch_cases_keep_composition_fixed_and_cycle_design_free(monkeypatch, tmp_path):
    calls = []

    def capture(**kwargs):
        calls.append(kwargs)
        return []

    monkeypatch.setattr(sweep, "optimise_recuperative_heat_pump_cases", capture)
    sweep.run_sweep(results_dir=tmp_path, n_cpu=2)
    [call] = calls
    cases = _build_optimisation_cases(param_kwargs=call["param_kwargs"], config=call["config"])
    assert len(cases) == len({case.signature for case in cases}) == 41
    assert call["results_dir"] == tmp_path
    assert call["n_cpu"] == 2
    for case in cases:
        assert case.base_case["T_air_tar"] == 200.0
        assert case.config.vary_composition is False
        assert case.config.preserve_incumbent is False
        assert case.config.vary_dT_rhx1_cold is True
        assert case.config.vary_dT_rhx2_cold is True
        assert case.config.vary_comp2_power_frac is True


@pytest.fixture
def completed_sweep(tmp_path):
    rows = []
    for index, refrigerant in enumerate(dict.fromkeys(p.refrigerant for p in sweep.composition_points())):
        spec = normalise_refrigerant_input(refrigerant)
        rows.append({
            "input_refrigerant": spec.format(),
            "T_air_tar": 200,
            "config_vary_composition": False,
            "batch_success": True,
            "best_mass_fractions": ";".join(str(value) for value in spec.fractions),
            "case_signature": str(index),
            "best_COP": 2.5 if index == 0 else 3.0,
            "best_SEC": 2 / 3 if index == 0 else 0.5,
            "best_Pmax": 4000,
            "best_dT_rhx1_cold": 20,
            "best_dT_rhx2_cold": 30,
            "best_comp2_power_frac": 0.4,
            "exact_evaluations": 80,
            "failed_evaluations": 0,
            "timed_out_evaluations": 0,
            "elapsed_seconds": 50,
        })

    def write(records=rows):
        with (tmp_path / f"{sweep.DEFAULT_FILENAME}.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(records)

    write()
    return rows, write


def test_analysis_converts_pressure_and_compares_power_at_equal_heat_duty(completed_sweep, tmp_path):
    rows = sweep.load_sweep_rows(tmp_path)
    assert len(rows) == 45
    assert rows[0]["Pmax_MPa"] == pytest.approx(4.0)
    assert rows[0]["power_reduction_vs_R600_percent"] == 0
    assert rows[2]["power_reduction_vs_R600_percent"] == pytest.approx(100 * (1 - 2.5 / 3))
    assert rows[0]["case_signature"] == rows[1]["case_signature"]


def test_analysis_rejects_changed_composition(completed_sweep, tmp_path):
    rows, write = completed_sweep
    rows[1]["best_mass_fractions"] = "0.72;0.28"
    write()
    with pytest.raises(ValueError, match="changed composition"):
        sweep.load_sweep_rows(tmp_path)


def test_analysis_rejects_incomplete_grid_and_preserves_failed_cases(completed_sweep, tmp_path):
    rows, write = completed_sweep
    write(rows[:-1])
    with pytest.raises(ValueError, match="incomplete"):
        sweep.load_sweep_rows(tmp_path)
    rows[1]["batch_success"] = False
    rows[1]["best_COP"] = ""
    write()
    failed = sweep.load_sweep_rows(tmp_path)[2]
    assert failed["success"] is False
    assert failed["COP"] is None
    assert failed["power_reduction_vs_R600_percent"] is None
