import csv

from HeatPumpLab.cycles.RecuperativeHeatPump.utils import cases as case_utils
from HeatPumpLab.cycles.RecuperativeHeatPump.utils.cases import (
    get_RHP_case_list_to_solve,
    get_unsolved_RHP_case_list_to_resolve,
)


def write_rows(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def default_case(**overrides):
    case = {
        "refrigerant": "Butane",
        "n_comp": 2,
        "dT_GC_min_tar": 20.0,
        "T_air_inlet": 15.0,
        "T_air_tar": 200.0,
        "T_source_in": 40.0,
        "Tref_evap_inlet": 10.0,
        "dTmin_rhx": 5.0,
        "dT_rhx1_cold": 0.0,
        "dT_rhx2_cold": 0.0,
        "dT_rhx3_cold": 0.0,
        "comp2_power_frac": 0.5,
        "comp3_power_frac": 0.0,
    }
    return case | overrides


def test_rhp_case_generation_uses_defaults_when_inputs_are_omitted(tmp_path):
    cases = get_RHP_case_list_to_solve(file_name="cases", results_dir=tmp_path)

    assert len(cases) == 1
    assert cases[0]["refrigerant"] == "Butane"
    assert cases[0]["sink"] == "Air"
    assert cases[0]["source"] == "water"


def test_rhp_case_generation_returns_requested_factorial_cases(tmp_path):
    cases = get_RHP_case_list_to_solve(
        file_name="cases",
        refrigerant_ls=["Butane", "Propane"],
        dT_rhx1_cold_ls=[0.0, 5.0],
        dT_rhx2_cold_ls=[0.0],
        comp2_power_frac_ls=[0.4],
        results_dir=tmp_path,
    )

    signatures = {
        (case["refrigerant"], case["dT_rhx1_cold"], case["comp2_power_frac"])
        for case in cases
    }
    assert signatures == {
        ("Butane", 0.0, 0.4),
        ("Butane", 5.0, 0.4),
        ("Propane", 0.0, 0.4),
        ("Propane", 5.0, 0.4),
    }


def test_rhp_case_generation_filters_completed_cases_from_results_dir(tmp_path):
    write_rows(tmp_path / "cases.csv", [default_case()])

    cases = get_RHP_case_list_to_solve(
        file_name="cases",
        refrigerant_ls=["Butane"],
        dT_rhx1_cold_ls=[0.0],
        dT_rhx2_cold_ls=[0.0],
        comp2_power_frac_ls=[0.5],
        results_dir=tmp_path,
    )

    assert cases == []


def test_rhp_case_generation_keeps_cases_without_completed_signature(tmp_path):
    write_rows(tmp_path / "cases.csv", [default_case(comp2_power_frac=0.5)])

    cases = get_RHP_case_list_to_solve(
        file_name="cases",
        refrigerant_ls=["Butane"],
        dT_rhx1_cold_ls=[0.0],
        dT_rhx2_cold_ls=[0.0],
        comp2_power_frac_ls=[0.6],
        results_dir=tmp_path,
    )

    assert len(cases) == 1
    assert cases[0]["comp2_power_frac"] == 0.6


def test_rhp_case_generation_ignores_malformed_previous_result_rows(tmp_path):
    write_rows(tmp_path / "cases.csv", [{"refrigerant": "Butane", "n_comp": "2"}])

    cases = get_RHP_case_list_to_solve(
        file_name="cases",
        refrigerant_ls=["Butane"],
        dT_rhx1_cold_ls=[0.0],
        dT_rhx2_cold_ls=[0.0],
        comp2_power_frac_ls=[0.5],
        results_dir=tmp_path,
    )

    assert len(cases) == 1
    assert cases[0]["refrigerant"] == "Butane"


def test_rhp_case_generation_ignores_previous_rows_without_case_identity(tmp_path):
    write_rows(tmp_path / "cases.csv", [{"refrigerant": "Butane", "n_comp": "bad"}])

    cases = get_RHP_case_list_to_solve(
        file_name="cases",
        refrigerant_ls=["Butane"],
        dT_rhx1_cold_ls=[0.0],
        dT_rhx2_cold_ls=[0.0],
        comp2_power_frac_ls=[0.5],
        results_dir=tmp_path,
    )

    assert len(cases) == 1
    assert cases[0]["refrigerant"] == "Butane"


def test_completed_case_filter_can_use_supplied_comparison_without_signatures(
    monkeypatch,
):
    monkeypatch.setattr(case_utils, "transfer_temp_files_to_results", lambda **_: None)
    monkeypatch.setattr(
        case_utils,
        "read_csv_to_results",
        lambda **_: [{"case": "done"}],
    )

    remaining = case_utils._remove_completed_cases_from_list(
        cases=[{"case": "done"}, {"case": "pending"}],
        file_name="cases",
        results_dir="unused",
        f_comparison=lambda result, case: result["case"] == case["case"],
        f_signature=None,
    )

    assert remaining == [{"case": "pending"}]


def test_completed_case_filter_tolerates_unusable_signatures(monkeypatch):
    monkeypatch.setattr(case_utils, "transfer_temp_files_to_results", lambda **_: None)
    monkeypatch.setattr(
        case_utils,
        "read_csv_to_results",
        lambda **_: [{"case": "done"}, {"case": "bad"}],
    )

    def signature(record):
        if record["case"] == "bad":
            raise ValueError("signature unavailable")
        return record["case"]

    remaining = case_utils._remove_completed_cases_from_list(
        cases=[{"case": "done"}, {"case": "pending"}, {"case": "bad"}],
        file_name="cases",
        results_dir="unused",
        f_comparison=lambda result, case: result["case"] == case["case"],
        f_signature=signature,
    )

    assert remaining == [{"case": "pending"}, {"case": "bad"}]


def test_unsolved_case_resolution_removes_cases_that_now_have_results(tmp_path):
    write_rows(
        tmp_path / "Cases-Unsolved.csv",
        [
            default_case(refrigerant="Butane"),
            default_case(refrigerant="Propane"),
        ],
    )
    write_rows(tmp_path / "Cases-Resolved.csv", [default_case(refrigerant="Butane")])

    remaining = get_unsolved_RHP_case_list_to_resolve(results_dir=tmp_path)

    assert len(remaining) == 1
    assert remaining[0]["refrigerant"] == "Propane"
