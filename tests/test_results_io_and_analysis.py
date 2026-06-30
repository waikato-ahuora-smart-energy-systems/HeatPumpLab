from types import SimpleNamespace

import numpy as np
import pytest

from HeatPumpLab.utils import analysis
from HeatPumpLab.utils.analysis import (
    compute_SEC_obj,
    get_fprime,
    get_local_optimum,
    get_parameter_sweep,
    get_results_report,
    get_unique_f_suffix,
    get_var_kwargs,
    print_results,
    solve_single_case_with_timer,
)
from HeatPumpLab.utils.readWriteOutput import (
    print_result_to_csv,
    print_result_to_terminal,
    read_csv_to_results,
    transfer_temp_files_to_results,
)


def test_csv_results_roundtrip_uses_requested_results_directory(tmp_path):
    print_result_to_csv(
        output={"case": "A", "COP": 3.2},
        file_name="summary",
        results_dir=tmp_path,
    )
    print_result_to_csv(
        output={"case": "B", "COP": 3.4},
        file_name="summary",
        results_dir=tmp_path,
    )

    rows = read_csv_to_results(file_name="summary", results_dir=tmp_path)

    assert rows == [{"case": "A", "COP": "3.2"}, {"case": "B", "COP": "3.4"}]


def test_result_writer_preserves_legacy_results_relative_file_names(tmp_path):
    print_result_to_csv(
        output={"case": "A"},
        file_name="results/summary",
        results_dir=tmp_path,
    )

    assert read_csv_to_results(file_name="summary", results_dir=tmp_path) == [
        {"case": "A"}
    ]
    assert not (tmp_path / "results").exists()


def test_missing_csv_returns_empty_result_list(tmp_path, capsys):
    rows = read_csv_to_results(file_name="missing", results_dir=tmp_path)

    assert rows == []
    assert "File not found" in capsys.readouterr().out


def test_terminal_result_printer_reports_dict_entries(capsys):
    print_result_to_terminal({"COP": 3.2, "SEC": 0.4})

    out = capsys.readouterr().out
    assert "COP 3.2" in out
    assert "SEC 0.4" in out


def test_terminal_result_printer_rejects_non_dict_output(capsys):
    print_result_to_terminal(["not", "a", "dict"])

    assert "wrong output format" in capsys.readouterr().out


def test_csv_writer_rejects_non_dict_output(tmp_path, capsys):
    print_result_to_csv(output=[], file_name="summary", results_dir=tmp_path)

    assert "wrong output format" in capsys.readouterr().out
    assert not (tmp_path / "summary.csv").exists()


def test_csv_writer_reports_file_system_failures(monkeypatch, tmp_path, capsys):
    def fail_open(*args, **kwargs):
        raise OSError("disk unavailable")

    monkeypatch.setattr("builtins.open", fail_open)

    print_result_to_csv(
        output={"case": "A"},
        file_name="summary",
        results_dir=tmp_path,
    )

    assert "Failed to append result to CSV file: disk unavailable" in (
        capsys.readouterr().out
    )


def test_csv_reader_reports_file_system_failures(monkeypatch, tmp_path, capsys):
    (tmp_path / "summary.csv").write_text("case\nA\n", encoding="utf-8")

    def fail_open(*args, **kwargs):
        raise OSError("disk unavailable")

    monkeypatch.setattr("builtins.open", fail_open)

    rows = read_csv_to_results(file_name="summary", results_dir=tmp_path)

    assert rows == []
    assert "Failed to read CSV file" in capsys.readouterr().out


def test_temporary_result_files_are_merged_into_requested_results_dir(tmp_path):
    temp_dir = tmp_path / "temp"
    temp_dir.mkdir()
    (temp_dir / "a.csv").write_text("case,COP\nA,3.2\n", encoding="utf-8")
    (temp_dir / "b.csv").write_text("case,COP\nB,3.4\n", encoding="utf-8")

    transfer_temp_files_to_results(file_name="merged", results_dir=tmp_path)

    rows = read_csv_to_results(file_name="merged", results_dir=tmp_path)
    assert {tuple(row.items()) for row in rows} == {
        (("case", "A"), ("COP", "3.2")),
        (("case", "B"), ("COP", "3.4")),
    }
    assert not list(temp_dir.glob("*.csv"))


def test_print_results_writes_to_temp_result_file(tmp_path):
    print_results(
        results={"case": "A", "COP": 3.2},
        fsuffix="-butane",
        print_to_terminal=False,
        results_dir=tmp_path,
    )

    assert read_csv_to_results(
        file_name="temp/R-butane", results_dir=tmp_path
    ) == [{"case": "A", "COP": "3.2"}]


def test_print_results_optionally_prints_and_persists(monkeypatch):
    calls = []

    monkeypatch.setattr(
        analysis,
        "print_result_to_terminal",
        lambda results: calls.append(("terminal", results)),
    )
    monkeypatch.setattr(
        analysis,
        "print_result_to_csv",
        lambda **kwargs: calls.append(("csv", kwargs)),
    )

    print_results(
        results={"case": "A"},
        fsuffix="-butane",
        print_to_terminal=True,
        results_dir="custom-results",
    )

    assert calls == [
        ("terminal", {"case": "A"}),
        (
            "csv",
            {
                "output": {"case": "A"},
                "suffix": "-butane",
                "file_name": "temp/R",
                "results_dir": "custom-results",
            },
        ),
    ]


def test_get_var_kwargs_maps_optimisation_vector_to_domain_names():
    param = {
        "additional_data": {
            "vars_input": {
                "dT_gc_out_min_offset": {"xi": 0.0},
                "comp2_power_frac": {"xi": 0.5},
            }
        }
    }

    assert get_var_kwargs([12.5, 0.4], param) == {
        "dT_gc_out_min_offset": 12.5,
        "comp2_power_frac": 0.4,
    }


def test_get_var_kwargs_accepts_scalar_for_single_variable():
    param = {"additional_data": {"vars_input": {"dT_gc_out_min_offset": {"xi": 0.0}}}}

    assert get_var_kwargs(8.0, param) == {"dT_gc_out_min_offset": 8.0}


def test_compute_sec_objective_reports_heat_pump_sec():
    assert compute_SEC_obj(SimpleNamespace(SEC=1.23), param={}) == pytest.approx(1.23)


def test_solve_single_case_with_timer_instantiates_and_solves_heat_pump():
    class FakeHP:
        def __init__(self, **args):
            self.args = args
            self.solved = False

        def solve(self):
            self.solved = True

    hp, dt = solve_single_case_with_timer(FakeHP, {"refrigerant": "Butane"})

    assert hp.args == {"refrigerant": "Butane"}
    assert hp.solved is True
    assert dt >= 0.0


def test_get_local_optimum_passes_domain_initial_values_and_bounds(monkeypatch):
    calls = []

    def fake_minimize(fun, x0, bounds, method, tol):
        calls.append({"x0": x0, "bounds": bounds, "method": method, "tol": tol})
        return SimpleNamespace(success=True)

    monkeypatch.setattr(analysis, "loc_min", fake_minimize)

    result, dt = get_local_optimum(
        HP_fn=object,
        param={},
        vars_input={
            "dT_gc_out_min_offset": {"xi": 1.0, "bnds": (0.0, 40.0)},
            "comp2_power_frac": {"xi": 0.5, "bnds": (0.0, 1.0)},
        },
        tol=1e-4,
        f_obj=lambda hp, param: hp.SEC,
        method="COBYQA",
    )

    assert result.success is True
    assert dt >= 0.0
    assert calls == [
        {
            "x0": [1.0, 0.5],
            "bounds": [(0.0, 40.0), (0.0, 1.0)],
            "method": "COBYQA",
            "tol": 1e-4,
        }
    ]


def test_get_local_optimum_evaluates_candidate_objective_with_domain_arguments(
    monkeypatch,
):
    solved_args = []

    def fake_solve_single_case_with_timer(HP_fn, args):
        solved_args.append((HP_fn, args))
        return SimpleNamespace(SEC=0.123), 0.0

    def fake_minimize(fun, x0, bounds, method, tol):
        return SimpleNamespace(fun=fun(np.array([0.75])))

    monkeypatch.setattr(
        analysis,
        "solve_single_case_with_timer",
        fake_solve_single_case_with_timer,
    )
    monkeypatch.setattr(analysis, "loc_min", fake_minimize)

    param = {
        "refrigerant": "Butane",
        "additional_data": {
            "vars_input": {
                "comp2_power_frac": {"xi": 0.5, "bnds": (0.0, 1.0)}
            }
        },
    }

    result, dt = get_local_optimum(
        HP_fn=object,
        param=param,
        vars_input=param["additional_data"]["vars_input"],
        tol=1e-4,
        f_obj=compute_SEC_obj,
        method="COBYQA",
    )

    assert result.fun == pytest.approx(123.0)
    assert dt >= 0.0
    assert solved_args == [
        (
            object,
            param | {"comp2_power_frac": 0.75},
        )
    ]


def test_get_fprime_returns_negative_objective_derivative(monkeypatch):
    monkeypatch.setattr(
        analysis,
        "approx_fprime",
        lambda xk, f, epsilon: np.array([2.5]),
    )

    fprime, dt = get_fprime(object, param={}, f_obj=lambda hp, param: hp.SEC)

    assert fprime == pytest.approx(np.array([-2.5]))
    assert dt >= 0.0


def test_get_results_report_adds_optimisation_metadata_and_weight(monkeypatch):
    reports = []

    class FakeHP:
        def __init__(self):
            self.output = {}

        def get_report(self):
            self.output = {"COP": 3.2}

    monkeypatch.setattr(analysis, "get_unique_f_suffix", lambda param: "-Butane")
    monkeypatch.setattr(
        analysis,
        "print_results",
        lambda results, fsuffix, print_to_terminal, results_dir: reports.append(
            (results, fsuffix, print_to_terminal, results_dir)
        ),
    )

    get_results_report(
        hp=FakeHP(),
        param={"additional_data": {"weight": 2.0, "SEC_0": 0.5, "P_0": 100.0}},
        result={"nfev": 7, "message": "ok", "success": True},
        obj=0.45,
        dt=1.2,
        print_to_terminal=True,
        results_dir="custom-results",
    )

    assert reports == [
        (
            {
                "COP": 3.2,
                "Weight": 2.0,
                "SEC_0": 0.5,
                "P_0": 100.0,
                "Obj_function": 0.45,
                "Function evaluations": 7,
                "Message": "ok",
                "Success": True,
                "Time_to_optimise": 1.2,
            },
            "-Butane",
            True,
            "custom-results",
        )
    ]


def test_get_unique_suffix_uses_refrigerant_name_for_pure_fluid(monkeypatch):
    monkeypatch.setattr(
        analysis,
        "PropertyPackage",
        lambda refrigerant: SimpleNamespace(
            fld_name=refrigerant,
            fld_composition=[1.0],
        ),
    )

    assert get_unique_f_suffix({"refrigerant": "Butane"}) == "-Butane"


def test_get_unique_suffix_includes_mixture_components_and_optional_parameters(
    monkeypatch,
):
    monkeypatch.setattr(
        analysis,
        "PropertyPackage",
        lambda refrigerant: SimpleNamespace(
            fld_name="Propane.FLD*Butane.FLD",
            fld_composition=[0.4, 0.6],
        ),
    )

    suffix = get_unique_f_suffix(
        {"refrigerant": "mix", "T_air_tar": 200.0, "additional_data": {}},
        fld_only=True,
    )

    assert suffix == "-Propane0_4*Butane0_6-200"


def test_parameter_sweep_generates_cases_and_solves_them_in_configured_results_dir(
    tmp_path,
):
    generated = []
    solved = []

    def case_list_fn(file_name, results_dir, refrigerant_ls=None):
        generated.append(
            {
                "file_name": file_name,
                "results_dir": results_dir,
                "refrigerant_ls": refrigerant_ls,
            }
        )
        return [{"case": "A"}, {"case": "B"}]

    def solve_fn(hp_cls, case):
        solved.append((hp_cls, case))

    fake_hp = SimpleNamespace(
        CASE_LIST_FN=case_list_fn,
        ANAYLSE_CASE_FN=solve_fn,
    )

    get_parameter_sweep(
        fake_hp,
        filename="sweep",
        param_kwargs={"refrigerant_ls": ["Butane"]},
        n_cpu=1,
        results_dir=tmp_path,
    )

    assert generated == [
        {
            "file_name": "sweep",
            "results_dir": tmp_path,
            "refrigerant_ls": ["Butane"],
        }
    ]
    assert solved == [(fake_hp, {"case": "A"}), (fake_hp, {"case": "B"})]
    assert (tmp_path / "sweep.csv").exists()


def test_parameter_sweep_accepts_omitted_parameter_kwargs(tmp_path):
    generated = []

    def case_list_fn(file_name, results_dir):
        generated.append({"file_name": file_name, "results_dir": results_dir})
        return []

    fake_hp = SimpleNamespace(
        CASE_LIST_FN=case_list_fn,
        ANAYLSE_CASE_FN=lambda hp_cls, case: None,
    )

    get_parameter_sweep(
        fake_hp,
        filename="sweep",
        param_kwargs=None,
        n_cpu=1,
        results_dir=tmp_path,
    )

    assert generated == [{"file_name": "sweep", "results_dir": tmp_path}]
    assert (tmp_path / "sweep.csv").exists()
