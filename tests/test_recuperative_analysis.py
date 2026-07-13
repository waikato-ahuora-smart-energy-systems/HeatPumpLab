from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.utils import analysis as rhp_analysis


class FakeHeatPump:
    def __init__(self, sec=1.5):
        self.SEC = sec
        self.args = {"T_air_tar": 200.0}
        self.T_air_tar = 200.0
        self.T_air_inlet = 15.0
        self.Tref_evap_outlet = 10.0
        self.COP = 2.0
        self.feasible = True
        self.converged = True
        self.output = {}
        self.comp_discharge = SimpleNamespace(pressure=1.2)
        self.gas_cooler1 = SimpleNamespace(
            source_out=SimpleNamespace(temperature=80.0)
        )
        self.dT_gc_out_min_offset = 0.0

    def get_report(self):
        self.output = {"SEC": self.SEC}
        return self.output


def test_exact_state_seed_is_extracted_from_successful_solved_hp():
    hp = FakeHeatPump()
    hp.comp_discharge.pressure = 4.2
    hp.Tref_evap_outlet = 21.5

    seed = rhp_analysis._get_exact_state_seed_from_hp(hp)

    assert seed == {
        "P_H_est": pytest.approx(4.2),
        "Tref_evap_outlet": pytest.approx(21.5),
    }


@pytest.mark.parametrize(
    "updates",
    [
        {"feasible": False},
        {"converged": False},
        {"comp_discharge": SimpleNamespace(pressure=0.0)},
        {"comp_discharge": SimpleNamespace(pressure=-1.0)},
        {"comp_discharge": SimpleNamespace(pressure=float("nan"))},
        {"comp_discharge": None},
        {"Tref_evap_outlet": float("inf")},
    ],
)
def test_exact_state_seed_refuses_invalid_hp_states(updates):
    hp = FakeHeatPump()
    for key, value in updates.items():
        setattr(hp, key, value)

    assert rhp_analysis._get_exact_state_seed_from_hp(hp) is None


def test_exact_state_seed_is_applied_without_mutating_original_args():
    args = {"refrigerant": "Butane"}
    seed = {"P_H_est": 4.2, "Tref_evap_outlet": 21.5}

    seeded_args = rhp_analysis._apply_exact_state_seed(args, seed)

    assert seeded_args == {
        "refrigerant": "Butane",
        "P_H_est": pytest.approx(4.2),
        "Tref_evap_outlet": pytest.approx(21.5),
    }
    assert args == {"refrigerant": "Butane"}


def test_exact_solve_args_apply_matching_exact_state_seed_and_strip_objective_flag():
    seed_hp = FakeHeatPump()
    seed_hp.comp_discharge.pressure = 4.2
    seed_hp.Tref_evap_outlet = 21.5
    seed_hp.dT_gc_out_min_offset = 6.5
    case = {
        "additional_data": {
            "vars_input": {"dT_gc_out_min_offset": {"xi": 0.0, "bnds": (0.0, 40.0)}}
        },
        "fixed_discharge_pressure_objective": True,
    }

    args = rhp_analysis._get_exact_solve_args(case, 6.5, seed_hp)

    assert args["dT_gc_out_min_offset"] == pytest.approx(6.5)
    assert args["P_H_est"] == pytest.approx(4.2)
    assert args["Tref_evap_outlet"] == pytest.approx(21.5)
    assert "fixed_discharge_pressure_objective" not in args


def test_exact_solve_args_do_not_apply_mismatched_exact_state_seed():
    seed_hp = FakeHeatPump()
    seed_hp.comp_discharge.pressure = 4.2
    seed_hp.Tref_evap_outlet = 21.5
    seed_hp.dT_gc_out_min_offset = 0.0
    case = {
        "additional_data": {
            "vars_input": {"dT_gc_out_min_offset": {"xi": 0.0, "bnds": (0.0, 40.0)}}
        },
        "fixed_discharge_pressure_objective": True,
    }

    args = rhp_analysis._get_exact_solve_args(case, 6.5, seed_hp)

    assert args["dT_gc_out_min_offset"] == pytest.approx(6.5)
    assert "P_H_est" not in args
    assert "Tref_evap_outlet" not in args
    assert "fixed_discharge_pressure_objective" not in args


def test_analyse_rhp_case_optimises_when_boundary_gradient_is_positive(monkeypatch):
    hp = FakeHeatPump(sec=1.7)
    solved_args = []
    reports = []

    monkeypatch.setattr(rhp_analysis, "get_fprime", lambda *args, **kwargs: (1.0, 0.2))
    monkeypatch.setattr(
        rhp_analysis,
        "get_local_optimum",
        lambda *args, **kwargs: (
            SimpleNamespace(x=[7.5], success=True, nfev=4, message="ok"),
            0.3,
        ),
    )
    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 9.0)

    def solve_single_case(hp_cls, args):
        solved_args.append(args)
        return hp, 0.4

    def report(**kwargs):
        reports.append(kwargs)

    monkeypatch.setattr(rhp_analysis, "solve_single_case_with_timer", solve_single_case)
    monkeypatch.setattr(rhp_analysis, "get_results_report", report)

    solved_hp, solve_time = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        results_dir="custom-results",
        use_fast_objective=False,
    )

    assert solved_hp is hp
    assert solve_time == pytest.approx(0.4)
    assert solved_args[0]["dT_gc_out_min_offset"] == pytest.approx(7.5)
    assert solved_args[0]["additional_data"]["vars_input"] == {
        "dT_gc_out_min_offset": {"xi": 0.0, "bnds": (0.0, 40.0)}
    }
    assert hp.dQdW == pytest.approx(9.0)
    assert reports[0]["obj"] == pytest.approx(1.7)
    assert reports[0]["dt"] == pytest.approx(0.9)
    assert reports[0]["results_dir"] == "custom-results"


def test_analyse_rhp_case_skips_optimisation_when_boundary_gradient_is_not_positive(
    monkeypatch,
):
    hp = FakeHeatPump(sec=1.3)
    solved_args = []

    monkeypatch.setattr(
        rhp_analysis, "get_fprime", lambda *args, **kwargs: (-0.1, 0.2)
    )
    monkeypatch.setattr(
        rhp_analysis,
        "get_local_optimum",
        lambda *args, **kwargs: pytest.fail("optimisation should not run"),
    )
    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 0.0)
    monkeypatch.setattr(rhp_analysis, "get_results_report", lambda **kwargs: None)

    def solve_single_case(hp_cls, args):
        solved_args.append(args)
        return hp, 0.4

    monkeypatch.setattr(rhp_analysis, "solve_single_case_with_timer", solve_single_case)

    solved_hp, _ = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane", "dT_gc_out_min_offset": 12.0},
        use_fast_objective=False,
    )

    assert solved_hp is hp
    assert solved_args[0]["dT_gc_out_min_offset"] == pytest.approx(0.0)


def test_analyse_rhp_case_can_skip_candidate_result_persistence(monkeypatch):
    hp = FakeHeatPump(sec=1.3)
    solved_args = []
    monkeypatch.setattr(
        rhp_analysis, "get_fprime", lambda *args, **kwargs: (-0.1, 0.2)
    )
    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 0.0)
    monkeypatch.setattr(
        rhp_analysis,
        "solve_single_case_with_timer",
        lambda hp_cls, args: solved_args.append(dict(args)) or (hp, 0.4),
    )
    monkeypatch.setattr(
        rhp_analysis,
        "get_results_report",
        lambda **kwargs: pytest.fail("candidate result should not be persisted"),
    )

    solved_hp, _ = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        use_fast_objective=False,
        record_result=False,
    )

    assert solved_hp is hp
    assert solved_args[0]["record_failed_solve"] is False


def test_analyse_rhp_case_keeps_failure_recording_for_standalone_solves(monkeypatch):
    hp = FakeHeatPump(sec=1.3)
    solved_args = []
    monkeypatch.setattr(
        rhp_analysis, "get_fprime", lambda *args, **kwargs: (-0.1, 0.2)
    )
    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 0.0)
    monkeypatch.setattr(
        rhp_analysis,
        "solve_single_case_with_timer",
        lambda hp_cls, args: solved_args.append(dict(args)) or (hp, 0.4),
    )
    monkeypatch.setattr(rhp_analysis, "get_results_report", lambda **kwargs: None)

    rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        use_fast_objective=False,
        record_result=True,
    )

    assert solved_args[0]["record_failed_solve"] is True


def test_analyse_rhp_case_skips_dqdw_follow_up_for_failed_exact_solve(monkeypatch):
    hp = FakeHeatPump(sec=float("inf"))
    hp.feasible = False
    hp.converged = False
    monkeypatch.setattr(
        rhp_analysis, "get_fprime", lambda *args, **kwargs: (-0.1, 0.2)
    )
    monkeypatch.setattr(
        rhp_analysis,
        "solve_single_case_with_timer",
        lambda hp_cls, args: (hp, 0.4),
    )
    monkeypatch.setattr(
        rhp_analysis,
        "_get_dQdW",
        lambda solved_hp: pytest.fail("failed exact solve must not launch dQdW solve"),
    )

    solved_hp, _ = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        use_fast_objective=False,
        record_result=False,
    )

    assert solved_hp is hp
    assert solved_hp.dQdW == pytest.approx(0.0)


def test_analyse_rhp_case_returns_unsolved_when_local_optimisation_fails(monkeypatch):
    monkeypatch.setattr(rhp_analysis, "get_fprime", lambda *args, **kwargs: (1.0, 0.2))
    monkeypatch.setattr(
        rhp_analysis,
        "get_local_optimum",
        lambda *args, **kwargs: (
            SimpleNamespace(x=[7.5], success=False, nfev=4, message="failed"),
            0.3,
        ),
    )
    monkeypatch.setattr(
        rhp_analysis,
        "solve_single_case_with_timer",
        lambda *args, **kwargs: pytest.fail("failed optimisation should not solve"),
    )

    solved_hp, solve_time = rhp_analysis.analyse_RHP_case(
        FakeHeatPump, {"refrigerant": "Butane"}, use_fast_objective=False
    )

    assert solved_hp is None
    assert solve_time == pytest.approx(0.3)


def test_analyse_rhp_case_returns_failed_exact_seed_for_diagnostics(monkeypatch):
    failed_hp = FakeHeatPump(sec=float("inf"))
    failed_hp.feasible = False
    failed_hp.converged = False
    failed_hp.constraint_errors = {"superheater_temperature_reversal": 4.0}
    monkeypatch.setattr(
        rhp_analysis,
        "_get_fast_objective_case",
        lambda hp_cls, case, offset: (case, 0.2, False, failed_hp),
    )
    monkeypatch.setattr(rhp_analysis, "get_fprime", lambda *args, **kwargs: (1.0, 0.2))
    monkeypatch.setattr(
        rhp_analysis,
        "get_local_optimum",
        lambda *args, **kwargs: (
            SimpleNamespace(x=[0.0], success=False, nfev=3, message="infeasible"),
            0.3,
        ),
    )

    solved_hp, solve_time = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        use_fast_objective=True,
        record_result=False,
    )

    assert solved_hp is failed_hp
    assert solve_time == pytest.approx(0.2)
    assert solved_hp.constraint_errors == {
        "superheater_temperature_reversal": pytest.approx(4.0)
    }


def test_analyse_rhp_case_returns_unsolved_when_boundary_gradient_is_nan(monkeypatch):
    monkeypatch.setattr(
        rhp_analysis, "get_fprime", lambda *args, **kwargs: (float("nan"), 0.2)
    )
    monkeypatch.setattr(
        rhp_analysis,
        "get_local_optimum",
        lambda *args, **kwargs: pytest.fail("NaN gradient should not optimise"),
    )
    monkeypatch.setattr(
        rhp_analysis,
        "solve_single_case_with_timer",
        lambda *args, **kwargs: pytest.fail("NaN gradient should not solve"),
    )

    solved_hp, solve_time = rhp_analysis.analyse_RHP_case(
        FakeHeatPump, {"refrigerant": "Butane"}, use_fast_objective=False
    )

    assert solved_hp is None
    assert solve_time == pytest.approx(0.2)


def test_analyse_rhp_case_still_optimises_when_gradient_estimate_fails(monkeypatch):
    hp = FakeHeatPump(sec=1.7)
    optimised = []

    def raise_gradient_error(*args, **kwargs):
        raise RuntimeError("no gradient")

    def optimise(*args, **kwargs):
        optimised.append(kwargs["method"])
        return SimpleNamespace(x=[4.0], success=True, nfev=3, message="ok"), 0.3

    monkeypatch.setattr(rhp_analysis, "get_fprime", raise_gradient_error)
    monkeypatch.setattr(rhp_analysis, "get_local_optimum", optimise)
    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 0.0)
    monkeypatch.setattr(rhp_analysis, "get_results_report", lambda **kwargs: None)
    monkeypatch.setattr(
        rhp_analysis,
        "solve_single_case_with_timer",
        lambda hp_cls, args: (hp, 0.4),
    )

    solved_hp, _ = rhp_analysis.analyse_RHP_case(
        FakeHeatPump, {"refrigerant": "Butane"}, use_fast_objective=False
    )

    assert solved_hp is hp
    assert optimised == ["COBYQA"]


def test_analyse_rhp_case_uses_fixed_pressure_objective_then_exact_final_solve(
    monkeypatch,
):
    seed_hp = FakeHeatPump(sec=9.0)
    seed_hp.comp_discharge.pressure = 4.2
    seed_hp.Tref_evap_outlet = 21.5
    final_hp = FakeHeatPump(sec=1.2)
    solved_args = []
    objective_params = []
    reports = []

    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 0.0)
    monkeypatch.setattr(rhp_analysis, "get_results_report", lambda **kwargs: reports.append(kwargs))

    def solve_single_case(hp_cls, args):
        solved_args.append(args)
        if len(solved_args) == 1:
            return seed_hp, 0.11
        return final_hp, 0.44

    def get_fprime(hp_cls, param, f_obj):
        objective_params.append(("fprime", param))
        return 1.0, 0.2

    def get_local_optimum(hp_cls, param, vars_input, tol, f_obj, method):
        objective_params.append(("optimise", param))
        return SimpleNamespace(x=[6.5], success=True, nfev=5, message="ok"), 0.3

    monkeypatch.setattr(rhp_analysis, "solve_single_case_with_timer", solve_single_case)
    monkeypatch.setattr(rhp_analysis, "get_fprime", get_fprime)
    monkeypatch.setattr(rhp_analysis, "get_local_optimum", get_local_optimum)

    solved_hp, solve_time = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        results_dir="custom-results",
    )

    assert solved_hp is final_hp
    assert solve_time == pytest.approx(0.44)
    assert solved_args[0]["dT_gc_out_min_offset"] == pytest.approx(0.0)
    assert "fixed_discharge_pressure_objective" not in solved_args[0]
    assert solved_args[1]["dT_gc_out_min_offset"] == pytest.approx(6.5)
    assert "fixed_discharge_pressure_objective" not in solved_args[1]
    assert "P_H_est" not in solved_args[1]
    assert "Tref_evap_outlet" not in solved_args[1]
    assert objective_params[0][1]["P_H_est"] == pytest.approx(4.2)
    assert objective_params[1][1]["P_H_est"] == pytest.approx(4.2)
    assert objective_params[0][1]["Tref_evap_outlet"] == pytest.approx(21.5)
    assert objective_params[1][1]["Tref_evap_outlet"] == pytest.approx(21.5)
    assert objective_params[0][1]["fixed_discharge_pressure_objective"] is True
    assert objective_params[1][1]["fixed_discharge_pressure_objective"] is True
    assert reports[0]["obj"] == pytest.approx(1.2)
    assert reports[0]["dt"] == pytest.approx(0.11 + 0.2 + 0.3 + 0.44)
    assert reports[0]["results_dir"] == "custom-results"


def test_analyse_rhp_case_reuses_exact_seed_when_selected_offset_is_unchanged(
    monkeypatch,
):
    seed_hp = FakeHeatPump(sec=1.3)
    seed_hp.comp_discharge.pressure = 4.2
    solved_args = []
    reports = []

    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 0.0)
    monkeypatch.setattr(
        rhp_analysis, "get_results_report", lambda **kwargs: reports.append(kwargs)
    )
    monkeypatch.setattr(rhp_analysis, "get_fprime", lambda *args, **kwargs: (-0.1, 0.2))
    monkeypatch.setattr(
        rhp_analysis,
        "get_local_optimum",
        lambda *args, **kwargs: pytest.fail("optimisation should not run"),
    )

    def solve_single_case(hp_cls, args):
        solved_args.append(args)
        return seed_hp, 0.11

    monkeypatch.setattr(rhp_analysis, "solve_single_case_with_timer", solve_single_case)

    solved_hp, solve_time = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        results_dir="custom-results",
    )

    assert solved_hp is seed_hp
    assert solve_time == pytest.approx(0.11)
    assert len(solved_args) == 1
    assert reports[0]["obj"] == pytest.approx(1.3)
    assert reports[0]["dt"] == pytest.approx(0.11 + 0.2)


def test_analyse_rhp_case_does_not_reuse_failed_seed_solve(monkeypatch):
    failed_seed = FakeHeatPump(sec=9.9)
    failed_seed.feasible = False
    failed_seed.converged = False
    final_hp = FakeHeatPump(sec=1.1)
    solved_args = []
    reports = []

    monkeypatch.setattr(rhp_analysis, "_get_dQdW", lambda solved_hp: 0.0)
    monkeypatch.setattr(
        rhp_analysis, "get_results_report", lambda **kwargs: reports.append(kwargs)
    )
    monkeypatch.setattr(rhp_analysis, "get_fprime", lambda *args, **kwargs: (-0.1, 0.2))
    monkeypatch.setattr(
        rhp_analysis,
        "get_local_optimum",
        lambda *args, **kwargs: pytest.fail("optimisation should not run"),
    )

    def solve_single_case(hp_cls, args):
        solved_args.append(args)
        if len(solved_args) == 1:
            return failed_seed, 0.11
        return final_hp, 0.44

    monkeypatch.setattr(rhp_analysis, "solve_single_case_with_timer", solve_single_case)

    solved_hp, solve_time = rhp_analysis.analyse_RHP_case(
        FakeHeatPump,
        {"refrigerant": "Butane"},
        results_dir="custom-results",
    )

    assert solved_hp is final_hp
    assert solve_time == pytest.approx(0.44)
    assert len(solved_args) == 2
    assert reports[0]["obj"] == pytest.approx(1.1)


def test_dqdw_uses_followup_case_to_estimate_heating_sensitivity(monkeypatch):
    hp0 = FakeHeatPump()
    hp0.Tref_evap_outlet = 21.5
    solved_args = []

    def solve_single_case(hp_cls, args):
        solved_args.append(args)
        hp1 = FakeHeatPump()
        hp1.T_air_tar = args["T_air_tar"]
        hp1.T_air_inlet = hp0.T_air_inlet
        hp1.COP = 1.9
        hp1.feasible = True
        return hp1, 0.1

    monkeypatch.setattr(rhp_analysis, "solve_single_case_with_timer", solve_single_case)

    dqdw = rhp_analysis._get_dQdW(hp0)

    base = (hp0.T_air_tar - hp0.T_air_inlet) / hp0.COP
    followup = (hp0.T_air_tar + 0.001 - hp0.T_air_inlet) / 1.9
    assert dqdw == pytest.approx(0.001 / (followup - base))
    assert solved_args[0]["P_H_est"] == pytest.approx(1.2)
    assert solved_args[0]["T_gc_out"] == pytest.approx(80.0)
    assert solved_args[0]["Tref_evap_outlet"] == pytest.approx(21.5)


def test_dqdw_is_zero_when_followup_case_is_infeasible(monkeypatch):
    hp0 = FakeHeatPump()

    def solve_single_case(hp_cls, args):
        hp1 = FakeHeatPump()
        hp1.feasible = False
        return hp1, 0.1

    monkeypatch.setattr(rhp_analysis, "solve_single_case_with_timer", solve_single_case)

    assert rhp_analysis._get_dQdW(hp0) == pytest.approx(0.0)


def test_analyse_unsolved_cases_retries_with_typed_args_and_records_resolved_results(
    monkeypatch,
    tmp_path,
):
    analysed = []
    printed = []
    transferred = []

    class RetryHeatPump(FakeHeatPump):
        @staticmethod
        def UNSOLVED_CASE_LIST_FN(results_dir):
            assert results_dir == tmp_path
            return [
                {
                    "refrigerant": "Butane",
                    "dT_GC_min_tar": "20",
                    "T_air_tar": "200",
                    "T_air_inlet": "15",
                    "T_source_in": "40",
                    "n_comp": "2",
                    "ihx_on": "1",
                    "dTmin_rhx": "5",
                    "dT_rhx1_cold": "0",
                    "dT_rhx2_cold": "0",
                    "dT_rhx3_cold": "0",
                    "comp2_power_frac": "0.5",
                    "comp3_power_frac": "0",
                    "comp_PR_uniform": "",
                    "comp_model_index": "1",
                    "motor_eff": "0.96",
                    "sink": "Air",
                    "source": "water",
                    "Tref_evap_inlet": "10",
                    "Pmax": "4200",
                    "Tref_evap_outlet": "21.5",
                }
            ]

    def analyse_case(hp_cls, args, results_dir):
        analysed.append((hp_cls, args, results_dir))
        return RetryHeatPump(), 0.1

    monkeypatch.setattr(rhp_analysis, "analyse_RHP_case", analyse_case)
    monkeypatch.setattr(rhp_analysis, "get_unique_f_suffix", lambda args: "-case")
    monkeypatch.setattr(
        rhp_analysis,
        "print_results",
        lambda output, fsuffix, results_dir: printed.append(
            (output, fsuffix, results_dir)
        ),
    )
    monkeypatch.setattr(
        rhp_analysis,
        "transfer_temp_files_to_results",
        lambda file_name, results_dir: transferred.append((file_name, results_dir)),
    )

    rhp_analysis.analyse_unsolved_RHP_cases(
        RetryHeatPump, debug=True, results_dir=tmp_path
    )

    assert analysed[0][0] is RetryHeatPump
    assert analysed[0][1]["n_comp"] == 2
    assert analysed[0][1]["T_air_tar"] == pytest.approx(200.0)
    assert "P_H_est" not in analysed[0][1]
    assert "Tref_evap_outlet" not in analysed[0][1]
    assert analysed[0][1]["debug"] is True
    assert analysed[0][2] == tmp_path
    assert printed == [({"SEC": 1.5}, "-case-RESOLVED", tmp_path)]
    assert transferred == [("Cases-Resolved", tmp_path)]
