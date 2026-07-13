from __future__ import annotations

import inspect
import math
from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump import optimisation as outer_opt
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    _objective as objective_helpers,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _batch as batch_helpers
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _cases as case_helpers
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import outer_design as outer_design_module
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _polish as polish_helpers
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    _results_io as results_io,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _types as opt_types


class FakeHeatPump:
    def __init__(self, sec: float, case: dict):
        self.SEC = sec
        self.COP = 2.5
        self.feasible = True
        self.converged = True
        self.args = dict(case)
        self.output = {
            "SEC": sec,
            "COP": self.COP,
            "Pmax": 4100.0,
            "converged": True,
            "Tref_evap_outlet": 21.0,
        }
        self.comp_discharge = SimpleNamespace(pressure=4.1)
        self.Tref_evap_outlet = 21.0
        self.dT_gc_out_min_offset = 0.0

    def get_report(self, dt=0):
        self.output["Solve time"] = dt
        return self.output


class SyntheticAnalyser:
    def __init__(
        self,
        *,
        dT1_target=20.0,
        dT2_target=35.0,
        comp2_target=0.45,
        composition_target=0.72,
        fail=False,
    ):
        self.calls = []
        self.hp_classes = []
        self.fail = fail
        self.dT1_target = dT1_target
        self.dT2_target = dT2_target
        self.comp2_target = comp2_target
        self.composition_target = composition_target

    def __call__(self, hp_cls, case, results_dir, use_fast_objective=True):
        del results_dir, use_fast_objective
        self.hp_classes.append(hp_cls)
        self.calls.append(dict(case))
        if self.fail:
            raise RuntimeError("synthetic solve failed")
        frac = _first_mass_fraction(case["refrigerant"])
        sec = (
            1.0
            + ((case["dT_rhx1_cold"] - self.dT1_target) / 50.0) ** 2
            + ((case["dT_rhx2_cold"] - self.dT2_target) / 50.0) ** 2
            + ((case["comp2_power_frac"] - self.comp2_target) / 0.3) ** 2
            + ((frac - self.composition_target) / 0.5) ** 2
        )
        return FakeHeatPump(sec=sec, case=case), 0.01


def test_public_api_exposes_only_the_rhp_optimiser_surface():
    assert outer_opt.__all__ == [
        "OPTIMISER_VERSION",
        "DesignEvaluation",
        "LineSearchResult",
        "LocalPolishResult",
        "OuterDesignOptimisationConfig",
        "OuterDesignOptimisationResult",
        "PhaseEvaluation",
        "optimise_recuperative_heat_pump_case",
        "optimise_recuperative_heat_pump_cases",
        "optimise_recuperative_heat_pump",
    ]
    assert hasattr(outer_opt, "optimise_recuperative_heat_pump")
    assert hasattr(outer_opt, "optimise_recuperative_heat_pump_case")
    assert hasattr(outer_opt, "optimise_recuperative_heat_pump_cases")
    assert not hasattr(outer_opt, "optimise_recuperative_outer_design")
    assert not hasattr(outer_opt, "OuterDesignObjective")
    assert not hasattr(outer_opt, "run_outer_design_optimisation")
    assert not hasattr(outer_opt, "run_local_polish")


def test_public_entry_point_has_no_test_only_dependency_parameters(
    monkeypatch,
    tmp_path,
):
    signature = inspect.signature(outer_opt.optimise_recuperative_heat_pump)
    assert "HP_fn" not in signature.parameters
    assert "analyser" not in signature.parameters
    assert "use_fast_objective" not in signature.parameters
    assert list(signature.parameters)[0] == "base_case"

    analyser = SyntheticAnalyser(dT1_target=0.0, dT2_target=35.0)
    objective_kwargs = []
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=8,
        line_budget=2,
        split_budget=0,
        polish_budget=0,
        use_fast_objective=False,
    )

    def build_synthetic_objective(**kwargs):
        objective_kwargs.append(dict(kwargs))
        return objective_helpers._OuterDesignObjective(
            **kwargs,
            analyser=analyser,
        )

    monkeypatch.setattr(
        outer_design_module,
        "_OuterDesignObjective",
        build_synthetic_objective,
    )
    result = outer_opt.optimise_recuperative_heat_pump(
        base_case={"T_air_tar": 180.0, "dT_rhx3_cold": 0.0},
        config=config,
        results_dir=tmp_path,
    )

    assert result.config is config
    assert result.exact_evaluations <= config.max_exact_evaluations
    assert result.best is not None
    assert objective_kwargs
    assert "analyser" not in objective_kwargs[0]
    assert objective_kwargs[0]["use_fast_objective"] is False
    assert analyser.calls


def test_private_objective_always_passes_recuperative_hp_to_analyser():
    analyser = SyntheticAnalyser()
    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=analyser,
    )

    record = objective.evaluate([20.0, 35.0, 0.45, 0.72])

    assert record.success is True
    assert {hp_cls.__name__ for hp_cls in analyser.hp_classes} == {"RecuperativeHP"}


def test_binary_refrigerant_string_uses_mass_basis_and_complement():
    refrigerant = objective_helpers._format_binary_refrigerant(
        "Propane",
        "Butane",
        0.568611,
    )

    assert refrigerant == "Propane;Butane|0.568611;0.431389|mass"


def test_objective_maps_vector_to_domain_case_without_mutating_base_case():
    base_case = {"T_air_tar": 180.0, "dT_rhx3_cold": 0.0}
    objective = objective_helpers._OuterDesignObjective(
        max_evals=5,
        base_case=base_case,
        fluids=("Propane", "Butane"),
        analyser=SyntheticAnalyser(),
    )

    case = objective.vector_to_case([12.0, 34.0, 0.42, 0.8])

    assert case["T_air_tar"] == pytest.approx(180.0)
    assert case["dT_rhx1_cold"] == pytest.approx(12.0)
    assert case["dT_rhx2_cold"] == pytest.approx(34.0)
    assert case["comp2_power_frac"] == pytest.approx(0.42)
    assert case["refrigerant"] == "Propane;Butane|0.8;0.2|mass"
    assert base_case == {"T_air_tar": 180.0, "dT_rhx3_cold": 0.0}


def test_objective_repairs_recuperator_lift_to_feasible_domain():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=5,
        base_case={"T_air_inlet": 10.0, "T_air_tar": 100.0, "dT_rhx3_cold": 0.0},
        fluids=("Propane", "Butane"),
        analyser=SyntheticAnalyser(),
    )

    case = objective.vector_to_case([60.0, 60.0, 0.42, 0.8])

    assert case["dT_rhx1_cold"] + case["dT_rhx2_cold"] == pytest.approx(60.0)
    assert case["dT_rhx1_cold"] == pytest.approx(30.0)
    assert case["dT_rhx2_cold"] == pytest.approx(30.0)


def test_objective_does_not_spend_exact_solve_on_setup_infeasible_lift():
    class RejectingAnalyser(SyntheticAnalyser):
        def __call__(self, hp_cls, case, results_dir, use_fast_objective=True):
            max_lift = (case["T_air_tar"] - case["T_air_inlet"]) * 2.0 / 3.0
            assert case["dT_rhx1_cold"] + case["dT_rhx2_cold"] <= max_lift
            return super().__call__(hp_cls, case, results_dir, use_fast_objective)

    analyser = RejectingAnalyser()
    objective = objective_helpers._OuterDesignObjective(
        max_evals=5,
        base_case={
            "T_air_inlet": 10.0,
            "T_air_tar": 100.0,
            "dT_rhx3_cold": 0.0,
        },
        fluids=("Propane", "Butane"),
        analyser=analyser,
    )

    record = objective.evaluate([80.0, 80.0, 0.42, 0.8])

    assert record.success is True
    assert record.x[0] + record.x[1] == pytest.approx(60.0)
    assert analyser.calls[0]["dT_rhx1_cold"] + analyser.calls[0]["dT_rhx2_cold"] <= 60.0


def test_objective_reuses_cached_exact_solution_for_duplicate_vectors():
    analyser = SyntheticAnalyser()
    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=analyser,
    )

    first = objective.evaluate([20.0, 35.0, 0.45, 0.72])
    second = objective.evaluate([20.0, 35.0, 0.45, 0.72])

    assert first.objective == pytest.approx(second.objective)
    assert objective.exact_evaluations == 1
    assert len(analyser.calls) == 1


def test_objective_returns_finite_penalty_for_failed_exact_solve():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=SyntheticAnalyser(fail=True),
    )

    record = objective.evaluate([20.0, 35.0, 0.45, 0.72])

    assert record.success is False
    assert math.isfinite(record.objective)
    assert record.objective == pytest.approx(opt_types.PENALTY_OBJECTIVE + 1.0)
    assert record.constraint_errors == {"solve_failure": pytest.approx(1.0)}
    assert "synthetic solve failed" in record.message
    assert objective.failed_evaluations == 1


def test_objective_penalty_uses_squared_physical_constraint_errors():
    def failed_analyser(error):
        def analyse(hp_cls, case, results_dir, use_fast_objective=True):
            del hp_cls, results_dir, use_fast_objective
            hp = FakeHeatPump(sec=1.0, case=case)
            hp.feasible = False
            hp.converged = False
            hp.constraint_errors = {"superheater_temperature_reversal": error}
            return hp, 0.01

        return analyse

    near = objective_helpers._OuterDesignObjective(
        max_evals=1,
        analyser=failed_analyser(2.0),
    ).evaluate([20.0, 35.0, 0.45, 0.72])
    far = objective_helpers._OuterDesignObjective(
        max_evals=1,
        analyser=failed_analyser(3.0),
    ).evaluate([20.0, 35.0, 0.45, 0.72])

    assert near.objective == pytest.approx(opt_types.PENALTY_OBJECTIVE + 4.0)
    assert far.objective == pytest.approx(opt_types.PENALTY_OBJECTIVE + 9.0)
    assert near.objective < far.objective
    assert near.constraint_errors == {
        "superheater_temperature_reversal": pytest.approx(2.0)
    }


def test_objective_does_not_exceed_exact_evaluation_budget():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        analyser=SyntheticAnalyser(),
    )

    first = objective.evaluate([20.0, 35.0, 0.45, 0.72])
    second = objective.evaluate([30.0, 35.0, 0.45, 0.72])

    assert first.success is True
    assert second.success is False
    assert second.message == "exact evaluation budget exhausted"
    assert objective.exact_evaluations == 1
    assert objective.budget_exhausted_calls == 1


def test_outer_design_optimisation_edge_a_holds_dT1_at_zero(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=0.0, dT2_target=42.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=12,
            split_budget=0,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.edge_a.success is True
    assert result.edge_a.best is not None
    assert result.edge_a.best.x[0] == pytest.approx(0.0)
    assert result.edge_a.best.x[1] == pytest.approx(42.0, abs=2.0)
    assert all(record.x[0] == pytest.approx(0.0) for record in result.edge_a.evaluations)


def test_outer_design_optimisation_edge_b_holds_dT2_at_zero(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=28.0, dT2_target=0.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=12,
            split_budget=0,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.edge_b.success is True
    assert result.edge_b.best is not None
    assert result.edge_b.best.x[0] == pytest.approx(28.0, abs=2.0)
    assert result.edge_b.best.x[1] == pytest.approx(0.0)
    assert all(record.x[1] == pytest.approx(0.0) for record in result.edge_b.evaluations)


def test_outer_design_optimisation_split_line_uses_sum_of_edge_optima(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=25.0, dT2_target=35.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=14,
            split_budget=10,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.split is not None
    assert result.edge_a.best is not None
    assert result.edge_b.best is not None
    expected_sum = result.edge_a.best.x[1] + result.edge_b.best.x[0]
    split_sums = [record.x[0] + record.x[1] for record in result.split.evaluations]
    assert split_sums
    assert all(total == pytest.approx(expected_sum, abs=1e-7) for total in split_sums)


def test_outer_design_optimisation_skips_split_when_both_edges_select_corner(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=0.0, dT2_target=0.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=10,
            split_budget=10,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.edge_a.best is not None
    assert result.edge_b.best is not None
    assert result.edge_a.best.x[:2] == pytest.approx((0.0, 0.0))
    assert result.edge_b.best.x[:2] == pytest.approx((0.0, 0.0))
    assert result.split is None


def test_outer_design_optimisation_fixed_dT1_searches_only_dT2(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=18.0, dT2_target=42.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            fixed_dT_rhx1_cold=18.0,
            vary_dT_rhx1_cold=False,
            line_budget=12,
            split_budget=8,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.edge_a.success is True
    assert result.edge_b.message == "dT_rhx1_cold is fixed"
    assert result.split is None
    assert result.edge_a.best is not None
    assert result.edge_a.best.x[0] == pytest.approx(18.0)
    assert result.edge_a.best.x[1] == pytest.approx(42.0, abs=2.0)
    assert all(record.x[0] == pytest.approx(18.0) for record in result.edge_a.evaluations)


def test_outer_design_optimisation_fixed_dT2_searches_only_dT1(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=27.0, dT2_target=16.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            fixed_dT_rhx2_cold=16.0,
            vary_dT_rhx2_cold=False,
            line_budget=12,
            split_budget=8,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.edge_a.message == "dT_rhx2_cold is fixed"
    assert result.edge_b.success is True
    assert result.split is None
    assert result.edge_b.best is not None
    assert result.edge_b.best.x[0] == pytest.approx(27.0, abs=2.0)
    assert result.edge_b.best.x[1] == pytest.approx(16.0)
    assert all(record.x[1] == pytest.approx(16.0) for record in result.edge_b.evaluations)


def test_outer_design_optimisation_fixed_dT_pair_skips_edge_and_split(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=11.0, dT2_target=22.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            fixed_dT_rhx1_cold=11.0,
            fixed_dT_rhx2_cold=22.0,
            vary_dT_rhx1_cold=False,
            vary_dT_rhx2_cold=False,
            line_budget=12,
            split_budget=8,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.corner.x[:2] == pytest.approx((11.0, 22.0))
    assert result.edge_a.evaluations == ()
    assert result.edge_b.evaluations == ()
    assert result.split is None
    assert result.best is not None
    assert result.best.x[:2] == pytest.approx((11.0, 22.0))


def test_outer_design_optimisation_polish_keeps_dT_values_fixed_and_varies_power_and_composition(
    monkeypatch,
):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(
            dT1_target=0.0,
            dT2_target=40.0,
            comp2_target=0.62,
            composition_target=0.8,
        ),
        config=outer_opt.OuterDesignOptimisationConfig(
            fixed_comp2_power_frac=0.4,
            fixed_first_component_mass_fraction=0.3,
            line_budget=10,
            split_budget=0,
            polish_budget=16,
            polish_top_n=1,
        ),
        monkeypatch=monkeypatch,
    )

    assert len(result.polish_results) == 1
    source = result.polish_results[0].source_x
    polish_records = result.polish_results[0].evaluations
    assert polish_records
    assert all(record.x[:2] == pytest.approx(source[:2]) for record in polish_records)
    assert any(
        record.x[2] != pytest.approx(source[2])
        or record.x[3] != pytest.approx(source[3])
        for record in polish_records
    )


def test_pattern_search_polish_improves_unsmooth_objective_without_leaving_bounds():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=20,
        analyser=SyntheticAnalyser(
            dT1_target=0.0,
            dT2_target=40.0,
            comp2_target=0.6,
            composition_target=0.6,
        ),
    )
    source = _phase_from_exact_record(
        objective,
        [0.0, 40.0, 0.5, 0.5],
        phase="source",
    )
    phase_records = []

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(budget=10),
        evaluate_phase=_phase_evaluator(objective, phase_records),
    )

    assert result.status == "ok"
    assert result.best is not None
    assert result.best.objective < source.objective
    assert opt_types.DEFAULT_BOUNDS[2, 0] <= result.best.x[2] <= opt_types.DEFAULT_BOUNDS[2, 1]
    assert opt_types.DEFAULT_BOUNDS[3, 0] <= result.best.x[3] <= opt_types.DEFAULT_BOUNDS[3, 1]


def test_pattern_search_uses_cache_for_duplicate_current_point():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=SyntheticAnalyser(),
    )
    source = _phase_from_exact_record(
        objective,
        [20.0, 35.0, 0.45, 0.72],
        phase="source",
    )
    exact_before = objective.exact_evaluations
    phase_records = []

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(budget=1),
        evaluate_phase=_phase_evaluator(objective, phase_records),
    )

    assert result.evaluations[0].x == pytest.approx(source.x)
    assert result.evaluations[0].cumulative_exact_evaluations == exact_before
    assert result.exact_evaluations == 1
    assert objective.exact_evaluations == exact_before + 1


def test_pattern_search_shrinks_steps_when_no_poll_improves():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=12,
        analyser=SyntheticAnalyser(
            dT1_target=20.0,
            dT2_target=35.0,
            comp2_target=0.5,
            composition_target=0.5,
        ),
    )
    source = _phase_from_exact_record(
        objective,
        [20.0, 35.0, 0.5, 0.5],
        phase="source",
    )
    phase_records = []

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            initial_steps=(0.05, 0.05),
            min_steps=(0.01, 0.01),
        ),
        evaluate_phase=_phase_evaluator(objective, phase_records),
    )

    comp2_offsets = sorted({round(abs(record.x[2] - 0.5), 3) for record in result.evaluations})
    assert 0.025 in comp2_offsets
    assert result.best is not None
    assert result.best.x[2:4] == pytest.approx((0.5, 0.5))


def test_pattern_search_keeps_dT_values_fixed():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=20,
        analyser=SyntheticAnalyser(comp2_target=0.6, composition_target=0.6),
    )
    source = _phase_from_exact_record(
        objective,
        [12.0, 44.0, 0.5, 0.5],
        phase="source",
    )
    phase_records = []

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(budget=8),
        evaluate_phase=_phase_evaluator(objective, phase_records),
    )

    assert result.evaluations
    assert all(record.x[:2] == pytest.approx((12.0, 44.0)) for record in result.evaluations)


def test_pattern_search_respects_fixed_comp2_power_fraction():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=12,
        analyser=SyntheticAnalyser(comp2_target=0.8, composition_target=0.7),
    )
    source = _phase_from_exact_record(
        objective,
        [12.0, 44.0, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            vary_comp2_power_frac=False,
            vary_first_component_mass_fraction=True,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.evaluations
    assert all(record.x[2] == pytest.approx(0.5) for record in result.evaluations)
    assert any(record.x[3] != pytest.approx(0.5) for record in result.evaluations)


def test_pattern_search_respects_fixed_first_component_mass_fraction():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=12,
        analyser=SyntheticAnalyser(comp2_target=0.7, composition_target=0.9),
    )
    source = _phase_from_exact_record(
        objective,
        [12.0, 44.0, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            vary_comp2_power_frac=True,
            vary_first_component_mass_fraction=False,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.evaluations
    assert any(record.x[2] != pytest.approx(0.5) for record in result.evaluations)
    assert all(record.x[3] == pytest.approx(0.5) for record in result.evaluations)


def test_pattern_search_skips_when_all_polish_variables_are_fixed():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=4,
        analyser=SyntheticAnalyser(),
    )
    source = _phase_from_exact_record(
        objective,
        [12.0, 44.0, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            vary_comp2_power_frac=False,
            vary_first_component_mass_fraction=False,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.status == "skipped"
    assert result.best is source
    assert result.evaluations == ()
    assert result.exact_evaluations == 0


def test_local_polish_is_pattern_search_only():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=10,
        analyser=SyntheticAnalyser(comp2_target=0.6, composition_target=0.6),
    )
    source = _phase_from_exact_record(
        objective,
        [20.0, 35.0, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(budget=5),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.method == "pattern_search"
    assert result.status == "ok"
    assert result.best is not None
    assert all(record.phase == "polish_pattern_search" for record in result.evaluations)


def test_outer_design_optimisation_shares_exact_cache_across_phases(monkeypatch):
    analyser = SyntheticAnalyser(dT1_target=0.0, dT2_target=0.0)

    result = _optimise_synthetic(
        analyser=analyser,
        config=outer_opt.OuterDesignOptimisationConfig(
            max_exact_evaluations=30,
            line_budget=3,
            split_budget=0,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    corner_visits = [
        record for record in result.evaluations if record.x[:2] == pytest.approx((0.0, 0.0))
    ]
    assert len(corner_visits) >= 3
    assert result.exact_evaluations < len(result.evaluations)
    assert len(analyser.calls) == result.exact_evaluations


def test_outer_design_optimisation_result_identifies_final_best_candidate(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(
            dT1_target=0.0,
            dT2_target=45.0,
            comp2_target=0.5,
            composition_target=0.5,
        ),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=14,
            split_budget=8,
            polish_budget=8,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.optimiser_version == outer_opt.OPTIMISER_VERSION
    assert result.best is not None
    assert result.best.objective == min(
        record.objective for record in result.evaluations if record.success
    )
    assert result.to_summary_dict()["optimiser_version"] == outer_opt.OPTIMISER_VERSION


def test_case_definition_loads_fixed_case_and_reference_from_results_csv(tmp_path):
    csv_path = tmp_path / "results.csv"
    csv_path.write_text(
        "\n".join(
            [
                "SEC,COP,Pmax,Time_to_optimise,refrigerant,T_source_in,Tref_evap_inlet,T_air_inlet,T_air_tar,dT_GC_min_tar,dTmin_rhx,ihx_on,n_comp,comp_PR_uniform,comp_model_index,motor_eff,sink,source,dT_rhx1_cold,dT_rhx2_cold,dT_rhx3_cold,comp2_power_frac,comp3_power_frac,converged",
                "SEC,COP,Pmax,Time_to_optimise,refrigerant,T_source_in,Tref_evap_inlet,T_air_inlet,T_air_tar,dT_GC_min_tar,dTmin_rhx,ihx_on,n_comp,comp_PR_uniform,comp_model_index,motor_eff,sink,source,dT_rhx1_cold,dT_rhx2_cold,dT_rhx3_cold,comp2_power_frac,comp3_power_frac,converged",
                "0.70,2.4,1000,1.0,Propane;Butane|0.5;0.5|mass,41,11,16,190,19,4,1,2,False,1,0.95,Air,water,5,40,0,0.45,0,True",
                "0.62,2.7,900,3.0,Propane;Butane|0.6;0.4|mass,42,12,17,195,18,3,0,2,False,1,0.96,Air,water,2,55,0,0.65,0,True",
            ]
        )
    )

    case_definition = results_io._load_case_definition_from_results_csv(csv_path)
    config = results_io._outer_design_config_from_case_definition(
        case_definition,
        line_budget=4,
    )

    assert case_definition.fluids == ("Propane", "Butane")
    assert case_definition.base_case["T_source_in"] == pytest.approx(42.0)
    assert case_definition.base_case["Tref_evap_inlet"] == pytest.approx(12.0)
    assert case_definition.base_case["T_air_tar"] == pytest.approx(195.0)
    assert case_definition.base_case["dT_GC_min_tar"] == pytest.approx(18.0)
    assert case_definition.base_case["dTmin_rhx"] == pytest.approx(3.0)
    assert case_definition.base_case["ihx_on"] is False
    assert case_definition.base_case["motor_eff"] == pytest.approx(0.96)
    assert case_definition.best_reference["SEC"] == pytest.approx(0.62)
    assert case_definition.best_reference["dT_rhx2_cold"] == pytest.approx(55.0)
    assert case_definition.best_reference["first_component_mass_fraction"] == pytest.approx(
        0.6
    )
    assert config.fixed_comp2_power_frac == pytest.approx(0.65)
    assert config.fixed_first_component_mass_fraction == pytest.approx(0.6)
    assert config.line_budget == 4


def test_write_outer_design_optimisation_results_includes_reference_metadata(
    monkeypatch,
    tmp_path,
):
    case_definition = opt_types.ResultsCsvCaseDefinition(
        base_case={"ihx_on": False},
        fluids=("Propane", "Butane"),
        selected_reference={"SEC": 0.7},
        best_reference={
            "SEC": 0.62,
            "dT_rhx1_cold": 2.5,
            "dT_rhx2_cold": 55.0,
            "comp2_power_frac": 0.65,
            "first_component_mass_fraction": 0.5,
        },
        source_csv="results.csv",
    )
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(dT1_target=0.0, dT2_target=0.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=3,
            split_budget=0,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    results_io._write_outer_design_optimisation_results(
        result,
        tmp_path,
        case_definition=case_definition,
    )

    summary = (tmp_path / "outer_design_optimisation_summary.json").read_text()
    evaluations = (tmp_path / "outer_design_optimisation_evaluations.csv").read_text()
    summary_csv = (tmp_path / "outer_design_optimisation_summary.csv").read_text()
    assert "reference_best_SEC" in summary
    assert "0.62" in summary
    assert "edge_a" in evaluations
    assert outer_opt.OPTIMISER_VERSION in summary_csv


def test_optimisation_case_generation_expands_fixed_variable_lists():
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={
            "fluids_ls": [("Propane", "Butane")],
            "T_air_tar_ls": [180.0, 200.0],
            "dT_rhx1_cold_ls": [0.0, 10.0],
            "dT_rhx2_cold_ls": [5.0],
            "comp2_power_frac_ls": [0.4, 0.6],
            "first_component_mass_fraction_ls": [0.3, 0.7],
        },
        config=outer_opt.OuterDesignOptimisationConfig(line_budget=2),
    )

    assert len(cases) == 16
    assert {case.base_case["T_air_tar"] for case in cases} == {180.0, 200.0}
    assert {case.fluids for case in cases} == {("Propane", "Butane")}
    assert {case.config.fixed_dT_rhx1_cold for case in cases} == {0.0, 10.0}
    assert {case.config.fixed_dT_rhx2_cold for case in cases} == {5.0}
    assert {case.config.fixed_comp2_power_frac for case in cases} == {0.4, 0.6}
    assert {case.config.fixed_first_component_mass_fraction for case in cases} == {
        0.3,
        0.7,
    }
    assert all(not case.config.vary_dT_rhx1_cold for case in cases)
    assert all(not case.config.vary_dT_rhx2_cold for case in cases)
    assert all(not case.config.vary_comp2_power_frac for case in cases)
    assert all(not case.config.vary_first_component_mass_fraction for case in cases)
    assert len({case.signature for case in cases}) == len(cases)


def test_optimisation_case_generation_omitted_variable_lists_keep_base_config():
    config = outer_opt.OuterDesignOptimisationConfig(
        fixed_comp2_power_frac=0.42,
        fixed_first_component_mass_fraction=0.64,
        vary_comp2_power_frac=True,
        vary_first_component_mass_fraction=True,
    )

    [case] = case_helpers._build_optimisation_cases(
        param_kwargs={"fluids": ("Propane", "Butane")},
        config=config,
    )

    assert case.config.fixed_comp2_power_frac == pytest.approx(0.42)
    assert case.config.fixed_first_component_mass_fraction == pytest.approx(0.64)
    assert case.config.vary_comp2_power_frac is True
    assert case.config.vary_first_component_mass_fraction is True


def test_optimisation_case_generation_uses_refrigerant_composition_as_start_only():
    [case] = case_helpers._build_optimisation_cases(
        param_kwargs={"refrigerant": "Propane;Butane|0.62;0.38|mass"},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )

    assert case.fluids == ("Propane", "Butane")
    assert case.config.fixed_first_component_mass_fraction == pytest.approx(0.62)
    assert case.config.vary_first_component_mass_fraction is True


def test_optimisation_case_filtering_skips_completed_signatures(tmp_path):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0, 200.0]},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )
    results_path = tmp_path / "RHPopt_results.csv"
    results_path.write_text(f"case_signature\n{cases[0].signature}\n")

    remaining = case_helpers._get_optimisation_case_list_to_solve(
        file_name="RHPopt_results",
        results_dir=tmp_path,
        param_kwargs={"T_air_tar_ls": [180.0, 200.0]},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )

    assert [case.signature for case in remaining] == [cases[1].signature]


def test_optimisation_case_filtering_retries_unsuccessful_result_rows(tmp_path):
    [case] = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0]},
    )
    results_path = tmp_path / "RHPopt_results.csv"
    results_path.write_text(
        f"case_signature,batch_success\n{case.signature},False\n"
    )

    remaining = case_helpers._get_optimisation_case_list_to_solve(
        file_name="RHPopt_results",
        results_dir=tmp_path,
        param_kwargs={"T_air_tar_ls": [180.0]},
    )

    assert remaining == [case]


def test_successful_retry_is_transferred_after_failed_master_row(tmp_path):
    [case] = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0]},
    )
    master_path = tmp_path / "RHPopt_results.csv"
    master_path.write_text(
        "case_signature,batch_success,best_objective\n"
        f"{case.signature},False,1e30\n"
    )
    batch_helpers.print_result_to_csv(
        output={
            "case_signature": case.signature,
            "batch_success": True,
            "best_objective": 0.5,
        },
        file_name="temp/optimisation/R",
        suffix=case.temp_suffix,
        results_dir=tmp_path,
    )

    batch_helpers._transfer_optimal_results(
        "RHPopt_results",
        results_dir=tmp_path,
    )

    rows = master_path.read_text().splitlines()
    assert len(rows) == 3
    assert rows[-1] == f"{case.signature},True,0.5"


def test_batch_runner_executes_cases_serially_and_writes_rows(monkeypatch, tmp_path):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0, 200.0]},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )
    calls = []
    written = []

    def fake_optimise_recuperative_heat_pump(**kwargs):
        calls.append(kwargs)
        return _fake_batch_result()

    def fake_print_result_to_csv(**kwargs):
        written.append(kwargs)

    monkeypatch.setattr(
        batch_helpers,
        "optimise_recuperative_heat_pump",
        fake_optimise_recuperative_heat_pump,
    )
    monkeypatch.setattr(batch_helpers, "print_result_to_csv", fake_print_result_to_csv)

    rows = batch_helpers._run_optimisation_cases(
        cases,
        results_dir=tmp_path,
        n_cpu=1,
    )

    assert len(rows) == 2
    assert [call["base_case"]["T_air_tar"] for call in calls] == [180.0, 200.0]
    assert [row["case_signature"] for row in rows] == [case.signature for case in cases]
    assert [write["suffix"] for write in written] == [
        case.temp_suffix for case in cases
    ]
    assert all(write["file_name"] == "temp/optimisation/R" for write in written)


def test_optimal_result_transfer_ignores_inner_candidate_files(tmp_path):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={
            "fluids": ("Propane", "Butane"),
            "T_air_tar_ls": [180.0, 200.0],
        },
    )
    optimal_row = {
        "case_signature": cases[0].signature,
        "best_objective": 0.5,
    }
    batch_helpers.print_result_to_csv(
        output=optimal_row,
        file_name="temp/optimisation/R",
        suffix=cases[0].temp_suffix,
        results_dir=tmp_path,
    )
    batch_helpers.print_result_to_csv(
        output={"case_signature": cases[1].signature, "best_objective": 0.6},
        file_name="temp/R",
        suffix=cases[1].temp_suffix,
        results_dir=tmp_path,
    )
    candidate_path = tmp_path / "temp" / "R-Propane0_5-Butane0_5.csv"
    candidate_path.parent.mkdir(parents=True, exist_ok=True)
    candidate_path.write_text("SEC,refrigerant\n0.6,Propane;Butane\n")

    batch_helpers._transfer_optimal_results(
        "RHPopt_results",
        results_dir=tmp_path,
    )

    master = (tmp_path / "RHPopt_results.csv").read_text()
    assert cases[0].signature in master
    assert cases[1].signature in master
    assert "refrigerant" not in master
    assert not (tmp_path / "temp" / "optimisation").exists()
    assert candidate_path.exists()


def test_single_case_runner_resumes_optimal_result_before_solving(monkeypatch, tmp_path):
    config = outer_opt.OuterDesignOptimisationConfig(max_exact_evaluations=3)
    case = case_helpers._build_optimisation_case(
        base_case={"T_air_tar": 180.0},
        fluids=("Propane", "Butane"),
        config=config,
    )
    batch_helpers.print_result_to_csv(
        output={"case_signature": case.signature, "best_objective": 0.5},
        file_name="temp/optimisation/R",
        suffix=case.temp_suffix,
        results_dir=tmp_path,
    )
    monkeypatch.setattr(
        batch_helpers,
        "optimise_recuperative_heat_pump",
        lambda **kwargs: pytest.fail("completed case should not be optimised again"),
    )

    row = outer_opt.optimise_recuperative_heat_pump_case(
        filename="RHPopt_results",
        base_case={"T_air_tar": 180.0},
        fluids=("Propane", "Butane"),
        config=config,
        results_dir=tmp_path,
    )

    assert row is not None
    assert row["case_signature"] == case.signature
    assert (tmp_path / "RHPopt_results.csv").exists()


def test_batch_runner_submits_parallel_cases_to_executor(tmp_path):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0, 200.0]},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )
    executors = []

    class FakeFuture:
        def __init__(self, case):
            self.case = case

        def result(self):
            return {"case_signature": self.case.signature}

    class FakeExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers
            self.submitted = []
            self.shutdown_called = False
            executors.append(self)

        def submit(self, fn, *args, **kwargs):
            self.submitted.append((fn, args, kwargs))
            return FakeFuture(args[0])

        def shutdown(self):
            self.shutdown_called = True

    rows = batch_helpers._run_optimisation_cases(
        cases,
        results_dir=tmp_path,
        n_cpu=2,
        executor_factory=FakeExecutor,
    )

    assert rows == [{"case_signature": case.signature} for case in cases]
    assert len(executors) == 1
    assert executors[0].max_workers == 2
    assert [submission[1][0] for submission in executors[0].submitted] == cases
    assert executors[0].shutdown_called is True


def test_batch_runner_surfaces_parallel_worker_failures(tmp_path):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0]},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )

    class FailingFuture:
        def result(self):
            raise RuntimeError("worker failed")

    class FakeExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers

        def submit(self, fn, *args, **kwargs):
            return FailingFuture()

        def shutdown(self):
            pass

    with pytest.raises(RuntimeError, match="worker failed"):
        batch_helpers._run_optimisation_cases(
            cases,
            results_dir=tmp_path,
            n_cpu=2,
            executor_factory=FakeExecutor,
        )


def test_public_batch_runner_transfers_optimal_files_before_and_after(
    monkeypatch,
    tmp_path,
):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0]},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )
    transfers = []
    events = []

    monkeypatch.setattr(
        batch_helpers,
        "_transfer_optimal_results",
        lambda filename, results_dir: (
            transfers.append((filename, results_dir)),
            events.append("transfer"),
        ),
    )

    def get_cases(**kwargs):
        events.append("filter")
        return cases

    monkeypatch.setattr(
        batch_helpers,
        "_get_optimisation_case_list_to_solve",
        get_cases,
    )

    def run_cases(cases, results_dir, n_cpu):
        events.append("optimise")
        return [{"case_signature": cases[0].signature}]

    monkeypatch.setattr(
        batch_helpers,
        "_run_optimisation_cases",
        run_cases,
    )

    rows = outer_opt.optimise_recuperative_heat_pump_cases(
        filename="RHPopt_results",
        param_kwargs={"T_air_tar_ls": [180.0]},
        config=outer_opt.OuterDesignOptimisationConfig(),
        n_cpu=3,
        results_dir=tmp_path,
    )

    assert rows == [{"case_signature": cases[0].signature}]
    assert transfers == [
        ("RHPopt_results", tmp_path),
        ("RHPopt_results", tmp_path),
    ]
    assert events == ["transfer", "filter", "optimise", "transfer"]


def test_real_outer_design_optimisation_smoke_when_refprop_is_available(hpl, tmp_path):
    del hpl
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=4,
        line_budget=1,
        split_budget=0,
        polish_budget=0,
    )

    result = outer_opt.optimise_recuperative_heat_pump(config=config, results_dir=tmp_path)

    assert result.exact_evaluations <= 4
    assert result.best is not None
    assert math.isfinite(result.best.objective)
    assert result.best.message


def test_real_local_pattern_search_polish_smoke_when_refprop_is_available(hpl, tmp_path):
    del hpl
    objective = objective_helpers._OuterDesignObjective(
        max_evals=2,
        results_dir=tmp_path,
    )
    source = _phase_from_exact_record(
        objective,
        [0.0, 0.0, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(budget=1),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.method == "pattern_search"
    assert result.best is not None
    assert math.isfinite(result.best.objective)


def _optimise_synthetic(
    *,
    analyser: SyntheticAnalyser,
    config: outer_opt.OuterDesignOptimisationConfig,
    monkeypatch,
):
    def build_synthetic_objective(**kwargs):
        return objective_helpers._OuterDesignObjective(
            **kwargs,
            analyser=analyser,
        )

    monkeypatch.setattr(
        outer_design_module,
        "_OuterDesignObjective",
        build_synthetic_objective,
    )
    return outer_opt.optimise_recuperative_heat_pump(
        config=config,
    )


def _first_mass_fraction(refrigerant: str) -> float:
    composition = refrigerant.split("|")[1]
    return float(composition.split(";")[0])


def _phase_from_exact_record(objective, x, *, phase):
    record = objective.evaluate(x)
    return outer_opt.PhaseEvaluation.from_record(
        phase=phase,
        record=record,
        cumulative_exact_evaluations=objective.exact_evaluations,
    )


def _phase_evaluator(objective, records):
    def evaluate_phase(phase, x):
        record = objective.evaluate(x)
        phase_record = outer_opt.PhaseEvaluation.from_record(
            phase=phase,
            record=record,
            cumulative_exact_evaluations=objective.exact_evaluations,
        )
        records.append(phase_record)
        return phase_record

    return evaluate_phase


def _fake_batch_result():
    best = SimpleNamespace(
        success=True,
        message="success",
        sec=0.5,
        cop=2.0,
        pmax=1234.0,
    )
    return SimpleNamespace(
        best=best,
        to_summary_dict=lambda: {
            "optimiser_version": outer_opt.OPTIMISER_VERSION,
            "best_objective": 0.5,
            "best_phase": "corner",
            "exact_evaluations": 1,
            "failed_evaluations": 0,
            "best_dT_rhx1_cold": 0.0,
            "best_dT_rhx2_cold": 0.0,
            "best_comp2_power_frac": 0.5,
            "best_first_component_mass_fraction": 0.5,
        },
    )
