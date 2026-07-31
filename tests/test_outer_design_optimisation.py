from __future__ import annotations

from dataclasses import replace
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
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    _composition_starts as composition_start_helpers,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    _isolated_evaluation as isolated_helpers,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    _multistart as multistart_helpers,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import outer_design as outer_design_module
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _polish as polish_helpers
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    _results_io as results_io,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    _search as search_helpers,
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
        "CompositionSearchResult",
        "DesignEvaluation",
        "HistoricalIncumbent",
        "IncumbentReplay",
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
    assert "fluids" not in signature.parameters
    assert "refrigerant" in signature.parameters
    assert "n_cpu" in signature.parameters
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

    record = objective.evaluate([20.0, 35.0, 0.45, 0.72, 0.28])

    assert record.success is True
    assert {hp_cls.__name__ for hp_cls in analyser.hp_classes} == {"RecuperativeHP"}


def test_objective_formats_all_composition_coordinates_in_mass_basis():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        refrigerant=["Propane", "Butane", "Pentane"],
        analyser=SyntheticAnalyser(),
    )

    case = objective.vector_to_case([0.0, 0.0, 0.5, 0.2, 0.3, 0.5])

    assert case["refrigerant"] == "Propane;Butane;Pentane|0.2;0.3;0.5|mass"


def test_objective_maps_vector_to_domain_case_without_mutating_base_case():
    base_case = {"T_air_tar": 180.0, "dT_rhx3_cold": 0.0}
    objective = objective_helpers._OuterDesignObjective(
        max_evals=5,
        base_case=base_case,
        refrigerant="Propane;Butane|0.5;0.5|mass",
        analyser=SyntheticAnalyser(),
    )

    case = objective.vector_to_case([12.0, 34.0, 0.42, 0.8, 0.2])

    assert case["T_air_tar"] == pytest.approx(180.0)
    assert case["dT_rhx1_cold"] == pytest.approx(12.0)
    assert case["dT_rhx2_cold"] == pytest.approx(34.0)
    assert case["comp2_power_frac"] == pytest.approx(0.42)
    assert case["refrigerant"] == "Propane;Butane|0.8;0.2|mass"
    assert base_case == {"T_air_tar": 180.0, "dT_rhx3_cold": 0.0}
    assert objective.base_case["T_source_in"] == pytest.approx(
        opt_types.DEFAULT_BASE_CASE["T_source_in"]
    )
    assert objective.base_case["T_air_tar"] == pytest.approx(180.0)


def test_objective_projects_ternary_composition_onto_simplex():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        refrigerant=["Propane", "Butane", "Pentane"],
        analyser=SyntheticAnalyser(),
    )

    repaired = objective.repair_vector([0.0, 0.0, 0.5, -0.2, 0.4, 1.8])

    assert all(value >= 0.0 for value in repaired[3:])
    assert sum(repaired[3:]) == pytest.approx(1.0)


def test_binary_composition_multistart_covers_boundaries_and_quarter_points():
    starts = composition_start_helpers._generate_composition_starts(
        (0.5, 0.5),
        vary_composition=True,
        max_starts=8,
    )

    assert starts == (
        (0.5, 0.5),
        (0.0, 1.0),
        (1.0, 0.0),
        (0.25, 0.75),
        (0.75, 0.25),
    )


def test_ternary_composition_multistart_covers_vertices_and_binary_edges():
    starts = composition_start_helpers._generate_composition_starts(
        (1 / 3, 1 / 3, 1 / 3),
        vary_composition=True,
        max_starts=8,
    )

    expected = {
        (1 / 3, 1 / 3, 1 / 3),
        (1.0, 0.0, 0.0),
        (0.0, 1.0, 0.0),
        (0.0, 0.0, 1.0),
        (0.5, 0.5, 0.0),
        (0.5, 0.0, 0.5),
        (0.0, 0.5, 0.5),
    }
    assert len(starts) == 7
    assert set(starts) == expected
    assert all(sum(start) == pytest.approx(1.0) for start in starts)


@pytest.mark.parametrize(
    ("initial", "vary_composition"),
    [((1.0,), True), ((0.2, 0.3, 0.5), False)],
)
def test_composition_multistart_uses_only_input_when_composition_cannot_vary(
    initial,
    vary_composition,
):
    assert composition_start_helpers._generate_composition_starts(
        initial,
        vary_composition=vary_composition,
        max_starts=8,
    ) == (initial,)


def test_composition_multistart_is_deterministically_capped():
    initial = (0.1, 0.2, 0.3, 0.4)

    first = composition_start_helpers._generate_composition_starts(
        initial,
        vary_composition=True,
        max_starts=5,
    )
    second = composition_start_helpers._generate_composition_starts(
        initial,
        vary_composition=True,
        max_starts=5,
    )

    assert first == second
    assert first[0] == initial
    assert len(first) == 5
    assert all(min(start) >= 0.0 and sum(start) == pytest.approx(1.0) for start in first)


def test_outer_design_config_defaults_to_robust_multistart_budget():
    config = outer_opt.OuterDesignOptimisationConfig()

    assert config.max_exact_evaluations == 250
    assert config.max_composition_starts == 8
    assert config.polish_top_n == 3
    assert config.polish_budget == 24
    assert config.polish_pure_recuperator_lifts is True
    assert config.preserve_incumbent is True
    assert config.incumbent_replay_limit == 3
    assert config.incumbent_polish_budget == 12
    assert config.exact_evaluation_timeout_seconds == pytest.approx(180.0)
    assert config.exact_evaluation_timeout_retries == 1
    assert config.exact_evaluation_worker_max_tasks == 20


def test_outer_design_config_rejects_invalid_exact_evaluation_limits():
    with pytest.raises(ValueError, match="timeout_seconds must be positive"):
        outer_opt.OuterDesignOptimisationConfig(
            exact_evaluation_timeout_seconds=0.0
        )
    with pytest.raises(ValueError, match="worker_max_tasks must be at least 1"):
        outer_opt.OuterDesignOptimisationConfig(
            exact_evaluation_worker_max_tasks=0
        )
    with pytest.raises(ValueError, match="timeout_retries must be nonnegative"):
        outer_opt.OuterDesignOptimisationConfig(
            exact_evaluation_timeout_retries=-1
        )


def test_v10_never_selects_worse_than_successful_current_solver_incumbent(
    monkeypatch,
):
    analyser = SyntheticAnalyser(
        dT1_target=55.0,
        dT2_target=25.0,
        comp2_target=0.3,
        composition_target=0.8,
    )

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
    incumbent = outer_opt.HistoricalIncumbent(
        source_version="edge-first-pattern-search-v9",
        source_file="prior.csv",
        stored_sec=0.99,
        x=(55.0, 25.0, 0.3, 0.8, 0.2),
    )
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=3,
        max_composition_starts=1,
        line_budget=1,
        split_budget=0,
        polish_budget=0,
        incumbent_polish_budget=0,
    )

    result = outer_opt.optimise_recuperative_heat_pump(
        config=config,
        incumbents=(incumbent,),
    )
    summary = result.to_summary_dict()

    assert len(result.incumbent_replays) == 1
    assert result.incumbent_replays[0].evaluation.success is True
    assert result.best is not None
    assert result.best.objective <= result.incumbent_replays[0].evaluation.objective
    assert result.best.phase == "incumbent_replay"
    assert summary["incumbent_source_version"] == incumbent.source_version
    assert summary["incumbent_stored_SEC"] == pytest.approx(0.99)
    assert summary["incumbent_replayed_SEC"] == pytest.approx(1.0)
    assert summary["selected_incumbent"] is True
    assert summary["no_worse_than_replayed_incumbent"] is True
    assert summary["SEC_change_vs_incumbent_pct"] <= 0.0
    assert summary["historical_model_drift_pct"] == pytest.approx(
        100.0 * (1.0 / 0.99 - 1.0)
    )


def test_historical_incumbents_are_attached_and_included_in_case_signature(
    tmp_path,
):
    [case] = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0]},
    )
    row = batch_helpers._result_row(case, _fake_batch_result())
    row.update(
        {
            "optimiser_version": "edge-first-pattern-search-v9",
            "best_SEC": 0.48,
            "best_dT_rhx1_cold": 12.0,
            "best_dT_rhx2_cold": 44.0,
            "best_comp2_power_frac": 0.4,
            "best_mass_fractions": "0.6;0.4",
        }
    )
    batch_helpers.print_result_to_csv(
        output=row,
        file_name="historical",
        results_dir=tmp_path,
    )

    [attached] = case_helpers._attach_incumbents_and_refresh_signatures(
        [case],
        results_dir=tmp_path,
    )

    assert len(attached.incumbents) == 1
    incumbent = attached.incumbents[0]
    assert incumbent.source_version == "edge-first-pattern-search-v9"
    assert incumbent.source_file == "historical.csv"
    assert incumbent.stored_sec == pytest.approx(0.48)
    assert incumbent.x == pytest.approx((12.0, 44.0, 0.4, 0.6, 0.4))
    assert attached.signature != case.signature
    assert attached.temp_suffix == f"-{attached.signature}"


def test_historical_incumbent_requires_all_extra_case_inputs_to_match(tmp_path):
    case = case_helpers._build_optimisation_case(
        base_case={"T_air_tar": 180.0},
        refrigerant="Propane;Butane|0.5;0.5|mass",
    )
    row = batch_helpers._result_row(case, _fake_batch_result())
    row["optimiser_version"] = "edge-first-pattern-search-v9"
    batch_helpers.print_result_to_csv(
        output=row,
        file_name="historical",
        results_dir=tmp_path,
    )
    case_with_extra_input = case_helpers._build_optimisation_case(
        base_case={"T_air_tar": 180.0, "experimental_solver_setting": 2.0},
        refrigerant="Propane;Butane|0.5;0.5|mass",
    )

    [attached] = case_helpers._attach_incumbents_and_refresh_signatures(
        [case_with_extra_input],
        results_dir=tmp_path,
    )

    assert attached.incumbents == ()


@pytest.mark.parametrize("polish_top_n", [0, 1, 4])
def test_outer_design_config_requires_two_or_three_polish_sources(polish_top_n):
    with pytest.raises(ValueError, match="polish_top_n"):
        outer_opt.OuterDesignOptimisationConfig(polish_top_n=polish_top_n)


@pytest.mark.parametrize(
    ("max_evaluations", "expected_starts"),
    [(100, 2), (250, 7)],
)
def test_composition_multistart_plan_reserves_polish_budget(
    max_evaluations,
    expected_starts,
):
    starts = tuple((index / 10, 1 - index / 10) for index in range(8))
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=max_evaluations,
    )

    selected, budget_per_start = composition_start_helpers._plan_composition_searches(
        starts,
        config,
    )

    assert len(selected) == expected_starts
    expected_budget = 26 if max_evaluations == 100 else 25
    assert budget_per_start == expected_budget
    assert (
        len(selected) * budget_per_start
        + min(3, len(selected)) * config.polish_budget
        <= max_evaluations
    )


def test_outer_optimiser_runs_line_searches_from_each_selected_composition(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(composition_target=0.1),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=3,
            split_budget=3,
            polish_budget=0,
        ),
        monkeypatch=monkeypatch,
    )

    assert result.requested_composition_starts == 5
    assert len(result.composition_searches) == 5
    assert result.corner.x[3:] == pytest.approx((0.5, 0.5))
    assert {search.composition for search in result.composition_searches} == {
        (0.5, 0.5),
        (0.0, 1.0),
        (1.0, 0.0),
        (0.25, 0.75),
        (0.75, 0.25),
    }


def test_multistart_finds_better_remote_composition_than_single_start(monkeypatch):
    common = {
        "line_budget": 3,
        "split_budget": 0,
        "polish_budget": 0,
    }
    single = _optimise_synthetic(
        analyser=SyntheticAnalyser(composition_target=0.0),
        config=outer_opt.OuterDesignOptimisationConfig(
            max_composition_starts=1,
            **common,
        ),
        monkeypatch=monkeypatch,
    )
    multiple = _optimise_synthetic(
        analyser=SyntheticAnalyser(composition_target=0.0),
        config=outer_opt.OuterDesignOptimisationConfig(**common),
        monkeypatch=monkeypatch,
    )

    assert single.best is not None
    assert multiple.best is not None
    assert multiple.best.objective < single.best.objective
    assert multiple.best.x[3] == pytest.approx(0.0)


def test_polish_sources_are_best_candidates_from_distinct_line_directions(monkeypatch):
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(composition_target=0.3),
        config=outer_opt.OuterDesignOptimisationConfig(
            line_budget=4,
            split_budget=4,
            polish_budget=2,
            polish_top_n=3,
        ),
        monkeypatch=monkeypatch,
    )

    direction_candidates = {
        candidate.x
        for search in result.composition_searches
        for candidate in (
            search.edge_a.best,
            search.edge_b.best,
            None if search.split is None else search.split.best,
        )
        if candidate is not None and candidate.success
    }
    sources = [polish.source_x for polish in result.polish_results]
    assert 1 <= len(sources) <= 3
    assert len(set(sources)) == len(sources)
    assert set(sources) <= direction_candidates


def test_multistart_never_exceeds_case_wide_exact_evaluation_budget(monkeypatch):
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=40,
        line_budget=5,
        split_budget=5,
        polish_budget=6,
    )
    result = _optimise_synthetic(
        analyser=SyntheticAnalyser(),
        config=config,
        monkeypatch=monkeypatch,
    )

    assert result.exact_evaluations <= config.max_exact_evaluations
    assert result.to_summary_dict()["attempted_composition_starts"] == len(
        result.composition_searches
    )


def test_parallel_composition_search_preserves_start_order_and_worker_failures():
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=2,
        line_budget=0,
        split_budget=0,
        polish_budget=0,
    )
    tasks = multistart_helpers._build_composition_search_tasks(
        ((0.5, 0.5), (0.0, 1.0)),
        budget_per_start=1,
        base_case=opt_types.DEFAULT_BASE_CASE,
        refrigerant="Propane;Butane|0.5;0.5|mass",
        config=config,
        results_dir="results",
    )
    executors = []
    analyser = SyntheticAnalyser()

    def build_synthetic_objective(**kwargs):
        return objective_helpers._OuterDesignObjective(
            **kwargs,
            analyser=analyser,
            molar_masses=(44.0, 58.0),
        )

    successful = multistart_helpers._run_composition_search(
        tasks[0],
        objective_factory=build_synthetic_objective,
    )

    class FakeFuture:
        def __init__(self, task):
            self.task = task

        def result(self):
            if self.task.start_index == 0:
                return successful
            raise RuntimeError(f"worker {self.task.start_index} failed")

    class FakeExecutor:
        def __init__(self, max_workers):
            self.max_workers = max_workers
            self.submitted = []
            self.shutdown_called = False
            executors.append(self)

        def submit(self, fn, task):
            self.submitted.append((fn, task))
            return FakeFuture(task)

        def shutdown(self):
            self.shutdown_called = True

    results = multistart_helpers._run_composition_searches(
        tasks,
        n_cpu=4,
        executor_factory=FakeExecutor,
    )

    assert [result.start_index for result in results] == [0, 1]
    assert [result.success for result in results] == [True, False]
    assert executors[0].max_workers == 2
    assert executors[0].shutdown_called is True


def test_spare_binary_search_budget_strengthens_line_and_split_coverage():
    config = outer_opt.OuterDesignOptimisationConfig()

    expanded = multistart_helpers._search_config_for_budget(
        config,
        max_evaluations=35,
    )

    assert expanded.line_budget == 12
    assert expanded.split_budget == 10
    assert composition_start_helpers._full_search_budget(expanded) == 35


def test_composition_search_falls_back_to_serial_when_pool_setup_fails():
    analyser = SyntheticAnalyser()
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=2,
        line_budget=0,
        split_budget=0,
        polish_budget=0,
    )
    tasks = multistart_helpers._build_composition_search_tasks(
        ((0.5, 0.5), (0.0, 1.0)),
        budget_per_start=1,
        base_case=opt_types.DEFAULT_BASE_CASE,
        refrigerant="Propane;Butane|0.5;0.5|mass",
        config=config,
        results_dir="results",
    )

    def build_synthetic_objective(**kwargs):
        return objective_helpers._OuterDesignObjective(
            **kwargs,
            analyser=analyser,
            molar_masses=(44.0, 58.0),
        )

    def unavailable_executor(*, max_workers):
        del max_workers
        raise OSError("processes unavailable")

    results = multistart_helpers._run_composition_searches(
        tasks,
        n_cpu=2,
        executor_factory=unavailable_executor,
        objective_factory=build_synthetic_objective,
    )

    assert all(result.success for result in results)
    assert len(analyser.calls) == 2


def test_mole_basis_objective_solves_in_mass_basis_and_reports_in_mole_basis():
    analyser = SyntheticAnalyser()
    objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        refrigerant="Propane;Butane|0.4;0.6|mol",
        analyser=analyser,
        molar_masses=(44.0, 58.0),
    )
    x = [0.0, 0.0, 0.5, *objective.initial_mass_fractions]

    record = objective.evaluate(x)

    assert analyser.calls[0]["refrigerant"].endswith("|mass")
    assert _composition(analyser.calls[0]["refrigerant"]) == pytest.approx(
        objective.initial_mass_fractions
    )
    assert record.refrigerant.endswith("|mol")
    assert _composition(record.refrigerant) == pytest.approx((0.4, 0.6))


def test_equivalent_mass_and_mole_inputs_share_internal_starting_composition():
    masses = (44.0, 58.0)
    mole_objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        refrigerant="Propane;Butane|0.4;0.6|mol",
        analyser=SyntheticAnalyser(),
        molar_masses=masses,
    )
    mass_objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        refrigerant="Propane;Butane|0.335877862595;0.664122137405|mass",
        analyser=SyntheticAnalyser(),
    )

    assert mass_objective.initial_mass_fractions == pytest.approx(
        mole_objective.initial_mass_fractions
    )


def test_objective_repairs_recuperator_lift_to_feasible_domain():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=5,
        base_case={"T_air_inlet": 10.0, "T_air_tar": 100.0, "dT_rhx3_cold": 0.0},
        refrigerant="Propane;Butane|0.5;0.5|mass",
        analyser=SyntheticAnalyser(),
    )

    case = objective.vector_to_case([60.0, 60.0, 0.42, 0.8, 0.2])

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
        refrigerant="Propane;Butane|0.5;0.5|mass",
        analyser=analyser,
    )

    record = objective.evaluate([80.0, 80.0, 0.42, 0.8, 0.2])

    assert record.success is True
    assert record.x[0] + record.x[1] == pytest.approx(60.0)
    assert analyser.calls[0]["dT_rhx1_cold"] + analyser.calls[0]["dT_rhx2_cold"] <= 60.0


def test_objective_reuses_cached_exact_solution_for_duplicate_vectors():
    analyser = SyntheticAnalyser()
    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=analyser,
    )

    first = objective.evaluate([20.0, 35.0, 0.45, 0.72, 0.28])
    second = objective.evaluate([20.0, 35.0, 0.45, 0.72, 0.28])

    assert first.objective == pytest.approx(second.objective)
    assert objective.exact_evaluations == 1
    assert len(analyser.calls) == 1


def test_objective_does_not_carry_inner_solver_state_between_designs():
    calls = []

    def analyser(hp_cls, case, results_dir, use_fast_objective=True):
        del hp_cls, results_dir, use_fast_objective
        calls.append(dict(case))
        return FakeHeatPump(sec=1.0, case=case), 0.01

    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=analyser,
        refrigerant="Pentane",
    )

    objective.evaluate([20.0, 35.0, 0.45, 1.0])
    objective.evaluate([21.0, 35.0, 0.45, 1.0])

    assert "P_H_est" not in calls[0]
    assert "Tref_evap_outlet" not in calls[0]
    assert "P_H_est" not in calls[1]
    assert "Tref_evap_outlet" not in calls[1]


def test_objective_returns_finite_penalty_for_failed_exact_solve():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=SyntheticAnalyser(fail=True),
    )

    record = objective.evaluate([20.0, 35.0, 0.45, 0.72, 0.28])

    assert record.success is False
    assert math.isfinite(record.objective)
    assert record.objective == pytest.approx(opt_types.PENALTY_OBJECTIVE + 1.0)
    assert record.constraint_errors == {"solve_failure": pytest.approx(1.0)}
    assert "synthetic solve failed" in record.message
    assert objective.failed_evaluations == 1


def test_objective_records_isolated_exact_evaluation_timeout(monkeypatch):
    class TimingOutEvaluator:
        def __init__(self, **kwargs):
            self.kwargs = kwargs

        def evaluate(self, request):
            raise isolated_helpers.ExactEvaluationTimeout(
                "exact thermodynamic evaluation timed out after 0.01 seconds"
            )

        def close(self):
            pass

    monkeypatch.setattr(
        objective_helpers,
        "IsolatedExactEvaluator",
        TimingOutEvaluator,
    )
    objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        refrigerant="Propane;CO2|0.5;0.5|mass",
        exact_evaluation_timeout_seconds=0.01,
    )

    record = objective.evaluate([0.0, 0.0, 0.5, 0.5, 0.5])

    assert record.success is False
    assert record.objective == pytest.approx(opt_types.PENALTY_OBJECTIVE + 1.0)
    assert record.constraint_errors == {"evaluation_timeout": 1.0}
    assert "timed out after 0.01 seconds" in record.message
    assert objective.failed_evaluations == 1


def test_isolated_evaluator_replaces_worker_after_timeout(monkeypatch):
    class FakeConnection:
        def __init__(self, response):
            self.response = response
            self.sent = []
            self.closed = False

        def send(self, payload):
            self.sent.append(payload)

        def poll(self, timeout):
            del timeout
            return self.response is not None

        def recv(self):
            return self.response

        def close(self):
            self.closed = True

    class FakeProcess:
        def __init__(self):
            self.alive = False
            self.terminated = False
            self.closed = False

        def start(self):
            self.alive = True

        def is_alive(self):
            return self.alive

        def join(self, timeout):
            del timeout

        def terminate(self):
            self.terminated = True
            self.alive = False

        def kill(self):
            self.alive = False

        def close(self):
            self.closed = True

    class FakeContext:
        def __init__(self):
            self.responses = [None, ("ok", "recovered")]
            self.processes = []

        def Pipe(self, duplex):
            assert duplex is True
            return FakeConnection(self.responses.pop(0)), FakeConnection(None)

        def Process(self, **kwargs):
            assert kwargs["target"] is isolated_helpers._exact_evaluation_worker
            process = FakeProcess()
            self.processes.append(process)
            return process

    context = FakeContext()
    monkeypatch.setattr(isolated_helpers, "get_context", lambda name: context)
    evaluator = isolated_helpers.IsolatedExactEvaluator(
        timeout_seconds=0.01,
        max_tasks=20,
    )

    with pytest.raises(isolated_helpers.ExactEvaluationTimeout):
        evaluator.evaluate({"candidate": 1})
    result = evaluator.evaluate({"candidate": 2})

    assert result == "recovered"
    assert len(context.processes) == 2
    assert context.processes[0].terminated is True
    evaluator.close()


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
    ).evaluate([20.0, 35.0, 0.45, 0.72, 0.28])
    far = objective_helpers._OuterDesignObjective(
        max_evals=1,
        analyser=failed_analyser(3.0),
    ).evaluate([20.0, 35.0, 0.45, 0.72, 0.28])

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

    first = objective.evaluate([20.0, 35.0, 0.45, 0.72, 0.28])
    second = objective.evaluate([30.0, 35.0, 0.45, 0.72, 0.28])

    assert first.success is True
    assert second.success is False
    assert second.message == "exact evaluation budget exhausted"
    assert objective.exact_evaluations == 1
    assert objective.budget_exhausted_calls == 1


def test_line_search_refines_the_best_observed_sample_bracket():
    sampled_values = []

    def evaluate_phase(phase, x):
        value = float(x[0])
        sampled_values.append(value)
        return SimpleNamespace(
            phase=phase,
            x=(value,),
            objective=(value - 75.0) ** 2,
            success=True,
        )

    result = search_helpers._line_search(
        phase="test",
        variable_name="x",
        lo=0.0,
        hi=100.0,
        budget=8,
        build_x=lambda value: [value],
        evaluate_phase=evaluate_phase,
    )

    assert sampled_values[:5] == pytest.approx((0.0, 25.0, 50.0, 75.0, 100.0))
    assert all(50.0 <= value <= 100.0 for value in sampled_values[5:])
    assert result.best is not None
    assert result.best.x[0] == pytest.approx(75.0)


def test_line_search_uses_observed_multibasin_values_to_choose_refinement_region():
    def evaluate_phase(phase, x):
        value = float(x[0])
        broad_well = 0.2 + (value - 20.0) ** 2 / 4000.0
        narrow_well = (value - 75.0) ** 2 / 20.0
        return SimpleNamespace(
            phase=phase,
            x=(value,),
            objective=min(broad_well, narrow_well),
            success=True,
        )

    result = search_helpers._line_search(
        phase="test",
        variable_name="x",
        lo=0.0,
        hi=100.0,
        budget=8,
        build_x=lambda value: [value],
        evaluate_phase=evaluate_phase,
    )

    assert result.best is not None
    assert result.best.x[0] == pytest.approx(75.0)
    assert result.best.objective == pytest.approx(0.0)


def test_line_search_refines_distinct_interior_basin_when_screen_best_is_boundary():
    sampled_values = []

    def evaluate_phase(phase, x):
        value = float(x[0])
        sampled_values.append(value)
        objective = (
            0.1 + (value - 56.0) ** 2 / 200.0
            if 35.0 <= value <= 70.0
            else 0.2 + value / 100.0
        )
        return SimpleNamespace(
            phase=phase,
            x=(value,),
            objective=objective,
            success=True,
        )

    result = search_helpers._line_search(
        phase="test",
        variable_name="x",
        lo=0.0,
        hi=100.0,
        budget=9,
        build_x=lambda value: [value],
        evaluate_phase=evaluate_phase,
    )

    assert sampled_values[:5] == pytest.approx((0.0, 25.0, 50.0, 75.0, 100.0))
    assert any(25.0 <= value <= 75.0 for value in sampled_values[5:])
    assert result.best is not None
    assert result.best.objective < 0.2
    assert result.best.x[0] == pytest.approx(56.0, abs=1.0)


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
            line_budget=10,
            split_budget=0,
            polish_budget=16,
            polish_top_n=2,
        ),
        monkeypatch=monkeypatch,
    )

    assert len(result.polish_results) == 2
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
        [0.0, 40.0, 0.5, 0.5, 0.5],
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
    assert all(0.0 <= fraction <= 1.0 for fraction in result.best.x[3:])
    assert sum(result.best.x[3:]) == pytest.approx(1.0)


def test_pattern_search_screens_remote_compressor_split_basin():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=6,
        analyser=SyntheticAnalyser(comp2_target=0.1),
    )
    source = _phase_from_exact_record(
        objective,
        [20.0, 35.0, 0.5, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=5,
            comp2_screen_points=5,
            vary_composition=False,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.best is not None
    assert result.best.x[2] == pytest.approx(0.1)
    assert result.best.objective < source.objective


def test_pattern_search_can_reopen_recuperator_lifts():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=8,
        analyser=SyntheticAnalyser(
            dT1_target=0.0,
            dT2_target=5.0,
            comp2_target=0.1,
            composition_target=0.5,
        ),
    )
    source = _phase_from_exact_record(
        objective,
        [0.0, 0.0, 0.5, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=7,
            vary_dT_rhx2_cold=True,
            vary_composition=False,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.best is not None
    assert result.best.x[1] == pytest.approx(5.0)
    assert result.best.x[2] == pytest.approx(0.1)
    assert result.best.objective < source.objective


def test_mixture_lift_feedback_relines_after_composition_polish():
    analyser = SyntheticAnalyser(
        dT1_target=85.0,
        dT2_target=0.0,
        comp2_target=0.9,
        composition_target=0.85,
    )
    source_objective = objective_helpers._OuterDesignObjective(
        max_evals=1,
        analyser=analyser,
    )
    source = _phase_from_exact_record(
        source_objective,
        [68.0, 0.0, 0.9, 0.85, 0.15],
        phase="polish_pattern_search",
    )
    task = polish_helpers._LiftFeedbackTask(
        source=source,
        budget=7,
        base_case=dict(opt_types.DEFAULT_BASE_CASE),
        refrigerant="Propane;Butane|0.5;0.5|mass",
        config=outer_opt.OuterDesignOptimisationConfig(),
        results_dir="results",
    )

    result = polish_helpers._run_lift_feedback(
        task,
        objective_factory=lambda **kwargs: objective_helpers._OuterDesignObjective(
            **kwargs,
            analyser=analyser,
        ),
    )

    assert result.success is True
    assert result.best is not None
    assert result.best.x[0] == pytest.approx(85.0, abs=3.0)
    assert result.best.objective < source.objective
    assert all(record.phase == "lift_feedback" for record in result.evaluations)


def test_pattern_search_transfers_mass_across_ternary_simplex():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=16,
        refrigerant="Propane;Butane;Pentane|0.2;0.3;0.5|mass",
        analyser=SyntheticAnalyser(composition_target=0.4),
    )
    source = _phase_from_exact_record(
        objective,
        [20.0, 35.0, 0.45, 0.2, 0.3, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=12,
            vary_comp2_power_frac=False,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.best is not None
    assert result.best.objective < source.objective
    assert all(sum(record.x[3:]) == pytest.approx(1.0) for record in result.evaluations)
    assert all(value >= 0.0 for record in result.evaluations for value in record.x[3:])


def test_binary_composition_poll_retains_positive_first_fraction_direction():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=2,
        analyser=SyntheticAnalyser(composition_target=0.6),
    )
    source = _phase_from_exact_record(
        objective,
        [20.0, 35.0, 0.45, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=1,
            vary_comp2_power_frac=False,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.best is not None
    assert result.best.x[3:] == pytest.approx((0.55, 0.45))


def test_pattern_search_uses_cache_for_duplicate_current_point():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=3,
        analyser=SyntheticAnalyser(),
    )
    source = _phase_from_exact_record(
        objective,
        [20.0, 35.0, 0.45, 0.72, 0.28],
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
        [20.0, 35.0, 0.5, 0.5, 0.5],
        phase="source",
    )
    phase_records = []

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            comp2_screen_points=0,
            initial_comp2_step=0.05,
            initial_composition_step=0.05,
            min_comp2_step=0.01,
            min_composition_step=0.01,
        ),
        evaluate_phase=_phase_evaluator(objective, phase_records),
    )

    comp2_offsets = sorted({round(abs(record.x[2] - 0.5), 3) for record in result.evaluations})
    assert 0.025 in comp2_offsets
    assert result.best is not None
    assert result.best.x[2:] == pytest.approx((0.5, 0.5, 0.5))


def test_pattern_search_keeps_dT_values_fixed():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=20,
        analyser=SyntheticAnalyser(comp2_target=0.6, composition_target=0.6),
    )
    source = _phase_from_exact_record(
        objective,
        [12.0, 44.0, 0.5, 0.5, 0.5],
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
        [12.0, 44.0, 0.5, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            vary_comp2_power_frac=False,
            vary_composition=True,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.evaluations
    assert all(record.x[2] == pytest.approx(0.5) for record in result.evaluations)
    assert any(record.x[3] != pytest.approx(0.5) for record in result.evaluations)


def test_pattern_search_respects_fixed_composition():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=12,
        analyser=SyntheticAnalyser(comp2_target=0.7, composition_target=0.9),
    )
    source = _phase_from_exact_record(
        objective,
        [12.0, 44.0, 0.5, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            vary_comp2_power_frac=True,
            vary_composition=False,
        ),
        evaluate_phase=_phase_evaluator(objective, []),
    )

    assert result.evaluations
    assert any(record.x[2] != pytest.approx(0.5) for record in result.evaluations)
    assert all(record.x[3:] == pytest.approx((0.5, 0.5)) for record in result.evaluations)


def test_pattern_search_skips_when_all_polish_variables_are_fixed():
    objective = objective_helpers._OuterDesignObjective(
        max_evals=4,
        analyser=SyntheticAnalyser(),
    )
    source = _phase_from_exact_record(
        objective,
        [12.0, 44.0, 0.5, 0.5, 0.5],
        phase="source",
    )

    result = polish_helpers._run_pattern_search_polish(
        source=source,
        objective=objective,
        config=opt_types.LocalPolishConfig(
            budget=8,
            vary_comp2_power_frac=False,
            vary_composition=False,
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
        [20.0, 35.0, 0.5, 0.5, 0.5],
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

    assert case_definition.refrigerant == "Propane;Butane|0.6;0.4|mass"
    assert case_definition.base_case["T_source_in"] == pytest.approx(42.0)
    assert case_definition.base_case["Tref_evap_inlet"] == pytest.approx(12.0)
    assert case_definition.base_case["T_air_tar"] == pytest.approx(195.0)
    assert case_definition.base_case["dT_GC_min_tar"] == pytest.approx(18.0)
    assert case_definition.base_case["dTmin_rhx"] == pytest.approx(3.0)
    assert case_definition.base_case["ihx_on"] is False
    assert case_definition.base_case["motor_eff"] == pytest.approx(0.96)
    assert case_definition.best_reference["SEC"] == pytest.approx(0.62)
    assert case_definition.best_reference["dT_rhx2_cold"] == pytest.approx(55.0)
    assert case_definition.best_reference["composition"] == "0.6;0.4"
    assert config.fixed_comp2_power_frac == pytest.approx(0.65)
    assert config.line_budget == 4


def test_write_outer_design_optimisation_results_includes_reference_metadata(
    monkeypatch,
    tmp_path,
):
    case_definition = opt_types.ResultsCsvCaseDefinition(
        base_case={"ihx_on": False},
        refrigerant="Propane;Butane|0.5;0.5|mass",
        selected_reference={"SEC": 0.7},
        best_reference={
            "SEC": 0.62,
            "dT_rhx1_cold": 2.5,
            "dT_rhx2_cold": 55.0,
            "comp2_power_frac": 0.65,
            "refrigerant": "Propane;Butane|0.5;0.5|mass",
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
            "refrigerant_ls": [
                "Propane;Butane|0.3;0.7|mass",
                "Propane;Butane|0.7;0.3|mass",
            ],
            "T_air_tar_ls": [180.0, 200.0],
            "dT_rhx1_cold_ls": [0.0, 10.0],
            "dT_rhx2_cold_ls": [5.0],
            "comp2_power_frac_ls": [0.4, 0.6],
        },
        config=outer_opt.OuterDesignOptimisationConfig(line_budget=2),
    )

    assert len(cases) == 16
    assert {case.base_case["T_air_tar"] for case in cases} == {180.0, 200.0}
    assert {case.refrigerant for case in cases} == {
        "Propane;Butane|0.3;0.7|mass",
        "Propane;Butane|0.7;0.3|mass",
    }
    assert {case.config.fixed_dT_rhx1_cold for case in cases} == {0.0, 10.0}
    assert {case.config.fixed_dT_rhx2_cold for case in cases} == {5.0}
    assert {case.config.fixed_comp2_power_frac for case in cases} == {0.4, 0.6}
    assert all(not case.config.vary_dT_rhx1_cold for case in cases)
    assert all(not case.config.vary_dT_rhx2_cold for case in cases)
    assert all(not case.config.vary_comp2_power_frac for case in cases)
    assert all(case.config.vary_composition for case in cases)
    assert len({case.signature for case in cases}) == len(cases)


def test_optimisation_case_generation_omitted_variable_lists_keep_base_config():
    config = outer_opt.OuterDesignOptimisationConfig(
        fixed_comp2_power_frac=0.42,
        vary_comp2_power_frac=True,
        vary_composition=False,
    )

    [case] = case_helpers._build_optimisation_cases(
        param_kwargs={"refrigerant": "Propane;Butane|0.64;0.36|mass"},
        config=config,
    )

    assert case.config.fixed_comp2_power_frac == pytest.approx(0.42)
    assert case.config.vary_comp2_power_frac is True
    assert case.config.vary_composition is False


def test_optimisation_case_generation_uses_refrigerant_composition_as_start_only():
    [case] = case_helpers._build_optimisation_cases(
        param_kwargs={"refrigerant": "Propane;Butane|0.62;0.38|mass"},
        config=outer_opt.OuterDesignOptimisationConfig(),
    )

    assert case.refrigerant == "Propane;Butane|0.62;0.38|mass"
    assert case.config.vary_composition is True


def test_optimisation_case_generation_keeps_nested_refrigerants_atomic():
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={
            "refrigerant_ls": [
                ["Propane", "Butane", "Pentane"],
                [("Propane", 0.2), ("Butane", 0.3), ("Pentane", 0.5)],
                "Propane;Butane;Pentane|0.3;0.4;0.3|mol",
            ]
        }
    )

    assert len(cases) == 3
    assert cases[0].refrigerant.endswith("|mass")
    assert _composition(cases[0].refrigerant) == pytest.approx((1 / 3, 1 / 3, 1 / 3))
    assert _composition(cases[1].refrigerant) == pytest.approx((0.2, 0.3, 0.5))
    assert cases[2].refrigerant.endswith("|mol")


def test_equivalent_list_and_mass_string_cases_share_signature():
    list_case = case_helpers._build_optimisation_case(
        base_case=None,
        refrigerant=["Propane", "Butane"],
    )
    string_case = case_helpers._build_optimisation_case(
        base_case=None,
        refrigerant="Propane;Butane|0.5;0.5|mass",
    )

    assert list_case.signature == string_case.signature


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
    assert all(call["n_cpu"] == 1 for call in calls)
    assert [row["case_signature"] for row in rows] == [case.signature for case in cases]
    assert [write["suffix"] for write in written] == [
        case.temp_suffix for case in cases
    ]
    assert all(write["file_name"] == "temp/optimisation/R" for write in written)


def test_batch_result_row_reports_generic_composition_fields():
    case = case_helpers._build_optimisation_case(
        base_case=None,
        refrigerant=["Propane", "Butane"],
    )

    row = batch_helpers._result_row(case, _fake_batch_result())

    assert row["input_refrigerant"] == "Propane;Butane|0.5;0.5|mass"
    assert row["refrigerant_components"] == "Propane;Butane"
    assert row["refrigerant_component_count"] == 2
    assert row["input_basis"] == "mass"
    assert row["optimisation_basis"] == "mass"
    assert row["initial_mass_fractions"] == "0.5;0.5"
    assert row["best_refrigerant"] == "Propane;Butane|0.5;0.5|mass"
    assert row["best_mass_fractions"] == "0.5;0.5"


def test_optimal_result_transfer_ignores_inner_candidate_files(tmp_path):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={
            "refrigerant": ["Propane", "Butane"],
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


def test_completed_result_transfer_only_merges_the_finished_case(tmp_path):
    cases = case_helpers._build_optimisation_cases(
        param_kwargs={"T_air_tar_ls": [180.0, 200.0]},
    )
    for index, case in enumerate(cases):
        batch_helpers.print_result_to_csv(
            output={
                "case_signature": case.signature,
                "batch_success": True,
                "best_objective": 0.5 + index,
            },
            file_name="temp/optimisation/R",
            suffix=case.temp_suffix,
            results_dir=tmp_path,
        )

    batch_helpers._transfer_completed_optimal_result(
        "RHPopt_results",
        cases[0],
        results_dir=tmp_path,
    )

    master = (tmp_path / "RHPopt_results.csv").read_text()
    assert cases[0].signature in master
    assert cases[1].signature not in master
    remaining_temp = (
        tmp_path
        / "temp"
        / "optimisation"
        / f"R{cases[1].temp_suffix}.csv"
    )
    assert remaining_temp.exists()


def test_single_case_runner_resumes_optimal_result_before_solving(monkeypatch, tmp_path):
    config = outer_opt.OuterDesignOptimisationConfig(max_exact_evaluations=3)
    case = case_helpers._build_optimisation_case(
        base_case={"T_air_tar": 180.0},
        refrigerant="Propane;Butane|0.5;0.5|mass",
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
        refrigerant="Propane;Butane|0.5;0.5|mass",
        config=config,
        results_dir=tmp_path,
    )

    assert row is not None
    assert row["case_signature"] == case.signature
    assert (tmp_path / "RHPopt_results.csv").exists()


def test_batch_worker_passes_attached_incumbents_to_outer_optimiser(
    monkeypatch,
    tmp_path,
):
    case = case_helpers._build_optimisation_case(
        base_case={"T_air_tar": 180.0},
        refrigerant="Propane;Butane|0.5;0.5|mass",
    )
    incumbent = outer_opt.HistoricalIncumbent(
        source_version="edge-first-pattern-search-v4",
        source_file="historical.csv",
        stored_sec=0.48,
        x=(12.0, 44.0, 0.4, 0.6, 0.4),
    )
    case = replace(case, incumbents=(incumbent,))
    calls = []

    def optimise(**kwargs):
        calls.append(kwargs)
        return _fake_batch_result()

    monkeypatch.setattr(batch_helpers, "optimise_recuperative_heat_pump", optimise)

    batch_helpers._solve_optimisation_case(
        case,
        results_dir=str(tmp_path),
        n_cpu=1,
    )

    assert calls[0]["incumbents"] == (incumbent,)


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
    completed_transfers = []
    events = []

    monkeypatch.setattr(
        batch_helpers,
        "_transfer_optimal_results",
        lambda filename, results_dir: (
            transfers.append((filename, results_dir)),
            events.append("transfer"),
        ),
    )
    monkeypatch.setattr(
        batch_helpers,
        "_transfer_completed_optimal_result",
        lambda filename, case, results_dir: (
            completed_transfers.append((filename, case.signature, results_dir)),
            events.append("completed-transfer"),
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

    def run_cases(cases, results_dir, n_cpu, result_callback):
        events.append("optimise")
        row = {"case_signature": cases[0].signature, "batch_success": True}
        result_callback(cases[0], row)
        return [row]

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

    assert rows == [
        {"case_signature": cases[0].signature, "batch_success": True}
    ]
    assert transfers == [
        ("RHPopt_results", tmp_path),
        ("RHPopt_results", tmp_path),
    ]
    assert completed_transfers == [
        ("RHPopt_results", cases[0].signature, tmp_path)
    ]
    assert events == [
        "transfer",
        "filter",
        "optimise",
        "completed-transfer",
        "transfer",
    ]


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


def test_real_serial_and_parallel_multistart_select_same_result(hpl, tmp_path):
    del hpl
    config = outer_opt.OuterDesignOptimisationConfig(
        max_exact_evaluations=2,
        max_composition_starts=2,
        vary_dT_rhx1_cold=False,
        vary_dT_rhx2_cold=False,
        vary_comp2_power_frac=False,
        line_budget=0,
        split_budget=0,
        polish_budget=0,
    )

    serial = outer_opt.optimise_recuperative_heat_pump(
        config=config,
        n_cpu=1,
        results_dir=tmp_path,
    )
    parallel = outer_opt.optimise_recuperative_heat_pump(
        config=config,
        n_cpu=2,
        results_dir=tmp_path,
    )

    assert serial.best is not None
    assert parallel.best is not None
    assert parallel.best.objective == pytest.approx(serial.best.objective)
    assert parallel.best.x == pytest.approx(serial.best.x)
    assert [search.composition for search in parallel.composition_searches] == [
        search.composition for search in serial.composition_searches
    ]


def test_real_local_pattern_search_polish_smoke_when_refprop_is_available(hpl, tmp_path):
    del hpl
    objective = objective_helpers._OuterDesignObjective(
        max_evals=2,
        results_dir=tmp_path,
    )
    source = _phase_from_exact_record(
        objective,
        [0.0, 0.0, 0.5, 0.5, 0.5],
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


def _composition(refrigerant: str) -> tuple[float, ...]:
    return tuple(float(value) for value in refrigerant.split("|")[1].split(";"))


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
        x=(0.0, 0.0, 0.5, 0.5, 0.5),
        refrigerant="Propane;Butane|0.5;0.5|mass",
    )
    return SimpleNamespace(
        best=best,
        corner=best,
        to_summary_dict=lambda: {
            "optimiser_version": outer_opt.OPTIMISER_VERSION,
            "best_objective": 0.5,
            "best_phase": "corner",
            "exact_evaluations": 1,
            "failed_evaluations": 0,
            "best_dT_rhx1_cold": 0.0,
            "best_dT_rhx2_cold": 0.0,
            "best_comp2_power_frac": 0.5,
            "best_mass_fractions": "0.5;0.5",
        },
    )
