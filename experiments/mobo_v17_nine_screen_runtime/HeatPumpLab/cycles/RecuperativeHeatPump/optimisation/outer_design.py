"""Public RHP outer-design optimiser entry point."""

from __future__ import annotations

from pathlib import Path
import time
from typing import Any

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import RefrigerantInput

from ._composition_starts import (
    _generate_composition_starts,
    _plan_composition_searches,
)
from ._incumbents import (
    best_successful_incumbent_replay,
    replay_historical_incumbents,
)
from ._multistart import (
    _best_direction_result,
    _best_evaluation,
    _build_composition_search_tasks,
    _run_composition_searches,
    _select_direction_polish_sources,
)
from ._objective import _OuterDesignObjective
from ._polish import (
    _build_lift_feedback_tasks,
    _build_polish_tasks,
    _run_lift_feedbacks,
    _run_polishes,
)
from ._types import (
    OPTIMISER_VERSION,
    HistoricalIncumbent,
    OuterDesignOptimisationConfig,
    OuterDesignOptimisationResult,
)

__all__ = ["optimise_recuperative_heat_pump"]


def optimise_recuperative_heat_pump(
    base_case: dict[str, Any] | None = None,
    *,
    refrigerant: RefrigerantInput = "Propane;Butane|0.5;0.5|mass",
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = 1,
    results_dir: str | Path = "results/outer_design_optimisation",
    incumbents: tuple[HistoricalIncumbent, ...] = (),
) -> OuterDesignOptimisationResult:
    """Optimise recuperative heat-pump outer design variables for one fixed case."""
    config = OuterDesignOptimisationConfig() if config is None else config
    start = time.perf_counter()

    problem = _OuterDesignObjective(
        max_evals=0,
        base_case=base_case,
        refrigerant=refrigerant,
        results_dir=results_dir,
        use_fast_objective=config.use_fast_objective,
    )
    requested_starts = _generate_composition_starts(
        problem.initial_mass_fractions,
        vary_composition=config.vary_composition,
        max_starts=config.max_composition_starts,
        cache_decimals=problem.cache_decimals,
    )
    selected_starts, search_budget = _plan_composition_searches(
        requested_starts,
        config,
    )
    search_tasks = _build_composition_search_tasks(
        selected_starts,
        budget_per_start=search_budget,
        base_case=problem.base_case,
        refrigerant=problem.input_refrigerant,
        config=config,
        results_dir=results_dir,
    )
    composition_searches = _run_composition_searches(
        search_tasks,
        n_cpu=n_cpu,
        objective_factory=_OuterDesignObjective,
    )

    corner = composition_searches[0].corner
    edge_a = _best_direction_result(composition_searches, "edge_a")
    edge_b = _best_direction_result(composition_searches, "edge_b")
    if edge_a is None or edge_b is None:
        raise RuntimeError("Composition search did not return both edge results.")
    split = _best_direction_result(composition_searches, "split")
    polish_sources = _select_direction_polish_sources(
        corner=corner,
        edge_a=edge_a,
        edge_b=edge_b,
        split=split,
        searches=composition_searches,
        top_n=config.polish_top_n,
    )
    search_exact = sum(search.exact_evaluations for search in composition_searches)
    polish_tasks = _build_polish_tasks(
        polish_sources,
        remaining_budget=max(0, config.max_exact_evaluations - search_exact),
        base_case=problem.base_case,
        refrigerant=problem.input_refrigerant,
        config=config,
        results_dir=results_dir,
    )
    polish_results = _run_polishes(
        polish_tasks,
        n_cpu=n_cpu,
        objective_factory=_OuterDesignObjective,
    )
    normal_polish_exact = sum(
        polish.exact_evaluations for polish in polish_results
    )
    lift_feedback_tasks = _build_lift_feedback_tasks(
        polish_results,
        remaining_budget=max(
            0,
            config.max_exact_evaluations
            - search_exact
            - normal_polish_exact,
        ),
        base_case=problem.base_case,
        refrigerant=problem.input_refrigerant,
        config=config,
        results_dir=results_dir,
    )
    lift_feedback_results = _run_lift_feedbacks(
        lift_feedback_tasks,
        n_cpu=n_cpu,
        objective_factory=_OuterDesignObjective,
    )

    incumbent_replays = replay_historical_incumbents(
        incumbents=incumbents,
        base_case=problem.base_case,
        refrigerant=problem.input_refrigerant,
        config=config,
        results_dir=results_dir,
        objective_factory=_OuterDesignObjective,
    )
    incumbent_source = best_successful_incumbent_replay(incumbent_replays)
    incumbent_polish_tasks = _build_polish_tasks(
        () if incumbent_source is None else (incumbent_source,),
        remaining_budget=config.incumbent_polish_budget,
        base_case=problem.base_case,
        refrigerant=problem.input_refrigerant,
        config=config,
        results_dir=results_dir,
        phase="incumbent_polish",
        per_source_budget=config.incumbent_polish_budget,
    )
    incumbent_polish_results = _run_polishes(
        incumbent_polish_tasks,
        n_cpu=n_cpu,
        objective_factory=_OuterDesignObjective,
    )

    evaluations = tuple(
        record for search in composition_searches for record in search.evaluations
    )
    evaluations += tuple(replay.evaluation for replay in incumbent_replays)
    evaluations += tuple(
        record for polish in polish_results for record in polish.evaluations
    )
    evaluations += tuple(
        record
        for result in lift_feedback_results
        for record in result.evaluations
    )
    evaluations += tuple(
        record
        for polish in incumbent_polish_results
        for record in polish.evaluations
    )
    feedback_exact = sum(
        len(result.evaluations) for result in lift_feedback_results
    )
    exact_evaluations = search_exact + normal_polish_exact + feedback_exact
    exact_evaluations += len(incumbent_replays) + sum(
        polish.exact_evaluations for polish in incumbent_polish_results
    )
    failed_evaluations = sum(
        search.failed_evaluations for search in composition_searches
    ) + sum(polish.failed_evaluations for polish in polish_results)
    failed_evaluations += sum(
        not record.success
        for result in lift_feedback_results
        for record in result.evaluations
    )
    failed_evaluations += sum(
        not replay.evaluation.success for replay in incumbent_replays
    ) + sum(
        polish.failed_evaluations for polish in incumbent_polish_results
    )
    additional_candidates = tuple(
        replay.evaluation for replay in incumbent_replays
    )
    additional_candidates += tuple(
        record
        for result in lift_feedback_results
        for record in result.evaluations
    )
    return OuterDesignOptimisationResult(
        optimiser_version=OPTIMISER_VERSION,
        config=config,
        corner=corner,
        edge_a=edge_a,
        edge_b=edge_b,
        split=split,
        composition_searches=composition_searches,
        requested_composition_starts=len(requested_starts),
        polish_results=polish_results,
        lift_feedback_results=lift_feedback_results,
        incumbent_replays=incumbent_replays,
        incumbent_polish_results=incumbent_polish_results,
        best=_best_evaluation(
            composition_searches,
            polish_results + incumbent_polish_results,
            additional_candidates,
        ),
        evaluations=evaluations,
        elapsed_seconds=time.perf_counter() - start,
        exact_evaluations=exact_evaluations,
        failed_evaluations=failed_evaluations,
    )
