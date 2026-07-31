"""Composition multistart orchestration for the RHP outer optimiser."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, replace
from pathlib import Path
import time
from typing import Any

import numpy as np

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import (
    normalise_refrigerant_input,
)

from ._composition_starts import _full_search_budget
from ._objective import _OuterDesignObjective
from ._search import (
    _best_phase_evaluation,
    _search_dT1_at_fixed_dT2,
    _search_dT2_at_fixed_dT1,
    _search_edge_a,
    _search_edge_b,
    _search_split_line,
    _select_polish_sources,
    _skipped_line_search,
)
from ._types import (
    PENALTY_OBJECTIVE,
    CompositionSearchResult,
    LineSearchResult,
    LocalPolishResult,
    OuterDesignOptimisationConfig,
    PhaseEvaluation,
)


@dataclass(frozen=True)
class _CompositionSearchTask:
    start_index: int
    composition: tuple[float, ...]
    max_evaluations: int
    base_case: dict[str, Any]
    refrigerant: str
    config: OuterDesignOptimisationConfig
    results_dir: str


def _build_composition_search_tasks(
    starts: tuple[tuple[float, ...], ...],
    *,
    budget_per_start: int,
    base_case: dict[str, Any],
    refrigerant: str,
    config: OuterDesignOptimisationConfig,
    results_dir: str | Path,
) -> tuple[_CompositionSearchTask, ...]:
    return tuple(
        _CompositionSearchTask(
            start_index=index,
            composition=composition,
            max_evaluations=budget_per_start,
            base_case=dict(base_case),
            refrigerant=refrigerant,
            config=config,
            results_dir=str(results_dir),
        )
        for index, composition in enumerate(starts)
    )


def _run_composition_searches(
    tasks: tuple[_CompositionSearchTask, ...],
    *,
    n_cpu: int,
    executor_factory=None,
    objective_factory=_OuterDesignObjective,
) -> tuple[CompositionSearchResult, ...]:
    if not tasks:
        return ()
    if max(1, int(n_cpu)) == 1:
        return tuple(
            _run_composition_search(task, objective_factory=objective_factory)
            for task in tasks
        )

    executor_factory = ProcessPoolExecutor if executor_factory is None else executor_factory
    try:
        executor = executor_factory(max_workers=min(max(1, int(n_cpu)), len(tasks)))
    except Exception:
        return tuple(
            _run_composition_search(task, objective_factory=objective_factory)
            for task in tasks
        )

    try:
        futures = [executor.submit(_run_composition_search_task, task) for task in tasks]
        results = []
        for task, future in zip(tasks, futures):
            try:
                results.append(future.result())
            except Exception as exc:
                results.append(
                    _failed_composition_search(
                        task,
                        f"{type(exc).__name__}: {exc}",
                    )
                )
        return tuple(results)
    finally:
        executor.shutdown()


def _run_composition_search_task(
    task: _CompositionSearchTask,
) -> CompositionSearchResult:
    return _run_composition_search(task, objective_factory=_OuterDesignObjective)


def _run_composition_search(
    task: _CompositionSearchTask,
    *,
    objective_factory,
) -> CompositionSearchResult:
    start = time.perf_counter()
    objective = objective_factory(
        max_evals=task.max_evaluations,
        base_case=task.base_case,
        refrigerant=task.refrigerant,
        results_dir=task.results_dir,
        use_fast_objective=task.config.use_fast_objective,
    )
    evaluations: list[PhaseEvaluation] = []

    def evaluate_phase(phase: str, x: np.ndarray | list[float]) -> PhaseEvaluation:
        record = objective.evaluate(x)
        phase_record = PhaseEvaluation.from_record(
            phase=phase,
            record=record,
            cumulative_exact_evaluations=objective.exact_evaluations,
        )
        evaluations.append(phase_record)
        return phase_record

    config = _search_config_for_budget(task.config, task.max_evaluations)
    fixed_comp2 = float(config.fixed_comp2_power_frac)
    fixed_dT1 = float(config.fixed_dT_rhx1_cold)
    fixed_dT2 = float(config.fixed_dT_rhx2_cold)
    corner = evaluate_phase(
        "corner",
        np.asarray(
            [
                fixed_dT1 if not config.vary_dT_rhx1_cold else 0.0,
                fixed_dT2 if not config.vary_dT_rhx2_cold else 0.0,
                fixed_comp2,
                *task.composition,
            ]
        ),
    )

    if config.vary_dT_rhx1_cold and config.vary_dT_rhx2_cold:
        edge_a = _search_edge_a(
            objective=objective,
            config=config,
            fixed_comp2=fixed_comp2,
            fixed_composition=task.composition,
            evaluate_phase=evaluate_phase,
        )
        edge_b = _search_edge_b(
            objective=objective,
            config=config,
            fixed_comp2=fixed_comp2,
            fixed_composition=task.composition,
            evaluate_phase=evaluate_phase,
        )
        split = _search_split_line(
            objective=objective,
            config=config,
            edge_a=edge_a,
            edge_b=edge_b,
            fixed_comp2=fixed_comp2,
            fixed_composition=task.composition,
            evaluate_phase=evaluate_phase,
        )
    elif config.vary_dT_rhx1_cold:
        edge_a = _skipped_line_search(
            "edge_a",
            "dT_rhx2_cold",
            "dT_rhx2_cold is fixed",
        )
        edge_b = _search_dT1_at_fixed_dT2(
            objective=objective,
            config=config,
            fixed_dT2=fixed_dT2,
            fixed_comp2=fixed_comp2,
            fixed_composition=task.composition,
            evaluate_phase=evaluate_phase,
        )
        split = None
    elif config.vary_dT_rhx2_cold:
        edge_a = _search_dT2_at_fixed_dT1(
            objective=objective,
            config=config,
            fixed_dT1=fixed_dT1,
            fixed_comp2=fixed_comp2,
            fixed_composition=task.composition,
            evaluate_phase=evaluate_phase,
        )
        edge_b = _skipped_line_search(
            "edge_b",
            "dT_rhx1_cold",
            "dT_rhx1_cold is fixed",
        )
        split = None
    else:
        edge_a = _skipped_line_search(
            "edge_a",
            "dT_rhx2_cold",
            "dT_rhx1_cold and dT_rhx2_cold are fixed",
        )
        edge_b = _skipped_line_search(
            "edge_b",
            "dT_rhx1_cold",
            "dT_rhx1_cold and dT_rhx2_cold are fixed",
        )
        split = None

    success = any(record.success for record in evaluations)
    return CompositionSearchResult(
        start_index=task.start_index,
        composition=task.composition,
        corner=corner,
        edge_a=edge_a,
        edge_b=edge_b,
        split=split,
        evaluations=tuple(evaluations),
        success=success,
        message="success" if success else "no successful evaluations",
        elapsed_seconds=time.perf_counter() - start,
        exact_evaluations=objective.exact_evaluations,
        failed_evaluations=objective.failed_evaluations,
    )


def _search_config_for_budget(
    config: OuterDesignOptimisationConfig,
    max_evaluations: int,
) -> OuterDesignOptimisationConfig:
    """Use spare per-start budget to strengthen line and split coverage."""
    extra = max(0, int(max_evaluations) - _full_search_budget(config))
    line_budget = int(config.line_budget)
    split_budget = int(config.split_budget)
    active_line_count = int(config.vary_dT_rhx1_cold) + int(
        config.vary_dT_rhx2_cold
    )
    target_line_budget = max(line_budget, 12)
    target_split_budget = max(split_budget, 12)

    while active_line_count and line_budget < target_line_budget:
        if extra < active_line_count:
            break
        line_budget += 1
        extra -= active_line_count
    if config.vary_dT_rhx1_cold and config.vary_dT_rhx2_cold:
        split_budget += min(extra, target_split_budget - split_budget)
    if (
        line_budget == config.line_budget
        and split_budget == config.split_budget
    ):
        return config
    return replace(
        config,
        line_budget=line_budget,
        split_budget=split_budget,
    )


def _failed_composition_search(
    task: _CompositionSearchTask,
    message: str,
) -> CompositionSearchResult:
    config = task.config
    x = (
        float(config.fixed_dT_rhx1_cold) if not config.vary_dT_rhx1_cold else 0.0,
        float(config.fixed_dT_rhx2_cold) if not config.vary_dT_rhx2_cold else 0.0,
        float(config.fixed_comp2_power_frac),
        *task.composition,
    )
    spec = normalise_refrigerant_input(task.refrigerant)
    corner = PhaseEvaluation(
        phase="corner",
        x=x,
        objective=PENALTY_OBJECTIVE + 1.0,
        success=False,
        message=message,
        sec=None,
        cop=None,
        pmax=None,
        elapsed_seconds=0.0,
        cumulative_exact_evaluations=0,
        constraint_errors={},
        refrigerant=task.refrigerant,
        mass_refrigerant=spec.format(task.composition, "mass"),
    )
    edge_a = _skipped_line_search("edge_a", "dT_rhx2_cold", message)
    edge_b = _skipped_line_search("edge_b", "dT_rhx1_cold", message)
    return CompositionSearchResult(
        start_index=task.start_index,
        composition=task.composition,
        corner=corner,
        edge_a=edge_a,
        edge_b=edge_b,
        split=None,
        evaluations=(corner,),
        success=False,
        message=message,
        elapsed_seconds=0.0,
        exact_evaluations=0,
        failed_evaluations=0,
    )


def _best_direction_result(
    searches: tuple[CompositionSearchResult, ...],
    direction: str,
) -> LineSearchResult | None:
    results = [
        result
        for search in searches
        if (result := getattr(search, direction)) is not None
    ]
    if not results:
        return None
    successful = [
        result for result in results if result.best is not None and result.best.success
    ]
    if not successful:
        return results[0]
    return min(successful, key=lambda result: result.best.objective)


def _select_direction_polish_sources(
    *,
    corner: PhaseEvaluation,
    edge_a: LineSearchResult,
    edge_b: LineSearchResult,
    split: LineSearchResult | None,
    searches: tuple[CompositionSearchResult, ...] = (),
    top_n: int,
) -> tuple[PhaseEvaluation, ...]:
    candidates = [
        candidate
        for search in searches
        for candidate in (
            search.edge_a.best,
            search.edge_b.best,
            None if search.split is None else search.split.best,
        )
    ]
    if not candidates:
        candidates = [edge_a.best, edge_b.best, None if split is None else split.best]
    sources = _select_polish_sources(candidates, top_n=top_n)
    if sources:
        return sources
    return _select_polish_sources([corner], top_n=1)


def _best_evaluation(
    searches: tuple[CompositionSearchResult, ...],
    polishes: tuple[LocalPolishResult, ...],
    additional: tuple[PhaseEvaluation, ...] = (),
) -> PhaseEvaluation | None:
    records = [record for search in searches for record in search.evaluations]
    records.extend(record for polish in polishes for record in polish.evaluations)
    records.extend(additional)
    return _best_phase_evaluation(records)
