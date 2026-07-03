"""Public RHP outer-design optimiser entry point."""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import numpy as np

from ._objective import _OuterDesignObjective
from ._polish import _run_pattern_search_polish
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
    OPTIMISER_VERSION,
    LocalPolishConfig,
    OuterDesignOptimisationConfig,
    OuterDesignOptimisationResult,
    PhaseEvaluation,
)

__all__ = ["optimise_recuperative_heat_pump"]


def optimise_recuperative_heat_pump(
    base_case: dict[str, Any] | None = None,
    *,
    fluids: tuple[str, str] = ("Propane", "Butane"),
    config: OuterDesignOptimisationConfig | None = None,
    results_dir: str | Path = "results/outer_design_optimisation",
) -> OuterDesignOptimisationResult:
    """Optimise recuperative heat-pump outer design variables for one fixed case."""
    config = OuterDesignOptimisationConfig() if config is None else config
    start = time.perf_counter()
    evaluations: list[PhaseEvaluation] = []

    objective = _OuterDesignObjective(
        max_evals=config.max_exact_evaluations,
        base_case=base_case,
        fluids=fluids,
        results_dir=results_dir,
        use_fast_objective=config.use_fast_objective,
    )

    def evaluate_phase(phase: str, x: np.ndarray | list[float]) -> PhaseEvaluation:
        record = objective.evaluate(x)
        phase_record = PhaseEvaluation.from_record(
            phase=phase,
            record=record,
            cumulative_exact_evaluations=objective.exact_evaluations,
        )
        evaluations.append(phase_record)
        return phase_record

    fixed_comp2 = float(config.fixed_comp2_power_frac)
    fixed_composition = float(config.fixed_first_component_mass_fraction)
    fixed_dT1 = float(config.fixed_dT_rhx1_cold)
    fixed_dT2 = float(config.fixed_dT_rhx2_cold)

    corner = evaluate_phase(
        "corner",
        [
            fixed_dT1 if not config.vary_dT_rhx1_cold else 0.0,
            fixed_dT2 if not config.vary_dT_rhx2_cold else 0.0,
            fixed_comp2,
            fixed_composition,
        ],
    )
    if config.vary_dT_rhx1_cold and config.vary_dT_rhx2_cold:
        edge_a = _search_edge_a(
            objective=objective,
            config=config,
            fixed_comp2=fixed_comp2,
            fixed_composition=fixed_composition,
            evaluate_phase=evaluate_phase,
        )
        edge_b = _search_edge_b(
            objective=objective,
            config=config,
            fixed_comp2=fixed_comp2,
            fixed_composition=fixed_composition,
            evaluate_phase=evaluate_phase,
        )
        split = _search_split_line(
            objective=objective,
            config=config,
            edge_a=edge_a,
            edge_b=edge_b,
            fixed_comp2=fixed_comp2,
            fixed_composition=fixed_composition,
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
            fixed_composition=fixed_composition,
            evaluate_phase=evaluate_phase,
        )
        split = None
    elif config.vary_dT_rhx2_cold:
        edge_a = _search_dT2_at_fixed_dT1(
            objective=objective,
            config=config,
            fixed_dT1=fixed_dT1,
            fixed_comp2=fixed_comp2,
            fixed_composition=fixed_composition,
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

    polish_sources = _select_polish_sources(
        [corner, edge_a.best, edge_b.best, None if split is None else split.best],
        top_n=config.polish_top_n,
    )
    polish_results = tuple(
        _run_pattern_search_polish(
            source=source,
            objective=objective,
            config=LocalPolishConfig(
                budget=config.polish_budget,
                vary_comp2_power_frac=config.vary_comp2_power_frac,
                vary_first_component_mass_fraction=(
                    config.vary_first_component_mass_fraction
                ),
            ),
            evaluate_phase=evaluate_phase,
        )
        for source in polish_sources
    )

    return OuterDesignOptimisationResult(
        optimiser_version=OPTIMISER_VERSION,
        config=config,
        corner=corner,
        edge_a=edge_a,
        edge_b=edge_b,
        split=split,
        polish_results=polish_results,
        best=_best_phase_evaluation(evaluations),
        evaluations=tuple(evaluations),
        elapsed_seconds=time.perf_counter() - start,
        exact_evaluations=objective.exact_evaluations,
        failed_evaluations=objective.failed_evaluations,
    )
