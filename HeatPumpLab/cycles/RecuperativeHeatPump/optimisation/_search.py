"""Edge and split-line searches for the RHP outer optimiser."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from scipy import optimize

from ._objective import _OuterDesignObjective
from ._types import (
    PENALTY_OBJECTIVE,
    LineSearchResult,
    OuterDesignOptimisationConfig,
    PhaseEvaluation,
)


def _search_edge_a(
    *,
    objective: _OuterDesignObjective,
    config: OuterDesignOptimisationConfig,
    fixed_comp2: float,
    fixed_composition: float,
    evaluate_phase: Callable[[str, np.ndarray | list[float]], PhaseEvaluation],
) -> LineSearchResult:
    return _line_search(
        phase="edge_a",
        variable_name="dT_rhx2_cold",
        lo=0.0,
        hi=_recuperator_edge_high(objective),
        budget=config.line_budget,
        build_x=lambda value: [0.0, value, fixed_comp2, fixed_composition],
        evaluate_phase=evaluate_phase,
    )


def _search_edge_b(
    *,
    objective: _OuterDesignObjective,
    config: OuterDesignOptimisationConfig,
    fixed_comp2: float,
    fixed_composition: float,
    evaluate_phase: Callable[[str, np.ndarray | list[float]], PhaseEvaluation],
) -> LineSearchResult:
    return _line_search(
        phase="edge_b",
        variable_name="dT_rhx1_cold",
        lo=0.0,
        hi=_recuperator_edge_high(objective),
        budget=config.line_budget,
        build_x=lambda value: [value, 0.0, fixed_comp2, fixed_composition],
        evaluate_phase=evaluate_phase,
    )


def _search_split_line(
    *,
    objective: _OuterDesignObjective,
    config: OuterDesignOptimisationConfig,
    edge_a: LineSearchResult,
    edge_b: LineSearchResult,
    fixed_comp2: float,
    fixed_composition: float,
    evaluate_phase: Callable[[str, np.ndarray | list[float]], PhaseEvaluation],
) -> LineSearchResult | None:
    dT2_star = 0.0 if edge_a.best is None else edge_a.best.x[1]
    dT1_star = 0.0 if edge_b.best is None else edge_b.best.x[0]
    if (
        abs(dT1_star) <= config.corner_tolerance
        and abs(dT2_star) <= config.corner_tolerance
    ):
        return None

    split_sum = min(objective.max_recuperator_lift, dT1_star + dT2_star)
    if split_sum <= config.corner_tolerance:
        return None

    return _line_search(
        phase="split",
        variable_name="s",
        lo=0.0,
        hi=1.0,
        budget=config.split_budget,
        build_x=lambda s: [
            split_sum * s,
            split_sum * (1.0 - s),
            fixed_comp2,
            fixed_composition,
        ],
        evaluate_phase=evaluate_phase,
    )


def _select_polish_sources(
    candidates: list[PhaseEvaluation | None],
    *,
    top_n: int,
) -> tuple[PhaseEvaluation, ...]:
    unique: dict[tuple[float, ...], PhaseEvaluation] = {}
    for candidate in candidates:
        if candidate is None or not candidate.success:
            continue
        key = tuple(round(value, 8) for value in candidate.x)
        if key not in unique or candidate.objective < unique[key].objective:
            unique[key] = candidate
    return tuple(
        sorted(unique.values(), key=lambda record: record.objective)[: max(0, int(top_n))]
    )


def _best_phase_evaluation(
    records: list[PhaseEvaluation] | tuple[PhaseEvaluation, ...],
) -> PhaseEvaluation | None:
    if not records:
        return None
    successful = [record for record in records if record.success]
    candidates = successful if successful else list(records)
    return min(candidates, key=lambda record: record.objective)


def _is_better_phase(
    candidate: PhaseEvaluation | None,
    incumbent: PhaseEvaluation | None,
    *,
    tolerance: float = 0.0,
) -> bool:
    if candidate is None:
        return False
    if incumbent is None:
        return True
    if candidate.success and not incumbent.success:
        return True
    if not candidate.success and incumbent.success:
        return False
    return candidate.objective < incumbent.objective - float(tolerance)


def _line_search(
    *,
    phase: str,
    variable_name: str,
    lo: float,
    hi: float,
    budget: int,
    build_x: Callable[[float], np.ndarray | list[float]],
    evaluate_phase: Callable[[str, np.ndarray | list[float]], PhaseEvaluation],
) -> LineSearchResult:
    budget = max(0, int(budget))
    if budget <= 0:
        return LineSearchResult(
            phase=phase,
            variable_name=variable_name,
            best=None,
            evaluations=(),
            success=False,
            message="line search budget was zero",
        )

    evaluations: list[PhaseEvaluation] = []
    sample_budget = min(budget, max(3, budget // 2))
    for value in _line_sample_points(lo, hi, sample_budget):
        evaluations.append(evaluate_phase(phase, build_x(float(value))))

    remaining = budget - len(evaluations)
    if remaining > 0 and hi > lo:
        calls = 0

        def scalar(value: float) -> float:
            nonlocal calls
            if calls >= remaining:
                return PENALTY_OBJECTIVE
            calls += 1
            record = evaluate_phase(phase, build_x(float(value)))
            evaluations.append(record)
            return record.objective

        optimize.minimize_scalar(
            scalar,
            bounds=(float(lo), float(hi)),
            method="bounded",
            options={"maxiter": remaining, "xatol": 1e-4},
        )

    best = _best_phase_evaluation(evaluations)
    return LineSearchResult(
        phase=phase,
        variable_name=variable_name,
        best=best,
        evaluations=tuple(evaluations),
        success=bool(best and best.success),
        message="success" if best and best.success else "no successful evaluations",
    )


def _line_sample_points(lo: float, hi: float, n_points: int) -> np.ndarray:
    n_points = max(0, int(n_points))
    if n_points == 0:
        return np.empty(0, dtype=float)
    if n_points == 1:
        return np.asarray([(lo + hi) / 2.0], dtype=float)
    return np.linspace(float(lo), float(hi), n_points)


def _recuperator_edge_high(objective: _OuterDesignObjective) -> float:
    return float(
        min(
            objective.bounds[0, 1],
            objective.bounds[1, 1],
            objective.max_recuperator_lift,
        )
    )
