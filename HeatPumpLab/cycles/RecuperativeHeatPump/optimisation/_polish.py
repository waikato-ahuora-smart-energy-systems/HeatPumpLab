"""Pattern-search polish for fixed recuperator lifts."""

from __future__ import annotations

import time
from collections.abc import Callable

import numpy as np

from ._objective import _OuterDesignObjective
from ._search import _best_phase_evaluation, _is_better_phase
from ._types import LocalPolishConfig, LocalPolishResult, PhaseEvaluation


def _run_pattern_search_polish(
    *,
    source: PhaseEvaluation,
    objective: _OuterDesignObjective,
    config: LocalPolishConfig | None = None,
    evaluate_phase: Callable[[str, np.ndarray | list[float]], PhaseEvaluation],
) -> LocalPolishResult:
    """Polish one fixed-recuperator source candidate with pattern search."""
    config = LocalPolishConfig() if config is None else config
    method = "pattern_search"
    start = time.perf_counter()
    start_exact = objective.exact_evaluations
    start_failed = objective.failed_evaluations

    if config.budget <= 0:
        return LocalPolishResult(
            method=method,
            status="skipped",
            source_x=source.x,
            best=None,
            evaluations=(),
            success=False,
            message="polish budget was zero",
            elapsed_seconds=time.perf_counter() - start,
            exact_evaluations=0,
            failed_evaluations=0,
        )

    evaluator = _LocalPolishEvaluator(
        source=source,
        objective=objective,
        config=config,
        evaluate_phase=evaluate_phase,
    )
    try:
        _poll_until_converged(evaluator)
        status = "ok"
        message = "success"
    except Exception as exc:
        status = "failed"
        message = f"{type(exc).__name__}: {exc}"

    best = _best_phase_evaluation(evaluator.evaluations)
    success = bool(best and best.success and status == "ok")
    if status == "ok" and not success:
        status = "failed"
        message = "no successful evaluations"
    return LocalPolishResult(
        method=method,
        status=status,
        source_x=source.x,
        best=best,
        evaluations=tuple(evaluator.evaluations),
        success=success,
        message=message,
        elapsed_seconds=time.perf_counter() - start,
        exact_evaluations=objective.exact_evaluations - start_exact,
        failed_evaluations=objective.failed_evaluations - start_failed,
    )


class _LocalPolishEvaluator:
    """Budget-aware fixed-dT evaluator for local polish methods."""

    def __init__(
        self,
        *,
        source: PhaseEvaluation,
        objective: _OuterDesignObjective,
        config: LocalPolishConfig,
        evaluate_phase: Callable[[str, np.ndarray | list[float]], PhaseEvaluation],
    ) -> None:
        self.source = source
        self.objective = objective
        self.config = config
        self.evaluate_phase = evaluate_phase
        self.phase = "polish_pattern_search"
        self.bounds = np.asarray(objective.bounds[2:4], dtype=float)
        self.fixed_dT = np.asarray(source.x[:2], dtype=float)
        self.start_comp = np.asarray(source.x[2:4], dtype=float)
        self.start_exact = objective.exact_evaluations
        self.evaluations: list[PhaseEvaluation] = []

    @property
    def exact_used(self) -> int:
        return self.objective.exact_evaluations - self.start_exact

    @property
    def remaining_exact(self) -> int:
        return max(0, int(self.config.budget) - self.exact_used)

    def build_x(self, comp_values: np.ndarray | list[float]) -> np.ndarray:
        comp = np.clip(
            np.asarray(comp_values, dtype=float),
            self.bounds[:, 0],
            self.bounds[:, 1],
        )
        return np.asarray([self.fixed_dT[0], self.fixed_dT[1], comp[0], comp[1]])

    def can_evaluate(self, comp_values: np.ndarray | list[float]) -> bool:
        x = self.objective.repair_vector(self.build_x(comp_values))
        cache_key = self.objective._cache_key(x)
        return cache_key in self.objective.cache or self.remaining_exact > 0

    def evaluate(self, comp_values: np.ndarray | list[float]) -> PhaseEvaluation | None:
        if not self.can_evaluate(comp_values):
            return None
        record = self.evaluate_phase(self.phase, self.build_x(comp_values))
        self.evaluations.append(record)
        return record


def _poll_until_converged(evaluator: _LocalPolishEvaluator) -> None:
    current = evaluator.evaluate(evaluator.start_comp)
    if current is None:
        return

    steps = np.asarray(evaluator.config.initial_steps, dtype=float)
    min_steps = np.asarray(evaluator.config.min_steps, dtype=float)
    shrink_factor = float(evaluator.config.shrink_factor)
    if not (0.0 < shrink_factor < 1.0):
        raise ValueError("shrink_factor must be between 0 and 1.")

    failed_shrinks = 0
    current_comp = np.asarray(current.x[2:4], dtype=float)
    while (
        evaluator.remaining_exact > 0
        and failed_shrinks <= evaluator.config.max_failed_shrinks
        and np.any(steps >= min_steps)
    ):
        best_poll = current
        best_comp = current_comp
        for idx in range(2):
            if steps[idx] < min_steps[idx]:
                continue
            for sign in (1.0, -1.0):
                candidate = current_comp.copy()
                candidate[idx] = np.clip(
                    candidate[idx] + sign * steps[idx],
                    evaluator.bounds[idx, 0],
                    evaluator.bounds[idx, 1],
                )
                if np.allclose(candidate, current_comp, atol=1e-14, rtol=0.0):
                    continue
                record = evaluator.evaluate(candidate)
                if record is None:
                    break
                if _is_better_phase(
                    record,
                    best_poll,
                    tolerance=evaluator.config.improvement_tolerance,
                ):
                    best_poll = record
                    best_comp = np.asarray(record.x[2:4], dtype=float)

        if best_poll is not current:
            current = best_poll
            current_comp = best_comp
            failed_shrinks = 0
        else:
            steps *= shrink_factor
            failed_shrinks += 1
