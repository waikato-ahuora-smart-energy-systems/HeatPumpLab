"""Simplex-aware pattern-search polish for fixed recuperator lifts."""

from __future__ import annotations

import time
from collections.abc import Callable, Iterable

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
    """Polish compressor split and composition without changing recuperator lifts."""
    config = LocalPolishConfig() if config is None else config
    start = time.perf_counter()
    start_exact = objective.exact_evaluations
    start_failed = objective.failed_evaluations
    can_vary_composition = config.vary_composition and len(source.x[3:]) > 1

    if config.budget <= 0:
        return _skipped_result(source, start, "polish budget was zero")
    if not config.vary_comp2_power_frac and not can_vary_composition:
        return _skipped_result(
            source,
            start,
            "no local polish variables enabled",
            best=source,
        )

    evaluator = _LocalPolishEvaluator(
        source=source,
        objective=objective,
        config=config,
        evaluate_phase=evaluate_phase,
        vary_composition=can_vary_composition,
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
        method="pattern_search",
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


def _skipped_result(
    source: PhaseEvaluation,
    start: float,
    message: str,
    *,
    best: PhaseEvaluation | None = None,
) -> LocalPolishResult:
    return LocalPolishResult(
        method="pattern_search",
        status="skipped",
        source_x=source.x,
        best=best,
        evaluations=(),
        success=bool(best and best.success),
        message=message,
        elapsed_seconds=time.perf_counter() - start,
        exact_evaluations=0,
        failed_evaluations=0,
    )


class _LocalPolishEvaluator:
    def __init__(
        self,
        *,
        source: PhaseEvaluation,
        objective: _OuterDesignObjective,
        config: LocalPolishConfig,
        evaluate_phase: Callable[[str, np.ndarray | list[float]], PhaseEvaluation],
        vary_composition: bool,
    ) -> None:
        self.source = source
        self.objective = objective
        self.config = config
        self.evaluate_phase = evaluate_phase
        self.vary_composition = vary_composition
        self.start_exact = objective.exact_evaluations
        self.evaluations: list[PhaseEvaluation] = []

    @property
    def exact_used(self) -> int:
        return self.objective.exact_evaluations - self.start_exact

    @property
    def remaining_exact(self) -> int:
        return max(0, int(self.config.budget) - self.exact_used)

    def can_evaluate(self, x: np.ndarray) -> bool:
        repaired = self.objective.repair_vector(x)
        return self.objective._cache_key(repaired) in self.objective.cache or self.remaining_exact > 0

    def evaluate(self, x: np.ndarray) -> PhaseEvaluation | None:
        if not self.can_evaluate(x):
            return None
        record = self.evaluate_phase("polish_pattern_search", x)
        self.evaluations.append(record)
        return record


def _poll_until_converged(evaluator: _LocalPolishEvaluator) -> None:
    current = evaluator.evaluate(np.asarray(evaluator.source.x, dtype=float))
    if current is None:
        return
    shrink_factor = float(evaluator.config.shrink_factor)
    if not 0.0 < shrink_factor < 1.0:
        raise ValueError("shrink_factor must be between 0 and 1.")

    comp2_step = float(evaluator.config.initial_comp2_step)
    composition_step = float(evaluator.config.initial_composition_step)
    failed_shrinks = 0
    while evaluator.remaining_exact > 0 and failed_shrinks <= evaluator.config.max_failed_shrinks:
        comp2_active = (
            evaluator.config.vary_comp2_power_frac
            and comp2_step >= evaluator.config.min_comp2_step
        )
        composition_active = (
            evaluator.vary_composition
            and composition_step >= evaluator.config.min_composition_step
        )
        if not comp2_active and not composition_active:
            break

        best_poll = current
        for candidate in _poll_candidates(
            np.asarray(current.x, dtype=float),
            objective=evaluator.objective,
            comp2_step=comp2_step if comp2_active else None,
            composition_step=composition_step if composition_active else None,
        ):
            record = evaluator.evaluate(candidate)
            if record is None:
                break
            if _is_better_phase(
                record,
                best_poll,
                tolerance=evaluator.config.improvement_tolerance,
            ):
                best_poll = record

        if best_poll is not current:
            current = best_poll
            failed_shrinks = 0
        else:
            if comp2_active:
                comp2_step *= shrink_factor
            if composition_active:
                composition_step *= shrink_factor
            failed_shrinks += 1


def _poll_candidates(
    current: np.ndarray,
    *,
    objective: _OuterDesignObjective,
    comp2_step: float | None,
    composition_step: float | None,
) -> Iterable[np.ndarray]:
    if comp2_step is not None:
        for sign in (1.0, -1.0):
            candidate = current.copy()
            candidate[2] = np.clip(
                candidate[2] + sign * comp2_step,
                objective.bounds[2, 0],
                objective.bounds[2, 1],
            )
            if not np.allclose(candidate, current, atol=1e-14, rtol=0.0):
                yield candidate

    if composition_step is None:
        return
    composition = current[3:]
    reference = int(np.argmax(composition))
    for index in range(len(composition)):
        if index == reference:
            continue
        for donor, receiver in ((index, reference), (reference, index)):
            amount = min(float(composition_step), float(composition[donor]))
            if amount <= 0.0:
                continue
            candidate = current.copy()
            candidate[3 + donor] -= amount
            candidate[3 + receiver] += amount
            yield candidate
