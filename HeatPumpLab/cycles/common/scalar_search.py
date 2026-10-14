"""Reusable scalar target searches for positive, logarithmic domains."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from math import exp, log

import numpy as np
from scipy.optimize import minimize_scalar, root_scalar

__all__ = [
    "ScalarSearchResult",
    "find_lowest_log_bracketed_target",
    "find_sampled_log_target",
]

ScalarEvaluator = Callable[[float], float]

# A sample that fails like this is unknown and skipped; any other exception is a defect and propagates.
# RefpropError is a ValueError.
_SAMPLE_FAILURES = (ValueError, RuntimeError, ArithmeticError)


@dataclass(frozen=True)
class ScalarSearchResult:
    """Result from a bounded scalar target search."""

    value: float
    residual: float
    success: bool


def find_lowest_log_bracketed_target(
    evaluate: ScalarEvaluator,
    bounds: tuple[float, float],
    residuals: tuple[float, float],
    *,
    target_tolerance: float,
    relative_tolerance: float,
    max_iterations: int = 50,
) -> ScalarSearchResult:
    """Find the lowest positive value entering the accepted residual band."""
    lower, upper = (float(value) for value in bounds)
    residual_lower, residual_upper = (float(value) for value in residuals)
    if abs(residual_lower) <= target_tolerance:
        return ScalarSearchResult(lower, abs(residual_lower), True)
    if not (residual_lower < 0.0 < residual_upper):
        return ScalarSearchResult(lower, abs(residual_lower), False)

    best_value = upper if abs(residual_upper) <= target_tolerance else None
    best_residual = residual_upper
    for _ in range(max_iterations):
        midpoint = exp((log(lower) + log(upper)) / 2.0)
        residual_midpoint = float(evaluate(midpoint))
        if residual_midpoint >= -target_tolerance:
            upper = midpoint
            residual_upper = residual_midpoint
            best_value = midpoint
            best_residual = residual_midpoint
        else:
            lower = midpoint
            residual_lower = residual_midpoint

        if (upper - lower) / upper <= max(relative_tolerance, 1e-7):
            break

    success = best_value is not None and abs(best_residual) <= target_tolerance
    return ScalarSearchResult(
        lower if best_value is None else best_value,
        abs(residual_lower if best_value is None else best_residual),
        success,
    )


def find_sampled_log_target(
    evaluate: ScalarEvaluator,
    bounds: tuple[float, float],
    initial: float,
    *,
    target_tolerance: float,
    x_tolerance: float,
    sample_count: int = 25,
) -> ScalarSearchResult:
    """Find the lowest accepted target in a potentially non-monotonic domain."""
    lower, upper = (float(value) for value in bounds)
    if lower <= 0.0 or upper <= lower:
        return ScalarSearchResult(float(initial), np.inf, False)

    values = np.geomspace(lower, upper, num=max(2, int(sample_count)))
    values = np.unique(np.append(values, np.clip(initial, lower, upper)))
    samples: list[tuple[float, float]] = []
    for value in values:
        try:
            residual = float(evaluate(float(value)))
        except _SAMPLE_FAILURES:
            continue
        if np.isfinite(residual):
            samples.append((float(value), residual))

    if not samples:
        return ScalarSearchResult(float(initial), np.inf, False)

    samples.sort(key=lambda sample: sample[0])
    candidates: list[tuple[float, float]] = []
    for index, (value, residual) in enumerate(samples):
        if abs(residual) <= target_tolerance:
            if index > 0 and abs(samples[index - 1][1]) > target_tolerance:
                value = _refine_tolerance_entry(
                    evaluate,
                    samples[index - 1][0],
                    value,
                    target_tolerance,
                    x_tolerance,
                )
                residual = float(evaluate(value))
            if abs(residual) <= target_tolerance:
                candidates.append((value, abs(residual)))

    for left, right in zip(samples, samples[1:]):
        if left[1] * right[1] >= 0.0:
            continue
        try:
            result = root_scalar(
                evaluate,
                method="brentq",
                bracket=(left[0], right[0]),
                xtol=x_tolerance,
                maxiter=50,
            )
        except (ValueError, RuntimeError):
            continue
        if result.converged:
            # a sign change across a discontinuity "converges" on the jump: keep only roots in the band
            residual = abs(float(evaluate(result.root)))
            if residual <= target_tolerance:
                candidates.append((float(result.root), residual))

    for index in range(1, len(samples) - 1):
        left, centre, right = samples[index - 1 : index + 2]
        if abs(centre[1]) > min(abs(left[1]), abs(right[1])):
            continue
        try:
            result = minimize_scalar(
                lambda value: abs(evaluate(value)),
                method="bounded",
                bounds=(left[0], right[0]),
                options={"xatol": x_tolerance},
            )
        except (ValueError, RuntimeError):
            continue
        if result.success and result.fun <= target_tolerance:
            value = float(result.x)
            if abs(left[1]) > target_tolerance:
                value = _refine_tolerance_entry(
                    evaluate,
                    left[0],
                    value,
                    target_tolerance,
                    x_tolerance,
                )
            residual = abs(float(evaluate(value)))
            if residual <= target_tolerance:
                candidates.append((value, residual))

    if candidates:
        value, residual = min(candidates, key=lambda candidate: candidate[0])
        return ScalarSearchResult(value, residual, True)

    best_residual = min(abs(residual) for _, residual in samples)
    return ScalarSearchResult(float(initial), best_residual, False)


def _refine_tolerance_entry(
    evaluate: ScalarEvaluator,
    outside: float,
    inside: float,
    target_tolerance: float,
    x_tolerance: float,
) -> float:
    try:
        result = root_scalar(
            lambda value: abs(evaluate(value)) - target_tolerance,
            method="brentq",
            bracket=(outside, inside),
            xtol=x_tolerance,
            maxiter=50,
        )
    except ValueError:
        return inside
    return float(result.root) if result.converged else inside
