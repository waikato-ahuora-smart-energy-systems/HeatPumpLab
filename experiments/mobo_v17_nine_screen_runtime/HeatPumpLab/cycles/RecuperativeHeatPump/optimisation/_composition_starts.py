"""Deterministic composition starts and case-wide search budgeting."""

from __future__ import annotations

from itertools import combinations

import numpy as np

from ._types import OuterDesignOptimisationConfig


def _generate_composition_starts(
    initial: tuple[float, ...],
    *,
    vary_composition: bool,
    max_starts: int,
    cache_decimals: int = 8,
) -> tuple[tuple[float, ...], ...]:
    """Return deterministic, well-separated mass-fraction starts."""
    initial_arr = _normalise_composition(initial)
    component_count = len(initial_arr)
    if component_count == 1 or not vary_composition or max_starts <= 1:
        return (tuple(float(value) for value in initial_arr),)

    candidates = [initial_arr, np.full(component_count, 1.0 / component_count)]
    candidates.extend(np.eye(component_count))
    for first, second in combinations(range(component_count), 2):
        candidate = np.zeros(component_count)
        candidate[first] = 0.5
        candidate[second] = 0.5
        candidates.append(candidate)
    if component_count == 2:
        candidates.extend((np.asarray((0.25, 0.75)), np.asarray((0.75, 0.25))))

    unique: dict[tuple[float, ...], np.ndarray] = {}
    for candidate in candidates:
        key = tuple(float(value) for value in np.round(candidate, cache_decimals))
        unique.setdefault(key, candidate)

    initial_key = tuple(float(value) for value in np.round(initial_arr, cache_decimals))
    selected = [initial_arr]
    remaining = [candidate for key, candidate in unique.items() if key != initial_key]
    while remaining and len(selected) < max(1, int(max_starts)):
        index = max(
            range(len(remaining)),
            key=lambda item: (
                min(
                    float(np.linalg.norm(remaining[item] - chosen))
                    for chosen in selected
                ),
                tuple(float(-value) for value in remaining[item]),
            ),
        )
        selected.append(remaining.pop(index))
    return tuple(tuple(float(value) for value in candidate) for candidate in selected)


def _plan_composition_searches(
    starts: tuple[tuple[float, ...], ...],
    config: OuterDesignOptimisationConfig,
) -> tuple[tuple[tuple[float, ...], ...], int]:
    """Select starts that fit while reserving the requested polish budget."""
    if not starts:
        raise ValueError("At least one composition start is required.")
    budget_per_start = _full_search_budget(config)
    if config.max_exact_evaluations < budget_per_start:
        return starts[:1], config.max_exact_evaluations

    polish_enabled = config.vary_comp2_power_frac or config.vary_composition
    maximum_starts = min(
        len(starts),
        max(1, config.max_exact_evaluations // budget_per_start),
    )
    for start_count in range(maximum_starts, 0, -1):
        polish_count = min(config.polish_top_n, start_count) if polish_enabled else 0
        requested = (
            start_count * budget_per_start
            + polish_count * config.polish_budget
        )
        if requested <= config.max_exact_evaluations:
            search_allowance = (
                config.max_exact_evaluations
                - polish_count * config.polish_budget
            )
            distributed_budget = min(
                budget_per_start + 12,
                search_allowance // start_count,
            )
            return starts[:start_count], distributed_budget
    return starts[:1], budget_per_start


def _full_search_budget(config: OuterDesignOptimisationConfig) -> int:
    budget = 1
    if config.vary_dT_rhx1_cold:
        budget += config.line_budget
    if config.vary_dT_rhx2_cold:
        budget += config.line_budget
    if config.vary_dT_rhx1_cold and config.vary_dT_rhx2_cold:
        budget += config.split_budget
    return budget


def _normalise_composition(values: tuple[float, ...]) -> np.ndarray:
    composition = np.asarray(values, dtype=float)
    if composition.ndim != 1 or composition.size == 0:
        raise ValueError("Composition must contain at least one fraction.")
    if np.any(~np.isfinite(composition)) or np.any(composition < 0.0):
        raise ValueError("Composition fractions must be finite and nonnegative.")
    total = float(composition.sum())
    if total <= 0.0:
        raise ValueError("Composition must contain at least one positive fraction.")
    return composition / total
