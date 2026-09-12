"""Cycle-agnostic solve status, convergence, and performance helpers."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Mapping

import numpy as np

if TYPE_CHECKING:
    from .heatPump import BaseHeatPump

__all__ = [
    "CycleFailure",
    "calculate_cop",
    "calculate_sec",
    "cycle_convergence_failure",
    "cycle_has_converged",
]


@dataclass(frozen=True)
class CycleFailure:
    """Describe an expected cycle-solve failure without using control-flow exceptions."""

    message: str
    constraint_errors: Mapping[str, float]

    def __post_init__(self) -> None:
        errors: dict[str, float] = {}
        for name, value in self.constraint_errors.items():
            numeric = float(value)
            if math.isfinite(numeric):
                errors[str(name)] = max(0.0, numeric)
        object.__setattr__(self, "constraint_errors", errors)


def cycle_has_converged(hp: BaseHeatPump) -> bool:
    """Return whether a cycle reached a feasible state with positive compressor work."""
    return bool(
        getattr(hp, "feasible", False)
        and getattr(hp, "iteration_converged", False)
        and all(
            getattr(operation, "converged", False)
            for operation in getattr(hp, "operations", ())
        )
        and sum(compressor.power for compressor in getattr(hp, "compressors", ()))
        > 0.0
    )


def cycle_convergence_failure(hp: BaseHeatPump) -> CycleFailure | None:
    """Return convergence diagnostics, or ``None`` for a valid cycle state."""
    if cycle_has_converged(hp):
        return None

    errors: dict[str, float] = {}
    if not getattr(hp, "iteration_converged", False):
        errors["cycle_iteration"] = 1.0

    unconverged = [
        operation
        for operation in getattr(hp, "operations", ())
        if not getattr(operation, "converged", False)
    ]
    temperature_error = max(
        (_heat_exchanger_temperature_error(operation) for operation in unconverged),
        default=0.0,
    )
    if temperature_error > 0.0:
        errors["heat_exchanger_temperature_reversal"] = temperature_error
    elif unconverged:
        errors["component_convergence"] = 1.0

    return CycleFailure(
        "Cycle operations failed to converge to a physically valid state.",
        errors or {"cycle_convergence": 1.0},
    )


def calculate_cop(
    heat_duty: float,
    compressor_power: float,
    *,
    converged: bool,
) -> float:
    """Return coefficient of performance for a converged heat-pump cycle."""
    return float(heat_duty / compressor_power) if converged else 0.0


def calculate_sec(
    heat_duty: float,
    compressor_power: float,
    *,
    converged: bool,
) -> float:
    """Return specific electricity consumption for a converged cycle."""
    return (
        float(compressor_power / (heat_duty - compressor_power))
        if converged
        else np.inf
    )


def _heat_exchanger_temperature_error(operation: object) -> float:
    try:
        source_in = operation.source_in.temperature
        source_out = operation.source_out.temperature
        sink_in = operation.sink_in.temperature
        sink_out = operation.sink_out.temperature
    except AttributeError:
        return 0.0
    return max(
        0.0,
        source_out - source_in,
        sink_in - sink_out,
        sink_in - source_out,
        sink_out - source_in,
    )
