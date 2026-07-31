"""Gas-cooler pinch calculations for recuperative heat-pump cycles."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

import numpy as np
from scipy.optimize import minimize_scalar

from ....utils.graph import plot_driving_force
from ....utils.solve_profile import profiled_section
from ._cycle_analysis import run_sequential_cycle_analysis

if TYPE_CHECKING:
    from ....components import GasCooler, Stream
    from ...common import BaseHeatPump

__all__ = [
    "find_gc_pinch_residual",
    "find_internal_gc_approach",
    "get_all_gc_pinch_candidates",
    "update_gc_pinch_properties",
]

PinchCandidate: TypeAlias = tuple[float, float, bool]

_INTERNAL_PINCH_RELATIVE_DUTY_TOLERANCE = 3e-3
_INTERNAL_PINCH_MAX_ITERATIONS = 40


@profiled_section("gc_dtmin_eval")
def find_gc_pinch_residual(
    discharge_pressure: float,
    hp: BaseHeatPump,
    *,
    return_abs: bool = False,
) -> float:
    """Return current minus target gas-cooler minimum temperature difference."""
    candidates = get_all_gc_pinch_candidates(discharge_pressure, hp)
    candidates.pop(-1)
    if candidates[-1][1] < 1e-3:
        candidates.pop(-1)
    residual = min(candidates)[0] - hp.dT_GC_min_tar
    return abs(residual) if return_abs else residual


def get_all_gc_pinch_candidates(
    discharge_pressure: float,
    hp: BaseHeatPump,
) -> list[PinchCandidate]:
    """Enumerate terminal and internal pinch candidates for every gas cooler."""
    run_sequential_cycle_analysis(discharge_pressure, hp)
    candidates: list[PinchCandidate] = []
    for gas_cooler in hp.gas_coolers:
        hot_approach = (
            gas_cooler.source_in.temperature - gas_cooler.sink_out.temperature
        )
        cold_approach = (
            gas_cooler.source_out.temperature - gas_cooler.sink_in.temperature
        )
        candidates.append((hot_approach, 1.0, True))
        candidates.append(find_internal_gc_approach(gas_cooler))
        candidates.append((cold_approach, 0.0, True))

    if hp.debug:
        plot_driving_force(hp.gas_coolers, dtmin=hp.dT_GC_min_tar)
    return candidates


def find_internal_gc_approach(gas_cooler: GasCooler) -> PinchCandidate:
    """Find the minimum internal temperature approach in one gas cooler."""
    gas_cooler.source_pinch.pressure = gas_cooler.source_out.pressure
    gas_cooler.sink_pinch.pressure = gas_cooler.sink_in.pressure
    if gas_cooler.duty <= 0.0:
        return (np.inf, 0.0, False)

    result = minimize_scalar(
        lambda duty: _pinch_temperature_difference(duty, gas_cooler),
        method="bounded",
        bounds=(0.001 * gas_cooler.duty, 0.999 * gas_cooler.duty),
        options={
            # The internal duty location only supports the explicitly checked
            # 0.05 K pressure-target band; cycle tolerance is unnecessarily
            # strict here and multiplies expensive mixture P-h flashes.
            "xatol": max(
                _INTERNAL_PINCH_RELATIVE_DUTY_TOLERANCE * gas_cooler.duty,
                1e-9,
            ),
            "maxiter": _INTERNAL_PINCH_MAX_ITERATIONS,
        },
    )
    return (result.fun, result.x / gas_cooler.duty, result.success)


def update_gc_pinch_properties(
    discharge_pressure: float,
    hp: BaseHeatPump,
) -> BaseHeatPump:
    """Refresh pinch stream properties after discharge pressure changes."""
    candidates = get_all_gc_pinch_candidates(discharge_pressure, hp)
    candidate_index = 1
    for gas_cooler in hp.gas_coolers:
        _set_pinch_properties(candidates[candidate_index][1], gas_cooler)
        candidate_index += 3
    return hp


def _pinch_temperature_difference(duty: float, gas_cooler: GasCooler) -> float:
    _set_pinch_properties(duty, gas_cooler)
    return _temperature_difference(
        gas_cooler.source_pinch,
        gas_cooler.sink_pinch,
    )


def _set_pinch_properties(duty: float, gas_cooler: GasCooler) -> GasCooler:
    gas_cooler.source_pinch.enthalpy = (
        gas_cooler.source_out.enthalpy
        + duty / gas_cooler.source_out.mass_flow
    )
    gas_cooler.source_pinch.props_in_mPH()
    gas_cooler.sink_pinch.enthalpy = (
        gas_cooler.sink_in.enthalpy + duty / gas_cooler.sink_in.mass_flow
    )
    gas_cooler.sink_pinch.props_in_mPH()
    return gas_cooler


def _temperature_difference(source: Stream, sink: Stream) -> float:
    return source.temperature - sink.temperature
