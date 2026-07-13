"""Evaporator and suction-pressure alignment for recuperative cycles."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ...common import CycleFailure
from ._discharge_pressure import solve_at_minimum_discharge_pressure

if TYPE_CHECKING:
    from ...common import BaseHeatPump

__all__ = [
    "align_evaporator_inlet_temperature",
    "evaporator_inlet_temperature_offset",
    "validate_initial_evaporator_outlet",
]


def evaporator_inlet_temperature_offset(hp: BaseHeatPump) -> float:
    """Return simulated minus maximum evaporator inlet temperature."""
    return hp.evaporator.sink_in.temperature - hp.Tref_evap_inlet


def validate_initial_evaporator_outlet(hp: BaseHeatPump) -> CycleFailure | None:
    """Validate an explicitly supplied evaporator outlet against the superheater."""
    return _superheater_reversal_failure(hp, hp.Tref_evap_outlet)


def align_evaporator_inlet_temperature(hp: BaseHeatPump) -> CycleFailure | None:
    """Maximise suction pressure without exceeding the evaporator inlet target."""
    tolerance = _alignment_tolerance(hp)
    maximum_temperature = _maximum_evaporator_outlet_temperature(hp)
    maximum_feasible_temperature = (
        maximum_temperature - tolerance
        if maximum_temperature is not None
        else None
    )
    for _ in range(25):
        offset = evaporator_inlet_temperature_offset(hp)
        if abs(offset) <= tolerance:
            hp.Tref_evap_outlet = hp.evaporator.sink_out.temperature
            return None
        if (
            offset < 0.0
            and maximum_feasible_temperature is not None
            and hp.evaporator.sink_out.temperature
            >= maximum_feasible_temperature - tolerance
        ):
            hp.Tref_evap_outlet = hp.evaporator.sink_out.temperature
            return None

        proposed_temperature = hp.evaporator.sink_out.temperature - offset
        if maximum_feasible_temperature is not None:
            proposed_temperature = min(
                proposed_temperature,
                maximum_feasible_temperature,
            )
        failure = _superheater_reversal_failure(hp, proposed_temperature)
        if failure is not None:
            return failure

        hp.evaporator.sink_out.props_in_ab("TQ", proposed_temperature, 1.0)
        hp.Tref_evap_outlet = hp.evaporator.sink_out.temperature
        hp.comp_discharge.pressure_bnds = None
        failure = solve_at_minimum_discharge_pressure(hp)
        if failure is not None:
            return failure

    return CycleFailure(
        "Evaporator inlet temperature correction failed to converge.",
        {},
    )


def _superheater_reversal_failure(
    hp: BaseHeatPump,
    evaporator_outlet_temperature: float,
) -> CycleFailure | None:
    maximum_temperature = _maximum_evaporator_outlet_temperature(hp)
    if maximum_temperature is None:
        return None
    violation = max(0.0, evaporator_outlet_temperature - maximum_temperature)
    if violation <= _alignment_tolerance(hp):
        return None
    return CycleFailure(
        "Evaporator outlet temperature would reverse superheater heat transfer.",
        {"superheater_temperature_reversal": violation},
    )


def _maximum_evaporator_outlet_temperature(hp: BaseHeatPump) -> float | None:
    try:
        return float(hp.evap_sup.sink_out.temperature)
    except AttributeError:
        return None


def _alignment_tolerance(hp: BaseHeatPump) -> float:
    return max(hp.tol, 1e-4)
