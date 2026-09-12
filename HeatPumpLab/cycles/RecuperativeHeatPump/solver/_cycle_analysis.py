"""Sequential component analysis for recuperative heat-pump cycles."""

from __future__ import annotations

from math import isclose, isfinite
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...common import BaseHeatPump

__all__ = ["run_sequential_cycle_analysis", "set_gas_cooler_sink_temperatures"]


class InvalidSinkStateError(ValueError):
    """A pressure probe cannot support a physical sink reconstruction."""

    def __init__(self, message, constraint_errors):
        super().__init__(message)
        self.constraint_errors = constraint_errors


def run_sequential_cycle_analysis(
    discharge_pressure: float,
    hp: BaseHeatPump,
) -> BaseHeatPump:
    """Sequentially solve each component until the coupled cycle converges."""
    hp.comp_discharge.pressure = discharge_pressure
    iteration = 0
    most_negative_residual = 0.0
    inlet_states_converged = False

    while (
        iteration == 0 or most_negative_residual < 0.0 or not inlet_states_converged
    ) and iteration < 10:
        inlet_states = _compressor_inlet_states(hp)
        enthalpy = hp.exp_valve.flow_in.enthalpy
        pressure_drop = hp.comp_discharge.pressure - hp.gas_cooler1.source_out.pressure
        compressor_power = sum(compressor.power for compressor in hp.compressors)
        gas_cooler_duty = sum(gas_cooler.duty for gas_cooler in hp.gas_coolers)
        iteration += 1

        for operation in hp.operations:
            operation.solve()

        # A recuperator can update a compressor inlet after that compressor has
        # run. Aggregate heat and power may still be unchanged in this sweep;
        # require the next sweep to consume the changed state before accepting it.
        inlet_states_converged = all(
            isclose(before, after, rel_tol=hp.tol, abs_tol=1e-12)
            for before, after in zip(inlet_states, _compressor_inlet_states(hp))
        )
        most_negative_residual = min(
            hp.tol * enthalpy - abs(hp.exp_valve.flow_in.enthalpy - enthalpy),
            hp.tol * pressure_drop
            - abs(
                hp.comp_discharge.pressure
                - hp.gas_cooler1.source_out.pressure
                - pressure_drop
            ),
            hp.tol * compressor_power
            - abs(
                sum(compressor.power for compressor in hp.compressors)
                - compressor_power
            ),
            hp.tol * gas_cooler_duty
            - abs(
                sum(gas_cooler.duty for gas_cooler in hp.gas_coolers) - gas_cooler_duty
            ),
        )

    hp.iteration_converged = most_negative_residual >= 0.0 and inlet_states_converged
    hp.converged = hp.iteration_converged and all(
        getattr(operation, "converged", False) for operation in hp.operations
    )
    set_gas_cooler_sink_temperatures(hp)
    if hp.debug:
        hp.get_cycle_plot(x_axis="s", y_axis="T")
        hp.get_cycle_plot(x_axis="h", y_axis="P")
        hp.get_cycle_plot(x_axis="h", y_axis="T")
    return hp


def _compressor_inlet_states(hp: BaseHeatPump) -> tuple[float, ...]:
    """Capture all stages' inlet states consumed by the compression solves."""
    return tuple(
        value
        for compressor in hp.compressors
        for value in (
            compressor.flow_in.pressure,
            compressor.flow_in.enthalpy,
            compressor.flow_in.mass_flow,
        )
    )


def set_gas_cooler_sink_temperatures(hp: BaseHeatPump) -> BaseHeatPump:
    """Back-calculate sink stream temperatures and mass flows for gas coolers."""
    duties = [gas_cooler.duty for gas_cooler in hp.gas_coolers]
    if not all(isfinite(duty) for duty in duties) or sum(duties) <= 0:
        raise InvalidSinkStateError("Nonpositive or nonfinite total gas-cooler duty.",
                                    {"invalid_sink_heat_duty": 1.0})
    reversal = max((-gc.duty for gc in hp.gas_coolers
                    if gc.duty < -getattr(gc, "tol", 1e-6)), default=0.)
    if reversal:
        raise InvalidSinkStateError("Negative gas-cooler duty before sink reconstruction.",
                                    {"gas_cooler_duty_reversal": reversal})
    hp.gas_coolers[0].sink_out.temperature = hp.T_air_tar
    hp.gas_coolers[0].sink_out.props_in_mPT()

    hot_enthalpy = hp.gas_coolers[0].sink_out.enthalpy
    cold_enthalpy = hp.gas_coolers[-1].sink_in.enthalpy
    enthalpy_rise = hot_enthalpy - cold_enthalpy
    if not isfinite(enthalpy_rise) or enthalpy_rise <= 0:
        raise InvalidSinkStateError("Nonpositive or nonfinite sink enthalpy rise.",
                                    {"invalid_sink_enthalpy_rise": 1.0})
    mass_flow = sum(duties) / enthalpy_rise
    if not isfinite(mass_flow) or mass_flow <= 0:
        raise InvalidSinkStateError("Nonpositive or nonfinite sink mass flow.",
                                    {"invalid_sink_mass_flow": 1.0})
    hp.gas_coolers[0].sink_out.mass_flow = mass_flow
    hp.gas_coolers[0].sink_out.heat_flow = (
        mass_flow * hp.gas_coolers[0].sink_out.enthalpy
    )

    for gas_cooler in hp.gas_coolers:
        gas_cooler.sink_pinch.mass_flow = mass_flow
        gas_cooler.sink_in.mass_flow = mass_flow
        gas_cooler.sink_in.enthalpy = (
            gas_cooler.sink_out.enthalpy - gas_cooler.duty / mass_flow
        )
        # Even a tolerated negative round-off duty can be amplified by a tiny
        # net flow. Do not send an unbounded intermediate enthalpy to REFPROP.
        enthalpy = gas_cooler.sink_in.enthalpy
        if not isfinite(enthalpy) or not cold_enthalpy-1e-6 <= enthalpy <= hot_enthalpy+1e-6:
            raise InvalidSinkStateError("Sink enthalpy outside its inlet/outlet bounds.",
                                        {"invalid_sink_enthalpy": 1.0})
        gas_cooler.sink_in.props_in_mPH()
    return hp
