"""Sequential component analysis for recuperative heat-pump cycles."""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ...common import BaseHeatPump

__all__ = ["run_sequential_cycle_analysis", "set_gas_cooler_sink_temperatures"]


def run_sequential_cycle_analysis(
    discharge_pressure: float,
    hp: BaseHeatPump,
) -> BaseHeatPump:
    """Sequentially solve each component until the coupled cycle converges."""
    hp.comp_discharge.pressure = discharge_pressure
    iteration = 0
    most_negative_residual = 0.0

    while (
        iteration == 0 or most_negative_residual < 0.0
    ) and iteration < 10:
        enthalpy = hp.exp_valve.flow_in.enthalpy
        pressure_drop = (
            hp.comp_discharge.pressure - hp.gas_cooler1.source_out.pressure
        )
        compressor_power = sum(compressor.power for compressor in hp.compressors)
        gas_cooler_duty = sum(gas_cooler.duty for gas_cooler in hp.gas_coolers)
        iteration += 1

        for operation in hp.operations:
            operation.solve()

        most_negative_residual = min(
            hp.tol * enthalpy
            - abs(hp.exp_valve.flow_in.enthalpy - enthalpy),
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
                sum(gas_cooler.duty for gas_cooler in hp.gas_coolers)
                - gas_cooler_duty
            ),
        )

    hp.iteration_converged = most_negative_residual >= 0.0
    hp.converged = hp.iteration_converged and all(
        getattr(operation, "converged", False) for operation in hp.operations
    )
    set_gas_cooler_sink_temperatures(hp)
    if hp.debug:
        hp.get_cycle_plot(x_axis="s", y_axis="T")
        hp.get_cycle_plot(x_axis="h", y_axis="P")
        hp.get_cycle_plot(x_axis="h", y_axis="T")
    return hp


def set_gas_cooler_sink_temperatures(hp: BaseHeatPump) -> BaseHeatPump:
    """Back-calculate sink stream temperatures and mass flows for gas coolers."""
    hp.gas_coolers[0].sink_out.temperature = hp.T_air_tar
    hp.gas_coolers[0].sink_out.props_in_mPT()

    mass_flow = sum(gas_cooler.duty for gas_cooler in hp.gas_coolers) / (
        hp.gas_coolers[0].sink_out.enthalpy
        - hp.gas_coolers[-1].sink_in.enthalpy
    )
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
        gas_cooler.sink_in.props_in_mPH()
    return hp
