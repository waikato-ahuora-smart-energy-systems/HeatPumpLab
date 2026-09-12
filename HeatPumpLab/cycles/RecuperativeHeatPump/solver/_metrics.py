"""Performance metrics for recuperative heat-pump cycles."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ...common import calculate_cop, calculate_sec, cycle_has_converged
from ....utils.solve_profile import profile_section
from ....utils.properties_package import RefpropError

if TYPE_CHECKING:
    from ...common import BaseHeatPump

__all__ = ["populate_full_metrics", "populate_objective_metrics"]


def populate_full_metrics(hp: BaseHeatPump) -> BaseHeatPump:
    """Populate all reported recuperative heat-pump metrics."""
    with profile_section("metrics"):
        hp.converged = cycle_has_converged(hp)
        _populate_headline_metrics(hp)
        # Total VHC uses the combined suction volume flow of active compressors.
        hp.VHC = (
            _total_heat_duty(hp) / sum(
                compressor.flow_in.vol_flow for compressor in hp.compressors
            )
            if hp.converged
            else 0.0
        )
        hot_exergy_temperature = _total_heat_duty(hp) / sum(
            gas_cooler.delta_entropy_src for gas_cooler in hp.gas_coolers
        )
        cold_exergy_temperature = sum(
            evaporator.duty for evaporator in hp.evaporators
        ) / sum(
            evaporator.delta_entropy_snk for evaporator in hp.evaporators
        )
        hp.eDTL = (
            hot_exergy_temperature - cold_exergy_temperature
            if hp.converged
            else 0.0
        )
        hp.SLE = hp.SEC / hp.eDTL * 1000 if hp.converged else 0.0
        hp.metric_errors = {}
        try:
            hp.LMTD_gcs = sum(
                gas_cooler.duty * gas_cooler.get_LMTD_ave()
                for gas_cooler in hp.gas_coolers
            ) / _total_heat_duty(hp)
        except RefpropError as exc:
            # Reporting samples do not define cycle feasibility. Keep the
            # unavailable metric explicit without reusing a failed flash.
            hp.LMTD_gcs = float("nan")
            hp.metric_errors["LMTD_gcs"] = str(exc)
        hp.cycle_type = (
            "Transcritical"
            if hp.comp_discharge.pressure > hp.comp_discharge.P_crit
            else "Subcritical"
        )
        return hp


def populate_objective_metrics(hp: BaseHeatPump) -> BaseHeatPump:
    """Populate only metrics required by objective evaluations."""
    with profile_section("objective_metrics"):
        hp.converged = cycle_has_converged(hp)
        _populate_headline_metrics(hp)
        return hp


def _populate_headline_metrics(hp: BaseHeatPump) -> None:
    heat_duty = _total_heat_duty(hp)
    compressor_power = _total_compressor_power(hp)
    hp.COP = calculate_cop(
        heat_duty,
        compressor_power,
        converged=hp.converged,
    )
    hp.SEC = calculate_sec(
        heat_duty,
        compressor_power,
        converged=hp.converged,
    )


def _total_heat_duty(hp: BaseHeatPump) -> float:
    return sum(gas_cooler.duty for gas_cooler in hp.gas_coolers)


def _total_compressor_power(hp: BaseHeatPump) -> float:
    return sum(compressor.power for compressor in hp.compressors)
