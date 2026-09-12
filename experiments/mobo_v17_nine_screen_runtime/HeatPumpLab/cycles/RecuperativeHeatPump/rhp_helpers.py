"""Helpers that assemble the Recuperative Heat Pump object graph."""

from __future__ import annotations
from typing import TYPE_CHECKING
import numpy as np
from scipy.optimize import root_scalar
from ...components import (
    Compressor,
    CustomRelationship,
    Evaporator,
    ExpansionValve,
    GasCooler,
    PlateHeatExchanger,
    SetRelationship,
    Stream,
)

if TYPE_CHECKING:
    from ...cycles import RecuperativeHP

__all__ = ["set_rhp_model"]

#######################################################################################################
# Public API
#######################################################################################################


def set_rhp_model(hp: RecuperativeHP) -> RecuperativeHP:
    """Populate the heat pump instance with components, streams, and logic."""
    hp = _get_RecuperativeHP_structure(hp)
    if _check_input_param_feasibility(hp):
        hp.feasible = True
        hp = _get_RecuperativeHP_components(hp)
        hp = _get_RecuperativeHP_streams(hp)
        hp = _get_RecuperativeHP_logic_blocks(hp)
        hp = _get_RecuperativeHP_solver_pathway(hp)
    else:
        hp.feasible = False
    return hp


#######################################################################################################
# Helper functions
#######################################################################################################


def _check_input_param_feasibility(hp: RecuperativeHP):
    """Return True when the provided design knobs represent a feasible case."""
    return (
        (hp.T_air_tar - hp.T_air_inlet) * 2 / 3
        >= hp.dT_rhx1_cold + hp.dT_rhx2_cold + hp.dT_rhx3_cold
        and 0.0 <= hp.dT_gc_out_min_offset <= 50.0
        and hp.comp2_power_frac + hp.comp3_power_frac <= 1.0
    )


def _get_RecuperativeHP_structure(hp: RecuperativeHP) -> RecuperativeHP:
    """Derive binary on/off flags for compressors and recuperators."""
    # Check the PHE combo
    hp.dT_rhx1_cold = hp.dT_rhx1_cold if hp.dT_rhx1_cold > 0.0 else 0.0
    hp.dT_rhx2_cold = (
        hp.dT_rhx2_cold if hp.n_comp > 1 and hp.dT_rhx2_cold > 0.0 else 0.0
    )
    hp.dT_rhx3_cold = (
        hp.dT_rhx3_cold if hp.n_comp > 2 and hp.dT_rhx3_cold > 0.0 else 0.0
    )

    # Check the compressor combo
    if hp.comp_PR_uniform or hp.comp2_power_frac + hp.comp3_power_frac >= 1.0:
        hp.comp2_power_frac = 1 / hp.n_comp if hp.n_comp > 1 else 0.0
        hp.comp3_power_frac = 1 / hp.n_comp if hp.n_comp > 2 else 0.0
    else:
        hp.comp2_power_frac = hp.comp2_power_frac if hp.n_comp > 1 else 0.0
        hp.comp3_power_frac = hp.comp3_power_frac if hp.n_comp > 2 else 0.0

    if hp.comp2_power_frac == 0.0:
        hp.comp3_power_frac = 0.0

    if hp.n_comp == 3 and hp.comp3_power_frac == 0.0:
        hp.n_comp = 2

    if hp.n_comp == 2 and hp.comp2_power_frac == 0.0:
        hp.n_comp = 1

    hp.comp1_power_frac = 1 - (hp.comp2_power_frac + hp.comp3_power_frac)

    # Binary variables for the heat pump structure
    hp.comp2_on = 1 if hp.n_comp > 1 else 0
    hp.comp3_on = 1 if hp.n_comp > 2 else 0
    hp.T_gc_out = hp.T_air_inlet + hp.dT_GC_min_tar + hp.dT_gc_out_min_offset
    hp.ihx_on = 1 if hp.T_source_in < hp.T_gc_out and hp.ihx_on == 1 else 0

    hp.rhx1_on = 1 if hp.dT_rhx1_cold > 0.0 else 0
    hp.rhx2_on = 1 if hp.dT_rhx2_cold > 0.0 and hp.comp2_on == 1 else 0
    hp.rhx3_on = 1 if hp.dT_rhx3_cold > 0.0 and hp.comp3_on == 1 else 0

    return hp


def _get_RecuperativeHP_components(hp: RecuperativeHP) -> RecuperativeHP:
    """Instantiate the component objects required for the configured structure."""
    # Initialisation of unit operations
    hp.comp1 = Compressor(
        name="Compressor",
        fluid=hp.refrigerant,
        motor_eff=hp.motor_eff,
        model_index=hp.comp_model_index,
    )
    if hp.comp2_on == 1:
        hp.comp2 = Compressor(
            name="Compressor",
            fluid=hp.refrigerant,
            motor_eff=hp.motor_eff,
            model_index=hp.comp_model_index,
        )
        hp.gas_cooler3 = GasCooler(
            name="Gas cooler",
            source_fluid=hp.refrigerant,
            sink_fluid=hp.sink,
            dTmin=hp.dT_GC_min_tar,
        )
        if hp.comp3_on == 1:
            hp.comp3 = Compressor(
                name="Compressor",
                fluid=hp.refrigerant,
                motor_eff=hp.motor_eff,
                model_index=hp.comp_model_index,
            )
            hp.gas_cooler4 = GasCooler(
                name="Gas cooler",
                source_fluid=hp.refrigerant,
                sink_fluid=hp.sink,
                dTmin=hp.dT_GC_min_tar,
            )

    hp.gas_cooler1 = GasCooler(
        name="Gas cooler",
        source_fluid=hp.refrigerant,
        sink_fluid=hp.sink,
        dTmin=hp.dT_GC_min_tar,
    )
    hp.gas_cooler2 = GasCooler(
        name="Gas cooler",
        source_fluid=hp.refrigerant,
        sink_fluid=hp.sink,
        dTmin=hp.dT_GC_min_tar,
    )

    hp.ihx = PlateHeatExchanger(
        name="Internal heat exchanger",
        source_fluid=hp.refrigerant,
        sink_fluid=hp.refrigerant,
        dTmin=hp.dTmin_rhx,
        model="dTmin_terminal",
    )
    if hp.rhx1_on == 1:
        hp.rhx1 = PlateHeatExchanger(
            name="Recuperator",
            source_fluid=hp.refrigerant,
            sink_fluid=hp.refrigerant,
            dT_c=hp.dT_rhx1_cold,
            model="fixed cold delta temperature",
        )
    if hp.rhx2_on == 1:
        hp.rhx2 = PlateHeatExchanger(
            name="Recuperator",
            source_fluid=hp.refrigerant,
            sink_fluid=hp.refrigerant,
            dT_c=hp.dT_rhx2_cold,
            model="fixed cold delta temperature",
        )
    if hp.rhx3_on == 1:
        hp.rhx3 = PlateHeatExchanger(
            name="Recuperator",
            source_fluid=hp.refrigerant,
            sink_fluid=hp.refrigerant,
            dT_c=hp.dT_rhx3_cold,
            model="fixed cold delta temperature",
        )

    hp.exp_valve = ExpansionValve(name="Expansion valve", fluid=hp.refrigerant)
    hp.evaporator = Evaporator(
        name="Evaporator",
        source_fluid=hp.source,
        sink_fluid=hp.refrigerant,
        model="saturated gas outlet evaporator",
    )
    hp.evap_sup = PlateHeatExchanger(
        name="Evaporator",
        source_fluid=hp.source,
        sink_fluid=hp.refrigerant,
        model="fixed cold temperature outlet",
    )

    return hp


def _get_RecuperativeHP_streams(hp: RecuperativeHP) -> RecuperativeHP:
    """Initialise thermodynamic streams that connect the unit operations."""
    temp_ref_out: Stream
    P_H: float
    P_L: float
    hp.gas_coolers = []
    hp.evaporators = []

    # Evaporator superheat section
    hp.evap_sup.source_in.props_in_ab(
        input_str="PT",
        prop1=0.5,
        prop2=hp.T_source_in,
        mass_flow=1e6,
        name="Source in",
        pp_instance="src",
    )
    hp.evap_sup.source_out.props_in_ab(
        input_str="PT",
        prop1=0.5,
        prop2=hp.T_source_in,
        mass_flow=1e6,
        name="Source out",
        pp_instance="src",
    )
    hp.Tref_evap_outlet = getattr(hp, "Tref_evap_outlet", hp.Tref_evap_inlet)
    hp.evap_sup.sink_in.props_in_ab(
        input_str="TQ",
        prop1=hp.Tref_evap_outlet,
        prop2=1.0,
        name="Evaporator sink outlet",
    )
    P_L = hp.evap_sup.sink_in.pressure
    hp.evap_sup.sink_out.props_in_ab(
        input_str="PT",
        prop1=P_L,
        prop2=(hp.T_source_in - hp.dTmin_rhx),
        name="Superheater sink outlet",
    )
    temp_ref_out = hp.evap_sup.sink_out

    hp.evaporators.append(hp.evap_sup)

    # Internal heat exchanger, if present
    if hp.ihx_on == 1:
        hp.ihx.sink_in = temp_ref_out
        hp.ihx.sink_out.props_in_ab(
            input_str="PT",
            prop1=P_L,
            prop2=(hp.T_gc_out - hp.dTmin_rhx),
            name="IHX sink outlet",
        )
        temp_ref_out = hp.ihx.sink_out

    # Plate heat exchanger 1 - cold-side, if present
    if hp.rhx1_on == 1:
        hp.rhx1.sink_in = temp_ref_out
        hp.rhx1.sink_out.props_in_ab(
            input_str="PT",
            prop1=P_L,
            prop2=(hp.rhx1.sink_in.temperature + hp.dT_rhx1_cold),
            name="rhx1 sink outlet",
        )
        temp_ref_out = hp.rhx1.sink_out

    def _get_discharge_pressure_ratio_estimate(
        P_H_est: float,
        P_L: float,
        T_eff: float,
        n_comp: int,
        T_tar: float,
        dT_min: float,
        fld: str,
    ) -> float:
        # Final discharge pressure
        if P_H_est > P_L:
            P_H = P_H_est
        else:
            P_H = P_L
            stream_in = Stream(fluid=fld, temperature=T_eff, pressure=P_L)
            stream_out = Stream(fluid=fld)
            stream_out_ideal = Stream(fluid=fld)
            T_out_min = hp.T_air_tar + hp.dT_GC_min_tar
            T_crit = temp_ref_out.T_crit

            try:
                # For critical point fluids, the cycle must be subcritical
                stream_pinch = Stream(fluid=fld)
                stream_cooled = Stream(fluid=fld)

                def _get_discharge_pressure_for_pinch_estimate(T):
                    if T > T_crit - 1:
                        stream_pinch.props_in_ab("TQ", T_crit - 10, 1)
                        D = stream_pinch.density
                        stream_pinch.props_in_ab("TD", T, D)
                    else:
                        stream_pinch.props_in_ab("TQ", T, 1)

                    stream_out_ideal.props_in_ab(
                        "PS", stream_pinch.pressure, stream_in.entropy
                    )
                    h_out = (
                        stream_in.enthalpy
                        + (stream_out_ideal.enthalpy - stream_in.enthalpy) / 0.7
                    )
                    stream_out.props_in_ab("PH", stream_pinch.pressure, h_out)
                    stream_cooled.props_in_ab("PT", stream_pinch.pressure, hp.T_gc_out)
                    err = (T - hp.dT_GC_min_tar - hp.T_gc_out) * (
                        stream_out.enthalpy - stream_cooled.enthalpy
                    ) - (stream_pinch.enthalpy - stream_cooled.enthalpy) * (
                        hp.T_air_tar - hp.T_gc_out
                    )
                    return err

                for multiplier in np.linspace(0, 10, 11):
                    try:
                        r = root_scalar(
                            f=_get_discharge_pressure_for_pinch_estimate,  # find the gas discharge pressure resulting in the target dtmin
                            method="brentq",
                            bracket=(
                                hp.T_gc_out,
                                max(
                                    temp_ref_out.T_crit - 1,
                                    hp.T_air_tar + hp.dT_GC_min_tar * multiplier,
                                ),
                            ),
                            x0=T_tar,
                            xtol=1e-3,
                        )
                        if r.converged:
                            P_H = stream_pinch.pressure
                            break
                    except Exception:
                        pass
                if P_H <= P_L:
                    P_H = min(20, P_L * 36)
            except Exception:
                try:
                    # Assume the high pressure = critical pressure
                    # Determine the minimum discharge temperature
                    P_H = temp_ref_out.P_crit
                    stream_out.props_in_ab("PS", P_H, stream_in.entropy)

                    if stream_out.temperature < T_out_min:
                        # For low critical point fluids, the end temperature is too low
                        # Find a pressure with sufficient temperature in the discharge
                        def _get_compressor_outlet_temperature(P_out):
                            stream_out.props_in_ab("PT", P_out, T_out_min)
                            stream_out_ideal.props_in_ab("PS", P_out, stream_in.entropy)
                            return (
                                stream_out_ideal.enthalpy - stream_in.enthalpy
                            ) / 0.7 - (stream_out.enthalpy - stream_in.enthalpy)

                        r = root_scalar(
                            f=_get_compressor_outlet_temperature,  # find the gas discharge pressure resulting in the target dtmin
                            method="brentq",
                            bracket=(P_L * 1.1, 20),
                            x0=P_H,
                            xtol=1e-3,
                        )
                        if r.converged:
                            P_H = r.root
                            pass
                        else:
                            P_H = 20
                except Exception:
                    P_H = 20

        return (P_H / P_L) ** (1 / n_comp)

    T_eff = temp_ref_out.temperature + sum(
        [hp.dT_rhx1_cold, hp.dT_rhx2_cold, hp.dT_rhx3_cold]
    )
    PR = _get_discharge_pressure_ratio_estimate(
        P_H_est=hp.P_H_est,
        P_L=P_L,
        T_eff=T_eff,
        n_comp=hp.n_comp,
        T_tar=hp.T_air_tar,
        dT_min=hp.dT_GC_min_tar,
        fld=hp.refrigerant,
    )

    # Compressor 1, always present
    hp.comp1.flow_in = temp_ref_out
    hp.comp1.flow_out.props_in_ab(
        input_str="PS",
        prop1=(P_L * PR),
        prop2=hp.comp1.flow_in.entropy,
        name="Compressor 1 outlet",
    )
    temp_ref_out = hp.comp1.flow_out
    hp.compressors = [hp.comp1]

    # Plate heat exchanger 2 - cold-side
    if hp.rhx2_on == 1:
        hp.rhx2.sink_in = temp_ref_out
        hp.rhx2.sink_out.props_in_ab(
            input_str="PT",
            prop1=hp.rhx2.sink_in.pressure,
            prop2=(hp.rhx2.sink_in.temperature + hp.dT_rhx2_cold),
            name="rhx2 sink outlet",
        )
        temp_ref_out = hp.rhx2.sink_out

    # Compressor 2
    if hp.comp2_on == 1:
        hp.comp2.flow_in = temp_ref_out
        hp.comp2.flow_out.props_in_ab(
            input_str="PS",
            prop1=(P_L * PR * PR),
            prop2=hp.comp2.flow_in.entropy,
            name="Compressor 2 outlet",
        )
        temp_ref_out = hp.comp2.flow_out
        hp.compressors.append(hp.comp2)

    # Plate heat exchanger 3 - cold-side
    if hp.rhx3_on == 1:
        hp.rhx3.sink_in = temp_ref_out
        hp.rhx3.sink_out.props_in_ab(
            input_str="PT",
            prop1=hp.rhx3.sink_in.pressure,
            prop2=(hp.rhx3.sink_in.temperature + hp.dT_rhx3_cold),
            name="rhx3 sink outlet",
        )
        temp_ref_out = hp.rhx3.sink_out

    # Compressor 3
    if hp.comp3_on == 1:
        hp.comp3.flow_in = temp_ref_out
        hp.comp3.flow_out.props_in_ab(
            input_str="PS",
            prop1=(P_L * PR * PR * PR),
            prop2=hp.comp3.flow_in.entropy,
            name="Compressor 3 outlet",
        )
        temp_ref_out = hp.comp3.flow_out
        hp.compressors.append(hp.comp3)

    hp.comp_discharge = temp_ref_out
    hp.comp_discharge.pressure_bnds = None
    P_H = temp_ref_out.pressure

    # Gas cooler 4 and Plate heat exchanger 3 - hot-sides
    hp.air_out = Stream(
        name="Air sink outlet",
        fluid=hp.sink,
        pp_instance="snk",
        pressure=0.101,
        temperature=hp.T_air_tar,
    )
    if hp.rhx3_on == 1:
        hp.gas_cooler4.sink_out = hp.air_out
        hp.gas_cooler4.sink_pinch.props_in_ab(
            input_str="PT",
            prop1=0.101,
            prop2=hp.T_air_tar,
            name="GC4 sink pinch point",
            pp_instance="snk",
        )
        hp.gas_cooler4.sink_in.props_in_ab(
            input_str="PT",
            prop1=0.101,
            prop2=hp.T_air_tar,
            name="GC4 sink inlet",
            pp_instance="snk",
        )
        temp_air_in = hp.gas_cooler4.sink_in

        hp.gas_coolers.append(hp.gas_cooler4)

        hp.gas_cooler4.source_in = temp_ref_out
        hp.gas_cooler4.source_pinch.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="GC4 source pinch point",
        )
        hp.gas_cooler4.source_out.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="GC4 source outlet",
        )

        hp.rhx3.source_in = hp.gas_cooler4.source_out
        hp.rhx3.source_out.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="rhx3 source outlet",
        )
        temp_ref_out = hp.rhx3.source_out
    else:
        temp_air_in = hp.air_out

    # Gas cooler 3 and Plate heat exchanger 2 - hot-sides
    if hp.rhx2_on == 1:
        hp.gas_cooler3.sink_out = temp_air_in
        hp.gas_cooler3.sink_pinch.props_in_ab(
            input_str="PT",
            prop1=0.101,
            prop2=hp.T_air_tar,
            name="GC3 sink pinch point",
            pp_instance="snk",
        )
        hp.gas_cooler3.sink_in.props_in_ab(
            input_str="PT",
            prop1=0.101,
            prop2=hp.T_air_tar,
            name="GC3 sink inlet",
            pp_instance="snk",
        )
        temp_air_in = hp.gas_cooler3.sink_in

        hp.gas_coolers.append(hp.gas_cooler3)

        hp.gas_cooler3.source_in = temp_ref_out
        hp.gas_cooler3.source_pinch.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="GC3 source pinch point",
        )
        hp.gas_cooler3.source_out.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="GC3 source outlet",
        )

        hp.rhx2.source_in = hp.gas_cooler3.source_out
        hp.rhx2.source_out.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="rhx2 source outlet",
        )
        temp_ref_out = hp.rhx2.source_out

    # Gas cooler 2 and Plate heat exchanger 1 - hot-sides
    if hp.rhx1_on == 1:
        hp.gas_cooler2.sink_out = temp_air_in
        hp.gas_cooler2.sink_pinch.props_in_ab(
            input_str="PT",
            prop1=0.101,
            prop2=hp.T_air_tar,
            name="GC2 sink pinch point",
            pp_instance="snk",
        )
        hp.gas_cooler2.sink_in.props_in_ab(
            input_str="PT",
            prop1=0.101,
            prop2=hp.T_air_tar,
            name="GC2 sink inlet",
            pp_instance="snk",
        )
        temp_air_in = hp.gas_cooler2.sink_in

        hp.gas_coolers.append(hp.gas_cooler2)

        hp.gas_cooler2.source_in = temp_ref_out
        hp.gas_cooler2.source_pinch.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="GC2 source pinch point",
        )
        hp.gas_cooler2.source_out.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="GC2 source outlet",
        )

        hp.rhx1.source_in = hp.gas_cooler2.source_out
        hp.rhx1.source_out.props_in_ab(
            input_str="PT",
            prop1=P_H,
            prop2=hp.comp_discharge.temperature,
            name="rhx1 source outlet",
        )
        temp_ref_out = hp.rhx1.source_out

    # Gas cooler 1
    hp.gas_cooler1.sink_out = temp_air_in
    hp.gas_cooler1.sink_pinch.props_in_ab(
        input_str="PT",
        prop1=0.101,
        prop2=hp.T_air_inlet,
        name="GC1 sink pinch point",
        pp_instance="snk",
    )
    hp.gas_cooler1.sink_in.props_in_ab(
        input_str="PT",
        prop1=0.101,
        prop2=hp.T_air_inlet,
        name="GC1 sink inlet",
        pp_instance="snk",
    )

    hp.gas_coolers.append(hp.gas_cooler1)

    hp.gas_cooler1.source_in = temp_ref_out
    hp.gas_cooler1.source_pinch.props_in_ab(
        input_str="PT", prop1=P_H, prop2=hp.T_gc_out, name="GC1 source pinch point"
    )
    hp.gas_cooler1.source_out.props_in_ab(
        input_str="PT", prop1=P_H, prop2=hp.T_gc_out, name="GC1 source outlet"
    )

    temp_ref_out = hp.gas_cooler1.source_out

    # Internal heat exchanger
    if hp.ihx_on == 1:
        hp.ihx.source_in = temp_ref_out
        hp.ihx.source_out.props_in_ab(
            input_str="PT", prop1=P_H, prop2=hp.T_gc_out, name="IHX source outlet"
        )
        temp_ref_out = hp.ihx.source_out

    # Expansion valve
    hp.exp_valve.flow_in = temp_ref_out
    hp.exp_valve.flow_out.props_in_ab(
        input_str="PH",
        prop1=P_L,
        prop2=hp.exp_valve.flow_in.enthalpy,
        name="Expansion valve outlet",
    )

    # Evaporator
    hp.evaporator.sink_in = hp.exp_valve.flow_out
    hp.evaporator.sink_out = hp.evap_sup.sink_in

    hp.evaporator.source_in.props_in_ab(
        input_str="PT",
        prop1=0.5,
        prop2=hp.T_source_in,
        mass_flow=1e6,
        name="Evaporator source inlet",
        pp_instance="src",
    )
    hp.evaporator.source_out.props_in_ab(
        input_str="PT",
        prop1=0.5,
        prop2=min(hp.evaporator.sink_in.temperature + hp.dTmin_rhx, hp.T_source_in),
        mass_flow=1e6,
        name="Evaporator source outlet",
        pp_instance="src",
    )

    hp.evaporators.append(hp.evaporator)

    return hp


def _get_RecuperativeHP_logic_blocks(hp: RecuperativeHP) -> RecuperativeHP:
    """Build algebraic relationships that enforce model constraints."""
    # Set special relationships between points in the cycle
    if hp.rhx3_on == 1:
        hp.set_rhx3_dT_cold = SetRelationship(
            source_obj=hp.rhx3.sink_in,
            target_obj=[hp.rhx3.sink_out],
            c=hp.dT_rhx3_cold,
            prop="temperature",
            prop2="pressure",
        )
        hp.set_gc4_hot_out = SetRelationship(
            source_obj=hp.rhx3.sink_out,
            target_obj=[hp.rhx3.source_in],
            c=hp.dTmin_rhx,
            prop="temperature",
            prop2="pressure",
        )

    if hp.rhx2_on == 1:
        hp.set_rhx2_dT_cold = SetRelationship(
            source_obj=hp.rhx2.sink_in,
            target_obj=[hp.rhx2.sink_out],
            c=hp.dT_rhx2_cold,
            prop="temperature",
            prop2="pressure",
        )
        hp.set_gc3_hot_out = SetRelationship(
            source_obj=hp.rhx2.sink_out,
            target_obj=[hp.rhx2.source_in],
            c=hp.dTmin_rhx,
            prop="temperature",
            prop2="pressure",
        )

    if hp.rhx1_on == 1:
        hp.set_rhx1_dT_cold = SetRelationship(
            source_obj=hp.rhx1.sink_in,
            target_obj=[hp.rhx1.sink_out],
            c=hp.dT_rhx1_cold,
            prop="temperature",
            prop2="pressure",
        )
        hp.set_gc2_hot_out = SetRelationship(
            source_obj=hp.rhx1.sink_out,
            target_obj=[hp.rhx1.source_in],
            c=hp.dTmin_rhx,
            prop="temperature",
            prop2="pressure",
        )

    if hp.n_comp == 2:
        hp.set_P_ints = CustomRelationship(
            source_obj=[
                hp.comp1.flow_in,
                hp.comp_discharge,
                hp.comp2_power_frac,
                hp.comp3_power_frac,
            ],
            target_obj=[hp.comp1.flow_out, hp.comp_discharge],
            relation=1,
            tol=hp.tol,
        )
    elif hp.n_comp == 3:
        hp.set_P_ints = CustomRelationship(
            source_obj=[
                hp.comp1.flow_in,
                hp.comp_discharge,
                hp.comp2_power_frac,
                hp.comp3_power_frac,
            ],
            target_obj=[hp.comp1.flow_out, hp.comp2.flow_out],
            relation=1,
            tol=hp.tol,
        )

    return hp


def _get_RecuperativeHP_solver_pathway(hp: RecuperativeHP) -> RecuperativeHP:
    """Persist an ordered sequence of operations for the cycle solver."""
    # Solution pathway set-up
    hp.operations = [hp.evap_sup]

    if hp.n_comp > 1:
        hp.operations.append(hp.set_P_ints)
    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)
    if hp.rhx1_on == 1:
        hp.operations.append(hp.set_rhx1_dT_cold)
        hp.operations.append(hp.set_gc2_hot_out)

    hp.operations.append(hp.comp1)

    if hp.rhx2_on == 1:
        hp.operations.append(hp.set_rhx2_dT_cold)
        hp.operations.append(hp.set_gc3_hot_out)
    if hp.comp2_on == 1:
        hp.operations.append(hp.comp2)
    if hp.rhx3_on == 1:
        hp.operations.append(hp.set_rhx3_dT_cold)
        hp.operations.append(hp.set_gc4_hot_out)
    if hp.comp3_on == 1:
        hp.operations.append(hp.comp3)

    if hp.rhx3_on == 1:
        hp.operations.append(hp.gas_cooler4)
        hp.operations.append(hp.rhx3)
    if hp.rhx2_on == 1:
        hp.operations.append(hp.gas_cooler3)
        hp.operations.append(hp.rhx2)
    if hp.rhx1_on == 1:
        hp.operations.append(hp.gas_cooler2)
        hp.operations.append(hp.rhx1)

    hp.operations.append(hp.gas_cooler1)

    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)

    hp.operations.append(hp.exp_valve)
    hp.operations.append(hp.evaporator)

    return hp
