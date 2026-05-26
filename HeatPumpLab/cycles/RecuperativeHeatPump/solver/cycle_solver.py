"""Numerical routines to solve Recuperative Heat Pump cycles."""

from __future__ import annotations

import numpy as np
from scipy.optimize import root_scalar, minimize, Bounds
from math import log, exp
from typing import TYPE_CHECKING

from ....utils.readWriteOutput import *
from ....utils.graph import *

if TYPE_CHECKING:
    from .. import BaseHeatPump
    from ....components import GasCooler
    from ....streams import Stream

__all__ = ["solve_heat_pump_cycle"]

#######################################################################################################
# Public API
#######################################################################################################


def solve_heat_pump_cycle(hp: BaseHeatPump) -> BaseHeatPump:
    """Solve a configured heat pump, mutating the model and returning it."""
    try:
        hp = _solve_case(hp)

    except:
        if hp.debug:
            hp.get_cycle_plot()
        _print_failed_to_solve_report(hp)
        print(
            "Heat pump solver failed with an exception. The failed case has been recorded."
        )
        hp.converged = False
    return hp


#######################################################################################################
# Helper functions
#######################################################################################################


def _solve_case(hp: BaseHeatPump) -> BaseHeatPump:
    """Drive the high-level solution loop for a heat pump instance."""
    hp = _find_hp_with_min_discharge_pressure(hp)
    if hp.debug:
        hp.get_cycle_plot()

    dt_offset = hp.evaporator.sink_in.temperature - hp.Tref_evap_inlet
    if abs(dt_offset) > 0.1 and _check_hp_model_convergance(hp):
        while (
            abs(
                _find_evap_inlet_temperature_offset(
                    hp.evaporator.sink_out.temperature - dt_offset, hp
                )
            )
            > 1e-6
        ):
            dt_offset = hp.evaporator.sink_in.temperature - hp.Tref_evap_inlet

        hp.Tref_evap_outlet = hp.Tref_evap_inlet + dt_offset
        hp.evaporator.sink_out.props_in_ab("TQ", hp.Tref_evap_outlet, 1.0)
        hp = _find_hp_with_min_discharge_pressure(hp)

    hp = _get_metrics(hp)
    return hp


def _check_hp_model_convergance(hp: BaseHeatPump) -> bool:
    """Return True when the cycle has converged and yielded positive work."""
    return (
        True
        if (
            sum([c.power for c in hp.compressors]) > 0.0
            and hp.converged == True
            and hp.feasible == True
        )
        else False
    )


def _get_metrics(hp: BaseHeatPump) -> BaseHeatPump:
    """Populate secondary performance metrics on the heat pump object."""

    def _get_COP(hp: BaseHeatPump) -> float:
        """Coefficient of performance of the solved case."""
        return (
            sum([g.duty for g in hp.gas_coolers])
            / sum([c.power for c in hp.compressors])
            if hp.converged
            else 0.0
        )

    def _get_SEC(hp: BaseHeatPump) -> float:
        """Specific electricity consumption for the solved cycle."""
        return (
            sum([c.power for c in hp.compressors])
            / (
                sum([g.duty for g in hp.gas_coolers])
                - sum([c.power for c in hp.compressors])
            )
            if hp.converged
            else np.inf
        )

    def _get_VHC(hp: BaseHeatPump) -> float:
        """Volumetric heating capacity delivered at the gas coolers."""
        return (
            sum([g.duty for g in hp.gas_coolers]) / hp.comp1.flow_in.vol_flow
            if hp.converged
            else 0.0
        )

    def _get_eDTL(hp: BaseHeatPump) -> float:
        """Exergy temperature lift between source and sink."""
        eT_h_ave = sum([g.duty for g in hp.gas_coolers]) / sum(
            [g.delta_entropy_src for g in hp.gas_coolers]
        )
        eT_c_ave = sum([e.duty for e in hp.evaporators]) / sum(
            [e.delta_entropy_snk for e in hp.evaporators]
        )
        return eT_h_ave - eT_c_ave if hp.converged else 0.0

    def _get_SLE(hp: BaseHeatPump) -> float:
        """Specific lift efficiency derived from SEC and eDTL."""
        return hp.SEC / hp.eDTL * 1000 if hp.converged else 0.0

    def _get_LMTD_gcs(hp: BaseHeatPump):
        """Overall log mean temperature difference across gas coolers."""
        return sum([g.duty * g.get_LMTD_ave() for g in hp.gas_coolers]) / sum(
            [g.duty for g in hp.gas_coolers]
        )

    def _get_hp_cycle_type(hp: BaseHeatPump) -> str:
        """Identify whether the cycle operated transcritically or subcritically."""
        return (
            "Transcritical"
            if hp.comp_discharge.pressure > hp.comp_discharge.P_crit
            else "Subcritical"
        )

    hp.converged = _check_hp_model_convergance(hp)
    hp.COP = _get_COP(hp)
    hp.SEC = _get_SEC(hp)
    hp.VHC = _get_VHC(hp)
    hp.eDTL = _get_eDTL(hp)
    hp.SLE = _get_SLE(hp)
    hp.LMTD_gcs = _get_LMTD_gcs(hp)
    hp.cycle_type = _get_hp_cycle_type(hp)
    return hp


def _find_evap_inlet_temperature_offset(
    Tref_evap_outlet: float, hp: BaseHeatPump
) -> float:
    """Return the mismatch between target and simulated evaporator inlet temperatures."""
    hp.evaporator.sink_out.props_in_ab("TQ", Tref_evap_outlet, 1.0)
    hp = _find_hp_with_min_discharge_pressure(hp)
    return hp.evaporator.sink_in.temperature - hp.Tref_evap_inlet


def _find_hp_with_min_discharge_pressure(hp: BaseHeatPump) -> BaseHeatPump:
    """Search for the discharge pressure that meets the gas-cooler pinch target."""
    bnds, x0, hp.feasible, _ = _get_discharge_pressure_bounds(hp)
    hp.comp_discharge.pressure_bnds = bnds
    if hp.feasible == True:
        r = root_scalar(
            f=_find_current_gc_dtmin_offset,  # find the gas discharge pressure resulting in the target dtmin
            method="brentq",
            bracket=bnds,
            x0=x0,
            xtol=hp.tol,
            args=(hp),
        )
        hp = _do_sequential_iterative_cycle_analysis(r.root, hp)
        hp = _update_pinch_props_for_gcs(P_comp_out=r.root, hp=hp)
        hp.converged = r.converged
    else:
        _print_failed_to_solve_report(hp)
        raise Exception(
            "Infeasible pressure bounds. The failed case has been recorded."
        )
    return hp


def _get_discharge_pressure_bounds(
    hp: BaseHeatPump,
) -> list[tuple[float, float], float, bool, tuple[float, float]]:
    """Determine feasible discharge pressure bounds and initial guess."""
    feasible: bool = False
    P_lb: float
    P_ub: float
    dt_neg: float
    dt_plus: float
    P_i: float

    def _recursive_find_bnds(
        P_suction: float,
        P_ub: float,
        dt_plus: float,
        P_lb: float,
        dt_neg: float,
        step: float,
        i: int,
    ):
        """Expand bounds until the pinch temperature difference changes sign."""
        if dt_plus < 0.0 and dt_neg < 0.0 and i < 1000:
            dt_neg = dt_plus
            P = P_lb = P_ub
            P_ub = max(P_suction + (P - P_suction) * step, P_suction)
            dt_plus = _find_current_gc_dtmin_offset(P_ub, hp)
            return _recursive_find_bnds(
                P_suction, P_ub, dt_plus, P_lb, dt_neg, step, i + 1
            )
        elif dt_plus > 0.0 and dt_neg > 0.0 and i < 1000:
            dt_plus = dt_neg
            P = P_ub = P_lb
            P_lb = max(P_suction + (P - P_suction) / step, P_suction)
            dt_neg = _find_current_gc_dtmin_offset(P_lb, hp)
            return _recursive_find_bnds(
                P_suction, P_ub, dt_plus, P_lb, dt_neg, step, i + 1
            )
        else:
            return P_ub, dt_plus, P_lb, dt_neg

    # Positive dtmin is P_ub and negative dtmin is P_lb
    if hp.comp_discharge.pressure_bnds is not None:
        if hp.compute_jac == False:
            P_lb, P_ub = hp.comp_discharge.pressure_bnds
        else:
            P_lb, P_ub = (
                hp.comp_discharge.pressure * 0.995,
                hp.comp_discharge.pressure * 1.005,
            )
        dt_neg, dt_plus = [_find_current_gc_dtmin_offset(P, hp) for P in (P_lb, P_ub)]
        pass
    else:
        dt_plus = dt_neg = _find_current_gc_dtmin_offset(hp.comp_discharge.pressure, hp)
        P_ub = P_lb = hp.comp_discharge.pressure
        P_ub, dt_plus, P_lb, dt_neg = _recursive_find_bnds(
            hp.comp1.flow_in.pressure, P_ub, dt_plus, P_lb, dt_neg, 1.1, 0
        )
        pass

    if hp.comp_discharge.pressure < P_lb or hp.comp_discharge.pressure > P_ub:
        P_i = exp(log(P_lb) - dt_neg / (dt_plus - dt_neg) * (log(P_ub) - log(P_lb)))
    else:
        P_i = hp.comp_discharge.pressure

    if abs(dt_plus * dt_neg) == -dt_plus * dt_neg and abs(dt_plus * dt_neg) > 0.0:
        feasible = True

    return (P_lb, P_ub), P_i, feasible, (dt_plus, dt_neg)


def _find_current_gc_dtmin_offset(
    P_comp_out: float, hp: BaseHeatPump, return_abs: bool = False
) -> float:
    """Compute the deviation between current and target gas-cooler pinch."""
    dtmin_ls: list[list[float, float]]

    dtmin_ls = _get_all_possible_gc_dtmin_pos(P_comp_out, hp)

    dtmin_ls.pop(-1)  # Ignore the cold-side dtmin of the inlet air to the gas cooler
    if (
        dtmin_ls[-1][1] < 1e-3
    ):  # Ignore the mid exchanger dtmin of the inlet air to the gas cooler if it is the cold-side
        dtmin_ls.pop(-1)

    if not return_abs:
        return min(dtmin_ls)[0] - hp.dT_GC_min_tar
    else:
        return abs(min(dtmin_ls)[0] - hp.dT_GC_min_tar)


def _get_all_possible_gc_dtmin_pos(
    P_comp_out: float, hp: BaseHeatPump
) -> list[list[float, float]]:
    """Enumerate pinch temperature differences for all gas cooler sections."""
    dtmin_ls: list[float]
    gc: GasCooler

    hp = _do_sequential_iterative_cycle_analysis(P_comp_out, hp)
    dtmin_ls = []
    for gc in hp.gas_coolers:
        dt_hot = gc.source_in.temperature - gc.sink_out.temperature
        dt_cold = gc.source_out.temperature - gc.sink_in.temperature
        dtmin_ls.append([dt_hot, 1.0, True])  # Pinch on the hot side of the gc
        dtmin_ls.append(_find_gc_dt_approaches(gc))  # Pinch in the middle of the gc
        dtmin_ls.append([dt_cold, 0.0, True])  # Pinch on the cold side of the gc

    if hp.debug:
        plot_driving_force(hp.gas_coolers, dtmin=hp.dT_GC_min_tar)
    return dtmin_ls


def _do_sequential_iterative_cycle_analysis(
    P_comp_out: float, hp: BaseHeatPump
) -> BaseHeatPump:
    """Sequentially solve each component until the coupled cycle converges."""
    hp.comp_discharge.pressure = P_comp_out
    i: int = 0
    h: float
    P: float
    W: float
    Q: float
    most_neg_residual: float = 0.0

    while (i == 0 or most_neg_residual < 0.0) and i < 10:  # Check for model convergance
        # Saves pressure and temperature properties to check for solution convergance
        h = hp.exp_valve.flow_in.enthalpy
        P = hp.comp_discharge.pressure - hp.gas_cooler1.source_out.pressure
        W = sum([c.power for c in hp.compressors])
        Q = sum([g.duty for g in hp.gas_coolers])
        i += 1

        # Solve the performance of each component and update each special relationship
        [obj.solve() for obj in hp.operations]

        most_neg_residual = min(
            hp.tol * h - abs(hp.exp_valve.flow_in.enthalpy - h),
            hp.tol * P
            - abs(hp.comp_discharge.pressure - hp.gas_cooler1.source_out.pressure - P),
            hp.tol * W - abs(sum([c.power for c in hp.compressors]) - W),
            hp.tol * Q - abs(sum([g.duty for g in hp.gas_coolers]) - Q),
        )
    hp = _set_gas_coolers_sink_temperatures(hp)
    if hp.debug:
        hp.get_cycle_plot(x_axis="s", y_axis="T")
        hp.get_cycle_plot(x_axis="h", y_axis="P")
        hp.get_cycle_plot(x_axis="h", y_axis="T")
        pass

    return hp


def _set_gas_coolers_sink_temperatures(hp: BaseHeatPump) -> BaseHeatPump:
    """Back-calculate sink stream temperatures and mass flows for gas coolers."""
    gc: GasCooler

    hp.gas_coolers[0].sink_out.temperature = hp.T_air_tar
    hp.gas_coolers[0].sink_out.props_in_mPT()

    m_air = sum([gc.duty for gc in hp.gas_coolers]) / (
        hp.gas_coolers[0].sink_out.enthalpy - hp.gas_coolers[-1].sink_in.enthalpy
    )
    hp.gas_coolers[0].sink_out.mass_flow = m_air
    hp.gas_coolers[0].sink_out.heat_flow = m_air * hp.gas_coolers[0].sink_out.enthalpy

    for gc in hp.gas_coolers:
        gc.sink_pinch.mass_flow = m_air
        gc.sink_in.mass_flow = m_air
        gc.sink_in.enthalpy = gc.sink_out.enthalpy - gc.duty / m_air
        gc.sink_in.props_in_mPH()

    return hp


def _find_gc_dt_approaches(gc: GasCooler) -> list[float, float]:
    """Optimise for the minimum pinch temperature difference inside a gas cooler."""
    gc.source_pinch.pressure = gc.source_out.pressure
    gc.sink_pinch.pressure = gc.sink_in.pressure

    bnds = Bounds(lb=0.0, ub=gc.duty)

    r0 = minimize(
        fun=_calc_pinch_dtmin_and_jac,
        args=(gc),
        method="L-BFGS-B",
        x0=[gc.duty * 0.25],
        bounds=bnds,
        jac=True,
    )  # TNC L-BFGS-B
    r1 = minimize(
        fun=_calc_pinch_dtmin_and_jac,
        args=(gc),
        method="L-BFGS-B",
        x0=[gc.duty * 0.50],
        bounds=bnds,
        jac=True,
    )  # TNC L-BFGS-B
    r2 = minimize(
        fun=_calc_pinch_dtmin_and_jac,
        args=(gc),
        method="L-BFGS-B",
        x0=[gc.duty * 0.75],
        bounds=bnds,
        jac=True,
    )  # TNC L-BFGS-B

    res = [r0, r1, r2]
    j = np.argmin([r["fun"] for r in res])
    if res[j]["x"] < 0.001 * gc.duty or res[j]["x"] > 0.999 * gc.duty:
        res.pop(j)
        j = np.argmin([r["fun"] for r in res])
        if res[j]["x"] < 0.001 * gc.duty or res[j]["x"] > 0.999 * gc.duty:
            res.pop(j)
            j = np.argmin([r["fun"] for r in res])

    r = res[j]

    return [r["fun"], r["x"][0] / gc.duty, r["success"]]


def _calc_pinch_dtmin_and_jac(Q_arr: np.ndarray, gc: GasCooler):
    """Return pinch delta-T and Jacobian for a given heat transfer split."""
    gc = _calc_pinch_props_at_Q(Q_arr[0], gc)
    return (
        _calc_dtmin_at_given_point(gc.source_pinch, gc.sink_pinch),  # objective
        _calc_jac_at_given_point(gc.source_pinch, gc.sink_pinch),  # Jacobian
    )


def _update_pinch_props_for_gcs(P_comp_out: float, hp: BaseHeatPump) -> BaseHeatPump:
    """Refresh pinch stream properties after the discharge pressure changes."""
    dtmin_ls = _get_all_possible_gc_dtmin_pos(P_comp_out, hp)
    i = 1
    for gc in hp.gas_coolers:
        gc = _calc_pinch_props_at_Q(Q=dtmin_ls[i][1], gc=gc)
        i += 3
    return hp


def _calc_pinch_props_at_Q(Q: float, gc: GasCooler) -> GasCooler:
    """Calculate pinch stream thermodynamic properties at a given load."""
    gc.source_pinch.enthalpy = gc.source_out.enthalpy + Q / gc.source_out.mass_flow
    gc.source_pinch.props_in_mPH()
    gc.sink_pinch.enthalpy = gc.sink_in.enthalpy + Q / gc.sink_in.mass_flow
    gc.sink_pinch.props_in_mPH()
    return gc


def _calc_dtmin_at_given_point(source: Stream, sink: Stream) -> float:
    """Delta-T between matched source and sink pinch points."""
    return source.temperature - sink.temperature


def _calc_jac_at_given_point(source: Stream, sink: Stream) -> float:
    """Sensitivity of delta-T to the heat load for the pinch match."""
    return 1 / (source.heat_capacity * source.mass_flow) - 1 / (
        sink.heat_capacity * sink.mass_flow
    )


def _print_failed_to_solve_report(
    hp: BaseHeatPump, filename: str = "Cases-Unsolved"
) -> None:
    """Persist the failed case to disk and emit a console warning."""
    print("Heat pump model failed to solve.")
    hp.get_report()
    print_result_to_csv(hp.output, file_name=filename)
