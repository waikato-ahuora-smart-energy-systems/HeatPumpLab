"""Numerical routines to solve Recuperative Heat Pump cycles."""

from __future__ import annotations

import numpy as np
from scipy.optimize import root_scalar, minimize_scalar
from math import log, exp
from typing import TYPE_CHECKING, TypeAlias

from ....utils.readWriteOutput import print_result_to_csv
from ....utils.graph import plot_driving_force
from ....utils.property_cache import property_call_cache
from ....utils.solve_profile import profiled_section, profile_section

if TYPE_CHECKING:
    from .. import BaseHeatPump
    from ....components import GasCooler, Stream

__all__ = ["solve_heat_pump_cycle"]

PinchCandidate: TypeAlias = tuple[float, float, bool]
_PRESSURE_FALLBACK_SAMPLES = 25


class _CycleConstraintError(Exception):
    """Expected cycle infeasibility with physical constraint errors."""

    def __init__(self, message: str, constraint_errors: dict[str, float]) -> None:
        super().__init__(message)
        self.constraint_errors = {
            str(name): max(0.0, float(error))
            for name, error in constraint_errors.items()
            if np.isfinite(float(error))
        }

#######################################################################################################
# Public API
#######################################################################################################


def solve_heat_pump_cycle(hp: BaseHeatPump) -> BaseHeatPump:
    """Solve a configured heat pump, mutating the model and returning it."""
    hp.constraint_errors = {}
    try:
        hp = _solve_case(hp)

    except Exception as exc:
        hp.converged = False
        hp.feasible = False
        hp.solve_error = str(exc)
        if isinstance(exc, _CycleConstraintError):
            hp.constraint_errors = dict(exc.constraint_errors)
        if hp.debug:
            hp.get_cycle_plot()
        if getattr(hp, "record_failed_solve", True):
            _print_failed_to_solve_report(hp)
            disposition = "The failed case has been recorded."
        else:
            disposition = "The failed case was not recorded."
        print(f"Heat pump solver failed with an exception. {disposition} Error: {exc}")
    return hp


#######################################################################################################
# Helper functions
#######################################################################################################


@profiled_section("solve_case")
def _solve_case(hp: BaseHeatPump) -> BaseHeatPump:
    """Drive the high-level solution loop for a heat pump instance."""
    with property_call_cache():
        if hasattr(hp, "Tref_evap_outlet"):
            _validate_superheater_outlet_temperature(hp, hp.Tref_evap_outlet)
        fixed_pressure_objective = _uses_fixed_discharge_pressure_objective(hp)
        if fixed_pressure_objective:
            hp = _solve_cycle_at_discharge_pressure(hp.P_H_est, hp)
        else:
            hp = _find_hp_with_min_discharge_pressure(hp)

        if hp.debug:
            hp.get_cycle_plot()

        if not fixed_pressure_objective:
            dt_offset = _get_evaporator_inlet_temperature_offset(hp)
            if abs(dt_offset) > 0.1 and _check_hp_model_convergence(hp):
                hp = _align_evaporator_inlet_temperature(hp)

        hp = (
            _get_objective_metrics(hp)
            if fixed_pressure_objective
            else _get_metrics(hp)
        )
        return hp


def _uses_fixed_discharge_pressure_objective(hp: BaseHeatPump) -> bool:
    """Return True when this solve should evaluate at the seeded pressure only."""
    return (
        bool(getattr(hp, "fixed_discharge_pressure_objective", False))
        and getattr(hp, "P_H_est", 0.0) > 0.0
    )


def _check_hp_model_convergence(hp: BaseHeatPump) -> bool:
    """Return True when the cycle has converged and yielded positive work."""
    return (
        sum([c.power for c in hp.compressors]) > 0.0
        and hp.converged
        and hp.feasible
    )


def _get_metrics(hp: BaseHeatPump) -> BaseHeatPump:
    """Populate secondary performance metrics on the heat pump object."""
    with profile_section("metrics"):
        return _get_metrics_impl(hp)


def _get_objective_metrics(hp: BaseHeatPump) -> BaseHeatPump:
    """Populate only the metrics needed by objective evaluations."""
    with profile_section("objective_metrics"):
        hp.converged = _check_hp_model_convergence(hp)
        hp.COP = _get_COP(hp)
        hp.SEC = _get_SEC(hp)
        return hp


def _get_COP(hp: BaseHeatPump) -> float:
    """Coefficient of performance of the solved case."""
    return (
        sum([g.duty for g in hp.gas_coolers]) / sum([c.power for c in hp.compressors])
        if hp.converged
        else 0.0
    )


def _get_SEC(hp: BaseHeatPump) -> float:
    """Specific electricity consumption for the solved cycle."""
    return (
        sum([c.power for c in hp.compressors])
        / (sum([g.duty for g in hp.gas_coolers]) - sum([c.power for c in hp.compressors]))
        if hp.converged
        else np.inf
    )


def _get_metrics_impl(hp: BaseHeatPump) -> BaseHeatPump:
    """Populate secondary performance metrics on the heat pump object."""

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

    hp.converged = _check_hp_model_convergence(hp)
    hp.COP = _get_COP(hp)
    hp.SEC = _get_SEC(hp)
    hp.VHC = _get_VHC(hp)
    hp.eDTL = _get_eDTL(hp)
    hp.SLE = _get_SLE(hp)
    hp.LMTD_gcs = _get_LMTD_gcs(hp)
    hp.cycle_type = _get_hp_cycle_type(hp)
    return hp


def _get_evaporator_inlet_temperature_offset(hp: BaseHeatPump) -> float:
    """Return simulated minus target evaporator inlet temperature."""
    return hp.evaporator.sink_in.temperature - hp.Tref_evap_inlet


def _get_evaporator_alignment_tolerance(hp: BaseHeatPump) -> float:
    """Return the practical temperature tolerance for evaporator alignment."""
    return max(hp.tol, 1e-4)


def _align_evaporator_inlet_temperature(hp: BaseHeatPump) -> BaseHeatPump:
    """Maximise suction pressure without exceeding the evaporator inlet target."""
    target_tol = _get_evaporator_alignment_tolerance(hp)
    maximum_outlet_temperature = _get_maximum_evaporator_outlet_temperature(hp)
    maximum_feasible_outlet = (
        maximum_outlet_temperature - target_tol
        if maximum_outlet_temperature is not None
        else None
    )
    for _ in range(25):
        dt_offset = _get_evaporator_inlet_temperature_offset(hp)
        if abs(dt_offset) <= target_tol:
            hp.Tref_evap_outlet = hp.evaporator.sink_out.temperature
            return hp

        if (
            dt_offset < 0.0
            and maximum_feasible_outlet is not None
            and hp.evaporator.sink_out.temperature
            >= maximum_feasible_outlet - target_tol
        ):
            hp.Tref_evap_outlet = hp.evaporator.sink_out.temperature
            return hp

        proposed_temperature = hp.evaporator.sink_out.temperature - dt_offset
        if maximum_feasible_outlet is not None:
            proposed_temperature = min(
                proposed_temperature,
                maximum_feasible_outlet,
            )
        _validate_superheater_outlet_temperature(hp, proposed_temperature)
        hp.evaporator.sink_out.props_in_ab("TQ", proposed_temperature, 1.0)
        hp.Tref_evap_outlet = hp.evaporator.sink_out.temperature
        _reset_discharge_pressure_bounds(hp)
        hp = _find_hp_with_min_discharge_pressure(hp)

    raise Exception("Evaporator inlet temperature correction failed to converge.")


def _validate_superheater_outlet_temperature(
    hp: BaseHeatPump, evaporator_outlet_temperature: float
) -> None:
    """Reject an evaporator outlet that would reverse superheater duty."""
    try:
        maximum_temperature = hp.evap_sup.sink_out.temperature
    except AttributeError:
        return
    violation = max(0.0, evaporator_outlet_temperature - maximum_temperature)
    if violation > _get_evaporator_alignment_tolerance(hp):
        raise _CycleConstraintError(
            "Evaporator outlet temperature would reverse superheater heat transfer.",
            {"superheater_temperature_reversal": violation},
        )


def _get_maximum_evaporator_outlet_temperature(hp: BaseHeatPump) -> float | None:
    """Return the superheater-limited evaporator dew-point temperature."""
    try:
        return float(hp.evap_sup.sink_out.temperature)
    except AttributeError:
        return None


def _reset_discharge_pressure_bounds(hp: BaseHeatPump) -> None:
    """Invalidate cached pressure bounds after thermodynamic state changes."""
    hp.comp_discharge.pressure_bnds = None


def _get_pressure_target_tolerance(hp: BaseHeatPump) -> float:
    """Return the absolute pinch residual tolerated for pressure targeting."""
    return max(10.0 * hp.tol, 0.05)


def _is_multi_component_working_fluid(hp: BaseHeatPump) -> bool:
    """Return True when the working-fluid string names a mixture."""
    return ";" in getattr(hp, "refrigerant", "")


def _bounds_cross_target(dt_bounds: tuple[float, float]) -> bool:
    """Return True when the pressure-bound residuals bracket the target."""
    dt_upper, dt_lower = dt_bounds
    return dt_upper * dt_lower < 0.0


def _select_pressure_bound_within_target(
    bnds: tuple[float, float], dt_bounds: tuple[float, float], hp: BaseHeatPump
) -> tuple[float, bool]:
    """Return the lowest pressure bound that already satisfies the pinch target."""
    dt_upper, dt_lower = dt_bounds
    target_tol = _get_pressure_target_tolerance(hp)
    if abs(dt_lower) <= target_tol:
        return bnds[0], True
    if abs(dt_upper) <= target_tol:
        return bnds[1], True
    return bnds[0], False


def _find_best_pressure_within_target(
    bnds: tuple[float, float], hp: BaseHeatPump
) -> tuple[float, bool]:
    """Minimise the pressure-target residual when a clean root is unavailable."""
    if bnds[0] == bnds[1]:
        return bnds[0], False

    r = minimize_scalar(
        fun=lambda pressure: _find_current_gc_dtmin_offset(
            pressure, hp, return_abs=True
        ),
        method="bounded",
        bounds=bnds,
        options={"xatol": hp.tol},
    )
    return r.x, r.success and r.fun <= _get_pressure_target_tolerance(hp)


def _find_pressure_from_bracket(
    bnds: tuple[float, float], x0: float, hp: BaseHeatPump
) -> tuple[float, bool]:
    """Find the gas-cooler pinch target inside a sign-changing pressure bracket."""
    r = root_scalar(
        f=_find_current_gc_dtmin_offset,
        method="brentq",
        bracket=bnds,
        x0=x0,
        xtol=hp.tol,
        maxiter=50,
        args=(hp,),
    )
    return r.root, r.converged


def _find_lowest_pressure_within_target(
    bnds: tuple[float, float], dt_bounds: tuple[float, float], hp: BaseHeatPump
) -> tuple[float, bool]:
    """Find the lowest bracketed pressure that satisfies the pinch tolerance."""
    p_low, p_high = bnds
    dt_high, dt_low = dt_bounds
    target_tol = _get_pressure_target_tolerance(hp)

    if abs(dt_low) <= target_tol:
        return p_low, True
    if not (dt_low < 0.0 < dt_high):
        return p_low, False

    best_pressure = p_high if abs(dt_high) <= target_tol else None
    best_dt = dt_high
    for _ in range(50):
        p_mid = exp((log(p_low) + log(p_high)) / 2.0)
        dt_mid = _find_current_gc_dtmin_offset(p_mid, hp)
        if dt_mid >= -target_tol:
            p_high = p_mid
            dt_high = dt_mid
            best_pressure = p_mid
            best_dt = dt_mid
        else:
            p_low = p_mid
            dt_low = dt_mid

        if (p_high - p_low) / p_high <= max(hp.tol, 1e-7):
            break

    if best_pressure is not None and abs(best_dt) <= target_tol:
        return best_pressure, True
    return p_low, False


def _find_pressure_for_gc_target(
    bnds: tuple[float, float],
    x0: float,
    dt_bounds: tuple[float, float],
    hp: BaseHeatPump,
) -> tuple[float, bool]:
    """Select a discharge pressure that satisfies the gas-cooler pinch target.

    Multi-component fluids can make the pressure residual flat or branch-sensitive
    near saturation glides. Prefer a true bracketed root, but keep bounded
    endpoint/minimisation fallbacks for those near-target cases.
    """
    endpoint_pressure, endpoint_converged = _select_pressure_bound_within_target(
        bnds, dt_bounds, hp
    )

    if _bounds_cross_target(dt_bounds):
        if _is_multi_component_working_fluid(hp):
            pressure, converged = _find_lowest_pressure_within_target(
                bnds, dt_bounds, hp
            )
            if converged:
                return pressure, True

        try:
            return _find_pressure_from_bracket(bnds, x0, hp)
        except ValueError:
            if endpoint_converged:
                return endpoint_pressure, True

    if endpoint_converged:
        return endpoint_pressure, True

    return _find_best_pressure_within_target(bnds, hp)


def _find_unbracketed_multicomponent_pressure(
    hp: BaseHeatPump,
    *,
    pressure_start: float,
) -> tuple[float, bool, float]:
    """Search a sampled pressure domain when a mixture target is not bracketed."""
    pressure_suction = float(hp.comp1.flow_in.pressure)
    pressure_lower = pressure_suction * (1.0 + max(float(hp.tol), 1e-9))
    pressure_upper = float(getattr(hp.comp_discharge, "P_max", pressure_start))
    if pressure_upper <= pressure_lower:
        return pressure_start, False, 1.0

    pressures = np.geomspace(
        pressure_lower,
        pressure_upper,
        num=_PRESSURE_FALLBACK_SAMPLES,
    )
    pressures = np.unique(
        np.append(pressures, np.clip(pressure_start, pressure_lower, pressure_upper))
    )
    samples: list[tuple[float, float]] = []
    for pressure in pressures:
        try:
            residual = float(_find_current_gc_dtmin_offset(float(pressure), hp))
        except Exception:
            continue
        if np.isfinite(residual):
            samples.append((float(pressure), residual))

    if not samples:
        return pressure_start, False, 1.0

    samples.sort(key=lambda sample: sample[0])
    target_tol = _get_pressure_target_tolerance(hp)
    candidates: list[float] = []

    for index, (pressure, residual) in enumerate(samples):
        if abs(residual) <= target_tol:
            if index > 0 and abs(samples[index - 1][1]) > target_tol:
                pressure = _refine_pressure_tolerance_entry(
                    samples[index - 1][0], pressure, hp
                )
            candidates.append(pressure)

    for left, right in zip(samples, samples[1:]):
        if left[1] * right[1] < 0.0:
            try:
                pressure, converged = _find_pressure_from_bracket(
                    (left[0], right[0]),
                    (left[0] + right[0]) / 2.0,
                    hp,
                )
            except (ValueError, RuntimeError):
                continue
            if converged:
                candidates.append(pressure)

    for index in range(1, len(samples) - 1):
        left, centre, right = samples[index - 1 : index + 2]
        if abs(centre[1]) > min(abs(left[1]), abs(right[1])):
            continue
        try:
            result = minimize_scalar(
                fun=lambda pressure: abs(_find_current_gc_dtmin_offset(pressure, hp)),
                method="bounded",
                bounds=(left[0], right[0]),
                options={"xatol": hp.tol},
            )
        except (ValueError, RuntimeError):
            continue
        if result.success and result.fun <= target_tol:
            pressure = float(result.x)
            if abs(left[1]) > target_tol:
                pressure = _refine_pressure_tolerance_entry(left[0], pressure, hp)
            candidates.append(pressure)

    if candidates:
        return min(candidates), True, 0.0

    best_residual = min(abs(residual) for _, residual in samples)
    shortfall = max(0.0, best_residual - target_tol)
    return pressure_start, False, shortfall


def _refine_pressure_tolerance_entry(
    pressure_outside: float,
    pressure_inside: float,
    hp: BaseHeatPump,
) -> float:
    """Return the lowest pressure entering the accepted pinch-error band."""
    target_tol = _get_pressure_target_tolerance(hp)
    try:
        result = root_scalar(
            lambda pressure: abs(_find_current_gc_dtmin_offset(pressure, hp))
            - target_tol,
            method="brentq",
            bracket=(pressure_outside, pressure_inside),
            xtol=hp.tol,
            maxiter=50,
        )
    except ValueError:
        return pressure_inside
    return float(result.root) if result.converged else pressure_inside


def _solve_cycle_at_discharge_pressure(
    pressure: float, hp: BaseHeatPump
) -> BaseHeatPump:
    """Run the cycle analysis at the chosen discharge pressure."""
    hp = _do_sequential_iterative_cycle_analysis(pressure, hp)
    hp = _update_pinch_props_for_gcs(P_comp_out=pressure, hp=hp)
    _validate_cycle_convergence(hp)
    return hp


def _validate_cycle_convergence(hp: BaseHeatPump) -> None:
    """Reject a final state when iteration or an active operation did not converge."""
    if _check_hp_model_convergence(hp):
        return
    raise _CycleConstraintError(
        "Cycle operations failed to converge to a physically valid state.",
        _get_cycle_constraint_errors(hp),
    )


def _get_cycle_constraint_errors(hp: BaseHeatPump) -> dict[str, float]:
    """Return physical errors for the final unconverged cycle state."""
    errors: dict[str, float] = {}
    if not getattr(hp, "iteration_converged", False):
        errors["cycle_iteration"] = 1.0

    unconverged = [
        operation
        for operation in getattr(hp, "operations", ())
        if not getattr(operation, "converged", False)
    ]
    temperature_error = max(
        (_get_heat_exchanger_temperature_error(operation) for operation in unconverged),
        default=0.0,
    )
    if temperature_error > 0.0:
        errors["heat_exchanger_temperature_reversal"] = temperature_error
    elif unconverged:
        errors["component_convergence"] = 1.0

    return errors or {"cycle_convergence": 1.0}


def _get_heat_exchanger_temperature_error(operation: object) -> float:
    """Return the largest temperature-direction violation for an operation."""
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


def _raise_infeasible_pressure_bounds(shortfall: float = 1.0) -> None:
    """Raise a consistent infeasible-pressure failure for the public boundary."""
    raise _CycleConstraintError(
        "Infeasible pressure bounds within property model limits.",
        {"gas_cooler_pinch_shortfall": max(0.0, float(shortfall))},
    )


@profiled_section("pressure_search")
def _find_hp_with_min_discharge_pressure(hp: BaseHeatPump) -> BaseHeatPump:
    """Search for the discharge pressure that meets the gas-cooler pinch target."""
    pressure_start = float(hp.comp_discharge.pressure)
    bnds, x0, hp.feasible, dt_bounds = _get_discharge_pressure_bounds(hp)
    hp.comp_discharge.pressure_bnds = bnds

    pressure = pressure_start
    converged = False
    if hp.feasible:
        pressure, converged = _find_pressure_for_gc_target(bnds, x0, dt_bounds, hp)

    shortfall = max(
        0.0,
        min(abs(residual) for residual in dt_bounds)
        - _get_pressure_target_tolerance(hp),
    )
    if not converged and _is_multi_component_working_fluid(hp):
        pressure, converged, shortfall = _find_unbracketed_multicomponent_pressure(
            hp,
            pressure_start=pressure_start,
        )
    if not converged:
        hp.feasible = False
        _raise_infeasible_pressure_bounds(shortfall)

    hp.feasible = True
    return _solve_cycle_at_discharge_pressure(pressure, hp)


def _expand_discharge_pressure_bounds(
    hp: BaseHeatPump,
    pressure_suction: float,
    pressure_start: float,
    dt_start: float,
    step: float = 1.1,
    max_iter: int = 1000,
) -> tuple[float, float, float, float]:
    """Expand pressure bounds until the pinch residual crosses or reaches target."""
    p_lower = p_upper = pressure_start
    dt_lower = dt_upper = dt_start
    pressure_limit = max(
        float(getattr(hp.comp_discharge, "P_max", np.inf)),
        pressure_start,
    )

    for _ in range(max_iter):
        if dt_upper < 0.0 and dt_lower < 0.0:
            next_pressure = min(
                max(
                    pressure_suction + (p_upper - pressure_suction) * step,
                    pressure_suction,
                ),
                pressure_limit,
            )
            if next_pressure <= p_upper:
                break
            dt_lower = dt_upper
            p_lower = p_upper
            p_upper = next_pressure
            dt_upper = _find_current_gc_dtmin_offset(p_upper, hp)
            if abs(dt_upper) <= _get_pressure_target_tolerance(hp):
                break
        elif dt_upper > 0.0 and dt_lower > 0.0:
            next_pressure = max(
                pressure_suction + (p_lower - pressure_suction) / step,
                pressure_suction,
            )
            if next_pressure >= p_lower:
                break
            dt_upper = dt_lower
            p_upper = p_lower
            p_lower = next_pressure
            dt_lower = _find_current_gc_dtmin_offset(p_lower, hp)
            if abs(dt_lower) <= _get_pressure_target_tolerance(hp):
                break
        else:
            break

    return p_upper, dt_upper, p_lower, dt_lower


@profiled_section("pressure_bounds")
def _get_discharge_pressure_bounds(
    hp: BaseHeatPump,
) -> tuple[tuple[float, float], float, bool, tuple[float, float]]:
    """Determine feasible discharge pressure bounds and initial guess."""
    feasible: bool = False
    P_lb: float
    P_ub: float
    dt_neg: float
    dt_plus: float
    P_i: float

    # Positive dtmin is P_ub and negative dtmin is P_lb
    if hp.comp_discharge.pressure_bnds is not None:
        if not hp.compute_jac:
            P_lb, P_ub = hp.comp_discharge.pressure_bnds
        else:
            P_lb, P_ub = (
                hp.comp_discharge.pressure * 0.995,
                hp.comp_discharge.pressure * 1.005,
            )
        dt_neg, dt_plus = [_find_current_gc_dtmin_offset(P, hp) for P in (P_lb, P_ub)]
    else:
        dt_plus = dt_neg = _find_current_gc_dtmin_offset(hp.comp_discharge.pressure, hp)
        P_ub = P_lb = hp.comp_discharge.pressure
        P_ub, dt_plus, P_lb, dt_neg = _expand_discharge_pressure_bounds(
            hp,
            hp.comp1.flow_in.pressure,
            hp.comp_discharge.pressure,
            dt_plus,
        )

    if P_lb == P_ub or dt_plus == dt_neg:
        P_i = P_lb
    elif hp.comp_discharge.pressure < P_lb or hp.comp_discharge.pressure > P_ub:
        P_i = exp(log(P_lb) - dt_neg / (dt_plus - dt_neg) * (log(P_ub) - log(P_lb)))
    else:
        P_i = hp.comp_discharge.pressure

    if _bounds_cross_target((dt_plus, dt_neg)):
        feasible = True
    elif _select_pressure_bound_within_target((P_lb, P_ub), (dt_plus, dt_neg), hp)[1]:
        feasible = True

    return (P_lb, P_ub), P_i, feasible, (dt_plus, dt_neg)


@profiled_section("gc_dtmin_eval")
def _find_current_gc_dtmin_offset(
    P_comp_out: float, hp: BaseHeatPump, return_abs: bool = False
) -> float:
    """Compute the deviation between current and target gas-cooler pinch."""
    dtmin_ls: list[PinchCandidate]

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
) -> list[PinchCandidate]:
    """Enumerate pinch temperature differences for all gas cooler sections."""
    dtmin_ls: list[PinchCandidate]
    gc: GasCooler

    hp = _do_sequential_iterative_cycle_analysis(P_comp_out, hp)
    dtmin_ls = []
    for gc in hp.gas_coolers:
        dt_hot = gc.source_in.temperature - gc.sink_out.temperature
        dt_cold = gc.source_out.temperature - gc.sink_in.temperature
        dtmin_ls.append((dt_hot, 1.0, True))  # Pinch on the hot side of the gc
        dtmin_ls.append(_find_gc_dt_approaches(gc))  # Pinch in the middle of the gc
        dtmin_ls.append((dt_cold, 0.0, True))  # Pinch on the cold side of the gc

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
    hp.iteration_converged = most_neg_residual >= 0.0
    hp.converged = hp.iteration_converged and all(
        getattr(operation, "converged", False) for operation in hp.operations
    )
    hp = _set_gas_coolers_sink_temperatures(hp)
    if hp.debug:
        hp.get_cycle_plot(x_axis="s", y_axis="T")
        hp.get_cycle_plot(x_axis="h", y_axis="P")
        hp.get_cycle_plot(x_axis="h", y_axis="T")

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


def _find_gc_dt_approaches(gc: GasCooler) -> PinchCandidate:
    """Optimise for the minimum pinch temperature difference inside a gas cooler."""
    gc.source_pinch.pressure = gc.source_out.pressure
    gc.sink_pinch.pressure = gc.sink_in.pressure

    if gc.duty <= 0.0:
        return (np.inf, 0.0, False)

    q_min = 0.001 * gc.duty
    q_max = 0.999 * gc.duty
    r = minimize_scalar(
        fun=lambda Q: _calc_pinch_dtmin(Q, gc),
        method="bounded",
        bounds=(q_min, q_max),
        options={"xatol": max(gc.tol * gc.duty, 1e-9)},
    )

    return (r.fun, r.x / gc.duty, r.success)


def _calc_pinch_dtmin(Q: float, gc: GasCooler) -> float:
    """Return pinch delta-T for a given heat transfer split."""
    gc = _calc_pinch_props_at_Q(Q, gc)
    return _calc_dtmin_at_given_point(gc.source_pinch, gc.sink_pinch)


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


def _print_failed_to_solve_report(
    hp: BaseHeatPump, filename: str = "Cases-Unsolved"
) -> None:
    """Persist the failed case to disk and emit a console warning."""
    print("Heat pump model failed to solve.")
    hp.get_report()
    print_result_to_csv(hp.output, file_name=filename)
