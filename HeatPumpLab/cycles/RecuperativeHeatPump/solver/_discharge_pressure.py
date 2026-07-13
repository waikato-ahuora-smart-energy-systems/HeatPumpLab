"""Discharge-pressure targeting for recuperative heat-pump cycles."""

from __future__ import annotations

from math import exp, log
from typing import TYPE_CHECKING

import numpy as np
from scipy.optimize import minimize_scalar, root_scalar

from ...common import (
    CycleFailure,
    cycle_convergence_failure,
    find_lowest_log_bracketed_target,
    find_sampled_log_target,
)
from ....utils.solve_profile import profiled_section
from ._cycle_analysis import run_sequential_cycle_analysis
from ._gas_cooler_pinch import find_gc_pinch_residual, update_gc_pinch_properties

if TYPE_CHECKING:
    from ...common import BaseHeatPump

__all__ = [
    "expand_discharge_pressure_bounds",
    "get_discharge_pressure_bounds",
    "select_gc_target_pressure",
    "solve_at_discharge_pressure",
    "solve_at_minimum_discharge_pressure",
]

_PRESSURE_FALLBACK_SAMPLES = 25


def pressure_target_tolerance(hp: BaseHeatPump) -> float:
    """Return the absolute pinch residual tolerated during pressure targeting."""
    return max(10.0 * hp.tol, 0.05)


def is_multicomponent_working_fluid(hp: BaseHeatPump) -> bool:
    """Return whether the working-fluid string names a mixture."""
    return ";" in getattr(hp, "refrigerant", "")


def select_gc_target_pressure(
    bounds: tuple[float, float],
    initial: float,
    residuals: tuple[float, float],
    hp: BaseHeatPump,
) -> tuple[float, bool]:
    """Select a pressure satisfying the gas-cooler pinch target."""
    endpoint_pressure, endpoint_success = _select_accepted_endpoint(
        bounds,
        residuals,
        hp,
    )
    if _bounds_cross_target(residuals):
        if is_multicomponent_working_fluid(hp):
            result = find_lowest_log_bracketed_target(
                lambda pressure: find_gc_pinch_residual(pressure, hp),
                bounds,
                (residuals[1], residuals[0]),
                target_tolerance=pressure_target_tolerance(hp),
                relative_tolerance=hp.tol,
            )
            if result.success:
                return result.value, True

        try:
            return _find_bracketed_pressure(bounds, initial, hp)
        except ValueError:
            if endpoint_success:
                return endpoint_pressure, True

    if endpoint_success:
        return endpoint_pressure, True
    return _find_minimum_residual_pressure(bounds, hp)


def solve_at_discharge_pressure(
    pressure: float,
    hp: BaseHeatPump,
) -> CycleFailure | None:
    """Solve and validate the cycle at one discharge pressure."""
    run_sequential_cycle_analysis(pressure, hp)
    update_gc_pinch_properties(pressure, hp)
    return cycle_convergence_failure(hp)


@profiled_section("pressure_search")
def solve_at_minimum_discharge_pressure(hp: BaseHeatPump) -> CycleFailure | None:
    """Solve at the lowest discharge pressure satisfying the pinch target."""
    pressure_start = float(hp.comp_discharge.pressure)
    bounds, initial, bounds_feasible, residuals = get_discharge_pressure_bounds(hp)
    hp.comp_discharge.pressure_bnds = bounds

    pressure = pressure_start
    converged = False
    if bounds_feasible:
        pressure, converged = select_gc_target_pressure(
            bounds,
            initial,
            residuals,
            hp,
        )

    shortfall = max(
        0.0,
        min(abs(residual) for residual in residuals)
        - pressure_target_tolerance(hp),
    )
    if not converged and is_multicomponent_working_fluid(hp):
        pressure, converged, shortfall = _find_unbracketed_mixture_pressure(
            hp,
            pressure_start=pressure_start,
        )
    if not converged:
        hp.feasible = False
        return CycleFailure(
            "Infeasible pressure bounds within property model limits.",
            {"gas_cooler_pinch_shortfall": shortfall},
        )

    hp.feasible = True
    return solve_at_discharge_pressure(pressure, hp)


def expand_discharge_pressure_bounds(
    hp: BaseHeatPump,
    suction_pressure: float,
    pressure_start: float,
    residual_start: float,
    step: float = 1.1,
    max_iterations: int = 1000,
) -> tuple[float, float, float, float]:
    """Expand pressure bounds until the pinch residual reaches or crosses target."""
    lower_pressure = upper_pressure = pressure_start
    lower_residual = upper_residual = residual_start
    pressure_limit = max(
        float(getattr(hp.comp_discharge, "P_max", np.inf)),
        pressure_start,
    )

    for _ in range(max_iterations):
        if upper_residual < 0.0 and lower_residual < 0.0:
            next_pressure = min(
                max(
                    suction_pressure
                    + (upper_pressure - suction_pressure) * step,
                    suction_pressure,
                ),
                pressure_limit,
            )
            if next_pressure <= upper_pressure:
                break
            lower_residual = upper_residual
            lower_pressure = upper_pressure
            upper_pressure = next_pressure
            upper_residual = find_gc_pinch_residual(upper_pressure, hp)
            if abs(upper_residual) <= pressure_target_tolerance(hp):
                break
        elif upper_residual > 0.0 and lower_residual > 0.0:
            next_pressure = max(
                suction_pressure
                + (lower_pressure - suction_pressure) / step,
                suction_pressure,
            )
            if next_pressure >= lower_pressure:
                break
            upper_residual = lower_residual
            upper_pressure = lower_pressure
            lower_pressure = next_pressure
            lower_residual = find_gc_pinch_residual(lower_pressure, hp)
            if abs(lower_residual) <= pressure_target_tolerance(hp):
                break
        else:
            break

    return upper_pressure, upper_residual, lower_pressure, lower_residual


@profiled_section("pressure_bounds")
def get_discharge_pressure_bounds(
    hp: BaseHeatPump,
) -> tuple[tuple[float, float], float, bool, tuple[float, float]]:
    """Determine pressure bounds and residuals around the pinch target."""
    if hp.comp_discharge.pressure_bnds is not None:
        if not hp.compute_jac:
            lower_pressure, upper_pressure = hp.comp_discharge.pressure_bnds
        else:
            lower_pressure, upper_pressure = (
                hp.comp_discharge.pressure * 0.995,
                hp.comp_discharge.pressure * 1.005,
            )
        lower_residual, upper_residual = [
            find_gc_pinch_residual(pressure, hp)
            for pressure in (lower_pressure, upper_pressure)
        ]
    else:
        upper_residual = lower_residual = find_gc_pinch_residual(
            hp.comp_discharge.pressure,
            hp,
        )
        upper_pressure = lower_pressure = hp.comp_discharge.pressure
        (
            upper_pressure,
            upper_residual,
            lower_pressure,
            lower_residual,
        ) = expand_discharge_pressure_bounds(
            hp,
            hp.comp1.flow_in.pressure,
            hp.comp_discharge.pressure,
            upper_residual,
        )

    if lower_pressure == upper_pressure or upper_residual == lower_residual:
        initial = lower_pressure
    elif not lower_pressure <= hp.comp_discharge.pressure <= upper_pressure:
        initial = exp(
            log(lower_pressure)
            - lower_residual
            / (upper_residual - lower_residual)
            * (log(upper_pressure) - log(lower_pressure))
        )
    else:
        initial = hp.comp_discharge.pressure

    residuals = (upper_residual, lower_residual)
    feasible = _bounds_cross_target(residuals) or _select_accepted_endpoint(
        (lower_pressure, upper_pressure),
        residuals,
        hp,
    )[1]
    return (
        (lower_pressure, upper_pressure),
        initial,
        feasible,
        residuals,
    )


def _bounds_cross_target(residuals: tuple[float, float]) -> bool:
    upper_residual, lower_residual = residuals
    return upper_residual * lower_residual < 0.0


def _select_accepted_endpoint(
    bounds: tuple[float, float],
    residuals: tuple[float, float],
    hp: BaseHeatPump,
) -> tuple[float, bool]:
    upper_residual, lower_residual = residuals
    tolerance = pressure_target_tolerance(hp)
    if abs(lower_residual) <= tolerance:
        return bounds[0], True
    if abs(upper_residual) <= tolerance:
        return bounds[1], True
    return bounds[0], False


def _find_minimum_residual_pressure(
    bounds: tuple[float, float],
    hp: BaseHeatPump,
) -> tuple[float, bool]:
    if bounds[0] == bounds[1]:
        return bounds[0], False
    result = minimize_scalar(
        lambda pressure: find_gc_pinch_residual(pressure, hp, return_abs=True),
        method="bounded",
        bounds=bounds,
        options={"xatol": hp.tol},
    )
    return (
        float(result.x),
        bool(result.success and result.fun <= pressure_target_tolerance(hp)),
    )


def _find_bracketed_pressure(
    bounds: tuple[float, float],
    initial: float,
    hp: BaseHeatPump,
) -> tuple[float, bool]:
    result = root_scalar(
        lambda pressure: find_gc_pinch_residual(pressure, hp),
        method="brentq",
        bracket=bounds,
        x0=initial,
        xtol=hp.tol,
        maxiter=50,
    )
    return float(result.root), bool(result.converged)


def _find_unbracketed_mixture_pressure(
    hp: BaseHeatPump,
    *,
    pressure_start: float,
) -> tuple[float, bool, float]:
    suction_pressure = float(hp.comp1.flow_in.pressure)
    lower_pressure = suction_pressure * (1.0 + max(float(hp.tol), 1e-9))
    upper_pressure = float(getattr(hp.comp_discharge, "P_max", pressure_start))
    if upper_pressure <= lower_pressure:
        return pressure_start, False, 1.0

    result = find_sampled_log_target(
        lambda pressure: find_gc_pinch_residual(pressure, hp),
        (lower_pressure, upper_pressure),
        pressure_start,
        target_tolerance=pressure_target_tolerance(hp),
        x_tolerance=hp.tol,
        sample_count=_PRESSURE_FALLBACK_SAMPLES,
    )
    shortfall = max(0.0, result.residual - pressure_target_tolerance(hp))
    return result.value, result.success, shortfall
