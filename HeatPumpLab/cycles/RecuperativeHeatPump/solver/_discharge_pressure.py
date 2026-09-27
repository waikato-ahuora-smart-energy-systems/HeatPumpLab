"""Discharge-pressure targeting for recuperative heat-pump cycles."""

from __future__ import annotations

from math import exp, isfinite, log
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
from ....utils.properties_package import RefpropError
from ._cycle_analysis import run_sequential_cycle_analysis
from ._fluid_policy import uses_band_pressure_targeting
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
_PRESSURE_BOUND_EXPANSION_MAX_ITERATIONS = 100
# Coarse bisection remains a fallback for failed interpolation. Its accepted
# pressure grid is too coarse for smooth local COP differences.
_MIXTURE_PRESSURE_RELATIVE_TOLERANCE = 1e-3
_MIXTURE_PRESSURE_LOG_TOLERANCE = 1e-9
_PRESSURE_BAND_INTERIOR_MARGIN = 1e-5  # K


def pressure_target_tolerance(hp: BaseHeatPump) -> float:
    """Return the absolute pinch residual tolerated during pressure targeting."""
    return max(10.0 * hp.tol, 0.05)


def is_multicomponent_working_fluid(hp: BaseHeatPump) -> bool:
    """Return whether the band-target pressure search applies (name kept for callers).

    Since v21 this is true for pure fluids as well; see ``uses_band_pressure_targeting``.
    """
    return uses_band_pressure_targeting(hp)


def _mixture_target(hp):
    tolerance = pressure_target_tolerance(hp)
    return -tolerance + min(_PRESSURE_BAND_INTERIOR_MARGIN, tolerance / 10)


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
    mixture = is_multicomponent_working_fluid(hp)
    crosses_target = (
        residuals[1] < _mixture_target(hp) < residuals[0]
        if mixture else _bounds_cross_target(residuals)
    )
    if crosses_target:
        if is_multicomponent_working_fluid(hp):
            interpolated = _interpolate_mixture_pressure(bounds, residuals, hp)
            if interpolated is not None:
                return interpolated, True
            result = find_lowest_log_bracketed_target(
                lambda pressure: find_gc_pinch_residual(pressure, hp),
                bounds,
                (residuals[1], residuals[0]),
                target_tolerance=pressure_target_tolerance(hp),
                relative_tolerance=max(
                    hp.tol,
                    _MIXTURE_PRESSURE_RELATIVE_TOLERANCE,
                ),
            )
            if result.success:
                return result.value, True

        if mixture:
            return bounds[0], False
        try:
            return _find_bracketed_pressure(bounds, initial, hp)
        except ValueError:
            if endpoint_success:
                return endpoint_pressure, True

    if mixture:
        # An arbitrary feasible endpoint does not establish minimum pressure.
        for pressure, residual in zip(bounds, reversed(residuals)):
            if residual == _mixture_target(hp):
                return pressure, True
        return bounds[0], False
    if endpoint_success:
        return endpoint_pressure, True
    return _find_minimum_residual_pressure(bounds, hp)


def _interpolate_mixture_pressure(bounds, residuals, hp):
    """Resolve entry into the existing pinch band without pressure quantization."""
    tolerance = pressure_target_tolerance(hp)
    threshold = _mixture_target(hp)
    upper_residual, lower_residual = residuals
    if not lower_residual < threshold < upper_residual:
        return None
    try:
        result = root_scalar(
            lambda value: find_gc_pinch_residual(exp(value), hp) - threshold,
            method="brentq",
            bracket=tuple(log(pressure) for pressure in bounds),
            xtol=_MIXTURE_PRESSURE_LOG_TOLERANCE,
            maxiter=50,
        )
        pressure = exp(result.root)
        residual = find_gc_pinch_residual(pressure, hp)
        if result.converged and isfinite(residual) and abs(residual) <= tolerance:
            return pressure
    except RefpropError:
        # A failed native probe is unknown, not an interpolatable residual.
        raise
    except (ValueError, RuntimeError):
        pass
    hp.pressure_interpolation_fallbacks = getattr(hp, "pressure_interpolation_fallbacks", 0) + 1
    return None


def solve_at_discharge_pressure(
    pressure: float,
    hp: BaseHeatPump,
    *,
    validate_pinch: bool = True,
) -> CycleFailure | None:
    """Solve and validate the cycle at one discharge pressure."""
    run_sequential_cycle_analysis(pressure, hp)
    update_gc_pinch_properties(pressure, hp)
    failure = cycle_convergence_failure(hp)
    if failure is not None or not validate_pinch:
        return failure
    # A cached search residual can belong to a different coupled cycle state.
    # Validate the refreshed final approach, not just component convergence.
    minimum_approach = hp.gc_minimum_approach
    if not isfinite(minimum_approach):
        return CycleFailure("Final gas-cooler pinch is not finite.", {})
    shortfall = hp.dT_GC_min_tar - minimum_approach - pressure_target_tolerance(hp)
    if minimum_approach < hp.dT_GC_min_tar - pressure_target_tolerance(hp):
        return CycleFailure(
            "Final gas-cooler pinch is below the required approach.",
            {"gas_cooler_pinch_shortfall": shortfall},
        )
    return None


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
        min(abs(residual) for residual in residuals) - pressure_target_tolerance(hp),
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
    max_iterations: int = _PRESSURE_BOUND_EXPANSION_MAX_ITERATIONS,
) -> tuple[float, float, float, float]:
    """Expand pressure bounds until the pinch residual reaches or crosses target."""
    lower_pressure = upper_pressure = pressure_start
    lower_residual = upper_residual = residual_start
    mixture = is_multicomponent_working_fluid(hp)
    target = _mixture_target(hp) if mixture else 0.0
    pressure_limit = max(
        float(getattr(hp.comp_discharge, "P_max", np.inf)),
        pressure_start,
    )

    for _ in range(max_iterations):
        if upper_residual < target and lower_residual < target:
            next_pressure = min(
                max(
                    suction_pressure + (upper_pressure - suction_pressure) * step,
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
            if not mixture and abs(upper_residual) <= pressure_target_tolerance(hp):
                break
        elif upper_residual > target and lower_residual > target:
            next_pressure = max(
                suction_pressure + (lower_pressure - suction_pressure) / step,
                suction_pressure * (1.0 + 1e-9) if mixture else suction_pressure,
            )
            if next_pressure >= lower_pressure:
                break
            upper_residual = lower_residual
            upper_pressure = lower_pressure
            lower_pressure = next_pressure
            lower_residual = find_gc_pinch_residual(lower_pressure, hp)
            if not mixture and abs(lower_residual) <= pressure_target_tolerance(hp):
                break
        else:
            break

    return upper_pressure, upper_residual, lower_pressure, lower_residual


@profiled_section("pressure_bounds")
def get_discharge_pressure_bounds(
    hp: BaseHeatPump,
) -> tuple[tuple[float, float], float, bool, tuple[float, float]]:
    """Determine pressure bounds and residuals around the pinch target."""
    mixture = is_multicomponent_working_fluid(hp)
    # Coupled evaporator and recuperator states can move the band boundary.
    # Rebuild mixture brackets instead of accepting stale feasible endpoints.
    if hp.comp_discharge.pressure_bnds is not None and not mixture:
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
    feasible = (
        lower_residual <= _mixture_target(hp) <= upper_residual
        if mixture else (
        _bounds_cross_target(residuals)
        or _select_accepted_endpoint(
            (lower_pressure, upper_pressure),
            residuals,
            hp,
        )[1]
        )
    )
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
