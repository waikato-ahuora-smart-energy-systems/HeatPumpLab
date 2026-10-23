"""Gas-cooler pinch calculations for recuperative heat-pump cycles."""

from __future__ import annotations

from typing import TYPE_CHECKING, TypeAlias

import numpy as np
from scipy.optimize import minimize_scalar

from ....utils.graph import plot_driving_force
from ....utils.solve_profile import profiled_section
from ._cycle_analysis import run_sequential_cycle_analysis
from ._fluid_policy import uses_band_pressure_targeting

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

# The profile is sampled first (a subcritical gas cooler with a long subcooled section can have two
# local minima, e.g. near the bubble and dew points), then refined around the smallest sample. At a
# kinked profile (dew point) the approach error is first order in the location error, so the
# tolerance is kept well inside the 0.05 K pressure-target band: about 2e-4 of the sink temperature rise.
_INTERNAL_PINCH_SAMPLES = 8
_INTERNAL_PINCH_RELATIVE_DUTY_TOLERANCE = 2e-4
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
    if not uses_band_pressure_targeting(hp):
        candidates.pop(-1)
        if candidates[-1][1] < 1e-3:
            candidates.pop(-1)
    # Mixtures target the lower tolerance boundary, so the fixed cold terminal
    # is not a spurious zero root. Use the same candidates as final validation.
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

    duty = gas_cooler.duty
    lo, hi = 0.001 * duty, 0.999 * duty
    grid = np.linspace(lo, hi, _INTERNAL_PINCH_SAMPLES + 2)
    samples = [_pinch_temperature_difference(x, gas_cooler) for x in grid]
    best = int(np.argmin(samples))
    result = minimize_scalar(
        lambda x: _pinch_temperature_difference(x, gas_cooler),
        method="bounded",
        bounds=(grid[max(best - 1, 0)], grid[min(best + 1, len(grid) - 1)]),
        options={
            # The internal duty location only supports the explicitly checked
            # 0.05 K pressure-target band; cycle tolerance is unnecessarily
            # strict here and multiplies expensive mixture P-h flashes.
            "xatol": max(_INTERNAL_PINCH_RELATIVE_DUTY_TOLERANCE * duty, 1e-9),
            "maxiter": _INTERNAL_PINCH_MAX_ITERATIONS,
        },
    )
    if samples[best] < result.fun:                 # a sample is a real state too: never report worse
        return (samples[best], grid[best] / duty, True)
    return (result.fun, result.x / duty, result.success)


def update_gc_pinch_properties(
    discharge_pressure: float,
    hp: BaseHeatPump,
) -> BaseHeatPump:
    """Refresh pinch stream properties after discharge pressure changes."""
    candidates = get_all_gc_pinch_candidates(discharge_pressure, hp)
    # Include every terminal, even the fixed cold terminal omitted from the
    # pressure-search residual, when validating the final physical state.
    hp.gc_minimum_approach = float(np.min([candidate[0] for candidate in candidates]))
    # Active pinch: which gas cooler (its gas_coolerN number) and which location, for surrogates.
    # hp.gas_coolers runs hot to cold, so the position is not the number.
    active = int(np.argmin([candidate[0] for candidate in candidates]))
    hp.gc_pinch_active = (f"gc{_gas_cooler_number(hp, hp.gas_coolers[active // 3])}_"
                          f"{('hot', 'internal', 'cold')[active % 3]}")
    hp.gc_approaches = gas_cooler_approaches(candidates)
    candidate_index = 1
    for gas_cooler in hp.gas_coolers:
        _set_pinch_properties(
            candidates[candidate_index][1] * gas_cooler.duty, gas_cooler
        )
        candidate_index += 3
    return hp


def _gas_cooler_number(hp: BaseHeatPump, gas_cooler: GasCooler) -> int:
    """N of the hp.gas_coolerN attribute holding this gas cooler."""
    for number in range(1, 5):
        if getattr(hp, f"gas_cooler{number}", None) is gas_cooler:
            return number
    raise ValueError("gas cooler is not one of hp.gas_cooler1-4")


def gas_cooler_approaches(candidates: list[PinchCandidate]) -> list[float]:
    """Minimum approach (K) in each gas cooler, in hp.gas_coolers order (hot to cold).

    The coldest gas cooler's cold terminal is fixed by the sink inlet and outlet
    specification, so it (and an internal minimum at that terminal) is excluded:
    it says nothing about how close that gas cooler is to pinching elsewhere.
    Two adjacent gas coolers with equal approaches mark a pinch-balance ridge.
    """
    approaches = []
    count = len(candidates) // 3
    for index in range(count):
        group = list(candidates[3 * index: 3 * index + 3])
        if index == count - 1:
            group = [c for k, c in enumerate(group)
                     if not (k == 2 or (k == 1 and c[1] < 1e-3))]
        approaches.append(float(min(c[0] for c in group)) if group else float("inf"))
    return approaches


def _pinch_temperature_difference(duty: float, gas_cooler: GasCooler) -> float:
    _set_pinch_properties(duty, gas_cooler)
    return _temperature_difference(
        gas_cooler.source_pinch,
        gas_cooler.sink_pinch,
    )


def _set_pinch_properties(duty: float, gas_cooler: GasCooler) -> GasCooler:
    gas_cooler.source_pinch.enthalpy = (
        gas_cooler.source_out.enthalpy + duty / gas_cooler.source_out.mass_flow
    )
    gas_cooler.source_pinch.props_in_mPH()
    gas_cooler.sink_pinch.enthalpy = (
        gas_cooler.sink_in.enthalpy + duty / gas_cooler.sink_in.mass_flow
    )
    gas_cooler.sink_pinch.props_in_mPH()
    return gas_cooler


def _temperature_difference(source: Stream, sink: Stream) -> float:
    return source.temperature - sink.temperature
