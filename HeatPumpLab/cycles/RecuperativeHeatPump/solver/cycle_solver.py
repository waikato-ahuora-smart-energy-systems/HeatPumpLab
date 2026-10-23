"""Public orchestration boundary for recuperative heat-pump cycle solves."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ....components.common import TemperatureCrossError
from ...common import CycleFailure, cycle_has_converged
from ....utils.property_cache import property_call_cache
from ....utils.readWriteOutput import print_result_to_csv
from ....utils.solve_profile import profiled_section
from ....utils.properties_package import RefpropError
from ._discharge_pressure import (
    solve_at_discharge_pressure,
    solve_at_minimum_discharge_pressure,
)
from ._evaporator_alignment import (
    align_evaporator_inlet_temperature,
    evaporator_inlet_temperature_offset,
    validate_initial_evaporator_outlet,
)
from ._metrics import populate_full_metrics, populate_objective_metrics
from ._cycle_analysis import InvalidSinkStateError
from ._fluid_policy import uses_band_pressure_targeting

if TYPE_CHECKING:
    from ...common import BaseHeatPump

__all__ = ["solve_heat_pump_cycle"]


def solve_heat_pump_cycle(hp: BaseHeatPump) -> BaseHeatPump:
    """Solve a configured heat pump, mutating the model and returning it."""
    hp.constraint_errors = {}
    hp.numerical_retries = list(getattr(hp, "numerical_retry_history", []))
    hp.solve_failure_kind = None
    hp.property_failure_context = None
    hp.solve_stage = "initialisation"
    hp.solve_error_stage = None
    unexpected = False
    try:
        try:
            failure = _solve_case(hp)
        except RefpropError as exc:
            failure = _retry_evaporator_seed(hp, exc)
        if (failure is not None
                and failure.message == "Final gas-cooler pinch is below the required approach."
                and 0 < failure.constraint_errors.get("gas_cooler_pinch_shortfall", 0) <= 0.05
                and uses_band_pressure_targeting(hp)
                and hasattr(hp, "args") and hp.tol > 1e-8):
            retries = hp.numerical_retries + [{"reason": failure.message, "tol": 1e-8}]
            hp.__init__(**dict(hp.args, tol=1e-8))
            hp.numerical_retries = retries
            hp.constraint_errors = {}
            # A tiny final shortfall can reflect incomplete coupled-state
            # convergence. Retry once; never relax the physical pinch limit.
            failure = _solve_case(hp)
    except InvalidSinkStateError as exc:
        hp.solve_failure_kind = "invalid_sink_state"
        failure = CycleFailure(str(exc), exc.constraint_errors)
    except RefpropError as exc:
        hp.solve_failure_kind = "native_property"
        hp.property_failure_context = getattr(exc, "property_context", None)
        failure = CycleFailure(str(exc), {"native_property_failure": 1.0})
        unexpected = True
    except TemperatureCrossError as exc:
        failure = CycleFailure(
            str(exc),
            {"heat_exchanger_temperature_cross": exc.approach_violation},
        )
    except Exception as exc:
        failure = CycleFailure(str(exc), {})
        unexpected = True

    if failure is not None:
        _apply_failure(hp, failure, unexpected=unexpected)
    return hp


def _retry_evaporator_seed(hp, original_error):
    """Retry one numerical failure with freshly constructed coupled state."""
    if (not uses_band_pressure_targeting(hp)
            or not hasattr(hp, "args")
            or _uses_fixed_discharge_pressure_objective(hp)):
        raise original_error
    seed = float(hp.evap_sup.sink_out.temperature) - max(hp.tol, 1e-4)
    initial_seed = hp.args.get("Tref_evap_outlet", hp.Tref_evap_inlet)
    if abs(seed - initial_seed) <= max(hp.tol, 1e-4):
        raise original_error
    args = dict(hp.args, Tref_evap_outlet=seed)
    diagnostics = hp.numerical_retries + [{"reason": str(original_error), "Tref_evap_outlet": seed,
                                          "property_context": getattr(original_error, "property_context", None)}]
    # Reinitialise the same object so relationships bind to its fresh streams.
    hp.__init__(**args)
    hp.numerical_retries = diagnostics
    hp.constraint_errors = {}
    return _solve_case(hp)


@profiled_section("solve_case")
def _solve_case(hp: BaseHeatPump) -> CycleFailure | None:
    """Run the explicit recuperative-cycle solve steps in order."""
    with property_call_cache():
        hp.solve_stage = "evaporator_seed_validation"
        if hasattr(hp, "Tref_evap_outlet"):
            failure = validate_initial_evaporator_outlet(hp)
            if failure is not None:
                return failure

        fixed_pressure_objective = _uses_fixed_discharge_pressure_objective(hp)
        hp.solve_stage = "discharge_pressure_search"
        failure = (
            # Offset-gradient probes intentionally keep the seed pressure fixed.
            # The selected design is always replayed through the exact path.
            solve_at_discharge_pressure(hp.P_H_est, hp, validate_pinch=False)
            if fixed_pressure_objective
            else solve_at_minimum_discharge_pressure(hp)
        )
        if failure is not None:
            return failure

        if hp.debug:
            hp.get_cycle_plot()

        if not fixed_pressure_objective:
            inlet_offset = evaporator_inlet_temperature_offset(hp)
            if abs(inlet_offset) > 0.1 and cycle_has_converged(hp):
                hp.solve_stage = "evaporator_realignment"
                failure = align_evaporator_inlet_temperature(hp)
                if failure is not None:
                    return failure

        hp.solve_stage = "performance_metrics"
        if fixed_pressure_objective:
            populate_objective_metrics(hp)
        else:
            populate_full_metrics(hp)
        return None


def _uses_fixed_discharge_pressure_objective(hp: BaseHeatPump) -> bool:
    return bool(
        getattr(hp, "fixed_discharge_pressure_objective", False)
        and getattr(hp, "P_H_est", 0.0) > 0.0
    )


def _apply_failure(
    hp: BaseHeatPump,
    failure: CycleFailure,
    *,
    unexpected: bool,
) -> None:
    hp.converged = False
    hp.feasible = False
    hp.solve_error = failure.message
    hp.solve_error_stage = getattr(hp, "solve_stage", None)
    hp.constraint_errors = dict(failure.constraint_errors)
    if hp.debug:
        hp.get_cycle_plot()

    if getattr(hp, "record_failed_solve", True):
        _record_failed_solve(hp)
        disposition = "The failed case has been recorded."
    else:
        disposition = "The failed case was not recorded."
    failure_kind = " with an unexpected exception" if unexpected else ""
    if unexpected or getattr(hp, "report_failed_solve", True):
        print(
            f"Heat pump solver failed{failure_kind}. "
            f"{disposition} Error: {failure.message}"
        )


def _record_failed_solve(
    hp: BaseHeatPump,
    filename: str = "Cases-Unsolved",
) -> None:
    print("Heat pump model failed to solve.")
    hp.get_report()
    print_result_to_csv(hp.output, file_name=filename)
