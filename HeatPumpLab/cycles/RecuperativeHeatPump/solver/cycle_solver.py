"""Public orchestration boundary for recuperative heat-pump cycle solves."""

from __future__ import annotations

from typing import TYPE_CHECKING

from ....components.common import TemperatureCrossError
from ...common import CycleFailure, cycle_has_converged
from ....utils.property_cache import property_call_cache
from ....utils.readWriteOutput import print_result_to_csv
from ....utils.solve_profile import profiled_section
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

if TYPE_CHECKING:
    from ...common import BaseHeatPump

__all__ = ["solve_heat_pump_cycle"]


def solve_heat_pump_cycle(hp: BaseHeatPump) -> BaseHeatPump:
    """Solve a configured heat pump, mutating the model and returning it."""
    hp.constraint_errors = {}
    unexpected = False
    try:
        failure = _solve_case(hp)
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


@profiled_section("solve_case")
def _solve_case(hp: BaseHeatPump) -> CycleFailure | None:
    """Run the explicit recuperative-cycle solve steps in order."""
    with property_call_cache():
        if hasattr(hp, "Tref_evap_outlet"):
            failure = validate_initial_evaporator_outlet(hp)
            if failure is not None:
                return failure

        fixed_pressure_objective = _uses_fixed_discharge_pressure_objective(hp)
        failure = (
            solve_at_discharge_pressure(hp.P_H_est, hp)
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
                failure = align_evaporator_inlet_temperature(hp)
                if failure is not None:
                    return failure

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
