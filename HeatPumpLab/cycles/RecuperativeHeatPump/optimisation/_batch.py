"""Parallel batch runner for RHP outer-design optimisation cases."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import os
from pathlib import Path
from typing import Any

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import (
    RefrigerantInput,
    normalise_refrigerant_input,
)
from HeatPumpLab.utils.readWriteOutput import print_result_to_csv

from ._cases import (
    _build_optimisation_case,
    _get_optimisation_case_list_to_solve,
)
from ._results_io import _completed_result_row, _transfer_optimal_results
from ._types import OptimisationCase, OuterDesignOptimisationConfig
from .outer_design import optimise_recuperative_heat_pump


def optimise_recuperative_heat_pump_case(
    filename: str = "RHPopt_ncomponent_results",
    base_case: dict[str, Any] | None = None,
    *,
    refrigerant: RefrigerantInput = "Propane;Butane|0.5;0.5|mass",
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = 1,
    results_dir: str | Path = "results/RHPopt",
) -> dict[str, Any]:
    """Optimise, persist, and resume one fixed RHP case."""
    config = OuterDesignOptimisationConfig() if config is None else config
    case = _build_optimisation_case(
        base_case=base_case,
        refrigerant=refrigerant,
        config=config,
    )
    _transfer_optimal_results(filename, results_dir=results_dir)
    if completed := _completed_result_row(
        case.signature,
        filename=filename,
        results_dir=results_dir,
    ):
        return completed

    row = _solve_optimisation_case(
        case,
        results_dir=str(results_dir),
        n_cpu=n_cpu,
    )
    _transfer_optimal_results(filename, results_dir=results_dir)
    return row


def optimise_recuperative_heat_pump_cases(
    filename: str = "RHPopt_ncomponent_results",
    param_kwargs: dict[str, Any] | None = None,
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = max((os.cpu_count() or 2) - 1, 1),
    results_dir: str | Path = "results/RHPopt",
) -> list[dict[str, Any]]:
    """Optimise a batch of fixed RHP cases and persist one row per case."""
    config = OuterDesignOptimisationConfig() if config is None else config
    _transfer_optimal_results(filename, results_dir=results_dir)
    cases = _get_optimisation_case_list_to_solve(
        file_name=filename,
        results_dir=results_dir,
        param_kwargs=param_kwargs,
        config=config,
    )
    rows = _run_optimisation_cases(
        cases,
        results_dir=results_dir,
        n_cpu=n_cpu,
    )
    _transfer_optimal_results(filename, results_dir=results_dir)
    return rows


def _run_optimisation_cases(
    cases: list[OptimisationCase],
    *,
    results_dir: str | Path,
    n_cpu: int,
    executor_factory=None,
) -> list[dict[str, Any]]:
    if not cases:
        return []

    worker_count = max(1, int(n_cpu))
    if worker_count == 1:
        return [
            _solve_optimisation_case(case, results_dir=str(results_dir), n_cpu=1)
            for case in cases
        ]

    if executor_factory is None:
        executor_factory = ProcessPoolExecutor

    executor = executor_factory(max_workers=worker_count)
    try:
        futures = [
            executor.submit(
                _solve_optimisation_case,
                case,
                results_dir=str(results_dir),
                n_cpu=1,
            )
            for case in cases
        ]
        return [future.result() for future in futures]
    finally:
        executor.shutdown()


def _solve_optimisation_case(
    case: OptimisationCase,
    *,
    results_dir: str,
    n_cpu: int = 1,
) -> dict[str, Any]:
    result = optimise_recuperative_heat_pump(
        base_case=case.base_case,
        refrigerant=case.refrigerant,
        config=case.config,
        n_cpu=n_cpu,
        results_dir=results_dir,
    )
    row = _result_row(case, result)
    print_result_to_csv(
        output=row,
        file_name="temp/optimisation/R",
        suffix=case.temp_suffix,
        results_dir=results_dir,
    )
    return row


def _result_row(case: OptimisationCase, result) -> dict[str, Any]:
    input_spec = normalise_refrigerant_input(case.refrigerant)
    row: dict[str, Any] = {
        "case_signature": case.signature,
        "input_refrigerant": case.refrigerant,
        "refrigerant_components": ";".join(input_spec.components),
        "refrigerant_component_count": input_spec.component_count,
        "input_basis": input_spec.basis,
        "optimisation_basis": "mass",
        "initial_composition": _format_fractions(input_spec.fractions),
        "initial_mass_fractions": _format_fractions(result.corner.x[3:]),
        **case.base_case,
        **{f"config_{key}": value for key, value in asdict(case.config).items()},
        **result.to_summary_dict(),
    }
    row.update(
        {
            "batch_success": bool(result.best and result.best.success),
            "batch_message": None if result.best is None else result.best.message,
        }
    )
    if result.best is not None:
        best_spec = normalise_refrigerant_input(result.best.refrigerant)
        row.update(
            {
                "best_SEC": result.best.sec,
                "best_COP": result.best.cop,
                "best_Pmax": result.best.pmax,
                "best_refrigerant": result.best.refrigerant,
                "best_composition": _format_fractions(best_spec.fractions),
                "best_mass_fractions": _format_fractions(result.best.x[3:]),
            }
        )
    return row


def _format_fractions(values) -> str:
    return ";".join(f"{float(value):.12g}" for value in values)
