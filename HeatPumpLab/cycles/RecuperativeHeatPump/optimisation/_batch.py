"""Parallel batch runner for RHP outer-design optimisation cases."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import os
from pathlib import Path
from typing import Any

from HeatPumpLab.utils.readWriteOutput import (
    print_result_to_csv,
    transfer_temp_files_to_results,
)

from ._cases import _get_optimisation_case_list_to_solve
from ._types import OptimisationCase, OuterDesignOptimisationConfig
from .outer_design import optimise_recuperative_heat_pump


def optimise_recuperative_heat_pump_cases(
    filename: str = "RHPopt_results",
    param_kwargs: dict[str, Any] | None = None,
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = max((os.cpu_count() or 2) - 1, 1),
    results_dir: str | Path = "results/RHPopt",
) -> list[dict[str, Any]]:
    """Optimise a batch of fixed RHP cases and persist one row per case."""
    config = OuterDesignOptimisationConfig() if config is None else config
    transfer_temp_files_to_results(filename, results_dir=results_dir)
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
    transfer_temp_files_to_results(filename, results_dir=results_dir)
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
        return [_solve_optimisation_case(case, results_dir=str(results_dir)) for case in cases]

    if executor_factory is None:
        executor_factory = ProcessPoolExecutor

    executor = executor_factory(max_workers=worker_count)
    try:
        futures = [
            executor.submit(_solve_optimisation_case, case, results_dir=str(results_dir))
            for case in cases
        ]
        return [future.result() for future in futures]
    finally:
        executor.shutdown()


def _solve_optimisation_case(
    case: OptimisationCase,
    *,
    results_dir: str,
) -> dict[str, Any]:
    result = optimise_recuperative_heat_pump(
        base_case=case.base_case,
        fluids=case.fluids,
        config=case.config,
        results_dir=results_dir,
    )
    row = _result_row(case, result)
    print_result_to_csv(
        output=row,
        file_name="temp/R",
        suffix=case.temp_suffix,
        results_dir=results_dir,
    )
    return row


def _result_row(case: OptimisationCase, result) -> dict[str, Any]:
    row: dict[str, Any] = {
        "case_signature": case.signature,
        "fluid_a": case.fluids[0],
        "fluid_b": case.fluids[1],
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
        row.update(
            {
                "best_SEC": result.best.sec,
                "best_COP": result.best.cop,
                "best_Pmax": result.best.pmax,
            }
        )
    return row
