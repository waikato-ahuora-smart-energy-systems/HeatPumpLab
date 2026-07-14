"""Run the default recuperative heat-pump outer-design optimisation."""

from __future__ import annotations

from itertools import combinations
import os
from pathlib import Path
from typing import Any

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import RefrigerantInput
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    OuterDesignOptimisationConfig,
    optimise_recuperative_heat_pump_case,
    optimise_recuperative_heat_pump_cases,
)

DEFAULT_REFRIGERANT = "Propane;Butane|0.5;0.5|mass"
DEFAULT_RESULTS_DIR = "results/RHPopt"
DEFAULT_FILENAME = "RHPopt_ncomponent_results"
_REFRIGERANTS = (
    "Propane",
    "Isobutane",
    "Butane",
    "Isopentane",
    "Neopentane",
    "Pentane",
    "CO2",
)
DEFAULT_PARAM_KWARGS = {
    "refrigerant_ls": [
        *combinations(_REFRIGERANTS, 2),
        *combinations(_REFRIGERANTS, 3),
        ("Water", "Ammonia"),
    ],
    "T_air_tar_ls": list(range(100, 221, 5)),
}


def run_one(
    base_case: dict[str, Any] | None = None,
    *,
    filename: str = DEFAULT_FILENAME,
    refrigerant: RefrigerantInput = DEFAULT_REFRIGERANT,
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = max((os.cpu_count() or 2) - 1, 1),
    results_dir: str | Path = DEFAULT_RESULTS_DIR,
) -> dict[str, Any]:
    """Optimise and persist one recuperative heat-pump case."""
    return optimise_recuperative_heat_pump_case(
        filename=filename,
        base_case=base_case,
        refrigerant=refrigerant,
        config=OuterDesignOptimisationConfig() if config is None else config,
        n_cpu=n_cpu,
        results_dir=results_dir,
    )


def run_many(
    *,
    filename: str = DEFAULT_FILENAME,
    param_kwargs: dict[str, Any] | None = None,
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = max((os.cpu_count() or 2) - 1, 1),
    results_dir: str | Path = DEFAULT_RESULTS_DIR,
) -> list[dict[str, Any]]:
    """Run the default RHP outer-design optimisation batch."""
    return optimise_recuperative_heat_pump_cases(
        filename=filename,
        param_kwargs=DEFAULT_PARAM_KWARGS if param_kwargs is None else param_kwargs,
        config=OuterDesignOptimisationConfig() if config is None else config,
        n_cpu=n_cpu,
        results_dir=results_dir,
    )


if __name__ == "__main__":
    run_many()
    # run_one()
