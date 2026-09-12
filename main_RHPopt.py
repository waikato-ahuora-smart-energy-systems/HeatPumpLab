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
DEFAULT_FILENAME = "RHPopt_ncomponent_results_bo_v18"
_REFRIGERANTS = (
    "Propane",
    "Isobutane",
    "Butane",
    "Isopentane",
    "Neopentane",
    "Pentane",
    "CO2",
)
_PURE_REFRIGERANTS = (*_REFRIGERANTS, "Water", "Ammonia")
DEFAULT_PARAM_KWARGS = {
    "refrigerant_ls": [
        *_PURE_REFRIGERANTS,
        *combinations(_REFRIGERANTS, 2),
        ("Water", "Ammonia"),
        *combinations(_REFRIGERANTS, 3),
    ],
    "T_air_tar_ls": list(range(100, 221, 5)),
}


def _reverse_case_order_param_kwargs(param_kwargs: dict[str, Any]) -> dict[str, Any]:
    reversed_kwargs = dict(param_kwargs)
    for key, value in param_kwargs.items():
        if key.endswith("_ls") and isinstance(value, (list, tuple)):
            reversed_kwargs[key] = list(reversed(value))
    return reversed_kwargs


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
    reverse: bool = False,
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = max((os.cpu_count() or 2) - 4, 1),
    results_dir: str | Path = DEFAULT_RESULTS_DIR,
) -> list[dict[str, Any]]:
    """Run the default RHP outer-design optimisation batch."""
    effective_param_kwargs = (
        DEFAULT_PARAM_KWARGS if param_kwargs is None else param_kwargs
    )
    if reverse:
        effective_param_kwargs = _reverse_case_order_param_kwargs(effective_param_kwargs)
    return optimise_recuperative_heat_pump_cases(
        filename=filename,
        param_kwargs=effective_param_kwargs,
        config=OuterDesignOptimisationConfig() if config is None else config,
        n_cpu=n_cpu,
        results_dir=results_dir,
    )


if __name__ == "__main__":
    run_many()
    # run_one()
