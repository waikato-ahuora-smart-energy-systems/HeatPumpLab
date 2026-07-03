"""Run the default recuperative heat-pump outer-design optimisation."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    OuterDesignOptimisationConfig,
    optimise_recuperative_heat_pump_cases,
)

DEFAULT_FLUIDS = ("Propane", "Butane")
DEFAULT_RESULTS_DIR = "results/RHPopt"
DEFAULT_FILENAME = "RHPopt_results"
DEFAULT_PARAM_KWARGS = {
    "fluids_ls": [
        ("Propane", "Isobutane"),
        ("Propane", "Butane"),
        ("Propane", "Isopentane"),
        ("Propane", "Neopentane"),
        ("Propane", "Pentane"),
        ("Isobutane", "Butane"),
        ("Isobutane", "Isopentane"),
        ("Isobutane", "Neopentane"),
        ("Isobutane", "Pentane"),
        ("Butane", "Isopentane"),
        ("Butane", "Neopentane"),
        ("Butane", "Pentane"),
        ("Isopentane", "Neopentane"),
        ("Isopentane", "Pentane"),
        ("Neopentane", "Pentane"),
    ],
}


def run(
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
    run()
