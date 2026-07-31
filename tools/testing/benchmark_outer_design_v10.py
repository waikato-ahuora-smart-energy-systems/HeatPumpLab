"""Run selected no-incumbent v10 robustness benchmarks."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (  # noqa: E402
    OuterDesignOptimisationConfig,
    optimise_recuperative_heat_pump,
)


def _run_case(refrigerant, target: float) -> None:
    config = OuterDesignOptimisationConfig(
        preserve_incumbent=False,
        incumbent_polish_budget=0,
    )
    result = optimise_recuperative_heat_pump(
        base_case={"T_air_tar": target},
        refrigerant=refrigerant,
        config=config,
        n_cpu=5,
        results_dir="/private/tmp/hpl-v10-global-prototype",
    )
    print(
        "best",
        refrigerant,
        target,
        result.best.phase,
        result.best.sec,
        result.best.x,
        "exact",
        result.exact_evaluations,
        "failed",
        result.failed_evaluations,
        "elapsed",
        result.elapsed_seconds,
    )
    for search in result.composition_searches:
        print(
            "search",
            search.composition,
            search.exact_evaluations,
            "edge_a",
            None
            if search.edge_a.best is None
            else (search.edge_a.best.sec, search.edge_a.best.x),
            "edge_b",
            None
            if search.edge_b.best is None
            else (search.edge_b.best.sec, search.edge_b.best.x),
            "split",
            None
            if search.split is None or search.split.best is None
            else (search.split.best.sec, search.split.best.x),
        )
    for polish in result.polish_results:
        print(
            "polish",
            polish.source_x,
            "=>",
            None if polish.best is None else (polish.best.sec, polish.best.x),
            polish.exact_evaluations,
        )
    for record in result.evaluations:
        if not record.success:
            print(
                "failed",
                record.phase,
                record.x,
                record.message,
                record.constraint_errors,
            )


_CASES = {
    "propane-isopentane-180": (("Propane", "Isopentane"), 180.0),
    "co2-135": ("CO2", 135.0),
    "propane-185": ("Propane", 185.0),
    "neopentane-180": ("Neopentane", 180.0),
    "propane-pentane-195": (("Propane", "Pentane"), 195.0),
    "isobutane-pentane-195": (("Isobutane", "Pentane"), 195.0),
    "butane-co2-210": (("Butane", "CO2"), 210.0),
}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case", choices=tuple(_CASES))
    args = parser.parse_args()
    _run_case(*_CASES[args.case])
