"""Run the default recuperative heat-pump outer-design optimisation."""

from __future__ import annotations

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    OuterDesignOptimisationConfig,
    OuterDesignOptimisationResult,
    optimise_recuperative_heat_pump,
)

DEFAULT_FLUIDS = ("Propane", "Butane")
DEFAULT_RESULTS_DIR = "results/RHPopt"


def run() -> OuterDesignOptimisationResult:
    param_kwargs = {
        "fluids": DEFAULT_FLUIDS,
        "config": OuterDesignOptimisationConfig(),
        "results_dir": DEFAULT_RESULTS_DIR,
    }
    """Run the default RHP outer-design optimisation."""
    
    return optimise_recuperative_heat_pump(**param_kwargs)


def print_summary(result: OuterDesignOptimisationResult) -> None:
    """Print a concise optimisation summary."""
    summary = result.to_summary_dict()
    print("RHP optimisation complete")
    print(f"Best SEC: {summary.get('best_objective')}")
    print(f"Best phase: {summary.get('best_phase')}")
    print(f"Exact evaluations: {result.exact_evaluations}")
    print(f"Failed evaluations: {result.failed_evaluations}")
    if result.best is not None:
        print(f"dT_rhx1_cold: {result.best.x[0]}")
        print(f"dT_rhx2_cold: {result.best.x[1]}")
        print(f"comp2_power_frac: {result.best.x[2]}")
        print(f"first_component_mass_fraction: {result.best.x[3]}")


if __name__ == "__main__":
    print_summary(run())
