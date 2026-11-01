"""Public joint Bayesian RHP outer-design optimiser."""

from pathlib import Path
import time

import numpy as np
from typing import Any

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import RefrigerantInput
from ._bayesian import BayesianDesignSearch, DesignDomain
from ._incumbents import _incumbent_is_compatible
from ._objective import _OuterDesignObjective
from ._types import (
    OPTIMISER_VERSION,
    HistoricalIncumbent,
    IncumbentReplay,
    OuterDesignOptimisationConfig,
    OuterDesignOptimisationResult,
    PhaseEvaluation,
)

__all__ = ["optimise_recuperative_heat_pump"]


def optimise_recuperative_heat_pump(
    base_case: dict[str, Any] | None = None,
    *,
    refrigerant: RefrigerantInput = "Propane;Butane|0.5;0.5|mass",
    config: OuterDesignOptimisationConfig | None = None,
    n_cpu: int = 1,
    results_dir: str | Path = "results/outer_design_optimisation",
    incumbents: tuple[HistoricalIncumbent, ...] = (),
) -> OuterDesignOptimisationResult:
    """Joint BO with one shared budget, including origin, failures and replays.

    Candidate selection is sequential. Batch cases provide process parallelism;
    n_cpu is retained on this public entry point for existing callers.
    """
    config = OuterDesignOptimisationConfig() if config is None else config
    start = time.perf_counter()
    problem = _OuterDesignObjective(
        max_evals=config.max_exact_evaluations,
        base_case=base_case,
        refrigerant=refrigerant,
        results_dir=results_dir,
        use_fast_objective=config.use_fast_objective,
        exact_evaluation_timeout_seconds=config.exact_evaluation_timeout_seconds,
        exact_evaluation_timeout_retries=config.exact_evaluation_timeout_retries,
        exact_evaluation_worker_max_tasks=config.exact_evaluation_worker_max_tasks,
        max_stage_pressure_ratio=config.max_stage_pressure_ratio,
        min_discharge_superheat_K=config.min_discharge_superheat_K,
        bounds=config.bounds,
    )
    evaluations, replays = [], []
    try:
        domain = DesignDomain(
            problem.initial_mass_fractions, problem.max_recuperator_lift, config
        )
        search = BayesianDesignSearch(domain, config)

        def evaluate(x, phase):
            count = problem.exact_evaluations
            record = problem.evaluate(x)
            if problem.exact_evaluations != count + 1:
                raise RuntimeError("BO proposed an already cached design")
            if not np.allclose(record.x, x, atol=1e-10, rtol=0):
                raise RuntimeError("BO design required material repair")
            search.tell(
                record.x, record.objective, record.success,
                (record.report or {}).get("constraint_margins"),
            )
            phase_record = PhaseEvaluation.from_record(
                phase=phase,
                record=record,
                cumulative_exact_evaluations=problem.exact_evaluations,
            )
            evaluations.append(phase_record)
            return phase_record

        if problem.remaining_budget:
            evaluate(domain.base, "origin")
        if config.preserve_incumbent:
            for incumbent in incumbents[: config.incumbent_replay_limit]:
                if not problem.remaining_budget:
                    break
                if not _incumbent_is_compatible(incumbent, problem, config):
                    continue
                encoded = domain.encode(incumbent.x)
                if (
                    domain.contains(encoded)
                    and np.allclose(
                        domain.decode(encoded), incumbent.x, atol=1e-10, rtol=0
                    )
                    and not search.observed(incumbent.x)
                ):
                    replays.append(
                        IncumbentReplay(
                            incumbent, evaluate(incumbent.x, "incumbent_replay")
                        )
                    )
        while problem.remaining_budget:
            x = search.ask()
            if x is None:
                break
            evaluate(x, "bayesian_optimisation")
        successful = [r for r in evaluations if r.success]
        return OuterDesignOptimisationResult(
            optimiser_version=OPTIMISER_VERSION,
            config=config,
            initial_x=tuple(domain.base),
            incumbent_replays=tuple(replays),
            best=min(successful, key=lambda r: r.objective) if successful else None,
            evaluations=tuple(evaluations),
            elapsed_seconds=time.perf_counter() - start,
            exact_evaluations=problem.exact_evaluations,
            failed_evaluations=problem.failed_evaluations,
            diagnostics=tuple(search.diagnostics),
        )
    finally:
        problem.close()
