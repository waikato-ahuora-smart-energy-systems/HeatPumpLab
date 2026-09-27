"""Shared types for the recuperative heat-pump outer optimiser."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

# BO replaces the staged outer search; invalidate prior completed signatures.
OPTIMISER_VERSION = "isolated-refprop-botorch-logei-v21-unified-pure-path"
DEFAULT_MAX_STAGE_PRESSURE_RATIO = 10.0

BASE_VARIABLE_NAMES = (
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "comp2_power_frac",
)

DEFAULT_BOUNDS = np.asarray(
    [
        (0.0, 100.0),
        (0.0, 100.0),
        (0.1, 0.9),
    ],
    dtype=float,
)

DEFAULT_BASE_CASE = {
    "T_source_in": 40.0,
    "Tref_evap_inlet": 10.0,
    "T_air_inlet": 15.0,
    "T_air_tar": 200.0,
    "dT_GC_min_tar": 20.0,
    "dTmin_rhx": 5.0,
    "dT_rhx1_cold": 0.0,
    "dT_rhx2_cold": 0.0,
    "dT_rhx3_cold": 0.0,
    "comp2_power_frac": 0.5,
    "comp3_power_frac": 0.0,
    "comp_PR_uniform": False,
    "comp_model_index": 1,
    "n_comp": 2,
    "ihx_on": True,
    "motor_eff": 0.96,
    "sink": "Air",
    "source": "water",
    "debug": False,
}

PENALTY_OBJECTIVE = 1e6


@dataclass(frozen=True)
class DesignEvaluation:
    """One unique exact design evaluation."""

    x: tuple[float, ...]
    objective: float
    success: bool
    message: str
    elapsed_seconds: float
    sec: float | None
    cop: float | None
    pmax: float | None
    refrigerant: str
    mass_refrigerant: str
    case: dict[str, Any]
    report: dict[str, Any]
    constraint_errors: dict[str, float]


@dataclass(frozen=True)
class PhaseEvaluation:
    """One optimisation-phase observation of an exact or cached evaluation."""

    phase: str
    x: tuple[float, ...]
    objective: float
    success: bool
    message: str
    sec: float | None
    cop: float | None
    pmax: float | None
    elapsed_seconds: float
    cumulative_exact_evaluations: int
    constraint_errors: dict[str, float]
    refrigerant: str
    mass_refrigerant: str

    @classmethod
    def from_record(
        cls,
        *,
        phase: str,
        record: DesignEvaluation,
        cumulative_exact_evaluations: int,
    ) -> PhaseEvaluation:
        return cls(
            phase=phase,
            x=record.x,
            objective=record.objective,
            success=record.success,
            message=record.message,
            sec=record.sec,
            cop=record.cop,
            pmax=record.pmax,
            elapsed_seconds=record.elapsed_seconds,
            cumulative_exact_evaluations=cumulative_exact_evaluations,
            constraint_errors=dict(record.constraint_errors),
            refrigerant=record.refrigerant,
            mass_refrigerant=record.mass_refrigerant,
        )

    def to_output_dict(self) -> dict[str, Any]:
        output = asdict(self)
        output.update(
            {name: value for name, value in zip(BASE_VARIABLE_NAMES, self.x[:3])}
        )
        output["mass_fractions"] = ";".join(f"{value:.12g}" for value in self.x[3:])
        return output


@dataclass(frozen=True)
class HistoricalIncumbent:
    """One prior successful design available for current-solver replay."""

    source_version: str
    source_file: str
    stored_sec: float
    x: tuple[float, ...]


@dataclass(frozen=True)
class IncumbentReplay:
    """Current-solver evaluation paired with its historical provenance."""

    incumbent: HistoricalIncumbent
    evaluation: PhaseEvaluation


@dataclass(frozen=True)
class OuterDesignOptimisationConfig:
    """Joint Bayesian optimisation; fixed coordinates are removed from the GP."""

    max_exact_evaluations: int = 200
    seed: int = 20260910
    initial_evaluations: int = 16
    candidate_pool_size: int = 512
    acquisition_restarts: int = 6
    exploration_interval: int = 8
    fixed_dT_rhx1_cold: float = 0.0
    fixed_dT_rhx2_cold: float = 0.0
    fixed_comp2_power_frac: float = 0.5
    vary_dT_rhx1_cold: bool = True
    vary_dT_rhx2_cold: bool = True
    vary_comp2_power_frac: bool = True
    vary_composition: bool = True
    use_fast_objective: bool = True
    preserve_incumbent: bool = True
    incumbent_replay_limit: int = 3
    exact_evaluation_timeout_seconds: float | None = 180.0
    exact_evaluation_timeout_retries: int = 1
    exact_evaluation_worker_max_tasks: int = 20
    max_stage_pressure_ratio: float = DEFAULT_MAX_STAGE_PRESSURE_RATIO

    def __post_init__(self) -> None:
        if not np.isfinite(self.max_stage_pressure_ratio) or self.max_stage_pressure_ratio <= 1:
            raise ValueError("max_stage_pressure_ratio must be finite and greater than one.")
        for name in (
            "max_exact_evaluations",
            "incumbent_replay_limit",
            "exact_evaluation_timeout_retries",
        ):
            value = getattr(self, name)
            if not isinstance(value, int) or value < 0:
                raise ValueError(f"{name} must be a nonnegative integer.")
        for name in (
            "initial_evaluations",
            "candidate_pool_size",
            "acquisition_restarts",
            "exploration_interval",
            "exact_evaluation_worker_max_tasks",
        ):
            value = getattr(self, name)
            if not isinstance(value, int) or value < 1:
                raise ValueError(f"{name} must be a positive integer.")
        if not isinstance(self.seed, int) or self.seed < 0:
            raise ValueError("seed must be a nonnegative integer.")
        if self.exact_evaluation_timeout_seconds is not None and (
            not np.isfinite(self.exact_evaluation_timeout_seconds)
            or self.exact_evaluation_timeout_seconds <= 0
        ):
            raise ValueError("exact_evaluation_timeout_seconds must be positive.")
        for value, (lower, upper) in zip(
            (
                self.fixed_dT_rhx1_cold,
                self.fixed_dT_rhx2_cold,
                self.fixed_comp2_power_frac,
            ),
            DEFAULT_BOUNDS,
        ):
            if not np.isfinite(value) or not lower <= value <= upper:
                raise ValueError(
                    "Fixed design values must satisfy the physical box bounds."
                )


@dataclass(frozen=True)
class OuterDesignOptimisationResult:
    """One budgeted joint BO trajectory and its best successful observation."""

    optimiser_version: str
    config: OuterDesignOptimisationConfig
    initial_x: tuple[float, ...]
    incumbent_replays: tuple[IncumbentReplay, ...]
    best: PhaseEvaluation | None
    evaluations: tuple[PhaseEvaluation, ...]
    elapsed_seconds: float
    exact_evaluations: int
    failed_evaluations: int
    diagnostics: tuple[dict[str, Any], ...] = ()

    def to_summary_dict(self) -> dict[str, Any]:
        output = {
            "optimiser_version": self.optimiser_version,
            "optimisation_method": "botorch_logei",
            "elapsed_seconds": self.elapsed_seconds,
            "exact_evaluations": self.exact_evaluations,
            "failed_evaluations": self.failed_evaluations,
            "timed_out_evaluations": sum(
                "evaluation_timeout" in r.constraint_errors for r in self.evaluations
            ),
            "isolated_worker_failures": sum(
                "worker_failure" in r.constraint_errors for r in self.evaluations
            ),
            "incumbent_candidate_count": len(self.incumbent_replays),
            "successful_incumbent_replays": sum(
                r.evaluation.success for r in self.incumbent_replays
            ),
            **{f"config_{key}": value for key, value in asdict(self.config).items()},
        }
        replays = [r for r in self.incumbent_replays if r.evaluation.success]
        if replays:
            replay = min(replays, key=lambda r: r.evaluation.objective)
            output.update(
                incumbent_source_version=replay.incumbent.source_version,
                incumbent_source_file=replay.incumbent.source_file,
                incumbent_stored_SEC=replay.incumbent.stored_sec,
                incumbent_replayed_SEC=replay.evaluation.sec,
                incumbent_replay_success=True,
                selected_incumbent=any(self.best == r.evaluation for r in replays),
                no_worse_than_replayed_incumbent=self.best.objective
                <= replay.evaluation.objective,
                historical_model_drift_pct=100
                * (replay.evaluation.objective / replay.incumbent.stored_sec - 1),
                SEC_change_vs_incumbent_pct=100
                * (self.best.objective / replay.evaluation.objective - 1),
            )
        else:
            output.update(
                incumbent_source_version=None,
                incumbent_source_file=None,
                incumbent_stored_SEC=None,
                incumbent_replayed_SEC=None,
                incumbent_replay_success=False,
                selected_incumbent=False,
                no_worse_than_replayed_incumbent=None,
                historical_model_drift_pct=None,
                SEC_change_vs_incumbent_pct=None,
            )
        if self.best is not None:
            output.update(
                best_objective=self.best.objective,
                best_phase=self.best.phase,
                best_success=self.best.success,
                best_message=self.best.message,
                best_mass_fractions=";".join(f"{v:.12g}" for v in self.best.x[3:]),
                best_refrigerant=self.best.refrigerant,
            )
            output.update(
                {f"best_{k}": v for k, v in zip(BASE_VARIABLE_NAMES, self.best.x[:3])}
            )
        return output


@dataclass(frozen=True)
class ResultsCsvCaseDefinition:
    """Fixed case and reference point loaded from an existing results CSV."""

    base_case: dict[str, Any]
    refrigerant: str
    selected_reference: dict[str, Any]
    best_reference: dict[str, Any]
    source_csv: str


@dataclass(frozen=True)
class OptimisationCase:
    """One batch case for the RHP outer optimiser."""

    base_case: dict[str, Any]
    refrigerant: str
    config: OuterDesignOptimisationConfig
    signature: str
    temp_suffix: str
    incumbents: tuple[HistoricalIncumbent, ...] = ()
