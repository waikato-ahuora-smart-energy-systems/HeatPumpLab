"""Shared types for the recuperative heat-pump outer optimiser."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

# V17 adds validated nearby-state recovery inside failed binary PH calls.
# Earlier completed evaluations must be replayed with the new recovery path.
# Changing the version invalidates completed signatures from the old objective.
OPTIMISER_VERSION = "isolated-refprop-pattern-search-v17"

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
        output.update({name: value for name, value in zip(BASE_VARIABLE_NAMES, self.x[:3])})
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
class LineSearchResult:
    """Result of one edge or split-line search."""

    phase: str
    variable_name: str
    best: PhaseEvaluation | None
    evaluations: tuple[PhaseEvaluation, ...]
    success: bool
    message: str


@dataclass(frozen=True)
class LocalPolishConfig:
    """Configuration for screened local pattern-search polishing."""

    budget: int = 12
    comp2_screen_points: int = 5
    initial_dT_step: float = 5.0
    min_dT_step: float = 0.25
    initial_comp2_step: float = 0.05
    initial_composition_step: float = 0.05
    min_comp2_step: float = 0.002
    min_composition_step: float = 0.002
    shrink_factor: float = 0.5
    max_failed_shrinks: int = 20
    improvement_tolerance: float = 0.0
    vary_dT_rhx1_cold: bool = False
    vary_dT_rhx2_cold: bool = False
    vary_comp2_power_frac: bool = True
    vary_composition: bool = True


@dataclass(frozen=True)
class LocalPolishResult:
    """Result of polishing one fixed-recuperator candidate."""

    method: str
    status: str
    source_x: tuple[float, ...]
    best: PhaseEvaluation | None
    evaluations: tuple[PhaseEvaluation, ...]
    success: bool
    message: str
    elapsed_seconds: float
    exact_evaluations: int
    failed_evaluations: int


@dataclass(frozen=True)
class CompositionSearchResult:
    """Edge-first search result for one fixed starting composition."""

    start_index: int
    composition: tuple[float, ...]
    corner: PhaseEvaluation
    edge_a: LineSearchResult
    edge_b: LineSearchResult
    split: LineSearchResult | None
    evaluations: tuple[PhaseEvaluation, ...]
    success: bool
    message: str
    elapsed_seconds: float
    exact_evaluations: int
    failed_evaluations: int


@dataclass(frozen=True)
class OuterDesignOptimisationConfig:
    """Configuration for the edge-first optimiser."""

    max_exact_evaluations: int = 250
    max_composition_starts: int = 8
    fixed_dT_rhx1_cold: float = 0.0
    fixed_dT_rhx2_cold: float = 0.0
    fixed_comp2_power_frac: float = 0.5
    vary_dT_rhx1_cold: bool = True
    vary_dT_rhx2_cold: bool = True
    vary_comp2_power_frac: bool = True
    vary_composition: bool = True
    line_budget: int = 8
    split_budget: int = 8
    polish_budget: int = 24
    polish_top_n: int = 3
    polish_pure_recuperator_lifts: bool = True
    mixture_lift_feedback_budget: int = 7
    use_fast_objective: bool = True
    corner_tolerance: float = 1e-8
    preserve_incumbent: bool = True
    incumbent_replay_limit: int = 3
    incumbent_polish_budget: int = 12
    exact_evaluation_timeout_seconds: float | None = 180.0
    exact_evaluation_timeout_retries: int = 1
    exact_evaluation_worker_max_tasks: int = 20

    def __post_init__(self) -> None:
        if self.polish_top_n not in {2, 3}:
            raise ValueError("polish_top_n must be either 2 or 3.")
        if self.max_composition_starts < 1:
            raise ValueError("max_composition_starts must be at least 1.")
        if self.incumbent_replay_limit < 0:
            raise ValueError("incumbent_replay_limit must be nonnegative.")
        if (
            self.exact_evaluation_timeout_seconds is not None
            and self.exact_evaluation_timeout_seconds <= 0.0
        ):
            raise ValueError("exact_evaluation_timeout_seconds must be positive.")
        if self.exact_evaluation_worker_max_tasks < 1:
            raise ValueError("exact_evaluation_worker_max_tasks must be at least 1.")
        if self.exact_evaluation_timeout_retries < 0:
            raise ValueError("exact_evaluation_timeout_retries must be nonnegative.")
        for name in (
            "max_exact_evaluations",
            "line_budget",
            "split_budget",
            "polish_budget",
            "mixture_lift_feedback_budget",
            "incumbent_polish_budget",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be nonnegative.")


@dataclass(frozen=True)
class OuterDesignOptimisationResult:
    """Complete result from the edge-first optimiser."""

    optimiser_version: str
    config: OuterDesignOptimisationConfig
    corner: PhaseEvaluation
    edge_a: LineSearchResult
    edge_b: LineSearchResult
    split: LineSearchResult | None
    composition_searches: tuple[CompositionSearchResult, ...]
    requested_composition_starts: int
    polish_results: tuple[LocalPolishResult, ...]
    lift_feedback_results: tuple[LineSearchResult, ...]
    incumbent_replays: tuple[IncumbentReplay, ...]
    incumbent_polish_results: tuple[LocalPolishResult, ...]
    best: PhaseEvaluation | None
    evaluations: tuple[PhaseEvaluation, ...]
    elapsed_seconds: float
    exact_evaluations: int
    failed_evaluations: int

    def to_summary_dict(self) -> dict[str, Any]:
        output: dict[str, Any] = {
            "optimiser_version": self.optimiser_version,
            "elapsed_seconds": self.elapsed_seconds,
            "exact_evaluations": self.exact_evaluations,
            "failed_evaluations": self.failed_evaluations,
            "timed_out_evaluations": sum(
                "evaluation_timeout" in record.constraint_errors
                for record in self.evaluations
            ),
            "isolated_worker_failures": sum(
                "worker_failure" in record.constraint_errors
                for record in self.evaluations
            ),
            "requested_composition_starts": self.requested_composition_starts,
            "attempted_composition_starts": len(self.composition_searches),
            "successful_composition_starts": sum(
                search.success for search in self.composition_searches
            ),
            "polished_direction_count": sum(
                polish.status != "skipped" for polish in self.polish_results
            ),
            "lift_feedback_count": sum(
                result.success for result in self.lift_feedback_results
            ),
            "incumbent_candidate_count": len(self.incumbent_replays),
            "successful_incumbent_replays": sum(
                replay.evaluation.success for replay in self.incumbent_replays
            ),
            "incumbent_polished_direction_count": sum(
                polish.status != "skipped"
                for polish in self.incumbent_polish_results
            ),
            **{f"config_{key}": value for key, value in asdict(self.config).items()},
        }
        successful_replays = [
            replay for replay in self.incumbent_replays if replay.evaluation.success
        ]
        if successful_replays:
            incumbent_best = min(
                successful_replays,
                key=lambda replay: replay.evaluation.objective,
            )
            incumbent_evaluations = tuple(
                replay.evaluation for replay in self.incumbent_replays
            ) + tuple(
                record
                for polish in self.incumbent_polish_results
                for record in polish.evaluations
            )
            output.update(
                {
                    "incumbent_source_version": (
                        incumbent_best.incumbent.source_version
                    ),
                    "incumbent_source_file": incumbent_best.incumbent.source_file,
                    "incumbent_stored_SEC": incumbent_best.incumbent.stored_sec,
                    "incumbent_replayed_SEC": incumbent_best.evaluation.sec,
                    "incumbent_replay_success": True,
                    "selected_incumbent": self.best in incumbent_evaluations,
                    "no_worse_than_replayed_incumbent": (
                        self.best is not None
                        and self.best.objective
                        <= incumbent_best.evaluation.objective
                    ),
                    "SEC_change_vs_incumbent_pct": (
                        None
                        if self.best is None
                        else 100.0
                        * (
                            self.best.objective
                            / incumbent_best.evaluation.objective
                            - 1.0
                        )
                    ),
                    "historical_model_drift_pct": (
                        100.0
                        * (
                            incumbent_best.evaluation.objective
                            / incumbent_best.incumbent.stored_sec
                            - 1.0
                        )
                    ),
                }
            )
        else:
            output.update(
                {
                    "incumbent_source_version": None,
                    "incumbent_source_file": None,
                    "incumbent_stored_SEC": None,
                    "incumbent_replayed_SEC": None,
                    "incumbent_replay_success": False,
                    "selected_incumbent": False,
                    "no_worse_than_replayed_incumbent": None,
                    "SEC_change_vs_incumbent_pct": None,
                    "historical_model_drift_pct": None,
                }
            )
        if self.best is not None:
            output.update(
                {
                    "best_objective": self.best.objective,
                    "best_phase": self.best.phase,
                    "best_success": self.best.success,
                    "best_message": self.best.message,
                    **{
                        f"best_{name}": value
                        for name, value in zip(BASE_VARIABLE_NAMES, self.best.x[:3])
                    },
                    "best_mass_fractions": ";".join(
                        f"{value:.12g}" for value in self.best.x[3:]
                    ),
                    "best_refrigerant": self.best.refrigerant,
                }
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
