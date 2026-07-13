"""Shared types for the recuperative heat-pump outer optimiser."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

import numpy as np

OPTIMISER_VERSION = "edge-first-pattern-search-v2"

VARIABLE_NAMES = (
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "comp2_power_frac",
    "first_component_mass_fraction",
)

DEFAULT_BOUNDS = np.asarray(
    [
        (0.0, 100.0),
        (0.0, 100.0),
        (0.1, 0.9),
        (0.0, 1.0),
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
        )

    def to_output_dict(self) -> dict[str, Any]:
        output = asdict(self)
        output.update({name: value for name, value in zip(VARIABLE_NAMES, self.x)})
        return output


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
    """Configuration for fixed-recuperator local polishing."""

    budget: int = 12
    initial_steps: tuple[float, float] = (0.05, 0.05)
    min_steps: tuple[float, float] = (0.002, 0.002)
    shrink_factor: float = 0.5
    max_failed_shrinks: int = 20
    improvement_tolerance: float = 0.0
    vary_comp2_power_frac: bool = True
    vary_first_component_mass_fraction: bool = True


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
class OuterDesignOptimisationConfig:
    """Configuration for the edge-first optimiser."""

    max_exact_evaluations: int = 100
    fixed_dT_rhx1_cold: float = 0.0
    fixed_dT_rhx2_cold: float = 0.0
    fixed_comp2_power_frac: float = 0.5
    fixed_first_component_mass_fraction: float = 0.5
    vary_dT_rhx1_cold: bool = True
    vary_dT_rhx2_cold: bool = True
    vary_comp2_power_frac: bool = True
    vary_first_component_mass_fraction: bool = True
    line_budget: int = 8
    split_budget: int = 8
    polish_budget: int = 12
    polish_top_n: int = 3
    use_fast_objective: bool = True
    corner_tolerance: float = 1e-8


@dataclass(frozen=True)
class OuterDesignOptimisationResult:
    """Complete result from the edge-first optimiser."""

    optimiser_version: str
    config: OuterDesignOptimisationConfig
    corner: PhaseEvaluation
    edge_a: LineSearchResult
    edge_b: LineSearchResult
    split: LineSearchResult | None
    polish_results: tuple[LocalPolishResult, ...]
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
            **{f"config_{key}": value for key, value in asdict(self.config).items()},
        }
        if self.best is not None:
            output.update(
                {
                    "best_objective": self.best.objective,
                    "best_phase": self.best.phase,
                    "best_success": self.best.success,
                    "best_message": self.best.message,
                    **{
                        f"best_{name}": value
                        for name, value in zip(VARIABLE_NAMES, self.best.x)
                    },
                }
            )
        return output


@dataclass(frozen=True)
class ResultsCsvCaseDefinition:
    """Fixed case and reference point loaded from an existing results CSV."""

    base_case: dict[str, Any]
    fluids: tuple[str, str]
    selected_reference: dict[str, Any]
    best_reference: dict[str, Any]
    source_csv: str


@dataclass(frozen=True)
class OptimisationCase:
    """One batch case for the RHP outer optimiser."""

    base_case: dict[str, Any]
    fluids: tuple[str, str]
    config: OuterDesignOptimisationConfig
    signature: str
    temp_suffix: str
