"""Exact objective evaluation for RHP outer-design optimisation."""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

from ._types import (
    DEFAULT_BASE_CASE,
    DEFAULT_BOUNDS,
    PENALTY_OBJECTIVE,
    DesignEvaluation,
)


class _OuterDesignObjective:
    """Budgeted exact objective for outer recuperative design optimisation."""

    def __init__(
        self,
        *,
        max_evals: int,
        base_case: dict[str, Any] | None = None,
        fluids: tuple[str, str] = ("Propane", "Butane"),
        bounds: np.ndarray = DEFAULT_BOUNDS,
        analyser: Callable | None = None,
        results_dir: str | Path = "results/outer_design_optimisation",
        use_fast_objective: bool = True,
        penalty: float = PENALTY_OBJECTIVE,
        cache_decimals: int = 8,
    ) -> None:
        self.max_evals = int(max_evals)
        self.base_case = dict(DEFAULT_BASE_CASE if base_case is None else base_case)
        self.fluids = fluids
        self.bounds = np.asarray(bounds, dtype=float)
        self.analyser = analyser
        self.results_dir = Path(results_dir)
        self.use_fast_objective = bool(use_fast_objective)
        self.penalty = float(penalty)
        self.cache_decimals = int(cache_decimals)
        self.cache: dict[tuple[float, ...], DesignEvaluation] = {}
        self.records: list[DesignEvaluation] = []
        self.budget_exhausted_calls = 0

    @property
    def exact_evaluations(self) -> int:
        """Return the number of unique exact solves attempted."""
        return len(self.records)

    @property
    def failed_evaluations(self) -> int:
        """Return the number of exact solves that failed or were infeasible."""
        return sum(not record.success for record in self.records)

    @property
    def remaining_budget(self) -> int:
        """Return the number of unique exact solves still available."""
        return max(0, self.max_evals - self.exact_evaluations)

    @property
    def max_recuperator_lift(self) -> float:
        """Return the maximum feasible cold-side recuperator lift."""
        return _max_recuperator_cold_lift(self.base_case)

    def evaluate(self, x: np.ndarray | list[float] | tuple[float, ...]) -> DesignEvaluation:
        """Evaluate one candidate, reusing cached exact results for duplicates."""
        x_arr = self.repair_vector(x)
        cache_key = self._cache_key(x_arr)
        if cache_key in self.cache:
            return self.cache[cache_key]
        if self.remaining_budget <= 0:
            self.budget_exhausted_calls += 1
            return self._budget_exhausted_evaluation(x_arr)

        case = self.vector_to_case(x_arr)
        start = time.perf_counter()
        try:
            hp, solve_time = self._analyse_case(case)
            record = self._record_from_hp(
                x=x_arr,
                case=case,
                hp=hp,
                elapsed_seconds=time.perf_counter() - start,
                fallback_solve_time=solve_time,
            )
        except Exception as exc:
            record = self._failure_record(
                x=x_arr,
                case=case,
                message=f"{type(exc).__name__}: {exc}",
                elapsed_seconds=time.perf_counter() - start,
            )

        self.cache[cache_key] = record
        self.records.append(record)
        return record

    def vector_to_case(self, x: np.ndarray | list[float] | tuple[float, ...]) -> dict[str, Any]:
        """Map a bounded optimiser vector to a recuperative HP case."""
        x_arr = self.repair_vector(x)
        case = dict(self.base_case)
        case.update(
            {
                "dT_rhx1_cold": float(x_arr[0]),
                "dT_rhx2_cold": float(x_arr[1]),
                "comp2_power_frac": float(x_arr[2]),
                "refrigerant": _format_binary_refrigerant(
                    self.fluids[0],
                    self.fluids[1],
                    float(x_arr[3]),
                ),
            }
        )
        return case

    def repair_vector(self, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
        """Project optimiser coordinates into the domain-feasible design space."""
        return _repair_design_vector(
            x,
            self.bounds,
            max_recuperator_lift=self.max_recuperator_lift,
        )

    def _analyse_case(self, case: dict[str, Any]):
        from HeatPumpLab.cycles.RecuperativeHeatPump import (
            RecuperativeHP,
            analyse_RHP_case,
        )

        if self.analyser is None:
            return analyse_RHP_case(
                RecuperativeHP,
                dict(case),
                results_dir=str(self.results_dir),
                use_fast_objective=self.use_fast_objective,
                record_result=False,
            )
        return self.analyser(
            RecuperativeHP,
            dict(case),
            results_dir=str(self.results_dir),
            use_fast_objective=self.use_fast_objective,
        )

    def _record_from_hp(
        self,
        *,
        x: np.ndarray,
        case: dict[str, Any],
        hp,
        elapsed_seconds: float,
        fallback_solve_time: float | None,
    ) -> DesignEvaluation:
        if hp is None:
            return self._failure_record(
                x=x,
                case=case,
                message="analysis returned no heat pump",
                elapsed_seconds=elapsed_seconds,
            )

        report = dict(getattr(hp, "output", {}) or {})
        if not report and hasattr(hp, "get_report"):
            report = dict(hp.get_report(dt=fallback_solve_time or elapsed_seconds))

        sec = _finite_or_none(getattr(hp, "SEC", report.get("SEC", None)))
        cop = _finite_or_none(getattr(hp, "COP", report.get("COP", None)))
        pmax = _finite_or_none(report.get("Pmax", None))
        success = bool(
            getattr(hp, "feasible", False)
            and getattr(hp, "converged", False)
            and sec is not None
        )
        if not success:
            return self._failure_record(
                x=x,
                case=case,
                message="heat pump was infeasible, unconverged, or missing SEC",
                elapsed_seconds=elapsed_seconds,
                report=report,
                constraint_errors=getattr(hp, "constraint_errors", None),
            )

        return DesignEvaluation(
            x=tuple(float(v) for v in x),
            objective=float(sec),
            success=True,
            message="success",
            elapsed_seconds=float(elapsed_seconds),
            sec=float(sec),
            cop=cop,
            pmax=pmax,
            refrigerant=str(case["refrigerant"]),
            case=dict(case),
            report=report,
            constraint_errors={},
        )

    def _failure_record(
        self,
        *,
        x: np.ndarray,
        case: dict[str, Any],
        message: str,
        elapsed_seconds: float,
        report: dict[str, Any] | None = None,
        constraint_errors: dict[str, float] | None = None,
    ) -> DesignEvaluation:
        errors = _validated_constraint_errors(constraint_errors)
        return DesignEvaluation(
            x=tuple(float(v) for v in x),
            objective=self.penalty + sum(error**2 for error in errors.values()),
            success=False,
            message=message,
            elapsed_seconds=float(elapsed_seconds),
            sec=None,
            cop=None,
            pmax=None,
            refrigerant=str(case["refrigerant"]),
            case=dict(case),
            report={} if report is None else dict(report),
            constraint_errors=errors,
        )

    def _budget_exhausted_evaluation(self, x: np.ndarray) -> DesignEvaluation:
        case = self.vector_to_case(x)
        return self._failure_record(
            x=x,
            case=case,
            message="exact evaluation budget exhausted",
            elapsed_seconds=0.0,
            constraint_errors={"budget_exhausted": 1e6},
        )

    def _cache_key(self, x: np.ndarray) -> tuple[float, ...]:
        return tuple(float(v) for v in np.round(x, self.cache_decimals))


def _format_binary_refrigerant(
    fluid_a: str,
    fluid_b: str,
    first_mass_fraction: float,
) -> str:
    """Return a HeatPumpLab binary mass-basis refrigerant string."""
    frac_a = float(np.clip(first_mass_fraction, 0.0, 1.0))
    frac_b = 1.0 - frac_a
    return f"{fluid_a};{fluid_b}|{frac_a:.12g};{frac_b:.12g}|mass"


def _clip_to_bounds(x: np.ndarray, bounds: np.ndarray) -> np.ndarray:
    """Clip a vector to optimizer bounds."""
    bounds = np.asarray(bounds, dtype=float)
    return np.clip(np.asarray(x, dtype=float), bounds[:, 0], bounds[:, 1])


def _max_recuperator_cold_lift(base_case: dict[str, Any]) -> float:
    """Return the model's feasible cold-side recuperator lift limit."""
    t_air_tar = float(base_case.get("T_air_tar", DEFAULT_BASE_CASE["T_air_tar"]))
    t_air_inlet = float(base_case.get("T_air_inlet", DEFAULT_BASE_CASE["T_air_inlet"]))
    dT_rhx3_cold = float(
        base_case.get("dT_rhx3_cold", DEFAULT_BASE_CASE["dT_rhx3_cold"])
    )
    return max(0.0, (t_air_tar - t_air_inlet) * 2.0 / 3.0 - dT_rhx3_cold)


def _repair_design_vector(
    x: np.ndarray | list[float] | tuple[float, ...],
    bounds: np.ndarray,
    *,
    max_recuperator_lift: float,
) -> np.ndarray:
    """Return a bounded vector that satisfies the recuperator lift constraint."""
    repaired = _clip_to_bounds(np.asarray(x, dtype=float), bounds)
    available_lift = max(0.0, float(max_recuperator_lift))
    recuperator_lift = repaired[0] + repaired[1]
    if recuperator_lift > available_lift:
        if recuperator_lift > 0.0:
            repaired[:2] *= available_lift / recuperator_lift
        else:
            repaired[:2] = 0.0
    return repaired


def _finite_or_none(value: Any) -> float | None:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return None
    return numeric if math.isfinite(numeric) else None


def _validated_constraint_errors(
    constraint_errors: dict[str, float] | None,
) -> dict[str, float]:
    """Return finite nonnegative errors, with a fallback for unknown failures."""
    validated: dict[str, float] = {}
    for name, value in (constraint_errors or {}).items():
        numeric = _finite_or_none(value)
        if numeric is not None:
            validated[str(name)] = max(0.0, numeric)
    return validated or {"solve_failure": 1.0}
