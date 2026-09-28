"""Exact objective evaluation for RHP outer-design optimisation."""

from __future__ import annotations

import math
import time
from collections.abc import Callable
from pathlib import Path
from typing import Any

import numpy as np

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import (
    RefrigerantInput,
    convert_composition,
    normalise_refrigerant_input,
)
from HeatPumpLab.utils.properties_package import (
    COMPRESSOR_GAS_PHASES, PropertyPackage, phase_from_refprop_quality,
)

from ._isolated_evaluation import (
    ExactEvaluationTimeout,
    ExactEvaluationWorkerError,
    IsolatedExactEvaluator,
)
from ._types import (
    DEFAULT_BASE_CASE,
    DEFAULT_BOUNDS,
    DEFAULT_MAX_STAGE_PRESSURE_RATIO,
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
        refrigerant: RefrigerantInput = "Propane;Butane|0.5;0.5|mass",
        bounds: np.ndarray = DEFAULT_BOUNDS,
        analyser: Callable | None = None,
        results_dir: str | Path = "results/outer_design_optimisation",
        use_fast_objective: bool = True,
        penalty: float = PENALTY_OBJECTIVE,
        cache_decimals: int = 8,
        molar_masses: tuple[float, ...] | None = None,
        exact_evaluation_timeout_seconds: float | None = None,
        exact_evaluation_timeout_retries: int = 1,
        exact_evaluation_worker_max_tasks: int = 20,
        max_stage_pressure_ratio: float = DEFAULT_MAX_STAGE_PRESSURE_RATIO,
        min_discharge_superheat_K: float = 0.0,
    ) -> None:
        if not np.isfinite(max_stage_pressure_ratio) or max_stage_pressure_ratio <= 1:
            raise ValueError("max_stage_pressure_ratio must be finite and greater than one.")
        if not np.isfinite(min_discharge_superheat_K) or min_discharge_superheat_K < 0:
            raise ValueError("min_discharge_superheat_K must be finite and nonnegative.")
        self.max_stage_pressure_ratio = float(max_stage_pressure_ratio)
        self.min_discharge_superheat_K = float(min_discharge_superheat_K)
        self.max_evals = int(max_evals)
        self.base_case = dict(DEFAULT_BASE_CASE)
        self.base_case.update(base_case or {})
        self.refrigerant_spec = normalise_refrigerant_input(refrigerant)
        self.input_refrigerant = self.refrigerant_spec.format()
        self.molar_masses = self._get_molar_masses(molar_masses)
        self.initial_mass_fractions = convert_composition(
            self.refrigerant_spec.fractions,
            from_basis=self.refrigerant_spec.basis,
            to_basis="mass",
            molar_masses=self.molar_masses,
        )
        self.bounds = np.asarray(bounds, dtype=float)
        self.analyser = analyser
        self.results_dir = Path(results_dir)
        self.use_fast_objective = bool(use_fast_objective)
        self.penalty = float(penalty)
        self.cache_decimals = int(cache_decimals)
        self.cache: dict[tuple[float, ...], DesignEvaluation] = {}
        self.records: list[DesignEvaluation] = []
        self.budget_exhausted_calls = 0
        self.exact_evaluation_timeout_seconds = (
            None if analyser is not None else exact_evaluation_timeout_seconds
        )
        self.exact_evaluation_worker_max_tasks = int(
            exact_evaluation_worker_max_tasks
        )
        self.exact_evaluation_timeout_retries = int(
            exact_evaluation_timeout_retries
        )
        self._isolated_evaluator: IsolatedExactEvaluator | None = None

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
            if self.exact_evaluation_timeout_seconds is None:
                hp, solve_time = self._analyse_case(case)
                record = self._record_from_hp(
                    x=x_arr,
                    case=case,
                    hp=hp,
                    elapsed_seconds=time.perf_counter() - start,
                    fallback_solve_time=solve_time,
                )
            else:
                record = self._evaluate_isolated(x_arr)
        except ExactEvaluationTimeout as exc:
            record = self._failure_record(
                x=x_arr,
                case=case,
                message=str(exc),
                elapsed_seconds=time.perf_counter() - start,
                constraint_errors={"evaluation_timeout": 1.0},
            )
        except ExactEvaluationWorkerError as exc:
            record = self._failure_record(
                x=x_arr,
                case=case,
                message=f"isolated exact-evaluation worker failed: {exc}",
                elapsed_seconds=time.perf_counter() - start,
                constraint_errors={"worker_failure": 1.0},
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

    def close(self) -> None:
        """Release the isolated exact-evaluation worker, if one was started."""
        if self._isolated_evaluator is not None:
            self._isolated_evaluator.close()
            self._isolated_evaluator = None

    def _evaluate_isolated(self, x: np.ndarray) -> DesignEvaluation:
        if self._isolated_evaluator is None:
            assert self.exact_evaluation_timeout_seconds is not None
            self._isolated_evaluator = IsolatedExactEvaluator(
                timeout_seconds=self.exact_evaluation_timeout_seconds,
                max_tasks=self.exact_evaluation_worker_max_tasks,
            )
        request = {
            "x": tuple(float(value) for value in x),
            "base_case": self.base_case,
            "refrigerant": self.input_refrigerant,
            "bounds": self.bounds,
            "results_dir": str(self.results_dir),
            "use_fast_objective": self.use_fast_objective,
            "penalty": self.penalty,
            "cache_decimals": self.cache_decimals,
            "molar_masses": self.molar_masses,
            "max_stage_pressure_ratio": self.max_stage_pressure_ratio,
            "min_discharge_superheat_K": self.min_discharge_superheat_K,
        }
        attempts = self.exact_evaluation_timeout_retries + 1
        for attempt in range(1, attempts + 1):
            try:
                return self._isolated_evaluator.evaluate(request)
            except (ExactEvaluationTimeout, ExactEvaluationWorkerError) as exc:
                if attempt == attempts:
                    if attempts == 1:
                        raise
                    message = f"{exc} after {attempts} isolated attempts"
                    raise type(exc)(message) from exc
        raise AssertionError("isolated exact-evaluation retry loop did not return")

    def __del__(self) -> None:
        try:
            self.close()
        except Exception:
            pass

    def vector_to_case(self, x: np.ndarray | list[float] | tuple[float, ...]) -> dict[str, Any]:
        """Map a bounded optimiser vector to a recuperative HP case."""
        x_arr = self.repair_vector(x)
        case = dict(self.base_case)
        case.update(
            {
                "dT_rhx1_cold": float(x_arr[0]),
                "dT_rhx2_cold": float(x_arr[1]),
                "comp2_power_frac": float(x_arr[2]),
                "refrigerant": self._format_mass_refrigerant(x_arr[3:]),
            }
        )
        return case

    def repair_vector(self, x: np.ndarray | list[float] | tuple[float, ...]) -> np.ndarray:
        """Project optimiser coordinates into the domain-feasible design space."""
        return _repair_design_vector(
            x,
            self.bounds,
            max_recuperator_lift=self.max_recuperator_lift,
            component_count=self.refrigerant_spec.component_count,
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
            report["solve_error"] = getattr(hp, "solve_error", None)
            report["solve_failure_kind"] = getattr(hp, "solve_failure_kind", None)
            report["solve_error_stage"] = getattr(hp, "solve_error_stage", None)
            report["property_failure_context"] = getattr(hp, "property_failure_context", None)
            report["numerical_retries"] = list(getattr(hp, "numerical_retries", []))
            return self._failure_record(
                x=x,
                case=case,
                message=(getattr(hp, "solve_error", None)
                         or "heat pump was infeasible, unconverged, or missing SEC"),
                elapsed_seconds=elapsed_seconds,
                report=report,
                constraint_errors=getattr(hp, "constraint_errors", None),
            )

        # A converged cycle can still lie outside the compressor operating limits.
        # Check actual active-stage ratios, not the overall ratio or stage split.
        errors = {}
        violations = []
        for stage in range(1, int(case['n_comp'])+1):
            ratio = _finite_or_none(report.get(f'PR_{stage}'))
            if ratio is None or ratio < 1:
                errors[f'stage_{stage}_pressure_ratio_invalid'] = 1.0
                violations.append(f'stage {stage}: missing or invalid pressure ratio')
            elif ratio > self.max_stage_pressure_ratio:
                errors[f'stage_{stage}_pressure_ratio_excess'] = ratio-self.max_stage_pressure_ratio
                violations.append(f'stage {stage}: PR={ratio:.8g} exceeds {self.max_stage_pressure_ratio:g}')
        report['max_stage_pressure_ratio_limit'] = self.max_stage_pressure_ratio
        phase_errors = {}
        phase_violations = []
        for stage in range(1, int(case['n_comp'])+1):
            for end in ('in', 'out'):
                quality = _finite_or_none(report.get(f'Comp_{stage}_q_refprop_{end}'))
                phase = phase_from_refprop_quality(quality)
                report[f'Comp_{stage}_phase_{end}'] = phase
                if phase not in COMPRESSOR_GAS_PHASES:
                    phase_errors[f'stage_{stage}_{end}_not_gas'] = 1.0
                    phase_violations.append(f'stage {stage} {end}: {phase} (REFPROP q={quality})')
        report['compressor_gas_only'] = True
        # Continuous margins (>= 0 satisfied) for constraint models in the search.
        margins = constraint_margins(
            report, int(case['n_comp']), self.max_stage_pressure_ratio,
            self.min_discharge_superheat_K)
        report['constraint_margins'] = margins
        report['min_discharge_superheat_K'] = self.min_discharge_superheat_K
        if self.min_discharge_superheat_K > 0:
            for stage in range(1, int(case['n_comp'])+1):
                if report.get(f'Comp_{stage}_phase_out') == 'supercritical':
                    continue
                margin = margins.get(f'stage_{stage}_discharge_superheat')
                if margin is None or margin < 0:
                    phase_errors[f'stage_{stage}_discharge_superheat_shortfall'] = (
                        1.0 if margin is None else -margin)
                    phase_violations.append(
                        f'stage {stage} out: superheat below {self.min_discharge_superheat_K:g} K')
        if errors:
            report['solve_failure_kind'] = 'constraint_violation'
            report['solve_error_stage'] = 'compressor_pressure_ratio'
            return self._failure_record(
                x=x, case=case, message='Compressor pressure-ratio constraint: '
                +'; '.join(violations+phase_violations),
                elapsed_seconds=elapsed_seconds, report=report,
                constraint_errors={**errors, **phase_errors})
        if phase_errors:
            report['solve_failure_kind'] = 'constraint_violation'
            report['solve_error_stage'] = 'compressor_phase'
            return self._failure_record(
                x=x, case=case,
                message='Compressor gas-only constraint: '+'; '.join(phase_violations),
                elapsed_seconds=elapsed_seconds, report=report, constraint_errors=phase_errors)

        return DesignEvaluation(
            x=tuple(float(v) for v in x),
            objective=float(sec),
            success=True,
            message="success",
            elapsed_seconds=float(elapsed_seconds),
            sec=float(sec),
            cop=cop,
            pmax=pmax,
            refrigerant=self._format_reported_refrigerant(x[3:]),
            mass_refrigerant=str(case["refrigerant"]),
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
            refrigerant=self._format_reported_refrigerant(x[3:]),
            mass_refrigerant=str(case["refrigerant"]),
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

    def _get_molar_masses(
        self, supplied: tuple[float, ...] | None
    ) -> tuple[float, ...]:
        if supplied is not None:
            masses = tuple(float(value) for value in supplied)
        elif self.refrigerant_spec.basis == "mol":
            masses = PropertyPackage(self.input_refrigerant).get_component_molar_masses()
        else:
            masses = tuple(1.0 for _ in self.refrigerant_spec.components)
        if len(masses) != self.refrigerant_spec.component_count:
            raise ValueError("Molar mass count does not match refrigerant components.")
        return masses

    def _format_mass_refrigerant(self, mass_fractions) -> str:
        return self.refrigerant_spec.format(tuple(mass_fractions), "mass")

    def _format_reported_refrigerant(self, mass_fractions) -> str:
        reported = convert_composition(
            tuple(mass_fractions),
            from_basis="mass",
            to_basis=self.refrigerant_spec.basis,
            molar_masses=self.molar_masses,
        )
        return self.refrigerant_spec.format(reported, self.refrigerant_spec.basis)


# Margin given to a supercritical compressor state, which the phase check accepts.
SUPERCRITICAL_MARGIN_K = 0.5


def constraint_margins(report, n_comp, max_stage_pressure_ratio, min_discharge_superheat_K=0.0):
    """Return continuous constraint margins from a solved report (>= 0 satisfied).

    Discharge-superheat margins subtract the required superheat. A supercritical
    discharge counts as satisfied only if its pseudo-saturation margin is not
    negative; the explicit phase check remains the authority for acceptance.
    """
    margins = {}
    for stage in range(1, int(n_comp)+1):
        ratio = _finite_or_none(report.get(f'PR_{stage}'))
        margins[f'stage_{stage}_pressure_ratio'] = (
            None if ratio is None else float(max_stage_pressure_ratio) - ratio)
        for end, name in (('out', 'discharge'), ('in', 'suction')):
            value = _finite_or_none(report.get(f'Comp_{stage}_dT_sh_{end}'))
            required = min_discharge_superheat_K if end == 'out' else 0.0
            margin = None if value is None else value - required
            if report.get(f'Comp_{stage}_phase_{end}') == 'supercritical':
                # Accepted regardless of superheat; keep the label consistent.
                margin = SUPERCRITICAL_MARGIN_K if margin is None else max(margin, SUPERCRITICAL_MARGIN_K)
            margins[f'stage_{stage}_{name}_superheat'] = margin
    return margins


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
    component_count: int,
) -> np.ndarray:
    """Return a bounded vector that satisfies the recuperator lift constraint."""
    values = np.asarray(x, dtype=float)
    if values.shape != (3 + int(component_count),):
        raise ValueError(
            f"Expected {3 + int(component_count)} design values, got {len(values)}."
        )
    repaired = np.concatenate(
        [_clip_to_bounds(values[:3], bounds), _project_to_simplex(values[3:])]
    )
    available_lift = max(0.0, float(max_recuperator_lift))
    recuperator_lift = repaired[0] + repaired[1]
    if recuperator_lift > available_lift:
        if recuperator_lift > 0.0:
            repaired[:2] *= available_lift / recuperator_lift
        else:
            repaired[:2] = 0.0
    return repaired


def _project_to_simplex(values: np.ndarray | list[float]) -> np.ndarray:
    """Return the Euclidean projection onto the nonnegative unit simplex."""
    vector = np.asarray(values, dtype=float)
    if vector.ndim != 1 or vector.size == 0 or np.any(~np.isfinite(vector)):
        raise ValueError("Composition coordinates must be a finite non-empty vector.")
    ordered = np.sort(vector)[::-1]
    cumulative = np.cumsum(ordered) - 1.0
    indices = np.arange(1, vector.size + 1)
    valid = ordered - cumulative / indices > 0.0
    rho = int(np.nonzero(valid)[0][-1])
    theta = cumulative[rho] / float(rho + 1)
    projected = np.maximum(vector - theta, 0.0)
    projected /= projected.sum()
    return projected


def _build_design_vector(
    dT1: float,
    dT2: float,
    comp2: float,
    mass_fractions: tuple[float, ...] | list[float] | np.ndarray,
) -> np.ndarray:
    return np.asarray([dT1, dT2, comp2, *mass_fractions], dtype=float)


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
