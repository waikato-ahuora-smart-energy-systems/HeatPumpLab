"""Historical incumbent discovery and current-solver replay."""

from __future__ import annotations

import csv
import math
from dataclasses import replace
from pathlib import Path

import numpy as np

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import (
    normalise_refrigerant_input,
)

from ._objective import _OuterDesignObjective
from ._types import (
    BASE_VARIABLE_NAMES,
    OPTIMISER_VERSION,
    HistoricalIncumbent,
    IncumbentReplay,
    OptimisationCase,
    OuterDesignOptimisationConfig,
    PhaseEvaluation,
)

_DESIGN_INPUT_KEYS = {
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "comp2_power_frac",
    "refrigerant",
}


def attach_historical_incumbents(
    cases: list[OptimisationCase],
    *,
    results_dir: str | Path,
) -> list[OptimisationCase]:
    """Attach distinct prior designs and refresh signatures for batch cases."""
    if not cases:
        return []
    rows = _historical_result_rows(results_dir)
    attached = []
    for case in cases:
        incumbents = (
            _incumbents_for_case(case, rows)
            if case.config.preserve_incumbent
            else ()
        )
        attached.append(replace(case, incumbents=incumbents))
    return attached


def replay_historical_incumbents(
    *,
    incumbents: tuple[HistoricalIncumbent, ...],
    base_case: dict,
    refrigerant: str,
    config: OuterDesignOptimisationConfig,
    results_dir: str | Path,
    objective_factory=_OuterDesignObjective,
) -> tuple[IncumbentReplay, ...]:
    """Cold-replay compatible historical designs with the current solver."""
    if not config.preserve_incumbent or not incumbents:
        return ()

    probe = objective_factory(
        max_evals=0,
        base_case=base_case,
        refrigerant=refrigerant,
        results_dir=results_dir,
        use_fast_objective=config.use_fast_objective,
        exact_evaluation_timeout_seconds=None,
        max_stage_pressure_ratio=config.max_stage_pressure_ratio,
        min_discharge_superheat_K=config.min_discharge_superheat_K,
    )
    compatible = tuple(
        incumbent
        for incumbent in incumbents[: config.incumbent_replay_limit]
        if _incumbent_is_compatible(incumbent, probe, config)
    )
    objective = objective_factory(
        max_evals=len(compatible),
        base_case=base_case,
        refrigerant=refrigerant,
        results_dir=results_dir,
        use_fast_objective=config.use_fast_objective,
        exact_evaluation_timeout_seconds=config.exact_evaluation_timeout_seconds,
        exact_evaluation_timeout_retries=config.exact_evaluation_timeout_retries,
        exact_evaluation_worker_max_tasks=config.exact_evaluation_worker_max_tasks,
        max_stage_pressure_ratio=config.max_stage_pressure_ratio,
        min_discharge_superheat_K=config.min_discharge_superheat_K,
    )
    replays = []
    for incumbent in compatible:
        record = objective.evaluate(incumbent.x)
        replays.append(
            IncumbentReplay(
                incumbent=incumbent,
                evaluation=PhaseEvaluation.from_record(
                    phase="incumbent_replay",
                    record=record,
                    cumulative_exact_evaluations=objective.exact_evaluations,
                ),
            )
        )
    return tuple(replays)


def best_successful_incumbent_replay(
    replays: tuple[IncumbentReplay, ...],
) -> PhaseEvaluation | None:
    """Return the lowest-objective successful incumbent replay."""
    successful = [replay.evaluation for replay in replays if replay.evaluation.success]
    return min(successful, key=lambda record: record.objective) if successful else None


def _historical_result_rows(results_dir: str | Path) -> list[tuple[Path, dict[str, str]]]:
    rows = []
    for path in sorted(Path(results_dir).glob("*.csv")):
        try:
            with path.open(newline="") as csv_file:
                reader = csv.DictReader(csv_file)
                if not _is_historical_result_schema(reader.fieldnames):
                    continue
                rows.extend((path, row) for row in reader)
        except (OSError, csv.Error):
            continue
    return rows


def _is_historical_result_schema(fieldnames: list[str] | None) -> bool:
    required = {
        "input_refrigerant",
        "optimiser_version",
        "batch_success",
        "best_SEC",
        "best_mass_fractions",
        *(f"best_{name}" for name in BASE_VARIABLE_NAMES),
    }
    return bool(fieldnames and required.issubset(fieldnames))


def _incumbents_for_case(
    case: OptimisationCase,
    rows: list[tuple[Path, dict[str, str]]],
) -> tuple[HistoricalIncumbent, ...]:
    best_by_version: dict[str, HistoricalIncumbent] = {}
    for path, row in rows:
        incumbent = _incumbent_from_row(case, path, row)
        if incumbent is None:
            continue
        current = best_by_version.get(incumbent.source_version)
        if current is None or incumbent.stored_sec < current.stored_sec:
            best_by_version[incumbent.source_version] = incumbent

    distinct: dict[tuple[float, ...], HistoricalIncumbent] = {}
    for incumbent in sorted(
        best_by_version.values(),
        key=lambda item: (item.stored_sec, item.source_version, item.source_file),
    ):
        key = tuple(round(value, 8) for value in incumbent.x)
        distinct.setdefault(key, incumbent)
    limit = max(0, int(case.config.incumbent_replay_limit))
    return tuple(distinct.values())[:limit]


def _incumbent_from_row(
    case: OptimisationCase,
    path: Path,
    row: dict[str, str],
) -> HistoricalIncumbent | None:
    version = row.get("optimiser_version", "")
    if (
        not version
        or version == OPTIMISER_VERSION
        or not _row_is_successful(row)
        or not _row_matches_case(row, case)
    ):
        return None
    try:
        stored_sec = float(row["best_SEC"])
        composition = tuple(
            float(value) for value in row["best_mass_fractions"].split(";")
        )
        x = tuple(float(row[name]) for name in (f"best_{n}" for n in BASE_VARIABLE_NAMES))
        x += composition
    except (KeyError, TypeError, ValueError):
        return None
    if (
        not math.isfinite(stored_sec)
        or stored_sec <= 0.0
        or len(composition)
        != normalise_refrigerant_input(case.refrigerant).component_count
        or any(not math.isfinite(value) for value in x)
    ):
        return None
    return HistoricalIncumbent(
        source_version=version,
        source_file=path.name,
        stored_sec=stored_sec,
        x=x,
    )


def _row_is_successful(row: dict[str, str]) -> bool:
    return str(row.get("batch_success", "")).strip().lower() in {
        "1",
        "true",
        "yes",
    }


def _row_matches_case(row: dict[str, str], case: OptimisationCase) -> bool:
    try:
        row_refrigerant = normalise_refrigerant_input(row["input_refrigerant"]).format()
    except (KeyError, TypeError, ValueError):
        return False
    if row_refrigerant != normalise_refrigerant_input(case.refrigerant).format():
        return False
    match_keys = tuple(
        key for key in case.base_case if key not in _DESIGN_INPUT_KEYS
    )
    return all(
        key in row and _values_match(row[key], case.base_case[key])
        for key in match_keys
    )


def _values_match(row_value: str, case_value: object) -> bool:
    if isinstance(case_value, bool):
        return row_value.strip().lower() == str(case_value).lower()
    if isinstance(case_value, (int, float)) and not isinstance(case_value, bool):
        try:
            return math.isclose(
                float(row_value),
                float(case_value),
                rel_tol=0.0,
                abs_tol=1e-9,
            )
        except ValueError:
            return False
    return str(row_value) == str(case_value)


def _incumbent_is_compatible(
    incumbent: HistoricalIncumbent,
    objective: _OuterDesignObjective,
    config: OuterDesignOptimisationConfig,
) -> bool:
    try:
        original = np.asarray(incumbent.x, dtype=float)
        repaired = objective.repair_vector(original)
    except (TypeError, ValueError):
        return False
    if not np.allclose(original, repaired, atol=1e-8, rtol=0.0):
        return False
    fixed_checks = (
        (config.vary_dT_rhx1_cold, original[0], config.fixed_dT_rhx1_cold),
        (config.vary_dT_rhx2_cold, original[1], config.fixed_dT_rhx2_cold),
        (
            config.vary_comp2_power_frac,
            original[2],
            config.fixed_comp2_power_frac,
        ),
    )
    if any(
        not varies and not math.isclose(value, fixed, abs_tol=1e-8)
        for varies, value, fixed in fixed_checks
    ):
        return False
    if not config.vary_composition and not np.allclose(
        original[3:],
        objective.initial_mass_fractions,
        atol=1e-8,
        rtol=0.0,
    ):
        return False
    return True
