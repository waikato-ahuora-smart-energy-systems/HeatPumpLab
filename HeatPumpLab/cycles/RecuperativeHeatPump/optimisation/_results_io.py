"""Private result CSV helpers for RHP outer-design optimisation."""

from __future__ import annotations

import csv
import json
import math
import re
from collections.abc import Callable
from dataclasses import asdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from ._types import (
    DEFAULT_BASE_CASE,
    VARIABLE_NAMES,
    OuterDesignOptimisationConfig,
    OuterDesignOptimisationResult,
    ResultsCsvCaseDefinition,
)

_LEGACY_OPTIMAL_RESULT = re.compile(r"R-[0-9a-f]{16}\.csv")


def _outer_design_config_from_case_definition(
    case_definition: ResultsCsvCaseDefinition,
    **overrides,
) -> OuterDesignOptimisationConfig:
    """Build a config seeded from a selected result row."""
    reference = case_definition.selected_reference
    values = {
        "fixed_comp2_power_frac": reference.get(
            "comp2_power_frac",
            OuterDesignOptimisationConfig.fixed_comp2_power_frac,
        ),
        "fixed_first_component_mass_fraction": reference.get(
            "first_component_mass_fraction",
            OuterDesignOptimisationConfig.fixed_first_component_mass_fraction,
        ),
    }
    values.update(overrides)
    return OuterDesignOptimisationConfig(**values)


def _load_case_definition_from_results_csv(
    path: str | Path,
    *,
    selector: str = "best-sec",
    case_index: int | None = None,
) -> ResultsCsvCaseDefinition:
    """Load a fixed outer-optimisation case from a prior results CSV."""
    source_csv = Path(path)
    df = _read_numeric_results_csv(source_csv)
    if df.empty:
        raise ValueError(f"No solved numeric results found in {source_csv}.")

    selected_row = _select_reference_row(df, selector=selector, case_index=case_index)
    best_row = df.loc[pd.to_numeric(df["SEC"], errors="coerce").idxmin()]
    base_case = _base_case_from_results_row(selected_row)
    fluids, _composition = _parse_binary_refrigerant(str(selected_row["refrigerant"]))
    return ResultsCsvCaseDefinition(
        base_case=base_case,
        fluids=fluids,
        selected_reference=_reference_from_results_row(selected_row),
        best_reference=_reference_from_results_row(best_row),
        source_csv=str(source_csv),
    )


def _write_outer_design_optimisation_results(
    result: OuterDesignOptimisationResult,
    output_dir: str | Path,
    *,
    case_definition: ResultsCsvCaseDefinition | None = None,
) -> None:
    """Write programmatic result summaries."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summary = result.to_summary_dict()
    if case_definition is not None:
        summary.update(
            {
                "reference_best_SEC": case_definition.best_reference.get("SEC"),
                "reference_best_dT_rhx1_cold": case_definition.best_reference.get(
                    "dT_rhx1_cold"
                ),
                "reference_best_dT_rhx2_cold": case_definition.best_reference.get(
                    "dT_rhx2_cold"
                ),
                "reference_best_comp2_power_frac": case_definition.best_reference.get(
                    "comp2_power_frac"
                ),
                "reference_best_first_component_mass_fraction": case_definition.best_reference.get(
                    "first_component_mass_fraction"
                ),
            }
        )

    payload = {
        "summary": summary,
        "evaluations": [record.to_output_dict() for record in result.evaluations],
        "case_definition": None
        if case_definition is None
        else asdict(case_definition),
    }
    (output_dir / "outer_design_optimisation_summary.json").write_text(
        json.dumps(payload, indent=2, default=_json_default)
    )
    _write_csv(output_dir / "outer_design_optimisation_evaluations.csv", payload["evaluations"])
    _write_csv(output_dir / "outer_design_optimisation_summary.csv", [summary])


def _transfer_optimal_results(
    filename: str,
    *,
    results_dir: str | Path,
) -> None:
    """Merge completed optimum rows into the master CSV and remove temp files."""
    results_path = Path(results_dir)
    temp_paths = _optimal_result_temp_paths(results_path)
    if not temp_paths:
        return

    output_path = _results_csv_path(file_name=filename, results_dir=results_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames, completed_signatures = _master_result_schema(output_path)

    for temp_path in temp_paths:
        with temp_path.open(newline="") as csv_file:
            reader = csv.DictReader(csv_file)
            temp_fields = reader.fieldnames
            if not temp_fields:
                temp_path.unlink()
                continue
            if fieldnames is None:
                fieldnames = temp_fields
            if temp_fields != fieldnames:
                raise ValueError(
                    f"Optimal result schema in {temp_path} does not match {output_path}."
                )
            rows = [
                row
                for row in reader
                if (signature := row.get("case_signature"))
                and signature != "case_signature"
                and signature not in completed_signatures
            ]

        if rows:
            write_header = not output_path.exists() or output_path.stat().st_size == 0
            with output_path.open("a", newline="") as output_file:
                writer = csv.DictWriter(output_file, fieldnames=fieldnames)
                if write_header:
                    writer.writeheader()
                writer.writerows(rows)
            completed_signatures.update(
                row["case_signature"] for row in rows if _result_row_is_complete(row)
            )
        temp_path.unlink()

    optimal_temp_dir = results_path / "temp" / "optimisation"
    if optimal_temp_dir.exists() and not any(optimal_temp_dir.iterdir()):
        optimal_temp_dir.rmdir()


def _completed_result_row(
    signature: str,
    *,
    filename: str,
    results_dir: str | Path,
) -> dict[str, Any] | None:
    path = _results_csv_path(file_name=filename, results_dir=results_dir)
    if not path.exists():
        return None
    with path.open(newline="") as csv_file:
        return next(
            (
                row
                for row in csv.DictReader(csv_file)
                if row.get("case_signature") == signature
                and _result_row_is_complete(row)
            ),
            None,
        )


def _results_csv_path(*, file_name: str, results_dir: str | Path) -> Path:
    relative_path = Path(file_name + ".csv")
    if relative_path.parts and relative_path.parts[0] == "results":
        relative_path = Path(*relative_path.parts[1:])
    return Path(results_dir) / relative_path


def _result_row_is_complete(row: dict[str, Any]) -> bool:
    """Return whether a persisted row contains a usable selected optimum."""
    status = row.get("batch_success", row.get("best_success"))
    if status in (None, ""):
        return True
    return str(status).strip().lower() in {"1", "true", "yes"}


def _optimal_result_temp_paths(results_path: Path) -> list[Path]:
    dedicated = sorted((results_path / "temp" / "optimisation").glob("*.csv"))
    legacy = sorted(
        path
        for path in (results_path / "temp").glob("*.csv")
        if _LEGACY_OPTIMAL_RESULT.fullmatch(path.name)
    )
    return dedicated + legacy


def _master_result_schema(path: Path) -> tuple[list[str] | None, set[str]]:
    if not path.exists() or path.stat().st_size == 0:
        return None, set()
    with path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        return reader.fieldnames, {
            signature
            for row in reader
            if (signature := row.get("case_signature"))
            and signature != "case_signature"
            and _result_row_is_complete(row)
        }


def _parse_binary_refrigerant(
    refrigerant: str,
) -> tuple[tuple[str, str], tuple[float, float]]:
    """Parse a HeatPumpLab binary mass-basis refrigerant string."""
    try:
        fluids_part, composition_part, basis = refrigerant.split("|")
        fluids = tuple(part.strip() for part in fluids_part.split(";"))
        composition = tuple(float(part) for part in composition_part.split(";"))
    except ValueError as exc:
        raise ValueError(f"Expected binary refrigerant string, got {refrigerant!r}.") from exc
    if len(fluids) != 2 or len(composition) != 2 or basis.strip().lower() != "mass":
        raise ValueError(f"Expected binary mass-basis refrigerant, got {refrigerant!r}.")
    return (fluids[0], fluids[1]), (composition[0], composition[1])


def _read_numeric_results_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required_columns = {"SEC", "refrigerant", *VARIABLE_NAMES[:3]}
    missing = sorted(required_columns - set(df.columns))
    if missing:
        raise ValueError(f"Results CSV {path} is missing required columns: {missing}.")
    numeric_sec = pd.to_numeric(df["SEC"], errors="coerce")
    df = df[numeric_sec.notna()].copy()
    df["SEC"] = numeric_sec[numeric_sec.notna()].astype(float)
    return df.reset_index(drop=True)


def _select_reference_row(
    df: pd.DataFrame,
    *,
    selector: str,
    case_index: int | None,
) -> pd.Series:
    if case_index is not None:
        if case_index < 0 or case_index >= len(df):
            raise IndexError(
                f"case_index {case_index} is outside the numeric result range 0-{len(df) - 1}."
            )
        return df.iloc[int(case_index)]
    if selector == "best-sec":
        return df.loc[pd.to_numeric(df["SEC"], errors="coerce").idxmin()]
    if selector == "slowest":
        if "Time_to_optimise" not in df.columns:
            raise ValueError("Cannot select slowest row without Time_to_optimise column.")
        times = pd.to_numeric(df["Time_to_optimise"], errors="coerce")
        return df.loc[times.idxmax()]
    if selector == "first":
        return df.iloc[0]
    raise ValueError(f"Unsupported case selector {selector!r}.")


def _base_case_from_results_row(row: pd.Series) -> dict[str, Any]:
    case = dict(DEFAULT_BASE_CASE)
    converters: dict[str, Callable[[Any], Any]] = {
        "T_source_in": float,
        "Tref_evap_inlet": float,
        "T_air_inlet": float,
        "T_air_tar": float,
        "dT_GC_min_tar": float,
        "dTmin_rhx": float,
        "dT_rhx3_cold": float,
        "comp3_power_frac": float,
        "comp_PR_uniform": _to_bool,
        "comp_model_index": lambda value: int(float(value)),
        "n_comp": lambda value: int(float(value)),
        "ihx_on": _to_bool,
        "motor_eff": float,
        "sink": str,
        "source": str,
    }
    for key, converter in converters.items():
        if key in row and not pd.isna(row[key]):
            case[key] = converter(row[key])
    return case


def _reference_from_results_row(row: pd.Series) -> dict[str, Any]:
    fluids, composition = _parse_binary_refrigerant(str(row["refrigerant"]))
    keys = [
        "SEC",
        "COP",
        "Pmax",
        "Time_to_optimise",
        "dT_rhx1_cold",
        "dT_rhx2_cold",
        "comp2_power_frac",
        "refrigerant",
        "converged",
    ]
    reference = {key: _clean_reference_value(row[key]) for key in keys if key in row}
    reference["fluid_a"] = fluids[0]
    reference["fluid_b"] = fluids[1]
    reference["first_component_mass_fraction"] = float(composition[0])
    return reference


def _clean_reference_value(value: Any) -> Any:
    if isinstance(value, np.generic):
        value = value.item()
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return value
    return numeric if math.isfinite(numeric) else value


def _to_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and math.isfinite(float(value)):
        return bool(int(value))
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Cannot parse boolean value {value!r}.")


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key, value in row.items():
            if key not in fieldnames and not isinstance(value, (dict, list, tuple)):
                fieldnames.append(key)
    with path.open("w", newline="") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})


def _json_default(value):
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    return str(value)
