"""Utility to remove duplicate Transcritical RecuperativeHP results.

Rules:
    * Within each operating-point signature, retain a single row if all rows are
      numerically identical within a relative tolerance (default 1e-4).
    * If any row in the signature deviates beyond tolerance, drop the entire
      group and report it to a conflicts CSV for manual inspection.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List

import numpy as np
import pandas as pd

# Core parameters that define a unique operating point.
SIGNATURE_FLOAT_COLUMNS: tuple[str, ...] = (
    "dT_GC_min_tar",
    "T_air_inlet",
    "T_air_tar",
    "T_source_in",
    "Tref_evap_inlet",
    "dTmin_rhx",
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "dT_rhx3_cold",
    "comp1_power_frac",
    "comp2_power_frac",
    "comp3_power_frac",
    "motor_eff",
    "COP",
    "SEC",
    "VHC",
    "eDTL",
    "SLE",
    "LMTD_gcs",
    "Pmax",
)

SIGNATURE_INT_COLUMNS: tuple[str, ...] = (
    "n_comp",
    "ihx_on",
    "comp_model_index",
    "comp_PR_uniform",
)

SIGNATURE_STR_COLUMNS: tuple[str, ...] = (
    "refrigerant",
    "sink",
    "source",
)

DEFAULT_RTL = 1e-4
DEFAULT_ATL = 1e-8


@dataclass(frozen=True)
class DeduplicationResult:
    cleaned: pd.DataFrame
    conflicts: pd.DataFrame
    conflict_signature_count: int


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    results_dir = repo_root / "results"
    parser = argparse.ArgumentParser(
        description="Remove duplicate rows from RecuperativeHP_results.csv according to tolerance rules."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=results_dir / "RecuperativeHP_results.csv",
        help="Path to the input CSV file containing raw results.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=results_dir / "RecuperativeHP_results_cleaned.csv",
        help="Where to write the deduplicated results.",
    )
    parser.add_argument(
        "--conflicts",
        type=Path,
        default=results_dir / "RecuperativeHP_results_conflicts.csv",
        help="Where to write rows from groups that disagreed beyond tolerance.",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=DEFAULT_RTL,
        help="Relative tolerance used for numeric comparison (default: 1e-4).",
    )
    parser.add_argument(
        "--atol",
        type=float,
        default=DEFAULT_ATL,
        help="Absolute tolerance used for numeric comparison (default: 1e-8).",
    )
    return parser.parse_args()


def _normalize_str(value: object) -> str:
    return "" if pd.isna(value) else str(value).strip().lower()


def _round_float(value: object, ndigits: int = 6) -> float:
    if pd.isna(value):
        return float("nan")
    try:
        return round(float(value), ndigits)
    except (TypeError, ValueError):
        return float("nan")


def _int_value(value: object) -> int | None:
    if pd.isna(value):
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def build_signature_column(df: pd.DataFrame) -> pd.Series:
    payload: List[tuple] = []
    for _, row in df.iterrows():
        floats = tuple(_round_float(row.get(col)) for col in SIGNATURE_FLOAT_COLUMNS)
        ints = tuple(_int_value(row.get(col)) for col in SIGNATURE_INT_COLUMNS)
        strs = tuple(_normalize_str(row.get(col)) for col in SIGNATURE_STR_COLUMNS)
        payload.append(floats + ints + strs)
    return pd.Series(payload, index=df.index, name="__signature__")


def _split_columns_by_dtype(df: pd.DataFrame) -> tuple[List[str], List[str]]:
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    other_cols = [col for col in df.columns if col not in numeric_cols and col != "__signature__"]
    return numeric_cols, other_cols


def deduplicate_results(
    df: pd.DataFrame, *, rtol: float, atol: float
) -> DeduplicationResult:
    if df.empty:
        return DeduplicationResult(cleaned=df.copy(), conflicts=df.copy())

    df = df.copy()
    df["__signature__"] = build_signature_column(df)
    numeric_cols, other_cols = _split_columns_by_dtype(df)

    cleaned_rows: List[pd.DataFrame] = []
    conflict_rows: List[pd.DataFrame] = []
    conflict_signatures: set = set()

    groups = df.groupby("__signature__", dropna=False, sort=False)
    for _, group in groups:
        if len(group) == 1:
            cleaned_rows.append(group.drop(columns="__signature__"))
            continue

        numeric_equal = True
        if numeric_cols:
            numeric_values = group[numeric_cols].to_numpy(dtype=float)
            reference = numeric_values[0]
            numeric_equal = np.all(
                np.isclose(numeric_values, reference, rtol=rtol, atol=atol, equal_nan=True)
            )

        other_equal = True
        if other_cols:
            reference_other = group.iloc[0][other_cols].fillna("<NA>").astype(str)
            comparison = (
                group[other_cols]
                .fillna("<NA>")
                .astype(str)
                .eq(reference_other.to_dict())
            )
            other_equal = comparison.all(axis=None)

        if numeric_equal and other_equal:
            cleaned_rows.append(group.iloc[[0]].drop(columns="__signature__"))
        else:
            conflict_rows.append(group.copy())
            conflict_signatures.add(group["__signature__"].iloc[0])

    cleaned_df = pd.concat(cleaned_rows, ignore_index=True) if cleaned_rows else df.iloc[0:0].drop(columns="__signature__")
    conflicts_df = pd.concat(conflict_rows, ignore_index=True) if conflict_rows else df.iloc[0:0]
    return DeduplicationResult(
        cleaned=cleaned_df,
        conflicts=conflicts_df,
        conflict_signature_count=len(conflict_signatures),
    )


def main() -> None:
    args = parse_args()
    if not args.input.exists():
        raise FileNotFoundError(f"Input CSV not found: {args.input}")

    df = pd.read_csv(args.input)
    result = deduplicate_results(df, rtol=args.rtol, atol=args.atol)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.conflicts.parent.mkdir(parents=True, exist_ok=True)

    result.cleaned.to_csv(args.output, index=False)
    print(f"Wrote cleaned results without duplicates to: {args.output} (rows: {len(result.cleaned)})")

    if not result.conflicts.empty:
        result.conflicts.to_csv(args.conflicts, index=False)
        print(
            f"Found {len(result.conflicts)} conflicting rows across {result.conflict_signature_count} signatures. "
            f"Details saved to: {args.conflicts}"
        )
    else:
        if args.conflicts.exists():
            args.conflicts.unlink()
        print("No conflicting duplicates detected.")


if __name__ == "__main__":
    main()
