"""Summarise Transcritical Recuperative Heat Pump results by refrigerant."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd


DEFAULT_METRICS: Tuple[Tuple[str, str], ...] = (
    ("COP", "max"),
    # ("SEC", "min"),
    ("VHC", "max"),
    ("eDTL", "min"),
    ("SLE", "min"),
    # ("dQ_per_dW", "max"),
    ("LMTD_gcs", "max"),
    ("Pmax", "min"),
    # ("COP_pct_increase_vs_base", "max"),
    # ("Pmax_reduction_vs_base", "max"),
)

DEFAULT_INCLUDE_COLUMNS: Sequence[str] = (
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "comp2_power_frac",
    "COP",
    "SEC",
    "VHC",
    "eDTL",
    "SLE",
    "dQ_per_dW",
    "LMTD_gcs",
    "Pmax",
    "COP_pct_increase_vs_base",
    "Pmax_reduction_vs_base",
)

BASECASE_DTRHX1_COLD = 0.0
BASECASE_DTRHX2_COLD = 0.0
BASECASE_COMP2_FRAC = 0.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Identify best-performing cases for each refrigerant based on selected metrics."
        )
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "RHP_results.csv",
        help="Path to the results CSV file.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "Summary.csv",
        help="Optional path to write the summary CSV.",
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=[f"{metric}:{mode}" for metric, mode in DEFAULT_METRICS],
        help=(
            "Metrics to evaluate. Use METRIC or METRIC:max/METRIC:min to control optimisation direction."
        ),
    )
    parser.add_argument(
        "--include-columns",
        nargs="+",
        # default=list(DEFAULT_INCLUDE_COLUMNS),
        help="Additional columns to include in the summary output.",
    )
    parser.add_argument(
        "--refrigerants",
        nargs="+",
        help="Restrict analysis to specific refrigerants.",
    )
    parser.add_argument(
        "--comp2-values",
        nargs="+",
        type=float,
        help="Restrict analysis to specific comp2_power_frac values.",
    )
    parser.add_argument(
        "--skip-success-filter",
        action="store_true",
        help="Include non-successful optimisation runs.",
    )
    parser.add_argument(
        "--skip-converged-filter",
        action="store_true",
        help="Include non-converged optimisation runs.",
    )
    parser.add_argument(
        "--float-precision",
        type=int,
        default=4,
        help="Number of decimal places to display for floating point values.",
    )
    return parser.parse_args()


def load_results(
    csv_path: Path,
    filter_success: bool,
    filter_converged: bool,
) -> pd.DataFrame:
    df = pd.read_csv(csv_path)
    if filter_success and "Success" in df.columns:
        df = df[df["Success"] == True]  # noqa: E712
    if filter_converged and "converged" in df.columns:
        df = df[df["converged"] == True]  # noqa: E712
    return df


def parse_metric_specs(metric_specs: Sequence[str]) -> List[Tuple[str, str]]:
    """Parse metric specifications into (name, mode) pairs."""
    default_modes = {name: mode for name, mode in DEFAULT_METRICS}
    parsed: List[Tuple[str, str]] = []
    for spec in metric_specs:
        parts = spec.split(":", 1)
        name = parts[0].strip()
        if not name:
            raise ValueError("Metric name cannot be empty.")
        mode = parts[1].strip() if len(parts) > 1 else default_modes.get(name, "max")
        mode = mode.lower()
        if mode not in {"min", "max"}:
            raise ValueError(
                f"Invalid optimisation mode '{mode}' for metric '{name}'. "
                "Use 'min' or 'max'."
            )
        parsed.append((name, mode))
    return parsed


def add_basecase_metrics(df: pd.DataFrame) -> pd.DataFrame:
    required_columns = {
        "refrigerant",
        "COP",
        "Pmax",
        "dT_rhx1_cold",
        "dT_rhx2_cold",
        "comp2_power_frac",
    }
    missing = required_columns.difference(df.columns)
    if missing:
        raise ValueError(
            f"Cannot compute base-case metrics; missing columns: {sorted(missing)}"
        )

    base_mask = (
        np.isclose(df["dT_rhx1_cold"].astype(float), BASECASE_DTRHX1_COLD)
        & np.isclose(df["dT_rhx2_cold"].astype(float), BASECASE_DTRHX2_COLD)
        & np.isclose(df["comp2_power_frac"].astype(float), BASECASE_COMP2_FRAC)
    )
    base_rows = df.loc[base_mask, ["refrigerant", "COP", "Pmax"]]

    if base_rows.empty:
        print("Warning: No base-case rows found; base-relative metrics will be NaN.")
        df = df.copy()
        df["COP_pct_increase_vs_base"] = np.nan
        df["Pmax_reduction_vs_base"] = np.nan
        return df

    base_lookup = base_rows.groupby("refrigerant")[["COP", "Pmax"]].mean()
    cop_base_map = base_lookup["COP"].to_dict()
    pmax_base_map = base_lookup["Pmax"].to_dict()

    df = df.copy()
    df["COP_base"] = df["refrigerant"].map(cop_base_map)
    df["Pmax_base"] = df["refrigerant"].map(pmax_base_map)

    with np.errstate(divide="ignore", invalid="ignore"):
        df["COP_pct_increase_vs_base"] = np.where(
            df["COP_base"] > 0,
            (df["COP"] - df["COP_base"]) / df["COP_base"] * 100.0,
            np.nan,
        )
    with np.errstate(divide="ignore", invalid="ignore"):
        df["Pmax_reduction_vs_base"] = np.where(
            df["Pmax_base"] > 0,
            (df["Pmax_base"] - df["Pmax"]) / df["Pmax_base"] * 100.0,
            np.nan,
        )

    df.drop(columns=["COP_base", "Pmax_base"], inplace=True)
    return df


def unique_sorted(values: Iterable[float]) -> List[float]:
    rounded = {}
    for value in values:
        if not isinstance(value, (int, float)) or math.isnan(value):
            continue
        key = round(float(value), 8)
        rounded.setdefault(key, float(value))
    return sorted(rounded.values())


def filter_dataframe(
    df: pd.DataFrame,
    refrigerants: Sequence[str] | None,
    comp2_values: Sequence[float] | None,
) -> pd.DataFrame:
    if refrigerants:
        available = set(df["refrigerant"].dropna().unique())
        chosen = [ref for ref in refrigerants if ref in available]
        if not chosen:
            raise ValueError("No matching refrigerants found in the dataset.")
        df = df[df["refrigerant"].isin(chosen)]

    if comp2_values:
        available_comp2 = (
            unique_sorted(df["comp2_power_frac"].dropna().unique())
            if "comp2_power_frac" in df.columns
            else []
        )
        if not available_comp2:
            raise ValueError("Dataset lacks valid 'comp2_power_frac' values.")
        filtered_values = [
            value
            for value in comp2_values
            if any(math.isclose(value, avail, rel_tol=1e-9) for avail in available_comp2)
        ]
        if not filtered_values:
            raise ValueError(
                "No matching comp2_power_frac values found after applying the filter."
            )
        df = df[
            df["comp2_power_frac"].astype(float).isin(filtered_values)
        ]
    return df


def summarise_best_cases(
    df: pd.DataFrame,
    metrics: Sequence[Tuple[str, str]],
    include_columns: Sequence[str],
) -> pd.DataFrame:
    records: List[dict[str, object]] = []
    for metric, mode in metrics:
        if metric not in df.columns:
            print(f"Skipping metric '{metric}' (not found in data).")
            continue
        metric_series = pd.to_numeric(df[metric], errors="coerce")
        df_metric = df.loc[metric_series.notna()].copy()
        df_metric[metric] = metric_series.loc[df_metric.index]
        if df_metric.empty:
            print(f"Skipping metric '{metric}' (no valid data).")
            continue

        for refrigerant, group in df_metric.groupby("refrigerant"):
            if group.empty:
                continue
            if mode == "max":
                idx = group[metric].idxmax()
            else:
                idx = group[metric].idxmin()
            best_row = group.loc[idx]
            record = {
                "metric": metric,
                "mode": mode,
                "refrigerant": refrigerant,
                "metric_value": best_row[metric],
            }
            include_columns = include_columns if include_columns is not None else best_row.index
            for column in include_columns:
                if column in best_row.index:
                    record[column] = best_row[column]                                       
            records.append(record)

    summary_df = pd.DataFrame.from_records(records)
    if summary_df.empty:
        return summary_df

    summary_df.sort_values(["metric", "refrigerant"], inplace=True)
    summary_df.reset_index(drop=True, inplace=True)
    return summary_df


def main() -> None:
    args = parse_args()

    df = load_results(
        args.csv,
        filter_success=not args.skip_success_filter,
        filter_converged=not args.skip_converged_filter,
    )
    df = add_basecase_metrics(df)
    args.comp2_values = [0.5]
    df = filter_dataframe(df, args.refrigerants, args.comp2_values)

    metric_specs = parse_metric_specs(args.metrics)
    summary = summarise_best_cases(df, metric_specs, args.include_columns)

    if summary.empty:
        print("No summary results generated with the provided filters.")
        return

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        summary.to_csv(args.output, index=False)
        print(f"Summary written to {args.output}")

    float_format = f"{{:.{args.float_precision}f}}".format
    with pd.option_context("display.float_format", float_format):
        print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
