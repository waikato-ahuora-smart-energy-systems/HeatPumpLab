"""Visualise how the comp2_power_frac split influences key performance metrics."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

__all__ = ["get_comp2_power_frac_comparison_analysis"]

DEFAULT_METRICS: Tuple[Tuple[str, str], ...] = (
    ("COP_pct_increase_vs_base", "max"),
    # ("VHC", "max"),
    ("Pmax_reduction_vs_base", "max"),
    # ("dQ_per_dW", "max"),
)

BASECASE_COMP2_FRAC = 0.5


METRIC_LABELS: dict[str, str] = {
    "COP_pct_increase_vs_base": r"COP Increase [%]",
    "VHC": "Volumetric heating capacity",
    "Pmax_reduction_vs_base": r"$P_{\max}$ Reduction [%]",
    "Pmax": "Maximum pressure [kPa]",
    "dQ_per_dW": "ΔQ/ΔW",
}

DEFAULT_REFRIGERANTS: Sequence[str] = (
    # "Ammonia",
    # "Propane",
    "Butane",
    # "Isobutane",
    # "Pentane",
    # "Isopentane",
    # "Neopentane",
    # "Cyclopentane",
    "Carbon dioxide",
    # "R1233zd(E)",
    # "R1336mzz(E)",
    # "R1336mzz(Z)",
    # "R1234ze(E)",
    "R1234ze(Z)",
    # "R1234yf",
    # "R1243zf",
    # "R1224yd(Z)",
    # "PROPYLEN",
    "DME",
    # "Methanol",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create per-refrigerant plots showing how performance metrics vary with "
            "the comp2_power_frac split."
        )
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "results"
        / "RecuperativeHP_results.csv",
        help="Path to the optimisation results CSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "results"
        / "figures"
        / "comp2_power_analysis",
        help="Directory where generated figures will be written.",
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=[f"{metric}:{mode}" for metric, mode in DEFAULT_METRICS],
        help=(
            "Metrics to include. Use METRIC or METRIC:max/METRIC:min to set the "
            "optimisation direction for identifying the best cases."
        ),
    )
    parser.add_argument(
        "--refrigerants",
        nargs="+",
        default=list(DEFAULT_REFRIGERANTS),
        help="Restrict the analysis to these refrigerants.",
    )
    parser.add_argument(
        "--comp2-values",
        nargs="+",
        type=float,
        help="Restrict the analysis to these exact comp2_power_frac values.",
    )
    parser.add_argument(
        "--skip-success-filter",
        action="store_true",
        help="Include optimisation runs flagged as unsuccessful.",
    )
    parser.add_argument(
        "--skip-converged-filter",
        action="store_true",
        help="Include optimisation runs flagged as non-converged.",
    )
    parser.add_argument(
        "--include-basecase-line",
        action="store_true",
        help="Draw a vertical reference line at the nominal 0.5 split.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution to use for saved figures.",
    )
    parser.add_argument(
        "--fig-width",
        type=float,
        default=4.5,
        help="Width (in inches) for the output figure.",
    )
    parser.add_argument(
        "--fig-height",
        type=float,
        default=3.75,
        help="Height (in inches) for the output figure.",
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
        np.isclose(df["dT_rhx1_cold"].astype(float), 0.0)
        & np.isclose(df["dT_rhx2_cold"].astype(float), 0.0)
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
        if not isinstance(value, (int, float)) or np.isnan(value):
            continue
        key = round(float(value), 8)
        rounded.setdefault(key, float(value))
    return sorted(rounded.values())


def sanitize_metric_name(metric: str) -> str:
    cleaned = metric.strip().replace("/", "_per_").replace(" ", "_")
    return cleaned or "metric"


def parse_metric_specs(metric_specs: Sequence[str]) -> List[Tuple[str, str]]:
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


def filter_dataframe(
    df: pd.DataFrame,
    refrigerants: Sequence[str] | None,
    comp2_values: Sequence[float] | None,
) -> pd.DataFrame:
    if refrigerants:
        df = df[df["refrigerant"].isin(refrigerants)]
    if comp2_values:
        allowed = {round(float(v), 8) for v in comp2_values}
        mask = df["comp2_power_frac"].apply(lambda v: round(float(v), 8) in allowed)
        df = df[mask]
    return df


def plot_metric_for_refrigerants(
    df: pd.DataFrame,
    metric: str,
    mode: str,
    include_basecase_line: bool,
    fig_width: float,
    fig_height: float,
    dpi: int,
    output_dir: Path,
) -> None:
    if metric not in df.columns:
        print(f"Skipping metric '{metric}': column not found in dataset.")
        return

    refrigerants = sorted(df["refrigerant"].dropna().unique())
    if not refrigerants:
        print("No refrigerant data available after filtering; nothing to plot.")
        return

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    best_direction = "max" if mode == "max" else "min"
    ylabel = METRIC_LABELS.get(metric, metric)
    plotted = False

    for refrigerant in refrigerants:
        subset = df[df["refrigerant"] == refrigerant]
        if subset.empty:
            continue

        subset = subset.copy()
        subset["comp2_power_frac"] = pd.to_numeric(
            subset["comp2_power_frac"],
            errors="coerce",
        )
        subset[metric] = pd.to_numeric(subset[metric], errors="coerce")
        subset = subset.dropna(subset=["comp2_power_frac", metric])
        if subset.empty:
            continue

        grouped = subset.groupby("comp2_power_frac")[metric]
        if mode == "max":
            aggregated = grouped.max().sort_index()
        else:
            aggregated = grouped.min().sort_index()

        if aggregated.empty:
            continue

        (line,) = ax.plot(
            aggregated.index,
            aggregated.values,
            marker="o",
            linewidth=1.5,
            label=refrigerant,
        )
        plotted = True
        line_color = line.get_color()

        if mode == "max":
            best_frac = aggregated.idxmax()
        else:
            best_frac = aggregated.idxmin()
        best_metric = aggregated.loc[best_frac]
        ax.scatter(
            [best_frac],
            [best_metric],
            color=line_color,
            s=60,
            zorder=5,
            marker="s",
            label=None,
        )
    if not plotted:
        plt.close(fig)
        print(f"No data available to plot for metric '{metric}'.")
        return

    if include_basecase_line:
        ax.axvline(
            BASECASE_COMP2_FRAC,
            color="0.35",
            linestyle="--",
            linewidth=1.0,
            label="Nominal split",
        )
    x_label = r"$\omega_{2}$ [-]"
    ax.set_xlabel(x_label)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.2)
    ax.legend(loc="best", fontsize=9)

    fig.tight_layout()

    output_dir.mkdir(parents=True, exist_ok=True)
    filename = f"{sanitize_metric_name(metric)}_best_vs_comp2_power_frac.png"
    path = output_dir / filename
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {path}")


def get_comp2_power_frac_comparison_analysis() -> None:
    args = parse_args()
    metrics = parse_metric_specs(args.metrics)
    df = load_results(
        args.csv,
        filter_success=not args.skip_success_filter,
        filter_converged=not args.skip_converged_filter,
    )

    basecase_columns = {"COP_pct_increase_vs_base", "Pmax_reduction_vs_base"}
    missing_basecase = [col for col in basecase_columns if col not in df.columns]
    if missing_basecase:
        df = add_basecase_metrics(df)

    required_columns = {"comp2_power_frac", "refrigerant"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(
            f"Dataset is missing required columns: {sorted(missing_columns)}"
        )

    df = filter_dataframe(df, args.refrigerants, args.comp2_values)
    if df.empty:
        raise ValueError("No data remains after applying requested filters.")

    df["comp2_power_frac"] = pd.to_numeric(
        df["comp2_power_frac"],
        errors="coerce",
    )
    unique_comp_fractions = unique_sorted(df["comp2_power_frac"].unique())
    print(
        "Analysing "
        f"{len(df)} cases across {len(unique_comp_fractions)} comp2_power_frac values."
    )

    for metric, mode in metrics:
        plot_metric_for_refrigerants(
            df=df,
            metric=metric,
            mode=mode,
            include_basecase_line=args.include_basecase_line,
            fig_width=args.fig_width,
            fig_height=args.fig_height,
            dpi=args.dpi,
            output_dir=args.output_dir,
        )


METRIC_LABELS: dict[str, str] = {
    "COP_pct_increase_vs_base": "COP increase vs base [%]",
    "VHC": "Volumetric heating capacity",
    "Pmax_reduction_vs_base": "Pmax reduction vs base [%]",
    "Pmax": "Maximum pressure [kPa]",
    "dQ_per_dW": "ΔQ/ΔW",
}
