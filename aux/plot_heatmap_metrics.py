"""Generate refrigerant heatmaps for Recuperative Heat Pump results."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DEFAULT_METRICS: Tuple[Tuple[str, str], ...] = (
    ("COP", "max"),
    ("SEC", "min"),
    ("VHC", "max"),
    ("eDTL", "min"),
    ("SLE", "min"),
    ("dQ_per_dW", "max"),
    ("LMTD_gcs", "max"),
    ("Pmax", "min"),
    ("COP_pct_increase_vs_base", "max"),
    ("Pmax_reduction_vs_base", "max"),
)

DEFAULT_REFRIGERANTS: Sequence[str] = (
    "Ammonia",
    "Propane",
    "Butane",
    "Isobutane",
    "Pentane",
    "Isopentane",
    "Neopentane",
    "Cyclopentane",
    "Carbon dioxide",
    "R1233zd(E)",
    "R1336mzz(E)",
    "R1336mzz(Z)",
    "R1234ze(E)",
    "R1234ze(Z)",
    "R1234yf",
    "R1243zf",
    "R1224yd(Z)",
    "PROPYLEN",
    "DME",
    "Methanol",
)


AXIS_LABELS: dict[str, str] = {
    "dT_rhx1_cold": r"$\Delta T_{\mathrm{RHX1,sink}}$ [K]",
    "dT_rhx2_cold": r"$\Delta T_{\mathrm{RHX2,sink}}$ [K]",
    "COP": "COP",
    "SEC": "SEC",
    "VHC": r"VHC [MJ/m$^3$]",
    "eDTL": r"Effective Temperature Lift [K]",
    "SLE": r"SLE [kJ$_e$/(kJ$_h$K)]",
    "dQ/dW": r"$\Delta Q/\Delta W$",
    "LMTD_gcs": r"LMTD$_{\mathrm{GC}}$ [K]",
    "Pmax": r"$P_{\max}$ [kPa]",
    "COP_pct_increase_vs_base": r"COP Increase [%]",
    "Pmax_reduction_vs_base": r"$P_{\max}$ Reduction [%]",
}



DEFAULT_COMP_FRAC_2: Sequence[float] = (0.3, 0.4, 0.5, 0.6, 0.7)

BASECASE_DTRHX1_COLD = 0.0
BASECASE_DTRHX2_COLD = 0.0
BASECASE_COMP2_FRAC = 0.5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create heatmap figures for metric performance versus "
            "dT_rhx1_cold and dT_rhx2_cold."
        )
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "RHP_results.csv",
        help="Path to the results CSV file.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "results" / "figures" / "heatmaps",
        help="Directory where plots will be saved.",
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=[f"{metric}:{mode}" for metric, mode in DEFAULT_METRICS],
        help=(
            "Metrics to include in the heatmaps. "
            "Use METRIC or METRIC:max/METRIC:min to set optimisation direction."
        ),
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
        "--comp2-values",
        nargs="+",
        type=float,
        default=list(DEFAULT_COMP_FRAC_2),
        help=(
            "Specific comp2_power_frac values to include. "
            "Defaults to a predefined set of representative points."
        ),
    )
    parser.add_argument(
        "--refrigerants",
        nargs="+",
        default=list(DEFAULT_REFRIGERANTS),
        help=(
            "Specific refrigerants to include. "
            "Defaults to a predefined list of refrigerants."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Output figure resolution.",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap to use for heatmaps.",
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


def unique_sorted(values: Iterable[float]) -> List[float]:
    """Return sorted unique floating point values with tolerance handling."""
    rounded = {}
    for value in values:
        if not isinstance(value, (int, float)) or math.isnan(value):
            continue
        key = round(float(value), 8)
        rounded.setdefault(key, float(value))
    return sorted(rounded.values())


def sanitize_float(value: float) -> str:
    """Sanitise a float for filesystem-safe labels."""
    return f"{value:.4f}".rstrip("0").rstrip(".").replace("-", "neg").replace(".", "p")


def sanitize_metric_name(metric: str) -> str:
    """Sanitise a metric name for filesystem paths."""
    cleaned = metric.strip().replace("/", "_per_").replace(" ", "_")
    return cleaned or "metric"


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


def format_axis_label(name: str) -> str:
    """Return a descriptive axis label for the given column name."""
    return AXIS_LABELS.get(name, name)


def add_basecase_metrics(df: pd.DataFrame) -> pd.DataFrame:
    """Add base-case relative metrics for each refrigerant."""
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
        print(
            "Warning: No base-case rows found; base-relative metrics will be NaN."
        )
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


def build_grid(
    df: pd.DataFrame,
    metric: str,
    x_col: str,
    y_col: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_vals = unique_sorted(df[x_col].dropna().unique())
    y_vals = unique_sorted(df[y_col].dropna().unique())
    if not x_vals or not y_vals:
        return np.array([]), np.array([]), np.array([[]])

    pivot = (
        df.pivot_table(
            index=y_col,
            columns=x_col,
            values=metric,
            aggfunc="mean",
        )
        .reindex(index=y_vals, columns=x_vals)
    )

    grid = pivot.to_numpy(dtype=float)
    mask = np.isnan(grid)
    if mask.all():
        return np.array(x_vals), np.array(y_vals), np.array([[]])

    masked_grid = np.ma.masked_array(grid, mask=mask)
    return np.array(x_vals), np.array(y_vals), masked_grid


def plot_metric_heatmaps(
    df: pd.DataFrame,
    metrics: Sequence[Tuple[str, str]],
    output_dir: Path,
    comp2_values: Sequence[float],
    refrigerants_filter: Sequence[str],
    cmap: str,
    dpi: int,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    refrigerants = sorted(
        df["refrigerant"].dropna().unique().tolist()
        if "refrigerant" in df.columns
        else []
    )
    if not refrigerants:
        raise ValueError("No refrigerant data available in the provided CSV.")
    refrigerant_filter_set = set(refrigerants_filter)
    refrigerants = [r for r in refrigerants if r in refrigerant_filter_set]
    if not refrigerants:
        raise ValueError("No matching refrigerants after applying the provided filter.")
    df = df[df["refrigerant"].isin(refrigerants)]

    available_comp2 = (
        unique_sorted(df["comp2_power_frac"].dropna().unique())
        if "comp2_power_frac" in df.columns
        else []
    )
    if not available_comp2:
        raise ValueError("CSV is missing valid 'comp2_power_frac' values.")
    comp2_values = [
        value for value in comp2_values if any(math.isclose(value, avail, rel_tol=1e-9) for avail in available_comp2)
    ]
    if not comp2_values:
        raise ValueError("No matching comp2_power_frac values after applying the provided filter.")

    for metric_name, mode in metrics:
        if metric_name not in df.columns:
            print(f"Skipping metric '{metric_name}' (not found in data).")
            continue

        metric_values = df[metric_name].astype(float)
        valid_metric_values = metric_values.replace([np.inf, -np.inf], np.nan).dropna()
        if valid_metric_values.empty:
            print(f"Skipping metric '{metric_name}' (no valid data).")
            continue

        vmin = valid_metric_values.min()
        vmax = valid_metric_values.max()

        for comp2_value in comp2_values:
            comp_mask = np.isclose(df["comp2_power_frac"].astype(float), comp2_value)
            subset = df.loc[comp_mask]
            if subset.empty:
                continue

            for refrigerant in refrigerants:
                refrigerant_slice = subset[subset["refrigerant"] == refrigerant]
                if refrigerant_slice.empty:
                    continue

                x_vals, y_vals, grid = build_grid(
                    refrigerant_slice, metric_name, "dT_rhx1_cold", "dT_rhx2_cold"
                )

                if grid.size == 0:
                    continue

                X, Y = np.meshgrid(x_vals, y_vals)
                fig, ax = plt.subplots(figsize=(4.5, 3.75))
                quad_mesh = ax.pcolormesh(
                    X,
                    Y,
                    grid,
                    shading="auto",
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                )
                ax.set_title(
                    f"{format_axis_label(metric_name)}: {refrigerant} | $\omega_2$ = {comp2_value:.1f}"
                )
                ax.set_xlabel(format_axis_label("dT_rhx1_cold"))
                ax.set_ylabel(format_axis_label("dT_rhx2_cold"))

                grid_filled = grid.filled(np.nan)
                if np.isfinite(grid_filled).any():
                    with np.errstate(all="ignore"):
                        if mode == "max":
                            best_index = np.nanargmax(grid_filled)
                            marker_color = "gold"
                        else:
                            best_index = np.nanargmin(grid_filled)
                            marker_color = "gold"
                    best_y_idx, best_x_idx = np.unravel_index(best_index, grid.shape)

                    x_star = x_vals[best_x_idx]
                    y_star = y_vals[best_y_idx]
                    ax.scatter(
                        x_star,
                        y_star,
                        marker="*",
                        s=140,
                        color=marker_color,
                        edgecolors="black",
                        linewidths=0.6,
                        zorder=5,
                    )
                    z_value = grid_filled[best_y_idx, best_x_idx]
                    ax.annotate(
                        f"{z_value:.2f}",
                        xy=(x_star, y_star),
                        xytext=(6, 6),
                        textcoords="offset points",
                        color="black",
                        fontsize=10,
                        fontweight="bold",
                        bbox=dict(
                            boxstyle="round,pad=0.2",
                            facecolor="white",
                            edgecolor="black",
                            linewidth=0.75,
                            alpha=0.85,
                        ),
                    )

                cbar = fig.colorbar(quad_mesh, ax=ax, shrink=0.9)
                cbar.set_label(format_axis_label(metric_name))
                fig.tight_layout()

                metric_dir = output_dir / sanitize_metric_name(metric_name)
                refrigerant_dir = metric_dir / refrigerant
                comp_dir = refrigerant_dir / f"comp2_{sanitize_float(comp2_value)}"
                comp_dir.mkdir(parents=True, exist_ok=True)

                output_path = comp_dir / "heatmap.png"
                fig.savefig(output_path, dpi=dpi)
                plt.close(fig)


def main() -> None:
    args = parse_args()
    df = load_results(
        args.csv,
        filter_success=not args.skip_success_filter,
        filter_converged=not args.skip_converged_filter,
    )
    df = add_basecase_metrics(df)
    metric_settings = parse_metric_specs(args.metrics)
    plot_metric_heatmaps(
        df,
        metrics=metric_settings,
        output_dir=args.output_dir,
        comp2_values=args.comp2_values,
        refrigerants_filter=args.refrigerants,
        cmap=args.cmap,
        dpi=args.dpi,
    )


if __name__ == "__main__":
    main()
