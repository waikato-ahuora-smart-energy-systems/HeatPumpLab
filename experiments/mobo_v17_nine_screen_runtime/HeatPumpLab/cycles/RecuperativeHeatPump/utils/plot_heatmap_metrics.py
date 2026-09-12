"""Generate refrigerant heatmaps for Recuperative Heat Pump results."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

__all__ = ["get_heatmap_results_comparison"]


DEFAULT_METRICS_AND_OBJ_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("COP", "max"),
    ("SEC", "min"),
    ("VHC", "max"),
    ("eDTL", "min"),
    ("SLE", "min"),
    ("LMTD_gcs", "min"),
    ("Pmax", "min"),
    ("COP_pct_increase_vs_base", "max"),
    ("Pmax_reduction_vs_base", "max"),
    ("dT_gc_out_min_offset", "max"),
)
MIN_REQ_METRIC_COLUMNS: set[str] = {
    "refrigerant",
    "dT_rhx1_cold",
    "dT_rhx2_cold",
    "comp2_power_frac",
    "COP",
    "Pmax",
}
DEFAULT_REFRIGERANTS: List[str] = (
    "Butane",
    "Carbon dioxide",
    "R1234ze(Z)",
    "DME",
)
AXIS_LABELS: dict[str, str] = {
    "dT_rhx1_cold": r"$\Delta T_{\mathrm{RHX1,sink}}$ [K]",
    "dT_rhx2_cold": r"$\Delta T_{\mathrm{RHX2,sink}}$ [K]",
    "COP": "COP",
    "SEC": "SEC",
    "VHC": r"VHC [MJ/m$^3$]",
    "eDTL": r"ETL [K]",
    "SLE": r"SLE [kJ$_e$/(kJ$_h$K)]",
    "dQ/dW": r"$\Delta Q/\Delta W$",
    "LMTD_gcs": r"MTD [K]",
    "Pmax": r"$P_{\max}$ [kPa]",
    "COP_pct_increase_vs_base": r"COP Increase [%]",
    "Pmax_reduction_vs_base": r"$P_{\max}$ Reduction [%]",
    "dT_gc_out_min_offset": r"Subcooling Offset [K]",
}
DEFAULT_COMP_FRAC_2: List[float] = [0.5]
BASECASE_DTRHX1_COLD = 0.0
BASECASE_DTRHX2_COLD = 0.0
BASECASE_COMP2_FRAC = 0.5


def _parse_args(
    refrigerant_ls: List[str] = DEFAULT_REFRIGERANTS,
    metrics_ls: List[Tuple[str, str]] = DEFAULT_METRICS_AND_OBJ_PAIRS,
    comp2_load_frac_ls: List[float] = DEFAULT_COMP_FRAC_2,
) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create heatmap figures for metric performance versus "
            "dT_rhx1_cold and dT_rhx2_cold."
        )
    )
    parser.add_argument(
        "--metrics",
        nargs="+",
        default=[f"{metric}:{mode}" for metric, mode in metrics_ls],
        help=(
            "Metrics to include in the heatmaps. "
            "Use METRIC or METRIC:max/METRIC:min to set optimisation direction."
        ),
    )
    parser.add_argument(
        "--comp2-values",
        nargs="+",
        type=float,
        default=comp2_load_frac_ls,
        help=(
            "Specific comp2_power_frac values to include. "
            "Defaults to a predefined set of representative points."
        ),
    )
    parser.add_argument(
        "--refrigerants",
        nargs="+",
        default=refrigerant_ls,
        help=(
            "Specific refrigerants to include. "
            "Defaults to a predefined list of refrigerants."
        ),
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=900,
        help="Output figure resolution.",
    )
    parser.add_argument(
        "--cmap",
        default="viridis",
        help="Matplotlib colormap to use for heatmaps.",
    )
    parser.add_argument(
        "--log-pmax",
        dest="log_pmax",
        action="store_true",
        default=False,
        help="Render Pmax heatmaps with a logarithmic colour scale (default: on).",
    )
    parser.add_argument(
        "--no-log-pmax",
        dest="log_pmax",
        action="store_false",
        help="Disable logarithmic colour scaling for Pmax heatmaps.",
    )
    return parser.parse_args()


def _load_results(
    csv_path: Path,
) -> pd.DataFrame:
    from io import StringIO

    # Load with low_memory=False, remove bad rows, then reload as StringIO with low_memory=True
    df = pd.read_csv(csv_path, low_memory=False)

    # Remove rows where problematic columns contain their own header strings
    problematic_cols = [
        "dT_rhx1_cold",
        "dT_rhx2_cold",
        "comp2_power_frac",
        "COP",
        "Pmax",
    ]
    for col in problematic_cols:
        if col in df.columns:
            df = df[df[col].astype(str).str.strip() != col]

    # Reload from cleaned CSV to let pandas infer types properly
    csv_str = StringIO(df.to_csv(index=False))
    df = pd.read_csv(csv_str, low_memory=True)

    # Apply filters
    df = df[df["Success"].astype(str).str.lower().isin(["true", "1"])]

    return df


def _unique_sorted(values: Iterable[float]) -> List[float]:
    """Return sorted unique floating point values with tolerance handling."""
    rounded = {}
    for value in values:
        if value is None:
            continue
        try:
            float_value = float(value)
        except (TypeError, ValueError):
            continue
        if math.isnan(float_value):
            continue
        key = round(float_value, 8)
        rounded.setdefault(key, float_value)
    return sorted(rounded.values())


def _sanitize_float(value: float) -> str:
    """Sanitise a float for filesystem-safe labels."""
    return f"{value:.4f}".rstrip("0").rstrip(".").replace("-", "neg").replace(".", "p")


def _sanitize_metric_name(
    metric: str,
) -> str:
    """Sanitise a metric name for filesystem paths."""
    cleaned = metric.strip().replace("/", "_per_").replace(" ", "_")
    return cleaned or "metric"


def _parse_metric_specs(
    metric_specs: Sequence[str],
) -> List[Tuple[str, str]]:
    """Parse metric specifications into (name, mode) pairs."""
    parsed: List[Tuple[str, str]] = []
    for spec in metric_specs:
        parts = spec.split(":", 1)
        name = parts[0].strip()
        if not name:
            raise ValueError("Metric name cannot be empty.")
        mode = (
            parts[1].strip().lower() if len(parts) > 1 and parts[1].strip() else "max"
        )
        if mode not in {"min", "max"}:
            raise ValueError(
                f"Invalid optimisation mode '{mode}' for metric '{name}'. "
                "Use 'min' or 'max'."
            )
        parsed.append((name, mode))
    return parsed


def _format_axis_label(
    name: str,
    axis_labels: Dict[str, str] = AXIS_LABELS,
) -> str:
    """Return a descriptive axis label for the given column name."""
    return axis_labels.get(name, name)


def _add_basecase_metrics(
    df: pd.DataFrame,
) -> pd.DataFrame:
    """Add base-case relative metrics for each refrigerant."""
    missing = MIN_REQ_METRIC_COLUMNS.difference(df.columns)
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


def _build_grid(
    df: pd.DataFrame,
    metric: str,
    x_col: str,
    y_col: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x_vals = _unique_sorted(pd.to_numeric(df[x_col], errors="coerce").dropna().unique())
    y_vals = _unique_sorted(pd.to_numeric(df[y_col], errors="coerce").dropna().unique())
    if not x_vals or not y_vals:
        return np.array([]), np.array([]), np.array([[]])

    pivot = df.pivot_table(
        index=y_col,
        columns=x_col,
        values=metric,
        aggfunc="mean",
    ).reindex(index=y_vals, columns=x_vals)

    grid = pivot.to_numpy(dtype=float)
    mask = np.isnan(grid)
    if mask.all():
        return np.array(x_vals), np.array(y_vals), np.array([[]])

    masked_grid = np.ma.masked_array(grid, mask=mask)
    return np.array(x_vals), np.array(y_vals), masked_grid


def _create_heatmaps(
    df: pd.DataFrame,
    metrics: Sequence[Tuple[str, str]],
    output_dir: Path,
    comp2_values: Sequence[float],
    refrigerants_filter: Sequence[str],
    cmap: str,
    dpi: int,
    use_log_pmax: bool,
) -> None:
    output_dir = Path(output_dir)
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
        _unique_sorted(df["comp2_power_frac"].dropna().unique())
        if "comp2_power_frac" in df.columns
        else []
    )
    if not available_comp2:
        raise ValueError("CSV is missing valid 'comp2_power_frac' values.")
    comp2_values = [
        value
        for value in comp2_values
        if any(math.isclose(value, avail, rel_tol=1e-9) for avail in available_comp2)
    ]
    if not comp2_values:
        raise ValueError(
            "No matching comp2_power_frac values after applying the provided filter."
        )

    for metric_name, mode in metrics:
        if metric_name not in df.columns:
            print(f"Skipping metric '{metric_name}' (not found in data).")
            continue

        metric_values = pd.to_numeric(df[metric_name], errors="coerce")
        valid_metric_values = metric_values.replace([np.inf, -np.inf], np.nan).dropna()
        if valid_metric_values.empty:
            print(f"Skipping metric '{metric_name}' (no valid data).")
            continue

        vmin = valid_metric_values.min()
        vmax = valid_metric_values.max()

        for comp2_value in comp2_values:
            comp_mask = np.isclose(
                pd.to_numeric(df["comp2_power_frac"], errors="coerce"),
                comp2_value,
            )
            subset = df.loc[comp_mask]

            for refrigerant in refrigerants:
                refrigerant_slice = subset[subset["refrigerant"] == refrigerant]
                if refrigerant_slice.empty:
                    continue

                x_vals, y_vals, grid = _build_grid(
                    refrigerant_slice, metric_name, "dT_rhx1_cold", "dT_rhx2_cold"
                )

                if grid.size == 0:
                    continue

                X, Y = np.meshgrid(x_vals, y_vals)
                fig, ax = plt.subplots(figsize=(4.5, 3.75))

                current_norm = None
                grid_for_plot = grid
                if metric_name.lower() == "pmax" and use_log_pmax:
                    positive_values = valid_metric_values[valid_metric_values > 0]
                    if not positive_values.empty:
                        log_vmin = float(positive_values.min())
                        log_vmax = float(positive_values.max())
                        log_vmin = max(log_vmin, 1e-6)
                        if log_vmax <= log_vmin:
                            log_vmax = log_vmin * 1.05
                        current_norm = LogNorm(vmin=log_vmin, vmax=log_vmax)
                        grid_for_plot = np.ma.masked_where(grid <= 0, grid)
                        if np.ma.getmaskarray(grid_for_plot).all():
                            current_norm = None
                            grid_for_plot = grid

                quad_mesh = ax.pcolormesh(
                    X,
                    Y,
                    grid_for_plot,
                    shading="auto",
                    cmap=cmap,
                    norm=current_norm,
                    vmin=None if current_norm else vmin,
                    vmax=None if current_norm else vmax,
                )
                ax.set_title(
                    rf"{_format_axis_label(metric_name)}: {refrigerant} | $\omega_2$ = {comp2_value:.1f}"
                )
                ax.set_xlabel(_format_axis_label("dT_rhx1_cold"))
                ax.set_ylabel(_format_axis_label("dT_rhx2_cold"))

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
                cbar.set_label(_format_axis_label(metric_name))
                fig.tight_layout()

                metric_dir = output_dir / _sanitize_metric_name(metric_name)
                refrigerant_dir = metric_dir / refrigerant
                comp_dir = refrigerant_dir / f"comp2_{_sanitize_float(comp2_value)}"
                comp_dir.mkdir(parents=True, exist_ok=True)

                output_path = comp_dir / "heatmap.png"
                fig.savefig(output_path, dpi=dpi)
                plt.close(fig)


def get_heatmap_results_comparison(
    filepath: Path = Path(__file__).resolve().parents[1]
    / "results"
    / "RecuperativeHP_results.csv",
    output_dir: Path = Path(__file__).resolve().parents[1]
    / "results"
    / "figures"
    / "heatmaps",
    metrics_ls: List[Tuple[str, str]] = DEFAULT_METRICS_AND_OBJ_PAIRS,
    refrigerant_ls: List[str] = DEFAULT_REFRIGERANTS,
    comp2_load_frac_ls: List[float] = DEFAULT_COMP_FRAC_2,
) -> None:
    args = _parse_args(
        refrigerant_ls=refrigerant_ls,
        metrics_ls=metrics_ls,
        comp2_load_frac_ls=comp2_load_frac_ls,
    )
    df = _load_results(
        csv_path=filepath,
    )
    df = _add_basecase_metrics(df)
    metric_settings = _parse_metric_specs(args.metrics)
    _create_heatmaps(
        df,
        metrics=metric_settings,
        output_dir=output_dir,
        comp2_values=args.comp2_values,
        refrigerants_filter=args.refrigerants,
        cmap=args.cmap,
        dpi=args.dpi,
        use_log_pmax=args.log_pmax,
    )
