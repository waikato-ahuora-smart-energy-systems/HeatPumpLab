"""Replay optimised designs and compare COP with suction volumetric capacity."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor, as_completed
import json
import math
from pathlib import Path

from ..optimisation._types import OPTIMISER_VERSION


def _replay_design(row: dict, results_dir: str) -> dict:
    from ..optimisation._objective import _OuterDesignObjective
    from ..optimisation._types import DEFAULT_BASE_CASE

    base_case = {
        key: (
            str(row[key]).lower() == "true"
            if isinstance(default, bool) else type(default)(row[key])
        )
        for key, default in DEFAULT_BASE_CASE.items()
    }
    objective = _OuterDesignObjective(
        max_evals=1,
        base_case=base_case,
        refrigerant=row["input_refrigerant"],
        results_dir=results_dir,
        exact_evaluation_timeout_seconds=180.0,
        exact_evaluation_timeout_retries=1,
        exact_evaluation_worker_max_tasks=1,
    )
    x = tuple(float(row[f"best_{key}"]) for key in (
        "dT_rhx1_cold", "dT_rhx2_cold", "comp2_power_frac",
    )) + tuple(float(value) for value in row["best_mass_fractions"].split(";"))
    try:
        record = objective.evaluate(x)
    finally:
        objective.close()
    vhc = float(record.report.get("VHC", math.nan)) / 1000.0
    if not record.success or not math.isfinite(vhc) or vhc <= 0:
        raise RuntimeError(f"VHC replay failed for {row['input_refrigerant']}: {record.message}")
    return {
        "case_signature": row["case_signature"],
        "optimiser_version": OPTIMISER_VERSION,
        "COP": record.cop,
        "VHC_MJ_per_m3": vhc,
        "SEC": record.sec,
        "Pmax_MPa": record.pmax / 1000.0,
        "COP_change_vs_optimisation_percent": 100 * (record.cop / float(row["best_COP"]) - 1),
        "T_comp_discharge_C": record.report.get("T_comp_discharge"),
    }


def replay_designs(rows: list[dict], output_dir: str | Path, n_cpu: int = 3) -> dict:
    """Persist replays by original signature, invalidating obsolete model metrics."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_path = output_dir / "design_replays.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    cache = {key: value for key, value in cache.items()
             if value.get("optimiser_version") == OPTIMISER_VERSION}
    pending = [row for row in rows if row["case_signature"] not in cache and str(row["batch_success"]).lower() == "true"]
    if not pending:
        return cache
    with ProcessPoolExecutor(max_workers=max(1, n_cpu)) as executor:
        futures = {executor.submit(_replay_design, row, str(output_dir)): row for row in pending}
        for future in as_completed(futures):
            row = futures[future]
            replay = future.result()
            cache[row["case_signature"]] = replay
            temporary = cache_path.with_suffix(".tmp")
            temporary.write_text(json.dumps(cache, indent=2, allow_nan=False) + "\n")
            temporary.replace(cache_path)
            print(f"Replayed COP/VHC: {row['input_refrigerant']}", flush=True)
    return cache


def nondominated_mask(points: list[tuple[float, float]]) -> list[bool]:
    """Maximise both COP and VHC; ties remain non-dominated, invalid points do not."""
    valid = [all(math.isfinite(value) and value > 0 for value in point) for point in points]
    return [
        valid[i] and not any(
            valid[j] and other[0] >= point[0] and other[1] >= point[1]
            and (other[0] > point[0] or other[1] > point[1])
            for j, other in enumerate(points)
        )
        for i, point in enumerate(points)
    ]


def add_pareto_metrics(rows: list[dict], replays: dict) -> None:
    """Attach paired replay metrics and the global sampled Pareto classification."""
    for row in rows:
        replay = replays.get(row["case_signature"])
        if row["success"] and replay is None:
            raise ValueError(f"Missing COP/VHC replay for {row['input_refrigerant']}")
        row["replay_COP"] = replay["COP"] if replay else None
        row["VHC_MJ_per_m3"] = replay["VHC_MJ_per_m3"] if replay else None
        row["COP_replay_change_percent"] = replay["COP_change_vs_optimisation_percent"] if replay else None
    mask = nondominated_mask([
        (row["replay_COP"] or math.nan, row["VHC_MJ_per_m3"] or math.nan)
        for row in rows
    ])
    for row, on_front in zip(rows, mask):
        row["non_dominated"] = on_front


def plot_cop_vhc(
    rows: list[dict], output_dir: str | Path, target_temperature: float,
    *, status_note: str = "",
) -> None:
    """Plot composition trajectories and the front among best-found COP designs."""
    from matplotlib import rc_context
    from matplotlib.figure import Figure

    output_dir = Path(output_dir)
    with rc_context({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False}):
        fig = Figure(figsize=(12, 8), layout="constrained")
        ax = fig.subplots()
        for family, color, marker, offset in zip(
            ("R290 / R600", "R290 / R601A", "R290 / R601", "R290 / R600A"),
            ("#176B87", "#B65320", "#6B4DA0", "#34804A"),
            ("o", "s", "^", "P"),
            ((0, -15), (0, 9), (0, -15), (0, 9)),
        ):
            selected = [row for row in rows if row["family"] == family]
            ax.plot(
                [row["VHC_MJ_per_m3"] for row in selected],
                [row["replay_COP"] for row in selected],
                color=color, marker=marker, markersize=6, linewidth=1.2, label=family,
            )
            for row in selected:
                percent = row["R290_mass_percent"]
                if row["success"] and 0 < percent < 100:
                    label_offset = {
                        ("R290 / R600", 50): (14, -12),
                        ("R290 / R600", 90): (-10, -14),
                        ("R290 / R601A", 80): (-23, -2),
                        ("R290 / R601", 70): (-6, 10),
                        ("R290 / R601", 80): (0, 10),
                        ("R290 / R601", 90): (12, 7),
                        ("R290 / R600A", 40): (-16, -20),
                        ("R290 / R600A", 80): (0, -17),
                        ("R290 / R600A", 90): (18, -15),
                    }.get((family, percent), offset)
                    ax.annotate(f"{percent}%", (row["VHC_MJ_per_m3"], row["replay_COP"]),
                                xytext=label_offset, textcoords="offset points", ha="center", fontsize=9, color=color,
                                bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 0.2})

        unique = {row["case_signature"]: row for row in rows}
        frontier = sorted((row for row in unique.values() if row["non_dominated"]), key=lambda row: row["VHC_MJ_per_m3"])
        ax.plot([row["VHC_MJ_per_m3"] for row in frontier], [row["replay_COP"] for row in frontier],
                color="#222222", linestyle="--", linewidth=1.3, marker="o", markersize=10,
                markerfacecolor="none", markeredgewidth=1.2, label="Non-dominated sampled designs", zorder=5)
        for fluid, label in (("Butane", "R600"), ("Propane", "R290 (100%)"), ("Isopentane", "R601A (0%)"), ("Pentane", "R601 (0%)"), ("Isobutane", "R600A (0%)")):
            row = next(row for row in unique.values() if row["input_refrigerant"] == fluid)
            if row["success"]:
                ax.scatter(row["VHC_MJ_per_m3"], row["replay_COP"], marker="D", s=48, color="#333333", zorder=6)
                ax.annotate(label, (row["VHC_MJ_per_m3"], row["replay_COP"]), xytext=(7, 8),
                            textcoords="offset points", fontsize=10, fontweight="bold")
        ax.set(xlabel="Suction volumetric heating capacity (MJ/m³)", ylabel="Heating COP",
               title=f"COP–VHC trade-off at {target_temperature:g}°C" + (f"\n{status_note}" if status_note else ""))
        ax.margins(x=0.12, y=0.13)
        ax.grid(alpha=0.2)
        ax.legend(loc="best", frameon=True, fontsize=10)
        fig.supxlabel("Labels show R290 mass percentage. Each composition is fixed during COP optimisation.\n"
                      "Dashed front compares best-found designs; connecting lines are visual guides.", fontsize=10)
        fig.savefig(output_dir / "cop_vhc_pareto.png", dpi=200)
        fig.savefig(output_dir / "cop_vhc_pareto.svg")
