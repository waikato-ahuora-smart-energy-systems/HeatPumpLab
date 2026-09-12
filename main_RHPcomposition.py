"""Optimise fixed mass-composition sweeps for selected hydrocarbon refrigerants.

Run with ``.venv/bin/python main_RHPcomposition.py``. Completed cases resume
from the results CSV. Shared pure-fluid endpoints are solved only once.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import math
import os
from pathlib import Path

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import normalise_refrigerant_input
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    OuterDesignOptimisationConfig,
    optimise_recuperative_heat_pump_cases,
)

DEFAULT_RESULTS_DIR = "results/RHPcomposition"
DEFAULT_FILENAME = "RHP_fixed_mass_composition_bo_v18"
FIXED_COMPOSITION_CONFIG = OuterDesignOptimisationConfig(
    vary_composition=False,
    preserve_incumbent=False,
)


@dataclass(frozen=True)
class CompositionPoint:
    family: str
    propane_mass_percent: int
    refrigerant: str


def composition_points() -> tuple[CompositionPoint, ...]:
    """Return the R600 baseline and all eleven points for each binary family."""
    points = [CompositionPoint("R600", 0, "Butane")]
    for label, partner in (
        ("R600", "Butane"),
        ("R601A", "Isopentane"),
        ("R601", "Pentane"),
        ("R600A", "Isobutane"),
    ):
        for percent in range(0, 101, 10):
            if percent == 0:
                refrigerant = partner
            elif percent == 100:
                refrigerant = "Propane"
            else:
                refrigerant = (
                    f"Propane;{partner}|{percent / 100:.1f};"
                    f"{(100 - percent) / 100:.1f}|mass"
                )
            points.append(CompositionPoint(f"R290 / {label}", percent, refrigerant))
    return tuple(points)


def run_sweep(
    *,
    target_temperature: float = 200.0,
    n_cpu: int = min(6, max((os.cpu_count() or 2) - 4, 1)),
    results_dir: str | Path = DEFAULT_RESULTS_DIR,
) -> list[dict]:
    """Optimise 41 unique fixed-composition cases at one sink target (degrees C)."""
    refrigerants = list(dict.fromkeys(p.refrigerant for p in composition_points()))
    return optimise_recuperative_heat_pump_cases(
        filename=DEFAULT_FILENAME,
        param_kwargs={
            "refrigerant_ls": refrigerants,
            "T_air_tar": target_temperature,
        },
        config=FIXED_COMPOSITION_CONFIG,
        n_cpu=n_cpu,
        results_dir=results_dir,
    )


def load_sweep_rows(
    results_dir: str | Path = DEFAULT_RESULTS_DIR,
    target_temperature: float = 200.0,
) -> list[dict]:
    """Validate the complete fixed-composition sweep and expand shared endpoints."""
    source = Path(results_dir) / f"{DEFAULT_FILENAME}.csv"
    with source.open(newline="") as stream:
        results = {}
        for row in csv.DictReader(stream):
            if float(row["T_air_tar"]) != target_temperature:
                continue
            if row["config_vary_composition"].lower() != "false":
                raise ValueError("Results contain a case with variable composition.")
            key = row["input_refrigerant"]
            if key in results:
                raise ValueError(f"Duplicate result for {key} at {target_temperature} C.")
            results[key] = row

    output = []
    for point in composition_points():
        spec = normalise_refrigerant_input(point.refrigerant)
        key = spec.format()
        if key not in results:
            raise ValueError(f"Sweep is incomplete: missing {key} at {target_temperature} C.")
        row = results[key]
        success = row["batch_success"].lower() == "true"
        if success:
            actual = tuple(float(value) for value in row["best_mass_fractions"].split(";"))
            if len(actual) != len(spec.fractions) or any(
                not math.isclose(a, b, abs_tol=1e-10, rel_tol=0)
                for a, b in zip(actual, spec.fractions)
            ):
                raise ValueError(f"Optimisation changed composition for {key}.")
        record = {
            "family": point.family,
            "R290_mass_percent": point.propane_mass_percent,
            "input_refrigerant": key,
            "T_air_tar_C": target_temperature,
            "success": success,
            "case_signature": row["case_signature"],
        }
        for field, source_field, scale in (
            ("COP", "best_COP", 1),
            ("SEC_W_per_Qhot_minus_W", "best_SEC", 1),
            ("Pmax_MPa", "best_Pmax", 0.001),
            ("dT_rhx1_K", "best_dT_rhx1_cold", 1),
            ("dT_rhx2_K", "best_dT_rhx2_cold", 1),
            ("comp2_power_fraction", "best_comp2_power_frac", 1),
        ):
            record[field] = float(row[source_field]) * scale if success else None
        for field in ("exact_evaluations", "failed_evaluations", "timed_out_evaluations"):
            record[field] = int(row[field])
        record["elapsed_seconds"] = float(row["elapsed_seconds"])
        output.append(record)
    baseline = output[0]["COP"]
    for row in output:
        row["power_reduction_vs_R600_percent"] = (
            100 * (1 - baseline / row["COP"])
            if baseline is not None and row["COP"] is not None else None
        )
    return output


def write_analysis(
    results_dir: str | Path = DEFAULT_RESULTS_DIR,
    target_temperature: float = 200.0,
) -> Path:
    """Write the validated comparison table and six-panel scientific figure."""
    from matplotlib import rc_context
    from matplotlib.figure import Figure
    from HeatPumpLab.cycles.RecuperativeHeatPump.utils.composition_pareto import (
        add_pareto_metrics, plot_cop_vhc, replay_designs,
    )

    rows = load_sweep_rows(results_dir, target_temperature)
    output_dir = Path(results_dir) / f"analysis_{target_temperature:g}C"
    output_dir.mkdir(parents=True, exist_ok=True)
    with (Path(results_dir) / f"{DEFAULT_FILENAME}.csv").open(newline="") as stream:
        raw_rows = [row for row in csv.DictReader(stream) if float(row["T_air_tar"]) == target_temperature]
    replays = replay_designs(raw_rows, output_dir)
    add_pareto_metrics(rows, replays)
    with (output_dir / "composition_comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    plot_cop_vhc(rows, output_dir, target_temperature)
    unique = list({row["case_signature"]: row for row in rows}.values())
    successful = [row for row in unique if row["success"]]
    report = [
        f"# Fixed-composition refrigerant comparison at {target_temperature:g}°C",
        "",
        "![COP versus VHC](cop_vhc_pareto.png)",
        "",
        "Compositions are mass fractions. Each composition is held fixed while "
        "recuperator temperature lifts and compressor power split are optimised for COP. "
        "The gas-cooler outlet offset is optimised by the inner analysis.",
        "",
        f"{len(successful)}/{len(unique)} unique designs solved successfully. "
        "Shared pure-fluid endpoints are reused across curves. COP and VHC on the "
        "Pareto plot come from the same fresh-process replay of each selected design.",
        "",
        "Total VHC is heating duty divided by the sum of active compressor suction volume flows, in MJ/m³. "
        "The non-dominated set maximises both COP and VHC among these sampled "
        "best-found COP designs. It is not a multi-objective search of all cycle designs, "
        "and connecting lines do not represent calculated intermediate compositions.",
        "",
        "Conditions: source inlet 40°C; evaporator refrigerant inlet 10°C; sink inlet "
        "15°C; gas-cooler minimum approach 20 K; recuperator minimum approach 5 K; "
        "two compressors; motor efficiency 0.96. Optimiser v13, default search "
        "budgets, fixed composition, historical incumbent replay disabled.",
        "",
        "## Highest COP sampled in each family",
        "",
        "| Family | R290 mass % | COP | VHC (MJ/m³) | Maximum pressure (MPa) |",
        "|---|---:|---:|---:|---:|",
    ]
    for family in dict.fromkeys(row["family"] for row in rows):
        candidates = [row for row in rows if row["family"] == family and row["success"]]
        if candidates:
            best = max(candidates, key=lambda row: row["replay_COP"])
            report.append(f"| {family} | {best['R290_mass_percent']} | {best['replay_COP']:.4f} | "
                          f"{best['VHC_MJ_per_m3']:.4f} | {best['Pmax_MPa']:.3f} |")
    report.extend([
        "", "## Non-dominated sampled designs", "",
        "| Refrigerant | R290 mass % | COP | VHC (MJ/m³) |",
        "|---|---:|---:|---:|",
    ])
    for row in sorted((row for row in unique if row["non_dominated"]), key=lambda row: row["VHC_MJ_per_m3"]):
        label = row["family"] if ";" in row["input_refrigerant"] else row["input_refrigerant"]
        report.append(f"| {label} | {row['R290_mass_percent']} | "
                      f"{row['replay_COP']:.4f} | {row['VHC_MJ_per_m3']:.4f} |")
    drift = max((abs(row["COP_replay_change_percent"]) for row in successful), default=0)
    report.extend([
        "", f"Largest absolute COP change on fresh replay: {drift:.6g}%.", "",
        "Power saving is calculated at equal delivered heat duty as "
        "`100 × (1 − COP_R600 / COP_case)`. The separate SEC column retains "
        "the model definition `W / (Q_hot − W)`.", "",
        "[All compositions and design variables](composition_comparison.csv)",
        "",
    ])
    (output_dir / "analysis.md").write_text("\n".join(report))

    panels = (
        ("COP", "Heating COP"),
        ("Pmax_MPa", "Maximum pressure (MPa)"),
        ("power_reduction_vs_R600_percent", "Power saving vs R600 at equal heat duty (%)"),
        ("dT_rhx1_K", "RHX1 cold-side temperature lift (K)"),
        ("dT_rhx2_K", "RHX2 cold-side temperature lift (K)"),
        ("comp2_power_fraction", "Compressor 2 power fraction"),
    )
    with rc_context({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}):
        fig = Figure(figsize=(14, 8.5), layout="constrained")
        axes = fig.subplots(2, 3)
        for ax, (field, label) in zip(axes.flat, panels):
            for family, color, marker in zip(
                ("R290 / R600", "R290 / R601A", "R290 / R601", "R290 / R600A"),
                ("#176B87", "#B65320", "#6B4DA0", "#34804A"),
                ("o", "s", "^", "P"),
            ):
                selected = [row for row in rows if row["family"] == family]
                ax.plot(
                    [row["R290_mass_percent"] for row in selected],
                    [row[field] if row[field] is not None else math.nan for row in selected],
                    color=color, marker=marker, markersize=4, label=family,
                )
            if rows[0][field] is not None:
                ax.axhline(rows[0][field], color="#444444", linestyle="--", linewidth=1, label="R600 baseline")
            ax.set(xlabel="R290 mass percentage (%)", ylabel=label, xlim=(0, 100))
            ax.set_xticks(range(0, 101, 20))
            ax.grid(alpha=0.2)
        handles, labels = axes[0, 0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="outside lower center", ncol=5, frameon=False)
        fig.suptitle(
            f"Fixed-composition refrigerant optimisation at {target_temperature:g}°C\n"
            "0–100% R290 by mass; each point optimised independently",
            fontsize=14,
        )
        fig.savefig(output_dir / "composition_comparison.png", dpi=180)
        fig.savefig(output_dir / "composition_comparison.svg")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target", type=float, default=200.0, help="Sink target in degrees C")
    parser.add_argument("--n-cpu", type=int, default=min(6, max((os.cpu_count() or 2) - 4, 1)))
    parser.add_argument("--results-dir", type=Path, default=Path(DEFAULT_RESULTS_DIR))
    parser.add_argument("--dry-run", action="store_true", help="List cases without solving")
    parser.add_argument("--report-only", action="store_true", help="Regenerate analysis from completed results")
    args = parser.parse_args()
    if args.dry_run:
        for point in composition_points():
            print(f"{point.family:14} {point.propane_mass_percent:3}% R290  {point.refrigerant}")
        print("41 unique optimisations; compositions fixed on a mass basis.")
        return
    if not args.report_only:
        run_sweep(target_temperature=args.target, n_cpu=args.n_cpu, results_dir=args.results_dir)
    print(f"Analysis written to {write_analysis(args.results_dir, args.target)}")


if __name__ == "__main__":
    main()
