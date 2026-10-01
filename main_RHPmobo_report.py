"""Export validated COP/total-VHC fronts and convergence figures."""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from HeatPumpLab.cycles.RecuperativeHeatPump.utils.composition_pareto import nondominated_mask
from main_RHPmobo import PARTNERS, OUTPUT, frontier_indices

LABELS = {"Butane": "R290 / R600", "Isopentane": "R290 / R601A",
          "Pentane": "R290 / R601", "Isobutane": "R290 / R600A"}
COLORS = ("#187a70", "#c86b25", "#5366b3", "#b14972")


def hypervolume(records):
    # Integrate the two-dimensional union in O(n log n), including duplicate
    # points only once. This keeps full 4,000-evaluation traces inexpensive.
    front = sorted((r for r in records if r["success"] and not r.get("validation_rejected")
                    and np.isfinite([r["COP"], r["VHC_MJ_per_m3"]]).all()
                    and min(r["COP"], r["VHC_MJ_per_m3"]) > 0), key=lambda r: -r["COP"])
    previous = 0.
    value = 0.
    for r in front:
        if r["VHC_MJ_per_m3"] > previous:
            value += (r["VHC_MJ_per_m3"]-previous)*r["COP"]
            previous = r["VHC_MJ_per_m3"]
    return value


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_report(output):
    output = Path(output)
    root = Path(__file__).resolve().parent
    manifest = json.loads((output/"source_manifest.json").read_text())
    for name, digest in manifest["files"].items():
        assert hashlib.sha256((root/name).read_bytes()).hexdigest() == digest, name
    datasets = {partner: json.loads((output/f"Propane_{partner}.json").read_text()) for partner in PARTNERS}
    budgets = {data["metadata"]["budget"] for data in datasets.values()}
    assert len(budgets) == 1
    budget = budgets.pop()
    assert budget >= 400
    rows, summaries, traces = [], [], {}
    for partner, data in datasets.items():
        assert data["complete"] and len(data["records"]) == data["metadata"]["budget"] == budget
        assert data["frontier_indices"] == frontier_indices(data["records"])
        selected = []
        for index in data["frontier_indices"]:
            original = data["records"][index]
            replay = data["frontier_replays"][str(index)]
            assert replay["success"] and np.allclose(replay["x"], original["x"], atol=1e-10, rtol=0)
            assert all(np.isclose(replay[key], original[key], rtol=1e-6, atol=1e-9)
                       for key in ("COP", "SEC", "VHC_MJ_per_m3", "Pmax_MPa"))
            selected.append(dict(family=LABELS[partner], partner=partner, evaluation=index+1,
                                 COP=replay["COP"], VHC_total_MJ_per_m3=replay["VHC_MJ_per_m3"],
                                 propane_mass_percent=100*replay["x"][3], RHX1_K=replay["x"][0],
                                 RHX2_K=replay["x"][1], compressor_2_power_fraction=replay["x"][2],
                                 discharge_pressure_MPa=replay["Pmax_MPa"]))
        keep = nondominated_mask([(r["COP"], r["VHC_total_MJ_per_m3"]) for r in selected])
        selected = [r for r, valid in zip(selected, keep) if valid]
        rows.extend(selected)
        history = [hypervolume(data["records"][:n]) for n in range(1, budget+1)]
        traces[partner] = history
        summaries.append(dict(partner=partner, family=LABELS[partner], evaluations=budget,
                              failed=sum(not r["success"] for r in data["records"]),
                              timeouts=sum("evaluation_timeout" in r["constraint_errors"] for r in data["records"]),
                              validated_front_points=len(selected), replay_evaluations=len(data["frontier_replays"]),
                              validation_rejections=sum(bool(r.get("validation_rejected")) for r in data["records"]),
                              highest_COP=max(selected, key=lambda r: r["COP"]),
                              highest_VHC=max(selected, key=lambda r: r["VHC_total_MJ_per_m3"]),
                              hypervolume=history[-1],
                              hypervolume_at_400=history[399],
                              hypervolume_gain_since_400_percent=100*(history[-1]/history[399]-1),
                              highest_COP_at_400=max(r["COP"] for r in data["records"][:400]
                                                    if r["success"] and not r.get("validation_rejected")),
                              last_100_hypervolume_gain_percent=100*(history[-1]/history[-101]-1),
                              evaluation_wall_seconds=sum(r["evaluation_wall_seconds"] for r in data["records"]),
                              proposal_seconds=sum(r["proposal_seconds"] for r in data["records"])))
    mask = nondominated_mask([(r["COP"], r["VHC_total_MJ_per_m3"]) for r in rows])
    for row, keep in zip(rows, mask):
        row["global_front"] = keep
    for summary in summaries:
        summary["combined_front_points"] = sum(r["global_front"] and r["partner"] == summary["partner"] for r in rows)
    rows.sort(key=lambda r: (r["family"], r["VHC_total_MJ_per_m3"]))
    write_csv(output/"pareto_fronts.csv", rows)
    write_csv(output/"global_pareto_front.csv", [r for r in rows if r["global_front"]])
    (output/"summary.json").write_text(json.dumps(dict(families=summaries, validated_fronts=rows), indent=2)+"\n")

    fig, ax = plt.subplots(figsize=(9, 6), constrained_layout=True)
    for partner, color in zip(PARTNERS, COLORS):
        samples = [r for r in datasets[partner]["records"] if r["success"]]
        ax.scatter([r["VHC_MJ_per_m3"] for r in samples], [r["COP"] for r in samples],
                   color=color, alpha=.10, s=9, linewidths=0)
        front = [r for r in rows if r["partner"] == partner]
        ax.plot([r["VHC_total_MJ_per_m3"] for r in front], [r["COP"] for r in front],
                "o-", color=color, markersize=4, linewidth=1.5, label=LABELS[partner])
    global_rows = [r for r in rows if r["global_front"]]
    ax.scatter([r["VHC_total_MJ_per_m3"] for r in global_rows], [r["COP"] for r in global_rows],
               facecolors="none", edgecolors="#222222", marker="D", s=50, linewidths=.8, label="Combined observed front")
    ax.set(xlabel="Total volumetric heating capacity (MJ/m³)", ylabel="Heating COP",
           title=f"COP–VHC_total Pareto estimates · 200°C sink\n{budget:,} BO evaluations per binary family")
    ax.grid(alpha=.18)
    ax.legend(fontsize=9)
    fig.savefig(output/"pareto_fronts.png", dpi=190)
    fig.savefig(output/"pareto_fronts.svg")
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for ax, partner in zip(axes.flat, PARTNERS):
        front = [r for r in rows if r["partner"] == partner]
        ax.plot([r["VHC_total_MJ_per_m3"] for r in front], [r["COP"] for r in front], color="#bbbbbb", linewidth=1)
        plotted = ax.scatter([r["VHC_total_MJ_per_m3"] for r in front], [r["COP"] for r in front],
                            c=[r["propane_mass_percent"] for r in front], cmap="viridis", vmin=0, vmax=100, s=28)
        ax.set(title=LABELS[partner], xlabel="VHC_total (MJ/m³)", ylabel="Heating COP")
        ax.grid(alpha=.18)
    fig.colorbar(plotted, ax=axes, label="Propane mass fraction (%)", shrink=.85)
    fig.suptitle("Validated sampled trade-offs and mixture composition")
    fig.savefig(output/"pareto_compositions.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 4), constrained_layout=True)
    for partner, color in zip(PARTNERS, COLORS):
        ax.plot(range(1, budget+1), traces[partner], color=color, label=LABELS[partner])
    ax.axvline(32, color="#777777", linestyle=":", linewidth=1)
    if budget > 400:
        ax.axvline(400, color="#777777", linestyle="--", linewidth=1, label="Previous budget")
    ax.set(xlabel="Candidate evaluations (failures included)", ylabel="Dominated hypervolume (reference: 0, 0)",
           title="Multi-objective search progress")
    ax.grid(alpha=.18)
    ax.legend(fontsize=9)
    fig.savefig(output/"hypervolume_progress.png", dpi=180)
    plt.close(fig)

    mixed_partners = [p for p in PARTNERS if any(r["global_front"] and r["partner"] == p
                     and 1e-8 < r["propane_mass_percent"] < 100-1e-8 for r in rows)]
    pure_only = [p for p in PARTNERS if p not in mixed_partners
                 and any(r["global_front"] and r["partner"] == p for r in rows)]
    pooled_finding = "Mixed-fluid designs on the pooled observed front come from " + ", ".join(
        f"propane/{p.lower()}" for p in mixed_partners) + "."
    if pure_only:
        pooled_finding += " The " + ", ".join(p.lower() for p in pure_only) + (
            " runs contribute only pure-fluid designs to that combined front; their mixed-fluid designs "
            "are dominated by sampled alternatives.")
    text = ["# COP–total-VHC multi-objective Bayesian optimisation", "",
            f"Four binary families were optimised with {budget:,} recorded candidate evaluations each ({4*budget:,} total), "
            "using the v17 internal PH fallback and the summed compressor-suction definition of total VHC. "
            "All reported frontier designs were checked in fresh native processes.", "",
            pooled_finding, "",
            "![Validated Pareto fronts](pareto_fronts.png)", "",
            "| Binary family | Family / combined front points | Failed calls | Highest COP (VHC; propane mass %) | Highest VHC (COP; propane mass %) |",
            "|---|---:|---:|---|---|"]
    for s in summaries:
        c, v = s["highest_COP"], s["highest_VHC"]
        text.append(f"| {s['family']} | {s['validated_front_points']} / {s['combined_front_points']} | {s['failed']} | "
                    f"{c['COP']:.6f} ({c['VHC_total_MJ_per_m3']:.4f} MJ/m³; {c['propane_mass_percent']:.2f}%) | "
                    f"{v['VHC_total_MJ_per_m3']:.4f} MJ/m³ ({v['COP']:.6f}; {v['propane_mass_percent']:.2f}%) |")
    text += ["", "## Domain and method", "",
             "Both COP and VHC_total are maximised. The four BO variables are RHX1 and RHX2 cold-side lifts "
             "(0–100 K each, with a combined limit of 123.333333 K), compressor-2 power fraction (0.4–0.6), "
             "and propane mass fraction (0–1). The sink target is 200°C, the source inlet 40°C, and the "
             "gas-cooler minimum approach target 20 K. All other settings are the existing DEFAULT_BASE_CASE. "
             "Discharge pressure and gas-cooler offset retain the existing inner cycle-solver policy; "
             "this is a Pareto estimate within that policy, not a search over arbitrary independent pressure settings.", "",
             "BoTorch qLogEHVI uses separate GP outputs with standardisation, a fixed (0, 0) reference point, "
             "64 quasi-Monte-Carlo samples and sequential q=1 proposals. The initial 32 evaluations combine "
             "composition endpoints, zero-lift anchors and a Sobol design. Hyperparameters are fitted every "
             "five evaluations; each eighth adaptive proposal uses maximin exploration. A distance-based "
             "feasibility heuristic downweights neighbourhoods of failed evaluations; failures are excluded "
             "from the thermodynamic GP outputs. This heuristic is not a calibrated probability of feasibility. "
             "Acquisition optimisation combines a global/boundary pool and continuous constrained restarts. "
             "See [BoTorch multi-objective documentation](https://botorch.org/docs/v0.17.1/multi_objective).", "",
             "## Validation and limits", "",
             f"Frontier verification used {sum(s['replay_evaluations'] for s in summaries)} additional evaluations outside "
             f"the search budgets. {sum(s['validation_rejections'] for s in summaries)} candidates were excluded after "
             "failing fresh replay or exceeding 10⁻⁶ relative / 10⁻⁹ absolute replay tolerance for COP, VHC, SEC or pressure. "
             "Any newly exposed frontier points were then replayed. The plotted/front-table values are from the "
             "successful fresh replays. Search checkpoints retain all observations, errors and rejected candidates.", "",
             "These are validated observed Pareto estimates from one seed per family, not proof of global optimality. "
             f"The search includes failed evaluations in each {budget:,}-call budget and uses a 180-second limit with no "
             "timeout retry per candidate. Lines connect sampled designs and do not validate intermediate points. "
             "Native property failures are not automatic evidence of physical infeasibility. "
             "Pure endpoints may appear in more than one family's run.", "",
             "The isopentane run resumed after an observation-validation check rejected a boundary point "
             "for a 1.42 × 10⁻¹⁴ K floating-point overshoot. The fix permits only encoding round-off "
             "(10⁻¹⁰ K at the joint limit); it does not move points or relax proposal screening. "
             "The 361 saved observations and GP parameters were preserved. The interrupted candidate was "
             "evaluated once again, adding one repeated oracle attempt outside the 1,600 saved search slots. "
             "Original checkpoints/source and migration hashes are archived. A separate 33-candidate pilot "
             "with four frontier replays tested the implementation and supplies no observations to these fronts.", "",
             "| Family | Final hypervolume | Gain during last 100 evaluations | Timeouts |",
             "|---|---:|---:|---:|"]
    for s in summaries:
        text.append(f"| {s['family']} | {s['hypervolume']:.6f} | {s['last_100_hypervolume_gain_percent']:.3f}% | {s['timeouts']} |")
    if budget > 400:
        text += ["", "## Improvement beyond 400 evaluations", "",
                 "Each run continues the original 400 observations and saved GP parameters with the same seed, "
                 "model and acquisition settings. Original results remain in the separate 400-evaluation directory. "
                 "Fresh frontier replays already completed there are retained; new frontier designs are replayed here. "
                 "The evaluation budget remains the stopping rule; small recent gains do not prove convergence.", "",
                 "| Family | COP at 400 | Final COP | COP gain | Hypervolume gain |",
                 "|---|---:|---:|---:|---:|"]
        for s in summaries:
            gain = 100*(s['highest_COP']['COP']/s['highest_COP_at_400']-1)
            text.append(f"| {s['family']} | {s['highest_COP_at_400']:.6f} | {s['highest_COP']['COP']:.6f} | "
                        f"{gain:.4f}% | {s['hypervolume_gain_since_400_percent']:.4f}% |")
    text += ["", "![Composition along each front](pareto_compositions.png)", "",
             "![Hypervolume progress](hypervolume_progress.png)", "", "## Data", "",
             "- [All family frontier designs](pareto_fronts.csv)",
             "- [Combined observed frontier](global_pareto_front.csv)",
             "- [Numerical summary](summary.json)",
             "- [Frozen model and search source](source_manifest.json)", "",
             f"Resume the saved checkpoints with `.venv/bin/python main_RHPmobo.py --budget {budget} --workers 4 "
             f"--output {output}`; rebuild this report with `.venv/bin/python main_RHPmobo_report.py --output {output}`."]
    (output/"REPORT.md").write_text("\n".join(text)+"\n")
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=OUTPUT)
    build_report(parser.parse_args().output)
