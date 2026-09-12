"""Summarise completed optimiser pilots without changing the composition sweep."""

from collections import defaultdict
import json
from pathlib import Path

from matplotlib.figure import Figure
import numpy as np


LABELS = {"rbf": "RBF", "gp_ei": "GP expected improvement", "space_filling": "Space filling"}
COLORS = {"rbf": "#0072B2", "gp_ei": "#D55E00", "space_filling": "#8A8A8A"}


def case_label(refrigerant):
    fluids, fractions, _ = refrigerant.split("|")
    other = {"Pentane": "R601", "Isobutane": "R600A", "Isopentane": "R601A", "Butane": "R600"}
    percentages = "/".join(f"{100 * float(f):g}" for f in fractions.split(";"))
    return f"R290/{other.get(fluids.split(';')[1], fluids.split(';')[1])} {percentages}"


def write_benchmark_report(output):
    output = Path(output)
    trials = []
    for path in sorted(output.glob("*.json")):
        trial = json.loads(path.read_text())
        if "metadata" in trial and trial.get("replay"):
            trials.append(trial)
    cold = [t for t in trials if not t["metadata"].get("warm_start")]
    if not cold:
        raise ValueError("No completed cold-start trials")
    groups = defaultdict(list)
    for trial in cold:
        groups[trial["metadata"]["refrigerant"]].append(trial)
    fig = Figure(figsize=(5 * len(groups), 4.7), layout="constrained")
    axes = np.atleast_1d(fig.subplots(1, len(groups)))
    for ax, (ref, runs) in zip(axes, sorted(groups.items())):
        visible_gains = []
        for run in runs:
            meta = run["metadata"]
            values = [r["COP"] if r["success"] else -np.inf for r in run["records"]]
            best = np.maximum.accumulate(values)
            best[~np.isfinite(best)] = np.nan
            gain = 100 * (best / meta["baseline_COP"] - 1)
            visible_gains.extend(gain[23:][np.isfinite(gain[23:])])
            ax.step(np.arange(1, len(best) + 1), gain, where="post",
                    color=COLORS[meta["method"]], label=f"{LABELS[meta['method']]} (seed {meta['seed']})")
        ax.axhline(0, color="#333333", linestyle="--", linewidth=1)
        ax.axvline(24, color="#AAAAAA", linestyle=":", linewidth=1)
        ax.set_title(case_label(ref) + " by mass", fontsize=11, weight="bold")
        ax.set_xlabel("Exact evaluations (initial 24 included)")
        ax.set_ylabel("Best COP change from v11 (%)")
        ax.set_xlim(24, max(r["metadata"]["budget"] for r in runs))
        if visible_gains:
            ax.set_ylim(min(min(visible_gains), 0) - 0.15, max(max(visible_gains), 0) + 0.15)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle("Fixed composition at 200°C: joint optimisation pilot\n"
                 "Adaptive phase after the common 24-point initial design", fontsize=14, weight="bold")
    fig.savefig(output / "optimiser_convergence.png", dpi=180)
    fig.savefig(output / "optimiser_convergence.svg")

    lines = [
        "# Fixed-composition optimiser pilot at 200°C", "",
        "The cold-start trials use the same 24 initial designs, seed, bounds and model "
        "settings. Each searches both recuperator lifts and compressor power split jointly; "
        "mass composition is held constant. Historical optima are comparison values only, "
        "never search seeds. An independent fresh-process replay follows each trial and "
        "is outside its evaluation budget.", "",
        "These are small pilot implementations using existing NumPy/SciPy dependencies, "
        "not BoTorch or a complete implementation of the published MSRBF algorithm. "
        "One random seed per method/case does not establish a statistically robust ranking.", "",
        "![Convergence](optimiser_convergence.png)", "",
        "## Cold-start results", "",
        "| Composition by mass | Method | Calls | v11 calls | v11 COP | Replayed COP | COP change | VHC (MJ/m³) | Reported solve / proposal seconds | Failed calls |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for trial in cold:
        m, r = trial["metadata"], trial["replay"]
        valid = r["success"]
        cop = r["COP"] if valid else float("nan")
        gain = 100 * (cop / m["baseline_COP"] - 1)
        solve = sum(item["elapsed_seconds"] for item in trial["records"])
        proposal = sum(item.get("proposal_seconds", 0) for item in trial["records"])
        failed = sum(not item["success"] for item in trial["records"])
        lines.append(f"| {case_label(m['refrigerant'])} | {LABELS[m['method']]} | {len(trial['records'])} "
                     f"| {m['baseline_evaluations']} | {m['baseline_COP']:.6f} | {cop:.6f} "
                     f"| {gain:+.3f}% | {r.get('VHC_MJ_per_m3', float('nan')):.4f} "
                     f"| {solve:.1f} / {proposal:.1f} | {failed} |")
    lines += ["", "Times sum the recorded inner solve and proposal durations under concurrent load; "
              "they exclude worker startup/communication and are not end-to-end wall times. "
              "A timed-out candidate can include "
              "a retry, so candidate counts are not always native solver invocation counts.", "",
              "## Warm-start diagnostic: reopen both mixture recuperator lifts", "",
              "This uses the existing pattern search from the saved v11 optimum, with "
              "composition fixed. Its calls are **additional** to the original v11 run and "
              "are not directly comparable with a cold start.", "",
              "| Composition by mass | Extra calls | v11 COP | Replayed COP | COP change |",
              "|---|---:|---:|---:|---:|"]
    for trial in trials:
        m, r = trial["metadata"], trial["replay"]
        if not m.get("warm_start"):
            continue
        cop = r["COP"] if r["success"] else float("nan")
        lines.append(f"| {case_label(m['refrigerant'])} | {len(trial['records'])} | {m['baseline_COP']:.6f} "
                     f"| {cop:.6f} | {100*(cop/m['baseline_COP']-1):+.3f}% |")
    lines += ["", "## Replayed designs", "",
              "| Composition by mass | Method | RHX1 lift (K) | RHX2 lift (K) | Compressor 2 power fraction | Pmax (MPa) | COP replay drift |",
              "|---|---|---:|---:|---:|---:|---:|"]
    for trial in trials:
        m, r = trial["metadata"], trial["replay"]
        if not r["success"]:
            continue
        original = min((item for item in trial["records"] if item["success"]), key=lambda item: item["SEC"])
        x = r["x"]
        lines.append(f"| {case_label(m['refrigerant'])} | {LABELS.get(m['method'], 'Warm joint polish')} "
                     f"| {x[0]:.6f} | {x[1]:.6f} | {x[2]:.6f} | {r['Pmax_MPa']:.4f} "
                     f"| {100*(r['COP']/original['COP']-1):+.6f}% |")
    lines += ["", "## Method details and limits", "",
              "- RBF: regularised cubic response surface; candidate merit trades predicted SEC "
              "against distance from sampled points, with weights cycling 0.3/0.5/0.8/0.95. "
              "This follows the value/distance selection idea of "
              "[Regis and Shoemaker](https://pubsonline.informs.org/doi/10.1287/ijoc.1060.0182).",
              "- GP-EI: standardised SEC, Matérn 5/2 covariance with three fitted length scales, "
              "fixed 1e-6 diagonal regularisation, and expected improvement. Hyperparameters "
              "are refitted periodically. Acquisition selection is over a finite candidate pool.",
              "- Both use global, boundary and local candidates; every eighth adaptive call "
              "is global maximin exploration. The space-filling control selects maximin on "
              "the same candidate-generation scheme, including incumbent-centred candidates.",
              "- Failed solves are excluded from the SEC model and downweighted by a separate "
              "distance heuristic. This is not a calibrated probability of feasibility. For "
              "production constrained BO, use an explicit outcome model as described in "
              "[BoTorch constraints](https://botorch.org/docs/constraints).",
              "- Bounds are 0–100 K for each lift, 0.1–0.9 for compressor split, and "
              "RHX1 + RHX2 ≤ 123⅓ K. Candidates obey the lift constraint before evaluation.",
              "- REFPROP isolation, 180-second timeout, one timeout/worker retry, and worker "
              "recycling after 20 calls match v11. Winner replays use a fresh worker.",
              "- The inner gas-cooler offset optimisation and cycle solver are unchanged. "
              "Reproducibility does not establish that the inner or outer solution is globally optimal.",
              "- This minimises SEC (equivalent to maximising COP here). VHC is measured at "
              "the selected design. A complete COP–VHC front requires multiobjective optimisation; "
              "[BoTorch supports hypervolume-based acquisition functions](https://botorch.org/docs/multi_objective).", ""]
    path = output / "benchmark_report.md"
    path.write_text("\n".join(lines))
    landscape = output / "pentane30_landscape.json"
    if landscape.exists():
        _plot_landscape(json.loads(landscape.read_text()), output)
    return path


def _plot_landscape(data, output):
    rows = [r for r in data["records"] if r["fast_inner_objective"] and r["success"]]
    fig = Figure(figsize=(10, 4), layout="constrained")
    axes = fig.subplots(1, 2)
    for ax, metric, label, color in zip(
        axes, ("COP", "Pmax_MPa"), ("COP", "Maximum pressure (MPa)"), ("#0072B2", "#D55E00"),
    ):
        ax.plot([r["interpolation_fraction"] for r in rows], [r[metric] for r in rows], "o-", color=color)
        ax.set_xlabel("Design interpolation: 0 = v11, 1 = GP-EI")
        ax.set_ylabel(label)
        ax.grid(alpha=0.2)
    fig.suptitle("R290/R601 30/70: a tested path between two designs at 200°C", weight="bold")
    fig.savefig(output / "pentane30_landscape.png", dpi=180)
    fig.savefig(output / "pentane30_landscape.svg")
