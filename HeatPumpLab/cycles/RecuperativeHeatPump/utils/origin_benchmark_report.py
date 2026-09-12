"""Aggregate repeated origin-start trials, with replay and budget validation."""

from collections import defaultdict
import json
from pathlib import Path

from matplotlib.figure import Figure
import numpy as np

from .optimiser_benchmark_report import case_label


LABELS = {"gp_ei": "BO (GP-EI)", "rbf": "RBF", "cma_es": "CMA-ES"}
COLORS = {"gp_ei": "#D55E00", "rbf": "#0072B2", "cma_es": "#009E73"}


def validate_trial(trial):
    meta, rows = trial["metadata"], trial["records"]
    if len(rows) != meta["budget"]:
        raise ValueError("Trial has not completed its budget")
    if not np.array_equal(rows[0]["x"][:3], meta["origin"]):
        raise ValueError("First evaluation is not the declared origin")
    fractions = rows[0]["x"][3:]
    spec = meta["refrigerant"].split("|")
    if len(spec) == 3 and spec[2] == "mass":
        requested = [float(value) for value in spec[1].split(";")]
        if not np.allclose(fractions, requested, atol=1e-12, rtol=0):
            raise ValueError("Evaluated composition differs from requested mass fractions")
    if len({tuple(r["x"]) for r in rows}) != len(rows):
        raise ValueError("Duplicate candidate evaluations")
    for row in rows:
        x = np.asarray(row["x"])
        if not np.allclose(x[3:], fractions, atol=1e-12, rtol=0):
            raise ValueError("Composition varied")
        if np.any(x[:3] < [0, 0, 0.1]) or np.any(x[:3] > [100, 100, 0.9]):
            raise ValueError("Candidate outside individual bounds")
        limit = (meta["base_case"]["T_air_tar"] - meta["base_case"]["T_air_inlet"]) * 2/3
        limit -= meta["base_case"]["dT_rhx3_cold"]
        if sum(x[:2]) > limit + 1e-8:
            raise ValueError("Candidate outside combined lift limit")
    replay = trial.get("replay", {})
    if not replay.get("success"):
        raise ValueError("No successful winner replay")
    best = min((r for r in rows if r["success"]), key=lambda r: r["SEC"])
    if not np.array_equal(best["x"], replay["x"]):
        raise ValueError("Replay used a different design")
    if abs(replay["COP"] / best["COP"] - 1) > 1e-8:
        raise ValueError("Winning COP failed the fresh-process replay")
    if replay["VHC_MJ_per_m3"] <= 0:
        raise ValueError("Missing positive VHC")


def write_origin_report(output, expected_tasks):
    output = Path(output)
    expected = {(ref, method, seed, budget) for ref, method, seed, budget, *_ in expected_tasks}
    found = {}
    for path in sorted(output.glob("*.json")):
        trial = json.loads(path.read_text())
        if "metadata" not in trial:
            continue
        m = trial["metadata"]
        key = (m["refrigerant"], m["method"], m["seed"], m["budget"])
        if key in expected:
            if key in found:
                raise ValueError(f"Duplicate trial: {key}")
            try:
                validate_trial(trial)
            except ValueError as exc:
                raise ValueError(f"{path.name}: {exc}") from exc
            found[key] = trial
    missing = expected - found.keys()
    if missing:
        raise ValueError(f"Missing {len(missing)} expected trials")
    trials = list(found.values())
    if len({t["metadata"]["search_source_sha256"] for t in trials}) != 1:
        raise ValueError("Trials used different search-code versions")
    if len({json.dumps(t["metadata"]["versions"], sort_keys=True) for t in trials}) != 1:
        raise ValueError("Trials used different dependency versions")
    groups = defaultdict(list)
    for trial in trials:
        groups[trial["metadata"]["refrigerant"]].append(trial)
    for runs in groups.values():
        if len({json.dumps(r["metadata"]["base_case"], sort_keys=True) for r in runs}) != 1:
            raise ValueError("Model settings differ within a composition")
        origins = [r["records"][0]["COP"] for r in runs]
        if not all(value is not None for value in origins) or not np.allclose(origins, origins[0], atol=1e-10, rtol=0):
            raise ValueError("Origin results differ across repeated trials")

    fig = Figure(figsize=(5.2 * len(groups), 4.8), layout="constrained")
    axes = np.atleast_1d(fig.subplots(1, len(groups)))
    summary = []
    for ax, (ref, runs) in zip(axes, sorted(groups.items())):
        reference = max(r["replay"]["COP"] for r in runs)
        for method in LABELS:
            selected = [r for r in runs if r["metadata"]["method"] == method]
            if not selected:
                continue
            traces = []
            for run in selected:
                values = [row["COP"] if row["success"] else np.nan for row in run["records"]]
                traces.append(np.fmax.accumulate(values))
            traces = np.asarray(traces)
            x = np.arange(1, traces.shape[1] + 1)
            ax.fill_between(x, traces.min(axis=0), traces.max(axis=0), color=COLORS[method], alpha=0.13)
            ax.plot(x, np.median(traces, axis=0), color=COLORS[method], label=LABELS[method], linewidth=1.8)
            cops = np.asarray([r["replay"]["COP"] for r in selected])
            wall = [sum(row["evaluation_wall_seconds"] + row["proposal_seconds"] for row in r["records"])
                    for r in selected]
            summary.append(dict(refrigerant=ref, method=method, runs=len(selected),
                                median_COP=float(np.median(cops)), min_COP=float(min(cops)), max_COP=float(max(cops)),
                                median_COP_at_64=float(np.median(traces[:, min(63, traces.shape[1]-1)])),
                                origin_COP=selected[0]["records"][0]["COP"],
                                within_0_1_percent_best_observed=int(sum(cops >= reference * 0.999)),
                                beats_v11=int(sum(cops > selected[0]["metadata"]["baseline_COP"])),
                                best_observed_COP=reference, median_wall_seconds=float(np.median(wall)),
                                failed_evaluations=sum(sum(not row["success"] for row in r["records"]) for r in selected),
                                native_restarts=sum(len(r.get("restarts", [])) for r in selected)))
        ax.axhline(runs[0]["metadata"]["baseline_COP"], linestyle="--", color="#555555", linewidth=1, label="Previous v11")
        ax.set(title=case_label(ref) + " by mass", xlabel="Candidate evaluations (origin included)", ylabel="Best nominal COP found")
        ax.set_xlim(1, max(r["metadata"]["budget"] for r in runs))
        ax.grid(alpha=0.2)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("Nominal optimisation from (0 K, 0 K, 50/50 split) at 200°C\n"
                 "Lines: median across seeds · shading: minimum–maximum", fontsize=14, weight="bold")
    fig.savefig(output / "origin_convergence.png", dpi=180)
    fig.savefig(output / "origin_convergence.svg")

    lines = ["# Origin-start optimisation benchmark", "",
             "All methods start from RHX1 = 0 K, RHX2 = 0 K and compressor 2 power fraction = 0.5. "
             "Composition is fixed by mass, with a 200°C sink target and the original model settings. "
             "Previous optima are comparators only. The complete declared budget includes the origin, "
             "initialisation, unsuccessful model evaluations and any algorithm restarts. Winner replays "
             "are additional validation evaluations outside that budget.", "",
             "These comparisons use the fixed default inner-solver initialisation. See the "
             "separate sensitivity checks below before interpreting the highest COP as a "
             "robust thermodynamic optimum.", "",
             "![Convergence](origin_convergence.png)", "",
             "## Consistency across seeds", "",
             "| Composition by mass | Method | Seeds | Median COP at 64 calls | Final median COP | Final min–max COP | Within 0.1% of best observed | Median evaluation + proposal time (s) | Failed evaluations |",
             "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for row in summary:
        lines.append(f"| {case_label(row['refrigerant'])} | {LABELS[row['method']]} | {row['runs']} "
                     f"| {row['median_COP_at_64']:.6f} | {row['median_COP']:.6f} | {row['min_COP']:.6f}–{row['max_COP']:.6f} "
                     f"| {row['within_0_1_percent_best_observed']}/{row['runs']} "
                     f"| {row['median_wall_seconds']:.1f} | {row['failed_evaluations']} |")
    lines += ["", "“Best observed” is the highest replayed COP among these methods and seeds for that "
              "composition, not a known global optimum. Durations sum proposal computation and model "
              "evaluation (including worker startup and communication), excluding optimiser updates, "
              "checkpoint writes, queueing and final replay. Concurrent load "
              "can affect timing. Native timeout retries are included in candidate durations; they do "
              "not create a new candidate-budget entry.", "",
              "## Every repeated trial", "",
              "| Composition | Method | Seed | Evaluations | COP | VHC (MJ/m³) | Pmax (MPa) | RHX1 (K) | RHX2 (K) | Power fraction | Restarts |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for trial in sorted(trials, key=lambda r: (r["metadata"]["refrigerant"], r["metadata"]["method"], r["metadata"]["seed"])):
        m, r = trial["metadata"], trial["replay"]
        lines.append(f"| {case_label(m['refrigerant'])} | {LABELS[m['method']]} | {m['seed']} | {len(trial['records'])} "
                     f"| {r['COP']:.6f} | {r['VHC_MJ_per_m3']:.4f} | {r['Pmax_MPa']:.4f} "
                     f"| {r['x'][0]:.5f} | {r['x'][1]:.5f} | {r['x'][2]:.6f} | {len(trial.get('restarts', []))} |")
    lines += ["", "## Algorithms", "",
              "BO and RBF reuse the prior GP expected-improvement and cubic-RBF implementations. "
              "Their global 24-point prescreen is replaced by five local observations: the origin, "
              "5 K changes to each recuperator separately, and compressor fractions 0.54 and 0.46. "
              "These points span the three input dimensions; subsequent proposals may explore the full "
              "feasible domain. No saved optima or historical evaluations train either surrogate.", "",
              "CMA-ES uses the official `cma` package, with the same origin as its initial distribution "
              "centre, a normalised standard deviation of 0.25, and population size seven. Box bounds "
              "use pycma's boundary transform. Candidates exceeding the known combined-lift limit "
              "are resampled before the expensive model call. Native stopping triggers a fresh "
              "origin-centred restart within the same total budget. The final incomplete population "
              "is evaluated and eligible to win, but does not update the distribution. See the "
              "[official CMA-ES implementation](https://cma-es.github.io/apidocs-pycma/cma.evolution_strategy.CMAEvolutionStrategy.html).", "",
              "Starting point and expensive-call budgets are matched; the methods' native "
              "initialisation distributions are different. Hyperparameters are fixed before the "
              "live trials. This compares these particular configurations, not every possible "
              "implementation or tuning of BO, RBF and CMA-ES.", "",
              "REFPROP is isolated with a 180-second timeout, one retry and worker recycling "
              "after 20 candidates. Each winner is replayed in a fresh process. Failed SEC "
              "observations are excluded from the surrogate response models, while CMA-ES "
              "ranks failures using a finite penalty. The nested inner optimisation is unchanged.", "",
              "These repeated seeds and compositions provide a repeatability check, not a general "
              "robustness guarantee over all refrigerants. COP is the optimisation target; VHC "
              "and pressure describe each selected design. This is not a multiobjective Pareto search.", ""]
    lines.extend(_diagnostics(output))
    (output / "origin_report.md").write_text("\n".join(lines))
    validation = dict(trials=len(trials), candidate_evaluations=sum(len(t["records"]) for t in trials),
                      fresh_replays=len(trials), budgets_complete=True, origins_verified=True,
                      compositions_fixed=True, duplicate_candidates=0, winner_replays_passed=True,
                      failed_evaluations=sum(r["failed_evaluations"] for r in summary))
    (output / "validation.json").write_text(json.dumps(validation, indent=2))
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    return output / "origin_report.md"


def _diagnostics(output):
    pressure_path = output / "pressure_seed_diagnostic.json"
    boundary_path = output / "boundary_diagnostic.json"
    lines = ["## Inner-solver sensitivity checks", ""]
    if not pressure_path.exists() or not boundary_path.exists():
        return []
    pressure_rows = json.loads(pressure_path.read_text())["records"]
    boundary = json.loads(boundary_path.read_text())
    lines += [
        "These diagnostic calls are additional to the optimiser budgets and were never "
        "fed back to the running searches. They test R290/R601 30/70 at 200°C.", "",
        "The same final design from seed 20260909 was evaluated with different numerical "
        "pressure guesses (`P_H_est`, MPa); zero selects the automatic initial estimate. "
        "All physical design variables remain fixed within each design group.", "",
        "| Design selected by | Initial pressure guess (MPa) | COP | Solved Pmax (MPa) |",
        "|---|---:|---:|---:|",
    ]
    for row in pressure_rows:
        if row["success"]:
            guess = "automatic" if row["pressure_seed_MPa"] == 0 else f"{row['pressure_seed_MPa']:g}"
            lines.append(f"| {LABELS[row['design']]} | {guess} | {row['COP']:.6f} | {row['Pmax_MPa']:.6f} |")
    lines += ["", "The CMA-ES design's reported COP depends strongly on the numerical pressure "
              "guess; matching a fresh replay with the default guess does not establish "
              "initialisation independence. The BO design tested here is much less sensitive. "
              "This does not establish that the higher-COP branch is physically invalid. "
              "It shows that the inner solve does not consistently recover the same result "
              "for an identical design, which confounds outer-optimiser comparisons.", "",
              "![Pressure-seed sensitivity](pressure_seed_sensitivity.png)", "",
              "A second diagnostic held RHX2 and compressor split at an earlier CMA-ES "
              f"candidate captured after {boundary['captured_after_evaluations']} evaluations, "
              "and perturbed RHX1 near zero:", "",
              "| RHX1 lift (K) | Fast inner objective | COP | Pmax (MPa) |",
              "|---:|---|---:|---:|"]
    for row in boundary["records"]:
        if row["success"]:
            lines.append(f"| {row['x'][0]:.9g} | {row['use_fast_objective']} | {row['COP']:.6f} | {row['Pmax_MPa']:.6f} |")
    lines += ["", "This is a narrow high-COP response region, not simply a discontinuity at "
              "exactly zero: very small positive lifts also recover the lower-COP solution. "
              "Exact-inner checks reproduce both the captured high result and the zero-lift "
              "result, so disabling the fast gas-cooler offset objective alone does not resolve it.", "",
              "The current GP length-scale lower bound is 0.02 in normalised coordinates "
              "(2 K for a lift), and the smallest local candidate standard deviation in both "
              "surrogate implementations is also 2 K. Those settings are poorly matched to "
              "features at hundredths or thousandths of a kelvin. RBF can still discover such "
              "regions by chance; CMA-ES adapts its sampling scale. This is a limitation of "
              "these tested configurations, not a general ranking of the algorithm families.", "",
              "A robust production solution should first make pressure/evaporator branch "
              "selection consistent and validate candidate optima against numerical starting "
              "guesses and small design perturbations. Then repeat the optimiser comparison "
              "on that controlled objective. The nominal results alone are insufficient "
              "grounds to promote the highest-scoring method as a robust design optimiser.", ""]
    fig = Figure(figsize=(8, 4.2), layout="constrained")
    ax = fig.subplots()
    guesses = sorted({r["pressure_seed_MPa"] for r in pressure_rows})
    for method in ("cma_es", "gp_ei"):
        rows = [r for r in pressure_rows if r["design"] == method and r["success"]]
        ax.plot([guesses.index(r["pressure_seed_MPa"]) for r in rows], [r["COP"] for r in rows],
                "o-", color=COLORS[method], label=f"Design selected by {LABELS[method]}")
    ax.set_xticks(range(len(guesses)), ["Automatic" if g == 0 else f"{g:g}" for g in guesses])
    ax.set(xlabel="Numerical pressure starting guess (MPa)", ylabel="COP",
           title="Same physical design, different inner-solver initialisation\nR290/R601 30/70 by mass at 200°C")
    ax.legend()
    ax.grid(alpha=0.2)
    fig.savefig(output / "pressure_seed_sensitivity.png", dpi=180)
    fig.savefig(output / "pressure_seed_sensitivity.svg")
    return lines
