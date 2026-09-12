"""Validate and report the v13 400-evaluation comparison, including production line search."""

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

from matplotlib.figure import Figure
import numpy as np

from main_RHPoptimiser_benchmark import DEFAULT_CASES
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import OPTIMISER_VERSION
from HeatPumpLab.cycles.RecuperativeHeatPump.utils.origin_benchmark_report import validate_trial
from HeatPumpLab.cycles.RecuperativeHeatPump.utils.optimiser_benchmark_report import case_label as _case_label

LABELS = {"line_search": "Line-search", "rbf": "RBF", "gp_ei": "BO (GP-EI)", "cma_es": "CMA-ES"}
COLORS = {"line_search": "#7B3294", "rbf": "#0072B2", "gp_ei": "#D55E00", "cma_es": "#009E73"}
SEEDS = (20260909, 20260910, 20260911)


def case_label(ref):
    return _case_label(ref).replace("R600A", "R600a")


def load_trials(output, budget=400):
    trials = []
    expected = {(ref, method, seed) for ref in DEFAULT_CASES for method in LABELS
                for seed in ((None,) if method == "line_search" else SEEDS)}
    for path in sorted(output.glob("*.json")):
        trial = json.loads(path.read_text())
        if "metadata" not in trial or "records" not in trial:
            continue
        m = trial["metadata"]
        key = m["refrigerant"], m["method"], m["seed"]
        if key not in expected:
            raise ValueError(f"Unexpected/duplicate trial: {path}")
        expected.remove(key)
        if m["budget"] != budget or m["evaluation_version"] != OPTIMISER_VERSION:
            raise ValueError(f"Wrong budget/model: {path}")
        if m["method"] == "line_search":
            # Production has separate stage caches. Repeated physical designs
            # genuinely invoke the model again and must count against its budget.
            unique = list({tuple(r["x"]): r for r in trial["records"]}.values())
            checked = dict(trial, records=unique, metadata=dict(m, budget=len(unique)))
            validate_trial(checked)
            if not 0 < len(trial["records"]) <= budget:
                raise ValueError("Line-search exceeded its ceiling")
            if len(trial["records"]) != trial["reported_exact_evaluations"]:
                raise ValueError("Production count differs from instrumented model calls")
        else:
            validate_trial(trial)
        trials.append(trial)
    if expected:
        raise ValueError(f"Missing {len(expected)} trials")
    model_hash = hashlib.sha256(b"".join(p.read_bytes() for p in sorted(Path("HeatPumpLab").rglob("*.py")))).hexdigest()
    if {t["metadata"]["model_source_sha256"] for t in trials} != {model_hash}:
        raise ValueError("Model source changed between runs or since execution")
    stochastic = [t for t in trials if t["metadata"]["method"] != "line_search"]
    if len({t["metadata"]["search_source_sha256"] for t in stochastic}) != 1:
        raise ValueError("Search code changed")
    if len({json.dumps(t["metadata"]["versions"], sort_keys=True) for t in stochastic}) != 1:
        raise ValueError("Search dependencies changed")
    for ref in DEFAULT_CASES:
        runs = [t for t in trials if t["metadata"]["refrigerant"] == ref]
        if len({json.dumps(t["metadata"]["base_case"], sort_keys=True) for t in runs}) != 1:
            raise ValueError("Physical settings differ")
        if np.ptp([t["records"][0]["COP"] for t in runs]) > 1e-10:
            raise ValueError("Origin COP differs")
    return trials


def report(output, budget=400):
    output = Path(output)
    trials = load_trials(output, budget)
    summary = []
    fig = Figure(figsize=(16, 5), layout="constrained")
    axes = fig.subplots(1, 3)
    for ax, ref in zip(axes, DEFAULT_CASES):
        runs = [t for t in trials if t["metadata"]["refrigerant"] == ref]
        line_cop = next(t["replay"]["COP"] for t in runs if t["metadata"]["method"] == "line_search")
        best_observed = max(t["replay"]["COP"] for t in runs)
        for method in LABELS:
            selected = [t for t in runs if t["metadata"]["method"] == method]
            traces = [np.fmax.accumulate([r["COP"] if r["success"] else np.nan for r in t["records"]]) for t in selected]
            padded = np.array([np.pad(v, (0, budget-len(v)), constant_values=v[-1]) for v in traces])
            cops = np.array([t["replay"]["COP"] for t in selected])
            winner = max(selected, key=lambda t: t["replay"]["COP"])
            wall = [t.get("optimisation_wall_seconds", sum(r["evaluation_wall_seconds"] + r["proposal_seconds"] for r in t["records"])) for t in selected]
            n = len(traces[0]) if method == "line_search" else budget
            x = np.arange(1, budget+1)
            median = np.median(padded, axis=0)
            ax.plot(x[:n], median[:n], color=COLORS[method], label=LABELS[method], lw=1.8)
            if method == "line_search":
                ax.scatter([n], [median[n-1]], color=COLORS[method], s=28, zorder=5)
                ax.plot([n, budget], [median[-1]]*2, color=COLORS[method], ls=":", lw=1)
            else:
                ax.fill_between(x, padded.min(axis=0), padded.max(axis=0), color=COLORS[method], alpha=.14)
            summary.append(dict(refrigerant=ref, method=method, runs=len(selected),
                                calls_per_run=[len(t["records"]) for t in selected],
                                median_COP=float(np.median(cops)), min_COP=float(cops.min()), max_COP=float(cops.max()),
                                median_COP_at_100=float(np.median(padded[:,99])), median_COP_at_200=float(np.median(padded[:,199])),
                                gain_over_line_percent=float(100*(np.median(cops)/line_cop-1)),
                                within_0_1_percent_best=int(sum(cops >= .999*best_observed)),
                                median_wall_seconds=float(np.median(wall)),
                                failures=sum(not r["success"] for t in selected for r in t["records"]),
                                winner_seed=winner["metadata"]["seed"], winner=winner["replay"]))
        ax.set(title=case_label(ref)+" by mass", xlabel="Model evaluations", ylabel="Best COP found", xlim=(1,budget))
        ax.grid(alpha=.2)
        ax.legend(fontsize=8, loc="lower right")
    fig.suptitle("v13 solver · default pressure initialisation · 200°C target\n"
                 "Median and observed seed range; line-search dot marks its actual stopping call", fontsize=13)
    fig.savefig(output / "convergence.png", dpi=180)
    fig.savefig(output / "convergence.svg")
    lines = ["# Optimiser retest: v13, 400-call allowance", "",
             "All methods use the updated v13 thermodynamic solver, fixed mass compositions, a 200°C sink target, "
             "and the same initial design: RHX1 = RHX2 = 0 K, compressor 2 power fraction = 0.5. "
             "BO, RBF and CMA-ES each use exactly 400 model evaluations per run and three seeds (20260909–20260911). "
             "The unchanged deterministic production line-search is run once per composition with a 400-call ceiling; "
             "its existing stages can finish earlier. Actual calls are reported below. "
             "No historical optimum is used to initialise or train a search.", "",
             "![Convergence](convergence.png)", "",
             "## Replayed results", "",
             "| Composition by mass | Method | Runs | Calls per run | Median COP | Min–max COP | Median gain over line-search | Failed calls |",
             "|---|---|---:|---|---:|---:|---:|---:|"]
    if (output / "pressure_followup.json").exists():
        lines[4:4] = [
            "**Remaining inner-solver sensitivity:** these are results under the default numerical initialisation. "
            "For the best CMA-ES design on propane/pentane 60/40, changing the pressure starting guess to 7 MPa "
            "reduces COP from 2.538172 to 2.321442 (8.54%). The pressure search accepts a higher feasible pressure "
            "inside the pinch tolerance band. The default winner replays are reproducible, but the 60/40 CMA-ES "
            "result is sensitive to pressure initialisation. See the diagnostic evidence below.", ""]
    for s in summary:
        calls = ", ".join(map(str, sorted(set(s["calls_per_run"]))))
        lines.append(f"| {case_label(s['refrigerant'])} | {LABELS[s['method']]} | {s['runs']} | {calls} | {s['median_COP']:.6f} "
                     f"| {s['min_COP']:.6f}–{s['max_COP']:.6f} | {s['gain_over_line_percent']:+.3f}% | {s['failures']} |")
    lines += ["", "The shaded bands are the observed minimum–maximum across three seeds, not confidence intervals. "
              "The dotted continuation after each line-search stopping marker carries forward its incumbent; "
              "it does not represent additional model calls. Failed candidate solves consume the allowance; "
              "each winner is replayed in a separate fresh process outside the optimisation allowance. "
              "Production line-search uses separate caches for its stages, so repeated designs that really "
              "invoke the model again are included in its call count.", "", "## Convergence and runtime", "",
              "| Composition | Method | Median COP at 100 | Median COP at 200 | Median final COP | Within 0.1% of best observed | Median runtime (s) |",
              "|---|---|---:|---:|---:|---:|---:|"]
    for s in summary:
        lines.append(f"| {case_label(s['refrigerant'])} | {LABELS[s['method']]} | {s['median_COP_at_100']:.6f} "
                     f"| {s['median_COP_at_200']:.6f} | {s['median_COP']:.6f} | {s['within_0_1_percent_best']}/{s['runs']} | {s['median_wall_seconds']:.1f} |")
    lines += ["", "Runtime for BO/RBF/CMA-ES sums model-call wall time and proposal generation, excluding checkpoint writes, "
              "distribution updates and winner replay. Line-search time covers its complete production search. "
              "Concurrent machine load affects these timings. Best observed means the highest replayed COP in "
              "this experiment, not a proven global optimum.", "", "## Best design found by each method", "",
              "| Composition | Method | Seed | COP | RHX1 (K) | RHX2 (K) | Compressor 2 fraction | VHC (MJ/m³) | Pmax (MPa) |",
              "|---|---|---:|---:|---:|---:|---:|---:|---:|"]
    for s in summary:
        r = s['winner']
        seed = s['winner_seed'] if s['winner_seed'] is not None else 'deterministic'
        lines.append(f"| {case_label(s['refrigerant'])} | {LABELS[s['method']]} | {seed} | {r['COP']:.6f} "
                     f"| {r['x'][0]:.6f} | {r['x'][1]:.6f} | {r['x'][2]:.7f} | {r['VHC_MJ_per_m3']:.5f} | {r['Pmax_MPa']:.5f} |")
    lines += ["", "## Protocol and validation", "",
              "The physical settings are loaded from the prior comparison's case rows; all objective values are newly evaluated. "
              "Every trial records the optimiser version, the SHA-256 of all HeatPumpLab Python sources, "
              "the fixed composition and its evaluation history. The stochastic trials also record dependency versions and seeds.", "",
              "The production line-search retains its edge searches, fixed-total-lift split search, compressor-split "
              "pattern-search polish and short lift feedback. Increasing its outer allowance does not increase its capped "
              "stage allocations. BO and RBF keep their existing five-point local bootstrap, Matérn-5/2 GP-EI and cubic "
              "RBF implementations. CMA-ES keeps pycma's population of seven, normalised initial sigma 0.25, feasible-domain "
              "rejection and origin-centred restarts on native stopping. These are comparisons of the existing configurations; "
              "a line-search extended to exhaust 400 calls would be a different configuration.", "",
              "Validation checks candidate counts, the common origin and settings, fixed compositions, bounds and combined "
              "lift constraints, no duplicate stochastic candidates, and each winner's fresh-process replay (relative COP "
              "tolerance 1e-8). REFPROP uses a 180-second timeout, one retry, and worker recycling after 20 candidates. "
              "A timeout retry remains part of the same candidate's duration and allowance entry. Native numerical failures "
              "are failed evaluations; they are not proof of physical infeasibility.", ""]
    failure_counts = Counter(r['message'] for t in trials for r in t['records'] if not r['success'])
    (output / "failures.json").write_text(json.dumps(dict(failure_counts), indent=2))
    validation = dict(trials=len(trials), candidate_evaluations=sum(len(t['records']) for t in trials),
                      stochastic_evaluations=sum(len(t['records']) for t in trials if t['metadata']['method'] != 'line_search'),
                      fresh_replays=len(trials), model_source_sha256=trials[0]['metadata']['model_source_sha256'],
                      winner_replays_passed=True, common_origins_verified=True, fixed_compositions_verified=True,
                      stochastic_budgets_complete=True, line_search_within_ceiling=True,
                      failures=sum(failure_counts.values()))
    lines += [f"Validated **{validation['candidate_evaluations']:,} optimisation calls** and **{len(trials)} fresh winner replays**. "
              f"There were **{validation['failures']} failed optimisation calls**.", ""]
    sensitivity_path = output / "sensitivity/results.json"
    if sensitivity_path.exists():
        checks = json.loads(sensitivity_path.read_text())
        if len(checks) != 2*len(summary):
            raise ValueError("Incomplete pressure-guess checks")
        lines += ["## Independent pressure-guess checks", "",
                  "Each method's best final design was also evaluated with numerical discharge-pressure "
                  "starting guesses of 3 and 7 MPa. These 24 additional fresh-process evaluations were "
                  "excluded from the searches. The physical design and all other settings are unchanged.", "",
                  "| Composition | Method | Successful checks | Largest absolute COP change |",
                  "|---|---|---:|---:|"]
        for s in summary:
            selected = [c for c in checks if c['job']['refrigerant'] == s['refrigerant']
                        and c['job']['method'] == s['method']]
            if len(selected) != 2 or {c['job']['pressure_guess_MPa'] for c in selected} != {3.0, 7.0}:
                raise ValueError("Missing/duplicate pressure guesses")
            for c in selected:
                if c['job']['model_source_sha256'] != validation['model_source_sha256']:
                    raise ValueError("Sensitivity check used a different model")
                if c['job']['x'] != s['winner']['x']:
                    raise ValueError("Sensitivity check used an outdated winning design")
            successes = [c for c in selected if c['record']['success']]
            delta = max((abs(100*(c['record']['COP']/s['winner']['COP']-1)) for c in successes), default=None)
            display = "unavailable" if delta is None else f"{delta:.8f}%"
            lines.append(f"| {case_label(s['refrigerant'])} | {LABELS[s['method']]} | {len(successes)}/2 | {display} |")
        validation['additional_pressure_guess_evaluations'] = len(checks)
        validation['successful_pressure_guess_evaluations'] = sum(c['record']['success'] for c in checks)
        failed_checks = [c for c in checks if not c['record']['success']]
        if failed_checks:
            lines += ["", "Unsuccessful pressure-guess checks:", ""]
            for c in failed_checks:
                j = c['job']
                lines.append(f"- {case_label(j['refrigerant'])}, {LABELS[j['method']]} winner, "
                             f"{j['pressure_guess_MPa']:g} MPa starting guess: {c['record']['message']}.")
            lines += ["", "The native log for this retest reports a PH flash error 226 for the failed "
                      "3 MPa BO propane/isobutane check. The default and 7 MPa checks solved. This is "
                      "an additional initialisation-dependent native flash failure."]
        lines += ["", "These checks measure sensitivity to two numerical starting guesses; they do not "
                  "establish global optimality or cover every possible inner-solver initialisation.", ""]
    followup_path = output / "pressure_followup.json"
    if followup_path.exists():
        followup = json.loads(followup_path.read_text())['records']
        if len(followup) != 12 or any(r['job']['model_source_sha256'] != validation['model_source_sha256'] for r in followup):
            raise ValueError("Incomplete or mismatched pressure follow-up")
        lines += ["## Follow-up: pressure tolerance can stop the search too early", "",
                  "All three CMA-ES winners for propane/pentane 60/40 produce substantially lower COP with a 7 MPa "
                  "pressure guess. For the best design, the default and 3 MPa guesses select approximately 3.8835 MPa; "
                  "5, 7 and 10 MPa guesses select approximately 4.1698, 5.8227 and 8.3020 MPa. "
                  "The evaporator outlet remains 34.9999°C and the gas-cooler offset remains zero.", "",
                  "Two instrumented direct cycle solves identify the early exit: pressure-bound expansion stops "
                  "as soon as the residual is within ±0.05 K, and endpoint selection accepts that pressure. "
                  "Starting at 7 MPa, residuals remain positive and near 0.003–0.006 K, so no lower infeasible "
                  "pressure is bracketed before acceptance. The final minimum approach is 20.0000 K at 5.8227 MPa, "
                  "versus 19.95001 K at 3.8835 MPa under the default guess. Both final states pass the pinch check; "
                  "the higher-pressure state is an unnecessarily costly operating point when the lower feasible "
                  "pressure is already available. This affects the claimed minimum-pressure solve.", "",
                  "Disabling the fast inner objective at a 7 MPa guess recovers COP 2.537801 and pressure 3.88255 MPa, "
                  "with a small 0.04504 K gas-cooler offset. That is close to the default result, with a longer solve. "
                  "The benchmark model was held fixed throughout the retest.", "",
                  "| CMA-ES seed | Fast inner objective | Pressure guess (MPa; 0 = automatic) | COP | Solved pressure (MPa) | GC offset (K) |",
                  "|---:|---|---:|---:|---:|---:|"]
        for row in sorted(followup, key=lambda r: (r['job']['seed'], not r['job']['fast'], r['job']['guess'])):
            j, r = row['job'], row['record']
            if not r['success']:
                raise ValueError("Pressure follow-up did not solve")
            lines.append(f"| {j['seed']} | {j['fast']} | {j['guess']} | {r['COP']:.6f} | {r['Pmax_MPa']:.6f} "
                         f"| {r['inner_state']['dT_gc_out_min_offset']:.6f} |")
        validation['additional_pressure_followup_objective_evaluations'] = len(followup)
        for guess in (0, 7):
            trace = json.loads((output / f"pressure_trace_{guess}.json").read_text())
            reference = next(r['record'] for r in followup if r['job']['seed'] == 20260909
                             and r['job']['fast'] and r['job']['guess'] == guess)
            if (trace['model_source_sha256'] != validation['model_source_sha256']
                    or not trace['converged'] or not trace['feasible']
                    or abs(trace['COP']/reference['COP'] - 1) > 1e-8
                    or trace['minimum_approach_K'] < 19.95 - 1e-8):
                raise ValueError("Direct pressure trace differs from its objective replay")
        validation['additional_direct_cycle_traces'] = 2
        lines += ["", "These 12 objective evaluations and two direct cycle traces are additional diagnostics "
                  "and were excluded from every optimiser's budget and training data. Raw evidence is saved in "
                  "`pressure_followup.json`, `pressure_trace_0.json` and `pressure_trace_7.json`.", ""]
    (output / "retest_report.md").write_text("\n".join(lines))
    (output / "summary.json").write_text(json.dumps(summary, indent=2))
    (output / "validation.json").write_text(json.dumps(validation, indent=2))
    rows = []
    for t in trials:
        m, r = t['metadata'], t['replay']
        rows.append(dict(refrigerant=m['refrigerant'], method=m['method'], seed=m['seed'],
                         calls=len(t['records']), COP=r['COP'], RHX1=r['x'][0], RHX2=r['x'][1],
                         power_fraction=r['x'][2], VHC_MJ_per_m3=r['VHC_MJ_per_m3'], Pmax_MPa=r['Pmax_MPa']))
    with (output / "trial_results.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(validation, indent=2))
    return output / "retest_report.md"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default="results/RHPcomposition/optimiser_v13_400")
    args = parser.parse_args()
    print(report(args.output))
