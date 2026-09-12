"""Compare new BoTorch LogEI trials with matched archived custom-BO trials."""

import argparse
from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from main_RHPoptimiser_benchmark import DEFAULT_CASES
from main_RHPbotorch_benchmark import OUTPUT

METHODS = ("botorch_log_ei", "gp_ei")


NAMES = {"botorch_log_ei": "BoTorch LogEI", "gp_ei": "Custom BO (archived)"}
LABELS = {DEFAULT_CASES[0]: "Propane/Pentane 30/70", DEFAULT_CASES[1]: "Propane/Pentane 60/40",
          DEFAULT_CASES[2]: "Propane/Isobutane 30/70"}
SEEDS = (20260909, 20260910, 20260911)


def failure_kind(record):
    if record["success"]:
        return "success"
    message = record["message"].lower()
    if "timeout" in message or "timed out" in message or "evaluation_timeout" in record["constraint_errors"]:
        return "timeout"
    if "226" in message:
        return "REFPROP 226"
    if "infeasible, unconverged" in message:
        return "cycle infeasible/unconverged"
    return "other solver/native error"


def summarise_run(path):
    data = json.loads(path.read_text())
    assert data.get("complete"), f"Incomplete: {path}"
    rows = data["records"]
    assert len(rows) == data["metadata"]["budget"] == 400
    assert len(set(tuple(np.round(r["x"], 8)) for r in rows)) == 400
    for row in rows:
        assert np.allclose(row["x"][3:], rows[0]["x"][3:], rtol=0, atol=1e-12)
        assert 0 <= row["x"][0] <= 100 and 0 <= row["x"][1] <= 100
        assert row["x"][0] + row["x"][1] <= 123.333333334
        assert 0.1 <= row["x"][2] <= 0.9
    incumbent = None
    curve, sec_curve = [], []
    for row in rows:
        if row["success"] and (incumbent is None or row["SEC"] < incumbent["SEC"]):
            incumbent = row
        curve.append(incumbent["COP"] if incumbent else None)
        sec_curve.append(incumbent["SEC"] if incumbent else None)
    assert incumbent is not None
    diagnostics = data.get("search_diagnostics", [])
    return dict(path=str(path), metadata=data["metadata"], best=incumbent,
                best_evaluation=rows.index(incumbent)+1,
                failed=sum(not r["success"] for r in rows),
                failure_types=dict(Counter(failure_kind(r) for r in rows if not r["success"])),
                proposal_seconds=sum(r["proposal_seconds"] for r in rows),
                evaluation_wall_seconds=sum(r["evaluation_wall_seconds"] for r in rows),
                failed_oracle_seconds=sum(r["evaluation_wall_seconds"] for r in rows if not r["success"]),
                COP_at={str(n): curve[n-1] for n in (100, 200, 400)},
                COP_curve=curve, SEC_curve=sec_curve,
                replay=data.get("replay"), replay_matches=bool(data.get("replay_matches", False)
                    and data.get("replay", {}).get("success") and all(np.isclose(
                        data["replay"][key], incumbent[key], rtol=1e-6, atol=1e-9)
                        for key in ("SEC", "COP", "Pmax_MPa", "VHC_MJ_per_m3"))),
                diagnostic_events=diagnostics if isinstance(diagnostics, list) else [])



def write_report(output):
    output = Path(output)
    runs = []
    for ref in DEFAULT_CASES:
        key = hashlib.sha256(ref.encode()).hexdigest()[:10]
        for method in METHODS:
            for seed in SEEDS:
                runs.append(summarise_run(output / f"{key}_{method}_{seed}.json"))
    assert len({r["metadata"]["model_source_sha256"] for r in runs}) == 1
    for method in METHODS:
        assert len({r["metadata"]["search_source_sha256"] for r in runs
                    if r["metadata"]["method"] == method}) == 1
    assert all(r["metadata"]["versions"]["scipy"] == "1.18.1" for r in runs)
    observations = defaultdict(list)
    for run in runs:
        for index, record in enumerate(json.loads(Path(run["path"]).read_text())["records"]):
            key = (run["metadata"]["refrigerant"], tuple(np.round(record["x"], 8)))
            observations[key].append((Path(run["path"]).name, index+1, record))
    repeated = [group for group in observations.values() if len(group) > 1]
    mixed, mismatches, max_spread = [], [], 0.
    for group in repeated:
        good = [entry for entry in group if entry[2]["success"]]
        if good and len(good) < len(group):
            mixed.append([(p, i, record["success"]) for p, i, record in group])
        if len(good) > 1:
            scores = np.asarray([entry[2]["SEC"] for entry in good])
            spread = float(np.ptp(scores)/scores.min())
            max_spread = max(max_spread, spread)
            if spread > 1e-6:
                mismatches.append([(p, i, record["SEC"]) for p, i, record in good])
    audit = dict(repeated_design_groups=len(repeated), repeated_evaluations=sum(map(len, repeated)),
                 max_relative_SEC_spread=max_spread, mixed_outcomes=mixed,
                 successful_score_mismatches=mismatches, complete=True)
    (output / "consistency_audit.json").write_text(json.dumps(audit, indent=2))
    groups = []
    for ref in DEFAULT_CASES:
        for method in METHODS:
            selected = [r for r in runs if r["metadata"]["refrigerant"] == ref
                        and r["metadata"]["method"] == method]
            cops = [r["best"]["COP"] for r in selected]
            best = min(selected, key=lambda r: r["best"]["SEC"])
            groups.append(dict(refrigerant=ref, method=method, median_COP=float(np.median(cops)),
                               min_COP=min(cops), max_COP=max(cops),
                               median_SEC=float(np.median([r["best"]["SEC"] for r in selected])),
                               failures=sum(r["failed"] for r in selected),
                               failed_oracle_seconds=sum(r["failed_oracle_seconds"] for r in selected),
                               median_proposal_seconds=float(np.median([r["proposal_seconds"] for r in selected])),
                               median_oracle_seconds=float(np.median([r["evaluation_wall_seconds"] for r in selected])),
                               median_COP_at={str(n): float(np.median([r["COP_at"][str(n)] for r in selected]))
                                              for n in (100, 200, 400)},
                               best_seed=best["metadata"]["seed"], best=best["best"],
                               diagnostic_event_counts=[len(r["diagnostic_events"]) for r in selected]))
    comparisons = []
    for ref in DEFAULT_CASES:
        new = next(g for g in groups if g["refrigerant"] == ref and g["method"] == "botorch_log_ei")
        old = next(g for g in groups if g["refrigerant"] == ref and g["method"] == "gp_ei")
        comparisons.append(dict(refrigerant=ref,
                                median_COP_change_percent=100*(new["median_COP"]/old["median_COP"]-1),
                                median_proposal_time_ratio=new["median_proposal_seconds"]/old["median_proposal_seconds"]))
    report = dict(groups=groups, comparisons=comparisons, runs=runs, total_search_evaluations=7200,
                  new_search_evaluations=3600, archived_search_evaluations=3600,
                  total_failures=sum(r["failed"] for r in runs),
                  total_oracle_seconds=sum(r["evaluation_wall_seconds"] for r in runs),
                  failed_oracle_seconds=sum(r["failed_oracle_seconds"] for r in runs),
                  failure_types=dict(sum((Counter(r["failure_types"]) for r in runs), Counter())),
                  replay_matches=sum(r["replay_matches"] for r in runs),
                  scipy_version="1.18.1")
    (output / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    lines = ["# BoTorch LogEI versus custom BO on v14", "",
             "400 unique candidate evaluations per run; three seeds per method and composition. "
             "All compositions are mass fractions, at a 200 °C sink target. "
             "Every run uses SciPy 1.18.1 and the same frozen v14 thermodynamic model.", "",
             "3,600 new BoTorch evaluations compared with 3,600 archived custom-BO evaluations; " + str(report["total_failures"]) + " failed candidate evaluations. "
             + str(report["replay_matches"]) + "/18 winners match their saved fresh-worker replays (nine new, nine archived) "
             "(SEC, COP, pressure and VHC relative tolerance 1e-6, absolute 1e-9). Replays are outside the search budget.", "",
             f"Failed candidates consumed {100*report['failed_oracle_seconds']/report['total_oracle_seconds']:.1f}% "
             "of accumulated oracle wall time. Failure classifications: "
             + ", ".join(f"{name}: {count}" for name, count in report["failure_types"].items()) + ".", "",
             f"Cross-run consistency: {audit['repeated_design_groups']} repeated-design groups, "
             f"covering {audit['repeated_evaluations']} evaluations; "
             f"{len(mixed)} groups with mixed success/failure and {len(mismatches)} with SEC "
             f"disagreement above 1e-6 relative. Maximum relative SEC spread: {max_spread:.3g}. "
             "[Consistency audit](consistency_audit.json).", "",
             "See [protocol](PROTOCOL.md), [source manifest](manifest.json), "
             "[full numerical summary](summary.json), [six implementation tests](search_tests.log), "
             "[validation](validation_summary.json), and [archived baseline provenance](baseline_provenance.json).", "",
             "## Final COP", "", "Incumbents are selected by minimum SEC. Ranges show the three seeds.", "",
             "| Composition | Method | Median COP | COP range | Failed / 1,200 |",
             "|---|---|---:|---:|---:|"]
    for g in groups:
        lines.append(f"| {LABELS[g['refrigerant']]} | {NAMES[g['method']]} | {g['median_COP']:.6f} | "
                     f"{g['min_COP']:.6f}–{g['max_COP']:.6f} | {g['failures']} |")
    lines += ["", "BoTorch's median COP change relative to archived custom BO: "
              + "; ".join(f"{LABELS[c['refrigerant']]} {c['median_COP_change_percent']:+.4f}%"
                          for c in comparisons) + "."]
    lines += ["", "![Convergence](convergence.png)", "",
              "The lower panels show percentage COP gap to the best design observed in this experiment "
              "for each composition, not a proven global optimum. Zero gaps are plotted at 0.00001% "
              "to permit a logarithmic axis. Lines show medians and shading spans the three seeds.", "",
              "## Progress and cost", "",
              "Median incumbent COP at each budget. Timing is median accumulated proposal/oracle "
              "wall time per run, measured under concurrent execution, excluding final replay. "
              "The custom runs are from the preceding experiment; these are not matched serial hardware timings.", "",
              "| Composition | Method | 100 | 200 | 400 | Proposal s | Oracle s | Diagnostic events by seed |",
              "|---|---|---:|---:|---:|---:|---:|---|"]
    for g in groups:
        lines.append(f"| {LABELS[g['refrigerant']]} | {NAMES[g['method']]} | "
                     + " | ".join(f"{g['median_COP_at'][str(n)]:.6f}" for n in (100, 200, 400))
                     + f" | {g['median_proposal_seconds']:.1f} | {g['median_oracle_seconds']:.1f} | "
                     + (str(g['diagnostic_event_counts']) if g['method'] == 'botorch_log_ei' else 'not recorded') + " |")
    lines += ["", "## Best design per method and composition", "",
              "| Composition | Method | Seed | COP | SEC | RHX1 K | RHX2 K | Compressor-2 fraction |",
              "|---|---|---:|---:|---:|---:|---:|---:|"]
    for g in groups:
        b = g["best"]
        lines.append(f"| {LABELS[g['refrigerant']]} | {NAMES[g['method']]} | {g['best_seed']} | "
                     f"{b['COP']:.8f} | {b['SEC']:.8f} | {b['x'][0]:.6f} | {b['x'][1]:.6f} | {b['x'][2]:.8f} |")
    lines += ["", "## Per-seed outcomes", "",
              "| Composition | Method | Seed | Final COP | Best evaluation | Failures | Replay match |",
              "|---|---|---:|---:|---:|---:|---|"]
    for r in runs:
        m = r["metadata"]
        lines.append(f"| {LABELS[m['refrigerant']]} | {NAMES[m['method']]} | {m['seed']} | "
                     f"{r['best']['COP']:.8f} | {r['best_evaluation']} | {r['failed']} | {r['replay_matches']} |")
    events = Counter(event["stage"] for run in runs for event in run["diagnostic_events"])
    lines += ["", "## BoTorch diagnostics", "",
              "Recorded search diagnostic events: " + (", ".join(
                  f"{stage}: {count}" for stage, count in sorted(events.items())) or "none") + ". "
              "Warnings do not necessarily imply a failed fit or fallback. Full messages and evaluation "
              "numbers are retained in each trial checkpoint and the numerical summary. The fitting "
              "stage tries two length-scale initialisations; a GP-fit event records an exception from "
              "one attempt. Acquisition fallbacks and nonfinite acquisition values have separate event types.", "",
              "## Limits of the comparison", "",
              "Three seeds describe variability but do not establish statistical superiority. "
              "This is a controlled BoTorch implementation: it retains the custom BO bootstrap, "
              "Matérn GP assumptions, failed-region heuristic, boundary candidate pool and periodic "
              "maximin exploration. BoTorch supplies the GP fit, LogEI and continuous acquisition "
              "optimisation. It does not use all of BoTorch's default GP settings. The experiment "
              "does not isolate package fitting differences from the added continuous acquisition "
              "optimisation. Archived baseline observations never enter the BoTorch search.", "",
              "## Reproduction", "",
              "Install the locked benchmark dependencies with `uv sync --group benchmark --group dev`. "
              "The native REFPROP installation and original baseline input are also required. "
              "Run `python main_RHPbotorch_benchmark.py --n-cpu 9` from the repository root, with "
              "OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and VECLIB_MAXIMUM_THREADS set to 1. "
              "Existing completed checkpoints are reused; use a separate output directory for new trials. "
              "The comparison folder must also contain the nine custom-BO checkpoints listed in "
              "baseline_provenance.json. Generate the report with `python main_RHPbotorch_report.py` "
              "and audit this frozen experiment with `python main_RHPbotorch_validate.py`.", ""]
    (output / "REPORT.md").write_text("\n".join(lines))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    colours = {"botorch_log_ei": "#007f73", "gp_ei": "#425dc7"}
    for column, ref in enumerate(DEFAULT_CASES):
        ax, gap_ax = axes[:, column]
        best_observed = max(r["best"]["COP"] for r in runs if r["metadata"]["refrigerant"] == ref)
        for method in METHODS:
            curves = np.asarray([r["COP_curve"] for r in runs if r["metadata"]["refrigerant"] == ref
                                 and r["metadata"]["method"] == method])
            x = np.arange(1, 401)
            ax.plot(x, np.median(curves, axis=0), color=colours[method], label=NAMES[method], lw=1.8)
            ax.fill_between(x, curves.min(axis=0), curves.max(axis=0), color=colours[method], alpha=0.13)
            gaps = np.maximum(100 * (best_observed - curves) / best_observed, 1e-5)
            gap_ax.plot(x, np.median(gaps, axis=0), color=colours[method], lw=1.8)
            gap_ax.fill_between(x, gaps.min(axis=0), gaps.max(axis=0), color=colours[method], alpha=0.13)
        ax.set(title=LABELS[ref], xlabel="Expensive candidate evaluations", ylabel="Incumbent COP", xlim=(1, 400))
        ax.grid(alpha=0.2)
        gap_ax.set(xlabel="Expensive candidate evaluations", ylabel="Gap to best observed COP (%)",
                   xlim=(1, 400), yscale="log", ylim=(1e-5, 20))
        gap_ax.grid(alpha=0.2)
    axes[0, 0].legend(loc="lower right", fontsize=8)
    fig.suptitle("400 evaluations · median and range over 3 seeds · frozen v14 · SciPy 1.18.1")
    fig.savefig(output / "convergence.png", dpi=180)
    fig.savefig(output / "convergence.svg")
    plt.close(fig)
    print(output / "REPORT.md")
    for g in groups:
        print(LABELS[g["refrigerant"]], NAMES[g["method"]], round(g["median_COP"], 7),
              "range", round(g["min_COP"], 7), round(g["max_COP"], 7), "failures", g["failures"])
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=OUTPUT)
    write_report(parser.parse_args().output)
