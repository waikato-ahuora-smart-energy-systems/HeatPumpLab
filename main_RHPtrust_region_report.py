"""Validate and report the completed fixed-composition trust-region benchmark."""

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
from main_RHPtrust_region_benchmark import METHODS, OUTPUT


NAMES = {"multistart_cobyqa": "Multistart COBYQA", "turbo_1": "TuRBO-1 (adapted)", "gp_ei": "BO (GP-EI)"}
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
    diagnostics = data.get("search_diagnostics", {})
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
                    and data.get("replay", {}).get("success") and np.isclose(
                        data["replay"]["COP"], incumbent["COP"], rtol=1e-6, atol=1e-9)),
                number_of_restarts=(len(diagnostics["starts"])-1 if "starts" in diagnostics
                                    else len(diagnostics.get("restarts", []))))


def write_report(output):
    output = Path(output)
    runs = []
    for ref in DEFAULT_CASES:
        key = hashlib.sha256(ref.encode()).hexdigest()[:10]
        for method in METHODS:
            for seed in SEEDS:
                runs.append(summarise_run(output / f"{key}_{method}_{seed}.json"))
    assert len({r["metadata"]["model_source_sha256"] for r in runs}) == 1
    assert len({r["metadata"]["search_source_sha256"] for r in runs}) == 1
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
                               restarts=[r["number_of_restarts"] for r in selected]))
    report = dict(groups=groups, runs=runs, total_search_evaluations=10800,
                  total_failures=sum(r["failed"] for r in runs),
                  total_oracle_seconds=sum(r["evaluation_wall_seconds"] for r in runs),
                  failed_oracle_seconds=sum(r["failed_oracle_seconds"] for r in runs),
                  failure_types=dict(sum((Counter(r["failure_types"]) for r in runs), Counter())),
                  replay_matches=sum(r["replay_matches"] for r in runs),
                  scipy_version="1.18.1")
    (output / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    lines = ["# COBYQA / TuRBO / BO comparison on v14", "",
             "400 unique candidate evaluations per run; three seeds per method and composition. "
             "All compositions are mass fractions, at a 200 °C sink target. "
             "Every run uses SciPy 1.18.1 and the same frozen v14 thermodynamic model.", "",
             "10,800 search evaluations; " + str(report["total_failures"]) + " failed candidate evaluations. "
             + str(report["replay_matches"]) + "/27 winners match in a fresh worker "
             "(SEC and COP relative tolerance 1e-6, absolute 1e-9). Replays are outside the search budget.", "",
             f"Failed candidates consumed {100*report['failed_oracle_seconds']/report['total_oracle_seconds']:.1f}% "
             "of accumulated oracle wall time. Failure classifications: "
             + ", ".join(f"{name}: {count}" for name, count in report["failure_types"].items()) + ".", "",
             f"Cross-run consistency: {audit['repeated_design_groups']} repeated-design groups, "
             f"covering {audit['repeated_evaluations']} evaluations; "
             f"{len(mixed)} groups with mixed success/failure and {len(mismatches)} with SEC "
             f"disagreement above 1e-6 relative. Maximum relative SEC spread: {max_spread:.3g}. "
             "[Consistency audit](consistency_audit.json).", "",
             "See [protocol](PROTOCOL.md), [source manifest](manifest.json), "
             "[full numerical summary](summary.json), [506-test regression log](pytest.log), "
             "and [three interruption/resume tests](resume_tests.log).", "",
             "## Final COP", "", "Incumbents are selected by minimum SEC. Ranges show the three seeds.", "",
             "| Composition | Method | Median COP | COP range | Failed / 1,200 |",
             "|---|---|---:|---:|---:|"]
    for g in groups:
        lines.append(f"| {LABELS[g['refrigerant']]} | {NAMES[g['method']]} | {g['median_COP']:.6f} | "
                     f"{g['min_COP']:.6f}–{g['max_COP']:.6f} | {g['failures']} |")
    lines += ["", "![Convergence](convergence.png)", "",
              "The lower panels show percentage COP gap to the best design observed in this experiment "
              "for each composition, not a proven global optimum. Zero gaps are plotted at 0.00001% "
              "to permit a logarithmic axis. Lines show medians and shading spans the three seeds.", "",
              "## Progress and cost", "",
              "Median incumbent COP at each budget. Timing is median accumulated proposal/oracle "
              "wall time per run, measured under concurrent execution, excluding final replay. "
              "It is not a serial hardware benchmark.", "",
              "| Composition | Method | 100 | 200 | 400 | Proposal s | Oracle s | Restarts by seed |",
              "|---|---|---:|---:|---:|---:|---:|---|"]
    for g in groups:
        lines.append(f"| {LABELS[g['refrigerant']]} | {NAMES[g['method']]} | "
                     + " | ".join(f"{g['median_COP_at'][str(n)]:.6f}" for n in (100, 200, 400))
                     + f" | {g['median_proposal_seconds']:.1f} | {g['median_oracle_seconds']:.1f} | {g['restarts']} |")
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
    lines += ["", "## Limits of the comparison", "",
              "Three seeds show variability but do not establish statistical superiority. "
              "The TuRBO run uses an explicitly documented NumPy/SciPy adaptation: constrained "
              "candidate rejection, a failed-region heuristic and L-BFGS-B GP fitting. It is not "
              "a run of the authors' GPyTorch package. BO retains its existing initialisation, "
              "boundary candidates and GP settings. COBYQA uses a finite penalty for failed solves. "
              "Results compare these complete implementations at the fixed budget.", ""]
    (output / "REPORT.md").write_text("\n".join(lines))
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    colours = dict(zip(METHODS, ("#007f73", "#d47c17", "#425dc7")))
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
