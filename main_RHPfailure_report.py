"""Summarise the completed diagnostic investigation without changing benchmarks."""

from collections import Counter, defaultdict
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from main_RHPfailure_investigation import OUTPUT, SOURCE, ROOT


def load(path):
    return json.loads(Path(path).read_text())


def kind(request):
    errors = request["original"]["constraint_errors"]
    if "evaluation_timeout" in errors:
        return "Native timeout"
    if "heat_exchanger_temperature_reversal" in errors:
        return "Tiny temperature reversal"
    if "gas_cooler_pinch_shortfall" in errors:
        return "Final pinch shortfall"
    return "REFPROP PH error 226"


def main():
    base = load(OUTPUT / "baseline_requests.json")
    entries, counts = [], defaultdict(Counter)
    all_successes = defaultdict(list)
    case_bests = defaultdict(float)
    for path in SOURCE.glob("*_*.json"):
        data = load(path)
        if "records" not in data:
            continue
        ref = data["metadata"]["refrigerant"]
        good = [r for r in data["records"] if r["success"]]
        all_successes[ref].extend(good)
        case_bests[ref] = max(case_bests[ref], max(r["COP"] for r in good))
    for request in base:
        stem = Path(request["prefix"]).name.removesuffix("_baseline")
        method = "BoTorch" if "botorch" in stem else "Custom BO"
        counts[method][kind(request)] += 1
        variants = {}
        for label in ("baseline", "pressure3", "pressure7", "tol8", "guard", "boundary_details"):
            prefix = OUTPUT / f"{stem}_{label}"
            status_path = prefix.with_suffix(".status.json")
            if not status_path.exists():
                continue
            status = load(status_path)
            result_path = prefix.with_suffix(".result.json")
            result = load(result_path) if result_path.exists() else None
            assert status["timed_out"] or status["returncode"] == 0
            variants[label] = dict(status=status, result=result)
        assert not variants["baseline"]["result"] or not variants["baseline"]["result"]["record"]["success"]
        assert variants["baseline"]["status"]["timed_out"] == (kind(request) == "Native timeout")
        recovered = [label for label, v in variants.items() if v["result"] and v["result"]["record"]["success"]]
        points = np.asarray([r["x"][:3] for r in all_successes[request["ref"]]])
        distance = np.linalg.norm((points-np.asarray(request["x"][:3]))/[100, 100, .8], axis=1)
        nearest = int(distance.argmin())
        for label in recovered:
            result = variants[label]["result"]
            assert result["state"]["gc_minimum_approach"] >= 19.95-1e-8
            assert result["record"]["COP"] < case_bests[request["ref"]]
        entries.append(dict(id=stem, method=method, ref=request["ref"], x=request["x"],
                            evaluation=request["evaluation"], failure_kind=kind(request),
                            recovered_variants=recovered, variants=variants,
                            nearest_success_distance=float(distance[nearest]),
                            nearest_success=all_successes[request["ref"]][nearest]))
    statuses = [load(f) for f in OUTPUT.glob("*.status.json")]
    assert len(statuses) == 155
    assert all(s["timed_out"] or s["returncode"] == 0 for s in statuses)
    guard_controls = []
    for f in OUTPUT.glob("*winner_guard.result.json"):
        data = load(f)
        assert data["record"]["success"]
        assert all(data["record"][k] == data["request"]["original"][k]
                   for k in ("SEC", "COP", "Pmax_MPa", "VHC_MJ_per_m3"))
        guard_controls.append(f.name)
    assert len(guard_controls) == 3
    frozen = load(SOURCE / "manifest.json")
    assert all(hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == sha for name, sha in frozen["files"].items())
    summary = dict(original_failures=44, original_outcomes_reproduced=44,
                   diagnostic_processes=155, diagnostic_timeouts=sum(s["timed_out"] for s in statuses),
                   counts_by_method={k: dict(v) for k, v in counts.items()},
                   recovered_candidates=sum(bool(e["recovered_variants"]) for e in entries),
                   recovered_pressure3=sum("pressure3" in e["recovered_variants"] for e in entries),
                   recovered_pressure7=sum("pressure7" in e["recovered_variants"] for e in entries),
                   recovered_tighter_tolerance=sum("tol8" in e["recovered_variants"] for e in entries),
                   guarded_winners_exact=guard_controls, frozen_source_files_verified=len(frozen["files"]),
                   entries=entries)
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2))
    names = {"Propane;Pentane|0.3;0.7|mass": "Propane/Pentane 30/70",
             "Propane;Pentane|0.6;0.4|mass": "Propane/Pentane 60/40",
             "Propane;Isobutane|0.3;0.7|mass": "Propane/Isobutane 30/70"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.7), constrained_layout=True)
    styles = {"REFPROP PH error 226": ("x", "#b12b42"), "Native timeout": ("^", "#dc8518"),
              "Tiny temperature reversal": ("s", "#315eb2"), "Final pinch shortfall": ("D", "#8a53a2")}
    for ax, (ref, label) in zip(axes, names.items()):
        points = np.asarray([r["x"][:3] for r in all_successes[ref]])
        ax.scatter(points[:, 0], points[:, 1], s=7, color="#b7bec5", alpha=.3, label="Successful benchmark calls")
        for category, (marker, colour) in styles.items():
            selected = [e for e in entries if e["ref"] == ref and e["failure_kind"] == category]
            if selected:
                points = np.asarray([e["x"][:3] for e in selected])
                ax.scatter(points[:, 0], points[:, 1], s=45, marker=marker, color=colour, label=category, zorder=3)
        recovered = [e for e in entries if e["ref"] == ref and e["recovered_variants"]]
        if recovered:
            points = np.asarray([e["x"][:3] for e in recovered])
            ax.scatter(points[:, 0], points[:, 1], s=110, facecolors="none", edgecolors="#007f65",
                       linewidths=1.3, label="Recovered in sensitivity replay", zorder=4)
        ax.set(title=label, xlabel="RHX1 lift (K)", ylabel="RHX2 lift (K)", xlim=(-2, 102), ylim=(-2, 102))
        ax.grid(alpha=.15)
    handles, labels = [], []
    for ax in axes:
        h, ls = ax.get_legend_handles_labels()
        for handle, label in zip(h, ls):
            if label not in labels:
                handles.append(handle)
                labels.append(label)
    fig.legend(handles, labels, loc="outside lower center", ncol=3, fontsize=8)
    fig.suptitle("44 failed candidates · mass compositions · all compressor fractions projected into two dimensions")
    fig.savefig(OUTPUT / "failure_locations.png", dpi=180)
    fig.savefig(OUTPUT / "failure_locations.svg")
    plt.close(fig)
    rows = ["| Candidate | RHX1 K | RHX2 K | Comp-2 fraction | Cause | P=3 MPa | P=7 MPa | tol=1e-8 |",
            "|---|---:|---:|---:|---|---|---|---|"]
    for e in entries:
        cells = []
        for label in ("pressure3", "pressure7", "tol8"):
            v = e["variants"].get(label)
            cells.append("—" if v is None else "timeout" if v["status"]["timed_out"] else
                         f"COP {v['result']['record']['COP']:.8f}" if v["result"]["record"]["success"] else "failed")
        rows.append(f"| {e['id']} | {e['x'][0]:.6f} | {e['x'][1]:.6f} | {e['x'][2]:.6f} | {e['failure_kind']} | " + " | ".join(cells) + " |")
    (OUTPUT / "candidate_table.md").write_text("\n".join(rows)+"\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "entries"}, indent=2))


if __name__ == "__main__":
    main()
