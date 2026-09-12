"""Validate and report the frozen v15 candidate, winner and near-zero replays."""

from collections import Counter
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / "results/RHPcomposition/v15_failure_fixes/verified"


def main():
    manifest = json.loads((OUTPUT / "source_manifest.json").read_text())
    for name, digest in manifest["files"].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
    requests = json.loads((OUTPUT / "requests.json").read_text())
    assert len(requests) == 82
    trials = []
    for request in requests:
        prefix = Path(request["prefix"])
        status = json.loads(prefix.with_suffix(".status.json").read_text())
        assert status["timed_out"] or status["returncode"] == 0
        path = prefix.with_suffix(".result.json")
        result = json.loads(path.read_text()) if path.exists() else None
        if result:
            assert result["request"] == request
            assert np.allclose(result["record"]["x"], request["x"], atol=1e-10, rtol=0)
            if result["record"]["success"]:
                assert result["state"]["gc_minimum_approach"] >= 19.95-1e-8
        role = request.get("role", "winner" if "winner" in prefix.name else "old_failure")
        trials.append(dict(id=prefix.name, role=role, request=request, status=status, result=result))
    groups = {role: [t for t in trials if t["role"] == role]
              for role in ("old_failure", "winner", "near_zero_probe", "random_control")}
    assert [len(groups[r]) for r in groups] == [44, 18, 8, 12]
    def success(trial):
        return bool(trial["result"] and trial["result"]["record"]["success"])
    metrics = ("SEC", "COP", "Pmax_MPa", "VHC_MJ_per_m3")
    winner_changes = {}
    for metric in metrics:
        assert all(success(t) for t in groups["winner"])
        winner_changes[metric] = max(abs(t["result"]["record"][metric]/t["request"]["original"][metric]-1)
                                     for t in groups["winner"])
        assert winner_changes[metric] <= 1e-6
    assert all(success(t) for t in groups["near_zero_probe"])
    assert all(success(t) for t in groups["random_control"])
    random_changes = {metric: max(abs(t["result"]["record"][metric]/t["request"]["original"][metric]-1)
                                 for t in groups["random_control"]) for metric in metrics}
    summary = dict(model_source_sha256=manifest["model_source_sha256"], evaluations=82,
                   results_by_group={role: dict(count=len(ts), successful=sum(success(t) for t in ts),
                                               timed_out=sum(t["status"]["timed_out"] for t in ts))
                                     for role, ts in groups.items()},
                   maximum_relative_winner_changes=winner_changes,
                   maximum_relative_random_control_changes=random_changes,
                   remaining_failure_messages=dict(Counter(t["result"]["record"]["message"] if t["result"]
                       else "diagnostic deadline" for t in groups["old_failure"] if not success(t))), trials=trials)
    (OUTPUT / "summary.json").write_text(json.dumps(summary, indent=2))
    rows = ["| Original candidate | v15 outcome | COP | Solve seconds |", "|---|---|---:|---:|"]
    for t in groups["old_failure"]:
        record = t["result"]["record"] if t["result"] else None
        message = "recovered" if success(t) else record["message"] if record else "diagnostic deadline"
        rows.append(f"| {t['id']} | {message} | " + (f"{record['COP']:.8f}" if success(t) else "—")
                    + f" | {t['status']['wall_seconds']:.2f} |")
    (OUTPUT / "candidate_results.md").write_text("\n".join(rows)+"\n")
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    for ax, seed in zip(axes, (20260909, 20260911)):
        probes = sorted((t for t in groups["near_zero_probe"] if str(seed) in t["id"]),
                        key=lambda t: t["request"]["x"][1])
        zero = probes[0]["result"]["record"]["COP"]
        x = [t["request"]["x"][1] for t in probes]
        y = [(t["result"]["record"]["COP"]-zero)*1e8 for t in probes]
        ax.plot(x, y, "o-", color="#007f73")
        ax.set(xscale="symlog", xlabel="Requested RHX2 lift (K)", ylabel="COP change from zero (×10⁻⁸)",
               title=f"Boundary candidate · seed {seed}")
        ax.set_xscale("symlog", linthresh=1e-9)
        ax.grid(alpha=.2)
    fig.suptitle("v15 near-zero checks · all eight probes feasible")
    fig.savefig(OUTPUT / "near_zero.png", dpi=180)
    plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k != "trials"}, indent=2))


if __name__ == "__main__":
    main()
