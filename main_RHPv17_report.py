"""Audit the production internal PH fallback against frozen v15 replays."""

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results/RHPcomposition/v17_internal_ph_fallback"


def main():
    manifest = json.loads((OUT / "source_manifest.json").read_text())
    for name, digest in manifest["files"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, name
    requests = json.loads((OUT / "requests.json").read_text())
    assert len(requests) == 82
    trials = []
    for request in requests:
        prefix = Path(request["prefix"])
        status = json.loads(prefix.with_suffix(".status.json").read_text())
        path = prefix.with_suffix(".result.json")
        result = json.loads(path.read_text()) if path.exists() else None
        assert status["timed_out"] or status["returncode"] == 0
        record = result["record"] if result else None
        if result:
            assert result["request"] == request
            assert all(abs(a-b) <= 1e-10 for a, b in zip(record["x"], request["x"], strict=True))
        success = bool(record and record["success"])
        changes = {}
        if success:
            state = result["state"]
            assert state["gc_minimum_approach"] >= 19.95-1e-8
            volume_flows = state["compressor_suction_volume_flows"]
            assert volume_flows and all(v > 0 for v in volume_flows)
            expected_vhc = state["total_heat_duty"] / sum(volume_flows) / 1000
            assert abs(record["VHC_MJ_per_m3"]/expected_vhc-1) <= 1e-12
            if request["v15_record"]["success"]:
                for metric in ("SEC", "COP", "Pmax_MPa"):
                    changes[metric] = abs(record[metric]/request["v15_record"][metric]-1)
                first_suction_vhc = state["total_heat_duty"]/volume_flows[0]/1000
                changes["VHC_first_suction_for_legacy_comparison"] = abs(
                    first_suction_vhc/request["v15_record"]["VHC_MJ_per_m3"]-1)
        trials.append(dict(id=prefix.name, role=request["role"], success=success,
                           v15_success=request["v15_record"]["success"], status=status,
                           message=record["message"] if record else "Diagnostic timeout",
                           COP=record["COP"] if record else None,
                           VHC_MJ_per_m3=record["VHC_MJ_per_m3"] if record else None,
                           nearby_recovery_attempts=(result["profile"]["counts"].get("refprop_ph_nearby_recovery", 0)
                                                     if result else None),
                           relative_changes=changes))
    groups = {}
    for role in ("old_failure", "winner", "near_zero_probe", "random_control"):
        selected = [t for t in trials if t["role"] == role]
        groups[role] = dict(count=len(selected), success=sum(t["success"] for t in selected),
                            timed_out=sum(t["status"]["timed_out"] for t in selected),
                            regressions=sum(t["v15_success"] and not t["success"] for t in selected))
    old_failures = [t for t in trials if t["role"] == "old_failure"]
    changes = {role: {metric: max((t["relative_changes"].get(metric, 0.) for t in trials
                                    if t["role"] == role and t["v15_success"]), default=0.)
                      for metric in ("SEC", "COP", "Pmax_MPa", "VHC_first_suction_for_legacy_comparison")}
               for role in groups}
    assert [groups[role]["count"] for role in groups] == [44, 18, 8, 12]
    assert all(t["success"] for t in trials if t["v15_success"])
    assert all(value <= 1e-6 for metrics in changes.values() for value in metrics.values())
    summary = dict(model_source_sha256=manifest["model_source_sha256"], groups=groups,
                   recovered_of_remaining_25=sum(t["success"] and not t["v15_success"] for t in old_failures),
                   retained_of_previous_19=sum(t["success"] and t["v15_success"] for t in old_failures),
                   maximum_relative_changes=changes, trials=trials)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    rows = ["| Candidate | v15 solved | v17 outcome | COP | Wall seconds |", "|---|---|---|---:|---:|"]
    for t in old_failures:
        cop = f"{t['COP']:.8f}" if t["success"] else "—"
        outcome = "success" if t["success"] else t["message"]
        rows.append(f"| {t['id']} | {t['v15_success']} | {outcome} | {cop} | {t['status']['wall_seconds']:.2f} |")
    (OUT / "candidate_outcomes.md").write_text("\n".join(rows)+"\n")
    print(json.dumps({k: v for k, v in summary.items() if k != "trials"}, indent=2))


if __name__ == "__main__":
    main()
