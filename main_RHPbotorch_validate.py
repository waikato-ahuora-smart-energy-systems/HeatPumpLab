"""Audit completed BoTorch benchmark against its frozen source and baseline."""

import hashlib
import json
from pathlib import Path
import zipfile

import numpy as np

from main_RHPbotorch_report import DEFAULT_CASES, METHODS, SEEDS, summarise_run
from main_RHPbotorch_benchmark import OUTPUT


def validate(output=OUTPUT):
    output = Path(output)
    manifest = json.loads((output / "manifest.json").read_text())
    with zipfile.ZipFile(output / "source_snapshot.zip") as archive:
        for name, expected in manifest["files"].items():
            assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == expected, name
            assert hashlib.sha256(archive.read(name)).hexdigest() == expected, name
    provenance = json.loads((output / "baseline_provenance.json").read_text())
    for entry in provenance["files"]:
        for field in ("source", "copy"):
            assert hashlib.sha256(Path(entry[field]).read_bytes()).hexdigest() == entry["sha256"]
    model_hash = hashlib.sha256(b"".join(
        p.read_bytes() for p in sorted(Path("HeatPumpLab").rglob("*.py")))).hexdigest()
    runs = []
    bootstrap_matches = 0
    maximum_replay_relative_difference = dict.fromkeys(
        ("SEC", "COP", "Pmax_MPa", "VHC_MJ_per_m3"), 0.)
    for ref in DEFAULT_CASES:
        key = hashlib.sha256(ref.encode()).hexdigest()[:10]
        for seed in SEEDS:
            pair = []
            for method in METHODS:
                path = output / f"{key}_{method}_{seed}.json"
                data = json.loads(path.read_text())
                summary = summarise_run(path)
                assert summary["replay_matches"], path
                meta = data["metadata"]
                assert meta["model_source_sha256"] == model_hash
                assert meta["refrigerant"] == ref and meta["seed"] == seed and meta["method"] == method
                for name, version in meta["versions"].items():
                    assert manifest["versions"][name] == version
                assert meta["timeout_seconds"] == 180 and meta["timeout_retries"] == 1
                assert meta["worker_max_tasks"] == 20
                fractions = [float(x) for x in ref.split("|")[1].split(";")]
                assert all(np.allclose(row["x"][3:], fractions, rtol=0, atol=1e-12)
                           for row in data["records"])
                assert np.allclose(data["replay"]["x"], summary["best"]["x"], rtol=0, atol=1e-12)
                for metric in maximum_replay_relative_difference:
                    difference = abs(data["replay"][metric] - summary["best"][metric])
                    relative = difference / abs(summary["best"][metric])
                    maximum_replay_relative_difference[metric] = max(
                        relative, maximum_replay_relative_difference[metric])
                pair.append(data)
                runs.append(summary)
            assert pair[0]["metadata"]["base_case"] == pair[1]["metadata"]["base_case"]
            for new, old in zip(pair[0]["records"][:5], pair[1]["records"][:5], strict=True):
                assert all(new[field] == old[field] for field in ("x", "success", "SEC", "COP"))
                bootstrap_matches += 1
    assert len(runs) == 18 and bootstrap_matches == 45
    summary = dict(passed=True, completed_runs=18, new_runs=9, archived_runs=9,
                   unique_evaluations_per_run=400, total_evaluations=7200,
                   new_evaluations=3600, archived_evaluations=3600,
                   winner_replays_matched=18, new_winner_replays=9, archived_winner_replays=9,
                   replay_rtol=1e-6, replay_atol=1e-9,
                   maximum_replay_relative_difference=maximum_replay_relative_difference,
                   exact_shared_bootstrap_matches=bootstrap_matches,
                   frozen_files_verified=len(manifest["files"]), archived_files_verified=9,
                   model_source_sha256=model_hash, versions=manifest["versions"])
    (output / "validation_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return summary


if __name__ == "__main__":
    validate()
