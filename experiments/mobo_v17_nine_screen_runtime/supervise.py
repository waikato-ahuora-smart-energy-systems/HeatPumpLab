"""Resumable four-slot queue for the nine new 800-evaluation screens."""

import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT.parents[1] / "results/RHPcomposition/mobo_v17_nine_screen_800"
ORDER = (
    "Propane_Cyclopentane",
    "Butane_Pentane",
    "Isobutane_Pentane",
    "Butane_Isopentane",
    "Isobutane_Isopentane",
    "Butane_Cyclopentane",
    "Isobutane_Cyclopentane",
    "Butane_Isobutane",
    "Pentane_Isopentane",
)


def save(path, data):
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def progress(state):
    rows = []
    for key in ORDER:
        path = OUTPUT / f"{key}.json"
        d = json.loads(path.read_text()) if path.exists() else dict(records=[])
        good = [
            r for r in d["records"] if r["success"] and not r.get("validation_rejected")
        ]
        rows.append(
            dict(
                pair=key,
                evaluations=len(d["records"]),
                budget=800,
                complete=d.get("complete", False),
                best_COP=max((r["COP"] for r in good), default=None),
                failures=sum(not r["success"] for r in d["records"]),
                frontier_replays=len(d.get("frontier_replays", {})),
            )
        )
    save(
        OUTPUT / "progress.json",
        dict(
            status=state["status"],
            updated_utc=time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            families=rows,
        ),
    )
    lines = [
        "# Nine additional binary-fluid screens",
        "",
        f"Status: **{state['status']}**. Four concurrent searches, one numerical thread each.",
        "",
        "| Pair | Evaluations | Best observed COP | Failed | Frontier replays |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in rows:
        cop = f"{r['best_COP']:.6f}" if r["best_COP"] is not None else "—"
        lines.append(
            f"| {r['pair'].replace('_', '/')} | {r['evaluations']} / 800 | {cop} | {r['failures']} | {r['frontier_replays']} |"
        )
    lines += ["", "[Protocol](PROTOCOL.md)"]
    if state["status"] == "complete":
        lines += ["", "[Final report](REPORT.md)"]
    (OUTPUT / "STATUS.tmp").write_text("\n".join(lines) + "\n")
    (OUTPUT / "STATUS.tmp").replace(OUTPUT / "STATUS.md")


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for name, digest in json.loads((ROOT / "source_manifest.json").read_text())[
            "files"
        ].items():
            assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == digest, (
                name
            )
        for name in (
            "source_manifest.json",
            "source_snapshot.zip",
            "PROTOCOL.md",
            "preflight.json",
        ):
            shutil.copy2(ROOT / name, OUTPUT / name)
        pending = [
            k
            for k in ORDER
            if not (OUTPUT / f"{k}.json").exists()
            or not json.loads((OUTPUT / f"{k}.json").read_text()).get("complete")
        ]
        state = dict(
            status="running",
            pid=os.getpid(),
            process_group=os.getpgrp(),
            started_utc=time.strftime("%Y-%m-%d %H:%M:%S UTC", time.gmtime()),
            runtime=str(ROOT),
            output=str(OUTPUT),
            concurrent_searches=4,
            numerical_threads_per_search=1,
            budget_per_pair=800,
            search_pids={},
            exit_codes={},
        )
        children = {}
        env = dict(
            os.environ,
            OPENBLAS_NUM_THREADS="1",
            OMP_NUM_THREADS="1",
            MKL_NUM_THREADS="1",
            VECLIB_MAXIMUM_THREADS="1",
            NUMEXPR_NUM_THREADS="1",
            OMP_THREAD_LIMIT="1",
            MPLCONFIGDIR="/tmp/heatpumplab-botorch-mpl",
            TORCH_EXTENSIONS_DIR="/tmp/heatpumplab-torch-extensions",
        )
        try:
            while pending or children:
                for key, child in list(children.items()):
                    if child.poll() is not None:
                        state["exit_codes"][key] = child.returncode
                        del children[key]
                while pending and len(children) < 4:
                    key = pending.pop(0)
                    with (OUTPUT / f"{key}.log").open("a") as log:
                        children[key] = subprocess.Popen(
                            [
                                sys.executable,
                                "-u",
                                str(ROOT / "main_RHPmobo.py"),
                                "--partner",
                                key,
                                "--budget",
                                "800",
                                "--workers",
                                "1",
                                "--output",
                                str(OUTPUT),
                            ],
                            cwd=ROOT,
                            env=env,
                            stdin=subprocess.DEVNULL,
                            stdout=log,
                            stderr=subprocess.STDOUT,
                        )
                    state["search_pids"][key] = children[key].pid
                state["active_pairs"] = list(children)
                state["pending_pairs"] = pending.copy()
                state["status"] = (
                    "running with failed family"
                    if any(state["exit_codes"].values())
                    else "running"
                )
                save(OUTPUT / "runner_state.json", state)
                progress(state)
                if children:
                    time.sleep(20)
            if any(state["exit_codes"].values()):
                raise RuntimeError(state["exit_codes"])
            state["status"] = "building report"
            save(OUTPUT / "runner_state.json", state)
            progress(state)
            subprocess.run(
                [sys.executable, str(ROOT / "report.py")], cwd=ROOT, env=env, check=True
            )
            state["status"] = "complete"
            state["finished_utc"] = time.strftime(
                "%Y-%m-%d %H:%M:%S UTC", time.gmtime()
            )
        except BaseException as exc:
            state.update(status="failed", error=repr(exc))
            for child in children.values():
                if child.poll() is None:
                    child.terminate()
            raise
        finally:
            save(OUTPUT / "runner_state.json", state)
            progress(state)


if __name__ == "__main__":
    main()
