"""Phase 1 check: does a pure fluid give the same result via the pure and mixture solver paths?

Each design is solved twice in fresh isolated workers:
  A) plain pure name, e.g. "Butane"            -> pure-fluid solver path
  B) zero-fraction binary, e.g. "Propane;Butane|0;1|mass" -> mixture solver path
at several discharge-pressure seeds (P_H_est, MPa; 0 = automatic).
Designs come from successful pure-vertex records of the v23 fixed-split runs.
Diagnostic only: no production source is changed.
"""
import json, os, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

V23 = Path("results/RHPcomposition/fixed200_balanced_v23_5000")
OUT = Path(sys.argv[sys.argv.index("--output") + 1]) if "--output" in sys.argv else Path("results/RHPcomposition/phase1_spelling_check")
SPELLINGS = ("pure", "mixture", "legacy") if "--with-legacy" in sys.argv else ("pure", "mixture")
SEEDS = (0.0, 3.0, 4.3, 5.5, 7.0)
PER_FLUID = 8


def select_designs():
    """Return {pure fluid: [(x3, partner_string, vertex_fraction_first)]}."""
    designs = {}
    rng = np.random.default_rng(20260928)
    for case_dir in sorted(p for p in V23.iterdir() if p.is_dir() and p.name != "snapshots"):
        first, partner = case_dir.name.split("_")
        data = json.loads((case_dir / "result.json").read_text())
        base_case = data["metadata"]["base_case"]
        for vertex, fluid in ((1.0, first), (0.0, partner)):
            good = [r for r in data["records"] if r["success"] and r["x"][3] == vertex]
            if not good:
                continue
            good.sort(key=lambda r: -r["COP"])
            picks = good[:4]
            rest = good[4:]
            if rest:
                picks += [rest[i] for i in rng.choice(len(rest), min(3, len(rest)), replace=False)]
            entry = designs.setdefault(fluid, {"base_case": base_case, "items": []})
            for r in picks:
                entry["items"].append(dict(x3=r["x"][:3], mixture=f"{first};{partner}",
                                           vertex=vertex, v23_COP=r["COP"], source=case_dir.name))
    out = {}
    for fluid, entry in designs.items():
        seen, items = set(), []
        # origin design first, then v23 designs, de-duplicated by cycle coordinates
        for item in [dict(x3=[0.0, 0.0, 0.5], mixture=entry["items"][0]["mixture"],
                          vertex=entry["items"][0]["vertex"], v23_COP=None, source="origin")] + entry["items"]:
            key = tuple(np.round(item["x3"], 9))
            if key not in seen:
                seen.add(key)
                items.append(item)
        out[fluid] = dict(base_case=entry["base_case"], items=items[:PER_FLUID])
    return out


def solve(task):
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
    base_case = dict(task["base_case"], P_H_est=task["seed"])
    if task["spelling"] == "legacy":
        base_case["legacy_pure_solver_path"] = True
    if task["spelling"] in ("pure", "legacy"):
        refrigerant, x = task["fluid"], list(task["x3"]) + [1.0]
    else:
        fr = f"{task['vertex']:g};{1 - task['vertex']:g}"
        refrigerant, x = f"{task['mixture']}|{fr}|mass", list(task["x3"]) + [task["vertex"], 1 - task["vertex"]]
    objective = _OuterDesignObjective(max_evals=1, base_case=base_case, refrigerant=refrigerant,
                                      results_dir=str(OUT / "scratch"), exact_evaluation_timeout_seconds=240.0,
                                      exact_evaluation_timeout_retries=0)
    start = time.perf_counter()
    try:
        record = objective.evaluate(x)
    finally:
        objective.close()
    rep = record.report or {}
    return dict(task, refrigerant=refrigerant, success=record.success, message=record.message,
                COP=rep.get("COP", record.cop), SEC=rep.get("SEC", record.sec), VHC=rep.get("VHC"),
                Pmax_kPa=rep.get("Pmax"), T_comp_discharge=rep.get("T_comp_discharge"),
                PR_1=rep.get("PR_1"), PR_2=rep.get("PR_2"),
                retries=len(rep.get("numerical_retries") or []), wall_s=time.perf_counter() - start)


def main():
    n_cpu = int(sys.argv[sys.argv.index("--n-cpu") + 1]) if "--n-cpu" in sys.argv else 6
    OUT.mkdir(parents=True, exist_ok=True)
    designs = select_designs()
    tasks = []
    for fluid, entry in designs.items():
        for d, item in enumerate(entry["items"]):
            for seed in SEEDS:
                for spelling in SPELLINGS:
                    tasks.append(dict(fluid=fluid, design=d, seed=seed, spelling=spelling,
                                      base_case=entry["base_case"], **item))
    summary = ", ".join(f"{f} ({len(e['items'])} designs)" for f, e in designs.items())
    print(f"{len(tasks)} solves: {summary}", flush=True)
    results_path = OUT / "results.jsonl"
    done = set()
    if results_path.exists():
        for line in results_path.read_text().splitlines():
            r = json.loads(line)
            done.add((r["fluid"], r["design"], r["seed"], r["spelling"]))
    todo = [t for t in tasks if (t["fluid"], t["design"], t["seed"], t["spelling"]) not in done]
    print(f"{len(done)} already done, {len(todo)} to run", flush=True)
    with ProcessPoolExecutor(n_cpu) as pool, results_path.open("a") as fh:
        futures = [pool.submit(solve, t) for t in todo]
        for i, fut in enumerate(as_completed(futures), 1):
            try:
                r = fut.result()
            except Exception as exc:
                print(f"task error: {exc}", flush=True)
                continue
            r.pop("base_case", None)
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            if i % 20 == 0 or i == len(futures):
                print(f"{i}/{len(futures)} done", flush=True)
    (OUT / "DONE").write_text(time.strftime("%Y-%m-%d %H:%M:%S"))
    print("finished", flush=True)


if __name__ == "__main__":
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    main()
