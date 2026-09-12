"""400-evaluation, three-seed COBYQA / TuRBO / BO comparison on frozen v14."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import time

import numpy as np

from main_RHPoptimiser_benchmark import DEFAULT_CASES, DEFAULT_SOURCE, load_case, serialise
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import OPTIMISER_VERSION
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.origin_search import OriginSurrogateSearch
from rhp_trust_region_search import TurboOne, multistart_cobyqa


METHODS = ("multistart_cobyqa", "turbo_1", "gp_ei")
OUTPUT = "results/RHPcomposition/trust_region_v14_400"
ROOT = Path(__file__).resolve().parent


def run_trial(task):
    ref, method, seed, budget, source, output = task
    _, base_case = load_case(ref, source)
    source_paths = [ROOT / name for name in (
        "rhp_trust_region_search.py", "main_RHPtrust_region_benchmark.py", "main_RHPoptimiser_benchmark.py")]
    source_paths += [ROOT / "HeatPumpLab/cycles/RecuperativeHeatPump/optimisation" / name
                     for name in ("origin_search.py", "surrogate_search.py")]
    metadata = dict(version=1, refrigerant=ref, method=method, seed=seed, budget=budget,
                    evaluation_version=OPTIMISER_VERSION, base_case=base_case,
                    model_source_sha256=hashlib.sha256(b"".join(
                        p.read_bytes() for p in sorted((ROOT / "HeatPumpLab").rglob("*.py")))).hexdigest(),
                    search_source_sha256=hashlib.sha256(b"".join(p.read_bytes() for p in source_paths)).hexdigest(),
                    versions={m: version(m) for m in ("numpy", "scipy", "ctREFPROP")},
                    origin=[0., 0., 0.5], timeout_seconds=180, timeout_retries=1, worker_max_tasks=20,
                    protocol="PROTOCOL.md")
    path = Path(output) / f"{hashlib.sha256(ref.encode()).hexdigest()[:10]}_{method}_{seed}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    result = json.loads(path.read_text()) if path.exists() else dict(metadata=metadata, records=[])
    if result["metadata"] != metadata:
        raise ValueError(f"Checkpoint source/settings differ: {path}")
    if result.get("complete"):
        return str(path)
    options = dict(refrigerant=ref, base_case=base_case, exact_evaluation_timeout_seconds=180,
                   exact_evaluation_timeout_retries=1, exact_evaluation_worker_max_tasks=20)
    objective = _OuterDesignObjective(max_evals=max(1, budget-len(result["records"])), **options)
    saved_count = len(result["records"])
    cursor = 0
    timer = time.perf_counter()

    def checkpoint():
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False))
        temporary.replace(path)

    def evaluate(proposal):
        nonlocal cursor, timer
        if cursor < saved_count:
            item = result["records"][cursor]
            if not np.allclose(proposal, item["x"][:3], atol=1e-7, rtol=0):
                raise ValueError(f"Checkpoint candidate {cursor} does not reproduce")
            cursor += 1
            timer = time.perf_counter()
            return item
        proposal_seconds = time.perf_counter() - timer
        if cursor >= budget:
            raise AssertionError("Expensive budget overrun")
        start = time.perf_counter()
        before = objective.exact_evaluations
        record = objective.evaluate([*proposal, *objective.initial_mass_fractions])
        if objective.exact_evaluations != before + 1:
            raise AssertionError("Duplicate or uncounted exact evaluation")
        if not np.allclose(record.x[:3], proposal, atol=1e-7, rtol=0):
            raise AssertionError("Candidate unexpectedly required repair")
        if not np.allclose(record.x[3:], objective.initial_mass_fractions, atol=1e-12, rtol=0):
            raise AssertionError("Composition moved")
        item = serialise(record)
        item.update(proposal_seconds=proposal_seconds, evaluation_wall_seconds=time.perf_counter()-start)
        result["records"].append(item)
        cursor += 1
        checkpoint()
        if cursor == 1 or cursor % 25 == 0:
            best = max((r["COP"] for r in result["records"] if r["success"]), default=0.)
            print(ref, method, seed, cursor, f"best COP {best:.7f}", flush=True)
        timer = time.perf_counter()
        return item

    try:
        if method == "multistart_cobyqa":
            result["search_diagnostics"] = multistart_cobyqa(
                evaluate, max_lift=objective.max_recuperator_lift, seed=seed, budget=budget)
        else:
            search = (TurboOne(max_lift=objective.max_recuperator_lift, seed=seed)
                      if method == "turbo_1" else OriginSurrogateSearch(
                          method="gp_ei", max_lift=objective.max_recuperator_lift, seed=seed, origin=(0., 0., 0.5)))
            for _ in range(budget):
                proposal = search.ask()
                item = evaluate(proposal)
                search.tell(proposal, item["SEC"] if item["success"] else 1e6, item["success"])
            if method == "turbo_1":
                result["search_diagnostics"] = dict(restarts=search.restarts, trust_region_events=search.events)
    finally:
        objective.close()
    assert len(result["records"]) == budget
    successes = [r for r in result["records"] if r["success"]]
    if successes:
        best = min(successes, key=lambda r: r["SEC"])
        replay = _OuterDesignObjective(max_evals=1, **dict(options, exact_evaluation_worker_max_tasks=1))
        try:
            result["replay"] = serialise(replay.evaluate(best["x"]))
        finally:
            replay.close()
        result["replay_matches"] = bool(result["replay"]["success"] and np.isclose(
            result["replay"]["SEC"], best["SEC"], rtol=1e-6, atol=1e-9))
    result["complete"] = True
    checkpoint()
    return str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=int, default=400)
    parser.add_argument("--n-cpu", type=int, default=6)
    parser.add_argument("--seed", type=int, action="append")
    parser.add_argument("--method", choices=METHODS, action="append")
    parser.add_argument("--refrigerant", action="append")
    parser.add_argument("--source", default=DEFAULT_SOURCE)
    parser.add_argument("--output", default=OUTPUT)
    args = parser.parse_args()
    if args.budget < 6:
        parser.error("budget must be at least six")
    tasks = [(ref, method, seed, args.budget, args.source, args.output)
             for seed in (args.seed or (20260909, 20260910, 20260911))
             for ref in (args.refrigerant or DEFAULT_CASES)
             for method in (args.method or METHODS)]
    with ProcessPoolExecutor(max_workers=max(1, args.n_cpu)) as pool:
        futures = [pool.submit(run_trial, task) for task in dict.fromkeys(tasks)]
        errors = []
        for future in as_completed(futures):
            try:
                print("Completed", future.result(), flush=True)
            except Exception as exc:
                errors.append(repr(exc))
                print("FAILED", repr(exc), flush=True)
        if errors:
            raise RuntimeError(errors)


if __name__ == "__main__":
    main()
