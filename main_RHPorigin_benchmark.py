"""Repeat BO, RBF and CMA-ES from the physical origin at fixed composition."""

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
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.origin_search import OriginCMAES, OriginSurrogateSearch


ORIGIN = (0.0, 0.0, 0.5)
METHODS = ("gp_ei", "rbf", "cma_es")
OUTPUT = "results/RHPcomposition/origin_200"


def make_search(method, max_lift, seed):
    if method == "cma_es":
        return OriginCMAES(max_lift=max_lift, origin=ORIGIN, seed=seed)
    return OriginSurrogateSearch(method=method, max_lift=max_lift, origin=ORIGIN, seed=seed)


def run_trial(task):
    ref, method, seed, budget, source, output = task
    baseline, base_case = load_case(ref, source)
    code_dir = Path(__file__).resolve().parent / "HeatPumpLab/cycles/RecuperativeHeatPump/optimisation"
    source_hash = hashlib.sha256(b"".join((code_dir / f).read_bytes() for f in
                                ("origin_search.py", "surrogate_search.py"))).hexdigest()
    model_root = Path(__file__).resolve().parent / "HeatPumpLab"
    model_hash = hashlib.sha256(b"".join(p.read_bytes() for p in sorted(model_root.rglob("*.py")))).hexdigest()
    metadata = dict(version=1, refrigerant=ref, method=method, seed=seed, budget=budget,
                    evaluation_version=OPTIMISER_VERSION,
                    model_source_sha256=model_hash,
                    origin=list(ORIGIN), base_case=base_case, search_source_sha256=source_hash,
                    versions={m: version(m) for m in ("numpy", "scipy", "cma", "ctREFPROP")},
                    baseline_COP=float(baseline["best_COP"]), baseline_signature=baseline["case_signature"],
                    baseline_evaluations=int(baseline["exact_evaluations"]),
                    initial_design="origin then four local probes" if method != "cma_es" else
                                   "origin then CMA-ES, sigma=0.25, population=7",
                    timeout_seconds=180, timeout_retries=1, worker_max_tasks=20)
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    key = hashlib.sha256(ref.encode()).hexdigest()[:10]
    path = output / f"{key}_{method}_{seed}.json"
    result = json.loads(path.read_text()) if path.exists() else dict(metadata=metadata, records=[])
    if result["metadata"] != metadata:
        raise ValueError(f"Checkpoint settings or code differ: {path}")
    if result.get("replay", {}).get("success"):
        return str(path)
    options = dict(refrigerant=ref, base_case=base_case, exact_evaluation_timeout_seconds=180,
                   exact_evaluation_timeout_retries=1, exact_evaluation_worker_max_tasks=20)
    objective = _OuterDesignObjective(max_evals=budget-len(result["records"]), **options)
    search = make_search(method, objective.max_recuperator_lift, seed)

    def checkpoint():
        if method == "cma_es":
            result["restarts"] = search.restart_reasons
            result["rejected_parameter_candidates"] = search.rejected_candidates
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False))
        temporary.replace(path)

    try:
        for previous in result["records"]:
            proposal = search.ask()
            if not np.allclose(proposal, previous["x"][:3], atol=1e-7, rtol=0):
                raise ValueError("Checkpoint candidate does not reproduce")
            search.tell(proposal, previous["SEC"] or 1e6, previous["success"])
        while len(result["records"]) < budget:
            start = time.perf_counter()
            proposal = search.ask()
            proposal_seconds = time.perf_counter() - start
            start = time.perf_counter()
            before = objective.exact_evaluations
            record = objective.evaluate([*proposal, *objective.initial_mass_fractions])
            evaluation_wall_seconds = time.perf_counter() - start
            if objective.exact_evaluations != before + 1:
                raise AssertionError("Duplicate or uncounted exact evaluation")
            if not np.allclose(record.x[:3], proposal, atol=1e-7, rtol=0):
                raise AssertionError("Candidate unexpectedly needed domain repair")
            if not np.allclose(record.x[3:], objective.initial_mass_fractions, atol=1e-12, rtol=0):
                raise AssertionError("Composition moved")
            search.tell(proposal, record.objective, record.success)
            item = serialise(record)
            item.update(proposal_seconds=proposal_seconds, evaluation_wall_seconds=evaluation_wall_seconds)
            result["records"].append(item)
            checkpoint()
            n = len(result["records"])
            if n == 1 or n % 25 == 0:
                best = max((r["COP"] for r in result["records"] if r["success"]), default=0)
                print(ref, method, seed, n, f"best COP {best:.7f}", flush=True)
    finally:
        objective.close()
    successes = [r for r in result["records"] if r["success"]]
    if successes:
        best = min(successes, key=lambda r: r["SEC"])
        replay = _OuterDesignObjective(max_evals=1, **dict(options, exact_evaluation_worker_max_tasks=1))
        try:
            result["replay"] = serialise(replay.evaluate(best["x"]))
        finally:
            replay.close()
    else:
        result["replay"] = dict(success=False, message="No successful design")
    checkpoint()
    return str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=int, default=200)
    parser.add_argument("--n-cpu", type=int, default=6)
    parser.add_argument("--seed", type=int, action="append")
    parser.add_argument("--method", choices=METHODS, action="append")
    parser.add_argument("--refrigerant", action="append")
    parser.add_argument("--source", default=DEFAULT_SOURCE)
    parser.add_argument("--output", default=OUTPUT)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    if args.budget < 5:
        parser.error("budget must include at least the five local bootstrap observations")
    tasks = [(ref, method, seed, args.budget, args.source, args.output)
             for seed in (args.seed or (20260909, 20260910, 20260911))
             for ref in (args.refrigerant or DEFAULT_CASES)
             for method in (args.method or METHODS)]
    tasks = list(dict.fromkeys(tasks))
    if args.report_only:
        from HeatPumpLab.cycles.RecuperativeHeatPump.utils.origin_benchmark_report import write_origin_report
        print(write_origin_report(args.output, tasks))
        return
    with ProcessPoolExecutor(max_workers=max(1, args.n_cpu)) as pool:
        futures = [pool.submit(run_trial, task) for task in tasks]
        for future in as_completed(futures):
            print("Completed", future.result(), flush=True)


if __name__ == "__main__":
    main()
