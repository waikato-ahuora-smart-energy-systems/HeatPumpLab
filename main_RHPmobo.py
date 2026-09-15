"""Run resumable COP/total-VHC multi-objective BO for four binary families."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import time

import numpy as np

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
    DEFAULT_BASE_CASE, DEFAULT_MAX_STAGE_PRESSURE_RATIO, OPTIMISER_VERSION)
from HeatPumpLab.cycles.RecuperativeHeatPump.utils.composition_pareto import nondominated_mask
from main_RHPoptimiser_benchmark import serialise

ROOT = Path(__file__).resolve().parent
PARTNERS = ("Butane", "Isopentane", "Pentane", "Isobutane")
OUTPUT = "results/RHPcomposition/mobo_v17_400"


def frontier_indices(records):
    values = [(r["COP"], r["VHC_MJ_per_m3"]) if r["success"] and not r.get("validation_rejected")
              else (float("nan"), float("nan"))
              for r in records]
    return [i for i, keep in enumerate(nondominated_mask(values)) if keep]


def metadata_for(partner, budget, seed):
    return dict(partner=partner, refrigerant=f"Propane;{partner}|0.5;0.5|mass", budget=budget,
                seed=seed, base_case=DEFAULT_BASE_CASE, evaluation_version=OPTIMISER_VERSION,
                max_stage_pressure_ratio=DEFAULT_MAX_STAGE_PRESSURE_RATIO,
                compressor_gas_only=True,
                algorithm="BoTorch qLogEHVI, q=1, distance-weighted feasibility; every eighth adaptive step maximin",
                objectives=["COP", "VHC_total_MJ_per_m3"], reference_point=[0., 0.],
                variables=["RHX1_K", "RHX2_K", "compressor_2_power_fraction", "propane_mass_fraction"],
                bounds=[[0., 100.], [0., 100.], [.1, .9], [0., 1.]],
                timeout_seconds=180, timeout_retries=0, worker_max_tasks=20,
                model_source_sha256=hashlib.sha256(b"".join(p.read_bytes() for p in
                    sorted((ROOT/"HeatPumpLab").rglob("*.py")))).hexdigest(),
                search_source_sha256=hashlib.sha256(b"".join((ROOT/p).read_bytes() for p in
                    ("rhp_mobo_search.py", "main_RHPmobo.py", "main_RHPoptimiser_benchmark.py"))).hexdigest(),
                versions={name: version(name) for name in ("numpy", "scipy", "botorch", "torch", "gpytorch", "ctREFPROP")})


def run_family(task):
    partner, budget, seed, output = task
    path = Path(output)/f"Propane_{partner}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = metadata_for(partner, budget, seed)
    result = json.loads(path.read_text()) if path.exists() else dict(metadata=metadata, records=[])
    if result["metadata"] != metadata:
        raise ValueError(f"Checkpoint source or settings mismatch: {path}")
    if result.get("complete"):
        return str(path)
    from rhp_mobo_search import ParetoBOSearch
    options = dict(refrigerant=metadata["refrigerant"], base_case=metadata["base_case"],
                   exact_evaluation_timeout_seconds=180, exact_evaluation_timeout_retries=0,
                   exact_evaluation_worker_max_tasks=20, results_dir=str(Path(output)/"oracle"))
    objective = _OuterDesignObjective(max_evals=max(1, budget-len(result["records"])), **options)
    search = ParetoBOSearch(objective.max_recuperator_lift, seed)
    for record in result["records"]:
        search.tell(record["x"], record["COP"], record["VHC_MJ_per_m3"], record["success"])
    search.kernel_state = result.get("kernel_state")
    search.mean_state = result.get("mean_state")
    search.diagnostics = result.get("search_diagnostics", [])

    def checkpoint():
        result["kernel_state"] = search.kernel_state
        result["mean_state"] = search.mean_state
        result["search_diagnostics"] = search.diagnostics
        result["frontier_indices"] = frontier_indices(result["records"])
        temp = path.with_suffix(".tmp")
        temp.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        temp.replace(path)

    try:
        while len(result["records"]) < budget:
            started = time.perf_counter()
            proposal = search.ask()
            proposal_seconds = time.perf_counter()-started
            before = objective.exact_evaluations
            started = time.perf_counter()
            record = objective.evaluate(proposal)
            assert objective.exact_evaluations == before+1
            assert np.allclose(record.x, proposal, rtol=0, atol=1e-7)
            item = serialise(record)
            if item["success"] and (not np.isfinite([item["COP"], item["VHC_MJ_per_m3"]]).all()
                                    or min(item["COP"], item["VHC_MJ_per_m3"]) <= 0):
                raise ValueError("A successful evaluation lacks paired positive COP/VHC")
            item.update(proposal_seconds=proposal_seconds, evaluation_wall_seconds=time.perf_counter()-started)
            search.tell(item["x"], item["COP"], item["VHC_MJ_per_m3"], item["success"])
            result["records"].append(item)
            checkpoint()
            count = len(result["records"])
            if count == 1 or count % 10 == 0:
                print(partner, count, "front", len(result["frontier_indices"]),
                      "failures", sum(not r["success"] for r in result["records"]), flush=True)
    finally:
        objective.close()
    # Validate the complete observed frontier with a fresh process per point.
    result.setdefault("frontier_replays", {})
    replay = _OuterDesignObjective(max_evals=budget,
                                   **dict(options, exact_evaluation_worker_max_tasks=1))
    try:
        while True:
            pending = [i for i in result["frontier_indices"] if str(i) not in result["frontier_replays"]]
            if not pending:
                break
            for index in pending:
                checked = serialise(replay.evaluate(result["records"][index]["x"]))
                result["frontier_replays"][str(index)] = checked
                original = result["records"][index]
                matches = checked["success"] and all(np.isclose(checked[key], original[key], rtol=1e-6, atol=1e-9)
                                                     for key in ("COP", "VHC_MJ_per_m3", "SEC", "Pmax_MPa"))
                if not matches:
                    original["validation_rejected"] = True
                checkpoint()
    finally:
        replay.close()
    result["complete"] = True
    checkpoint()
    return str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=int, default=400)
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--partner", choices=PARTNERS, action="append")
    parser.add_argument("--output", default=OUTPUT)
    args = parser.parse_args()
    if args.budget < 33:
        parser.error("At least 33 evaluations are required to include a BO proposal")
    errors = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = [pool.submit(run_family, (partner, args.budget, args.seed, args.output))
                   for partner in dict.fromkeys(args.partner or PARTNERS)]
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
