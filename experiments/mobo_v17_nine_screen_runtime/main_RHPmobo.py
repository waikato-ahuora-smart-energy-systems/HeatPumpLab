"""Run resumable COP/total-VHC multi-objective BO for four binary families."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
from importlib.metadata import version
import json
from pathlib import Path
import time

import numpy as np

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import (
    _OuterDesignObjective,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
    DEFAULT_BASE_CASE,
    OPTIMISER_VERSION,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.utils.composition_pareto import (
    nondominated_mask,
)
from main_RHPoptimiser_benchmark import serialise as basic_serialise
import torch

ROOT = Path(__file__).resolve().parent
PAIRS = {
    "Propane_Cyclopentane": ("Propane", "Cyclopentane"),
    "Butane_Pentane": ("Butane", "Pentane"),
    "Isobutane_Pentane": ("Isobutane", "Pentane"),
    "Butane_Isopentane": ("Butane", "Isopentane"),
    "Isobutane_Isopentane": ("Isobutane", "Isopentane"),
    "Butane_Cyclopentane": ("Butane", "Cyclopentane"),
    "Isobutane_Cyclopentane": ("Isobutane", "Cyclopentane"),
    "Butane_Isobutane": ("Butane", "Isobutane"),
    "Pentane_Isopentane": ("Pentane", "Isopentane"),
}
PARTNERS = tuple(PAIRS)
OUTPUT = str(ROOT.parents[1] / "results/RHPcomposition/mobo_v17_nine_screen_800")


def serialise(record):
    item = basic_serialise(record)
    report = record.report
    names = (
        "Comp_1_T_out",
        "Comp_2_T_out",
        "T_comp_discharge",
        "T_gc_out",
        "Tref_evap_outlet",
        "PR_1",
        "PR_2",
        "Overall_PR",
    )
    for name in names:
        value = report.get(name)
        item[name] = float(value) if value is not None and np.isfinite(value) else None
    temperatures = [item[name] for name in names[:5] if item[name] is not None]
    item["max_reported_refrigerant_C"] = max(temperatures) if temperatures else None
    discharge = [item[name] for name in names[:3] if item[name] is not None]
    item["max_compressor_discharge_C"] = max(discharge) if discharge else None
    if record.success and any(
        item[name] is None for name in ("Comp_1_T_out", "Comp_2_T_out", "PR_1", "PR_2")
    ):
        raise ValueError(
            "Successful design is missing compressor temperature or pressure-ratio diagnostics"
        )
    return item


def frontier_indices(records):
    values = [
        (r["COP"], r["VHC_MJ_per_m3"])
        if r["success"] and not r.get("validation_rejected")
        else (float("nan"), float("nan"))
        for r in records
    ]
    return [i for i, keep in enumerate(nondominated_mask(values)) if keep]


def metadata_for(partner, budget, seed):
    return dict(
        partner=partner,
        components=list(PAIRS[partner]),
        property_model_status=next(
            r
            for r in json.loads((ROOT / "preflight.json").read_text())
            if r["pair"] == partner
        ),
        refrigerant=";".join(PAIRS[partner]) + "|0.5;0.5|mass",
        budget=budget,
        seed=seed,
        base_case=DEFAULT_BASE_CASE,
        evaluation_version=OPTIMISER_VERSION,
        algorithm="global BoTorch qLogEHVI, q=1, distance-weighted feasibility; every eighth adaptive step maximin",
        scoring_batch=8,
        restart_batch=1,
        init_batch=8,
        objectives=["COP", "VHC_total_MJ_per_m3"],
        reference_point=[0.0, 0.0],
        variables=[
            "RHX1_K",
            "RHX2_K",
            "compressor_2_power_fraction",
            "first_component_mass_fraction",
        ],
        bounds=[[0.0, 100.0], [0.0, 100.0], [0.1, 0.9], [0.0, 1.0]],
        timeout_seconds=180,
        timeout_retries=0,
        worker_max_tasks=1,
        model_source_sha256=hashlib.sha256(
            b"".join(
                p.read_bytes() for p in sorted((ROOT / "HeatPumpLab").rglob("*.py"))
            )
        ).hexdigest(),
        search_source_sha256=hashlib.sha256(
            b"".join(
                (ROOT / p).read_bytes()
                for p in (
                    "rhp_mobo_search.py",
                    "main_RHPmobo.py",
                    "main_RHPoptimiser_benchmark.py",
                )
            )
        ).hexdigest(),
        versions={
            name: version(name)
            for name in ("numpy", "scipy", "botorch", "torch", "gpytorch", "ctREFPROP")
        },
    )


def run_family(task):
    partner, budget, seed, output = task
    path = Path(output) / f"{partner}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    metadata = metadata_for(partner, budget, seed)
    result = (
        json.loads(path.read_text())
        if path.exists()
        else dict(metadata=metadata, records=[])
    )
    if result["metadata"] != metadata:
        raise ValueError(f"Checkpoint source or settings mismatch: {path}")
    if result.get("complete"):
        return str(path)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    from rhp_mobo_search import ParetoBOSearch

    options = dict(
        refrigerant=metadata["refrigerant"],
        base_case=metadata["base_case"],
        exact_evaluation_timeout_seconds=180,
        exact_evaluation_timeout_retries=0,
        exact_evaluation_worker_max_tasks=1,
        results_dir=str(Path(output) / "oracle"),
    )
    objective = _OuterDesignObjective(
        max_evals=max(1, budget - len(result["records"])), **options
    )
    search = ParetoBOSearch(objective.max_recuperator_lift, seed)
    for record in result["records"]:
        search.tell(
            record["x"],
            record["COP"],
            record["VHC_MJ_per_m3"],
            record["success"] and not record.get("validation_rejected"),
        )
    search.kernel_state = result.get("kernel_state")
    search.mean_state = result.get("mean_state")
    search.diagnostics = result.get("search_diagnostics", [])

    def checkpoint():
        result["kernel_state"] = search.kernel_state
        result["mean_state"] = search.mean_state
        result["search_diagnostics"] = search.diagnostics
        result["frontier_indices"] = frontier_indices(result["records"])
        temp = path.with_suffix(".tmp")
        temp.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        temp.replace(path)

    try:
        while len(result["records"]) < budget:
            started = time.perf_counter()
            proposal = search.ask()
            proposal_seconds = time.perf_counter() - started
            before = objective.exact_evaluations
            started = time.perf_counter()
            record = objective.evaluate(proposal)
            assert objective.exact_evaluations == before + 1
            assert np.allclose(record.x, proposal, rtol=0, atol=1e-7)
            item = serialise(record)
            if "worker_failure" in item["constraint_errors"]:
                raise RuntimeError(
                    f"Infrastructure failure; search paused without training on this candidate: {item['message']}"
                )
            if item["success"] and (
                not np.isfinite([item["COP"], item["VHC_MJ_per_m3"]]).all()
                or min(item["COP"], item["VHC_MJ_per_m3"]) <= 0
            ):
                raise ValueError(
                    "A successful evaluation lacks paired positive COP/VHC"
                )
            item.update(
                proposal_seconds=proposal_seconds,
                evaluation_wall_seconds=time.perf_counter() - started,
            )
            search.tell(item["x"], item["COP"], item["VHC_MJ_per_m3"], item["success"])
            result["records"].append(item)
            checkpoint()
            count = len(result["records"])
            if count == 1 or count % 10 == 0:
                print(
                    partner,
                    count,
                    "front",
                    len(result["frontier_indices"]),
                    "failures",
                    sum(not r["success"] for r in result["records"]),
                    flush=True,
                )
    finally:
        objective.close()
    # Validate the complete observed frontier with a fresh process per point.
    result.setdefault("frontier_replays", {})
    replay = _OuterDesignObjective(
        max_evals=budget, **dict(options, exact_evaluation_worker_max_tasks=1)
    )
    try:
        while True:
            pending = [
                i
                for i in result["frontier_indices"]
                if str(i) not in result["frontier_replays"]
            ]
            if not pending:
                break
            for index in pending:
                checked = serialise(replay.evaluate(result["records"][index]["x"]))
                result["frontier_replays"][str(index)] = checked
                original = result["records"][index]
                matches = checked["success"] and all(
                    np.isclose(checked[key], original[key], rtol=1e-6, atol=1e-9)
                    for key in (
                        "COP",
                        "VHC_MJ_per_m3",
                        "SEC",
                        "Pmax_MPa",
                        "max_reported_refrigerant_C",
                        "PR_1",
                        "PR_2",
                    )
                )
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
    parser.add_argument("--budget", type=int, default=800)
    parser.add_argument("--seed", type=int, default=20260912)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--partner", choices=PARTNERS, action="append")
    parser.add_argument("--output", default=OUTPUT)
    args = parser.parse_args()
    if args.budget < 33:
        parser.error("At least 33 evaluations are required to include a BO proposal")
    if args.workers == 1:
        for partner in dict.fromkeys(args.partner or PARTNERS):
            print(
                "Completed",
                run_family((partner, args.budget, args.seed, args.output)),
                flush=True,
            )
        return
    errors = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as pool:
        futures = [
            pool.submit(run_family, (partner, args.budget, args.seed, args.output))
            for partner in dict.fromkeys(args.partner or PARTNERS)
        ]
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
