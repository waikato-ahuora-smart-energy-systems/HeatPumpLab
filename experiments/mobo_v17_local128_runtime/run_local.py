"""Run resumable COP/total-VHC multi-objective BO for four binary families."""

import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import torch
from local_search import LocalParetoSearch
from main_RHPmobo import metadata_for

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
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


REFRESH_EVERY = 48


def restore_search(result):
    phase = result['local128_phase']
    start = phase['start_evaluations']
    epoch = start + ((len(result['records']) - start) // REFRESH_EVERY) * REFRESH_EVERY
    baseline = dict(records=result['records'][:epoch], kernel_state=result.get('kernel_state'),
                    mean_state=result.get('mean_state'))
    search = LocalParetoSearch(baseline, 128, seed=phase['seed'] + 104729 * (epoch-start))
    for r in result['records'][epoch:]:
        search.tell(r['x'], r['COP'], r['VHC_MJ_per_m3'], r['success'] and not r.get('validation_rejected'))
    saved = result.get('local128_state', {})
    if saved.get('epoch_start') == epoch:
        search.states = {int(k): v for k, v in saved['model_states'].items()}
    return search, epoch


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
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    phase = result.setdefault('local128_phase', dict(
        start_evaluations=len(result['records']), subset_size=128, seed=20260912,
        radius=.05, anchor_count=6, refresh_every=REFRESH_EVERY,
        algorithm='nearest-success local qLogEHVI; full observed Pareto baseline',
        runtime_source_sha256=hashlib.sha256(b''.join(p.read_bytes() for p in sorted(ROOT.glob('*.py')))).hexdigest()))
    expected = hashlib.sha256(b''.join(p.read_bytes() for p in sorted(ROOT.glob('*.py')))).hexdigest()
    assert phase['runtime_source_sha256'] == expected, 'Local runtime changed since checkpoint'

    options = dict(refrigerant=metadata["refrigerant"], base_case=metadata["base_case"],
                   exact_evaluation_timeout_seconds=180, exact_evaluation_timeout_retries=0,
                   exact_evaluation_worker_max_tasks=1, results_dir=str(Path(output)/"oracle"))
    objective = _OuterDesignObjective(max_evals=max(1, budget-len(result["records"])), **options)
    search, epoch = restore_search(result)
    search.max_lift = objective.max_recuperator_lift

    def checkpoint():
        result['local128_state'] = dict(epoch_start=epoch, model_states=search.states,
                                        anchor_evaluations=[i+1 for i in search.anchor_indices])
        result["frontier_indices"] = frontier_indices(result["records"])
        temp = path.with_suffix(".tmp")
        temp.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n")
        temp.replace(path)

    checkpoint()
    try:
        while len(result["records"]) < budget:
            target_epoch = phase['start_evaluations'] + ((len(result['records'])-phase['start_evaluations']) // REFRESH_EVERY) * REFRESH_EVERY
            if target_epoch != epoch:
                search, epoch = restore_search(result)
                search.max_lift = objective.max_recuperator_lift
            started = time.perf_counter()
            proposal = search.ask()
            proposal_seconds = time.perf_counter()-started
            before = objective.exact_evaluations
            started = time.perf_counter()
            record = objective.evaluate(proposal)
            assert objective.exact_evaluations == before+1
            assert np.allclose(record.x, proposal, rtol=0, atol=1e-7)
            item = serialise(record)
            if "worker_failure" in item["constraint_errors"]:
                raise RuntimeError(f"Infrastructure failure; search paused without training on this candidate: {item['message']}")
            if item["success"] and (not np.isfinite([item["COP"], item["VHC_MJ_per_m3"]]).all()
                                    or min(item["COP"], item["VHC_MJ_per_m3"]) <= 0):
                raise ValueError("A successful evaluation lacks paired positive COP/VHC")
            item.update(proposal_seconds=proposal_seconds, evaluation_wall_seconds=time.perf_counter()-started,
                        search_method="local128", search_diagnostic=search.last_diagnostic)
            search.tell(item["x"], item["COP"], item["VHC_MJ_per_m3"], item["success"])
            result["records"].append(item)
            checkpoint()
            count = len(result["records"])
            if count == phase["start_evaluations"]+1 or count % 10 == 0:
                print(partner, count, "front", len(result["frontier_indices"]),
                      "failures", sum(not r["success"] for r in result["records"]),
                      "proposal_seconds", round(proposal_seconds, 3), "training_points", search.last_diagnostic.get("training_points"), flush=True)
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
    parser.add_argument('--partner', choices=PARTNERS, required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    print('Completed', run_family((args.partner, 4000, 20260910, args.output)), flush=True)


if __name__ == '__main__':
    main()
