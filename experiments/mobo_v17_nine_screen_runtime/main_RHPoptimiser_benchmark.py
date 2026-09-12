"""Cold-start RBF / GP-EI / maximin pilot on fixed-composition RHP cases.

Run from the repository root. Every method receives the same initial design,
budget, model settings and seed. Existing v11 optima are comparators only.
Each evaluation is checkpointed; completed runs can be safely resumed.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
from pathlib import Path
import time

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import DEFAULT_BASE_CASE, OPTIMISER_VERSION
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation.surrogate_search import JointSurrogateSearch


DEFAULT_CASES = (
    "Propane;Pentane|0.3;0.7|mass",
    "Propane;Pentane|0.6;0.4|mass",
    "Propane;Isobutane|0.3;0.7|mass",
)
DEFAULT_SOURCE = "results/RHPcomposition/RHP_fixed_mass_composition_v11.csv"
DEFAULT_OUTPUT = "results/RHPcomposition/optimiser_investigation"


def load_case(refrigerant, source):
    with Path(source).open() as stream:
        rows = [row for row in csv.DictReader(stream)
                if row["input_refrigerant"] == refrigerant and float(row["T_air_tar"]) == 200]
    if len(rows) != 1:
        raise ValueError(f"Expected one 200 C baseline for {refrigerant}")
    row = rows[0]
    base_case = {
        key: (str(row[key]).lower() == "true" if isinstance(default, bool)
              else type(default)(row[key]))
        for key, default in DEFAULT_BASE_CASE.items()
    }
    return row, base_case


def serialise(record):
    return dict(x=list(record.x), success=record.success, SEC=record.sec,
                COP=record.cop, VHC_MJ_per_m3=record.report.get("VHC", 0) / 1000,
                Pmax_MPa=record.pmax / 1000 if record.pmax is not None else None,
                elapsed_seconds=record.elapsed_seconds, message=record.message,
                constraint_errors=record.constraint_errors,
                solve_failure_kind=record.report.get("solve_failure_kind"),
                solve_error_stage=record.report.get("solve_error_stage"),
                property_failure_context=record.report.get("property_failure_context"),
                numerical_retries=record.report.get("numerical_retries", []))


def run_trial(task):
    refrigerant, method, budget, seed, source, output = task
    baseline, base_case = load_case(refrigerant, source)
    objective_options = dict(
        refrigerant=refrigerant, base_case=base_case,
        exact_evaluation_timeout_seconds=180,
        exact_evaluation_timeout_retries=1,
        exact_evaluation_worker_max_tasks=20,
    )
    metadata = dict(version=1, refrigerant=refrigerant, method=method,
                    evaluation_version=OPTIMISER_VERSION,
                    budget=budget, seed=seed, base_case=base_case,
                    objective_options=objective_options,
                    baseline_COP=float(baseline["best_COP"]),
                    baseline_evaluations=int(baseline["exact_evaluations"]),
                    baseline_signature=baseline["case_signature"])
    key = hashlib.sha256(refrigerant.encode()).hexdigest()[:10]
    path = Path(output) / f"{key}_{method}_{seed}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    saved = json.loads(path.read_text()) if path.exists() else None
    if saved and saved["metadata"] != metadata:
        raise ValueError(f"Checkpoint settings differ: {path}")
    if saved and saved.get("replay"):
        return str(path)
    records = saved["records"] if saved else []
    objective = _OuterDesignObjective(max_evals=budget, **objective_options)
    search = JointSurrogateSearch(method=method, seed=seed, max_lift=objective.max_recuperator_lift)
    fractions = objective.initial_mass_fractions
    result = dict(metadata=metadata, records=records)

    def checkpoint():
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False))
        temporary.replace(path)

    try:
        for previous in records:
            # Replay ask calls to restore the exact random-stream state.
            proposed = search.ask()
            if max(abs(proposed[i] - previous["x"][i]) for i in range(3)) > 1e-7:
                raise ValueError("Checkpoint proposal does not reproduce")
            search.tell(previous["x"][:3], previous["SEC"] or 1e6, previous["success"])
        while len(records) < budget:
            start = time.perf_counter()
            proposed = search.ask()
            proposal_seconds = time.perf_counter() - start
            record = objective.evaluate([*proposed, *fractions])
            if max(abs(record.x[i] - proposed[i]) for i in range(3)) > 1e-7:
                raise AssertionError("Candidate unexpectedly needed repair")
            if max(abs(a-b) for a, b in zip(record.x[3:], fractions)) > 1e-10:
                raise AssertionError("Composition changed during search")
            search.tell(record.x[:3], record.objective, record.success)
            item = serialise(record)
            item["proposal_seconds"] = proposal_seconds
            records.append(item)
            checkpoint()
            if len(records) % 8 == 0:
                best = max((r["COP"] for r in records if r["success"]), default=0)
                print(refrigerant, method, len(records), "best COP", best, flush=True)
    finally:
        objective.close()
    successful = [r for r in records if r["success"]]
    if successful:
        best = min(successful, key=lambda r: r["SEC"])
        options = dict(objective_options, exact_evaluation_worker_max_tasks=1)
        replay = _OuterDesignObjective(max_evals=1, **options)
        try:
            result["replay"] = serialise(replay.evaluate(best["x"]))
        finally:
            replay.close()
    else:
        result["replay"] = dict(success=False, message="No feasible design")
    checkpoint()
    return str(path)


def run_joint_polish(task):
    """Diagnose frozen mixture lifts using the existing local search, warm-started."""
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._polish import _run_pattern_search_polish
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import LocalPolishConfig, PhaseEvaluation

    refrigerant, _, budget, seed, source, output = task
    baseline, base_case = load_case(refrigerant, source)
    options = dict(refrigerant=refrigerant, base_case=base_case,
                   exact_evaluation_timeout_seconds=180,
                   exact_evaluation_timeout_retries=1,
                   exact_evaluation_worker_max_tasks=20)
    key = hashlib.sha256(refrigerant.encode()).hexdigest()[:10]
    path = Path(output) / f"{key}_joint_polish_{seed}.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    result = dict(metadata=dict(
        version=1, method="joint_polish", refrigerant=refrigerant, budget=budget,
        evaluation_version=OPTIMISER_VERSION,
        seed=seed, warm_start=True, base_case=base_case, objective_options=options,
        baseline_COP=float(baseline["best_COP"]),
        baseline_evaluations=int(baseline["exact_evaluations"]),
        baseline_signature=baseline["case_signature"],
    ), records=[])
    if path.exists():
        saved = json.loads(path.read_text())
        if saved["metadata"] != result["metadata"]:
            raise ValueError(f"Checkpoint settings differ: {path}")
        if saved.get("replay"):
            return str(path)
        raise ValueError(f"Incomplete local polish cannot resume: {path}")
    objective = _OuterDesignObjective(max_evals=budget, **options)

    def checkpoint():
        temporary = path.with_suffix(".tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False))
        temporary.replace(path)

    def evaluate_phase(phase, x):
        record = objective.evaluate(x)
        result["records"] = [serialise(r) for r in objective.records]
        checkpoint()
        return PhaseEvaluation.from_record(
            phase=phase, record=record,
            cumulative_exact_evaluations=objective.exact_evaluations,
        )

    try:
        x = [float(baseline[f"best_{k}"]) for k in
             ("dT_rhx1_cold", "dT_rhx2_cold", "comp2_power_frac")]
        start = evaluate_phase("baseline", [*x, *objective.initial_mass_fractions])
        polish = _run_pattern_search_polish(
            source=start, objective=objective,
            config=LocalPolishConfig(budget=budget-1, vary_dT_rhx1_cold=True,
                                     vary_dT_rhx2_cold=True, vary_composition=False),
            evaluate_phase=evaluate_phase,
        )
        result["polish_status"] = polish.status
        result["polish_message"] = polish.message
        best = min((r for r in objective.records if r.success), key=lambda r: r.sec)
    finally:
        objective.close()
    replay = _OuterDesignObjective(max_evals=1, **dict(options, exact_evaluation_worker_max_tasks=1))
    try:
        result["replay"] = serialise(replay.evaluate(best.x))
    finally:
        replay.close()
    checkpoint()
    return str(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--budget", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--n-cpu", type=int, default=6)
    parser.add_argument("--source", default=DEFAULT_SOURCE)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--refrigerant", action="append")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--method", action="append", choices=(*JointSurrogateSearch.METHODS, "joint_polish"))
    args = parser.parse_args()
    if args.report_only:
        from HeatPumpLab.cycles.RecuperativeHeatPump.utils.optimiser_benchmark_report import write_benchmark_report
        print(write_benchmark_report(args.output))
        return
    if args.budget < 24:
        parser.error("budget must cover the 24-point initial design")
    tasks = [(ref, method, args.budget, args.seed, args.source, args.output)
             for ref in (args.refrigerant or DEFAULT_CASES)
             for method in (args.method or JointSurrogateSearch.METHODS)]
    with ProcessPoolExecutor(max_workers=max(1, args.n_cpu)) as pool:
        futures = [pool.submit(run_joint_polish if task[1] == "joint_polish" else run_trial, task)
                   for task in tasks]
        for future in as_completed(futures):
            print("Completed", future.result(), flush=True)


if __name__ == "__main__":
    main()
