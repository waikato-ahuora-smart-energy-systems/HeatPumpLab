"""One independent native integration probe per binary pair, outside search budgets."""

import json
from pathlib import Path
from main_RHPmobo import PAIRS, serialise, DEFAULT_BASE_CASE
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import (
    _OuterDesignObjective,
)


def main():
    records = []
    for key, pair in PAIRS.items():
        obj = _OuterDesignObjective(
            refrigerant=";".join(pair) + "|0.5;0.5|mass",
            base_case=DEFAULT_BASE_CASE,
            max_evals=1,
            exact_evaluation_timeout_seconds=180,
            exact_evaluation_timeout_retries=0,
            exact_evaluation_worker_max_tasks=1,
            results_dir="/tmp/heatpumplab-nine-screen-smoke",
        )
        try:
            record = serialise(obj.evaluate([0.0, 20.0, 0.3, 0.5, 0.5]))
            if "worker_failure" in record["constraint_errors"]:
                raise RuntimeError(record["message"])
            records.append(dict(pair=key, record=record))
            print(
                key,
                record["success"],
                record["COP"],
                record["max_reported_refrigerant_C"],
                record["message"],
                flush=True,
            )
        finally:
            obj.close()
        (Path(__file__).resolve().parent / "smoke_results.json").write_text(
            json.dumps(records, indent=2) + "\n"
        )


if __name__ == "__main__":
    main()
