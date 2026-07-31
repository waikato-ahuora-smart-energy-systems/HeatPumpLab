"""Probe cold repeatability of one outer-design vector in isolated processes."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))


def _solve(index: int):
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import (
        _OuterDesignObjective,
    )

    x = (0.0, 55.28599659661655, 0.45, 0.225, 0.775)
    objective = _OuterDesignObjective(
        max_evals=1,
        base_case={"T_air_tar": 180},
        refrigerant=("Propane", "Isopentane"),
        results_dir=f"/private/tmp/hpl-repeat-{index}",
    )
    record = objective.evaluate(x)
    report_keys = (
        "SEC",
        "COP",
        "Pmax",
        "Tref_evap_outlet",
        "dT_gc_out_min_offset",
        "converged",
    )
    return (
        index,
        record.success,
        record.sec,
        record.pmax,
        record.message,
        record.constraint_errors,
        {key: record.report.get(key) for key in report_keys},
    )


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=6) as pool:
        for answer in pool.map(_solve, range(12)):
            print(answer)
