"""Diagnostic-only, fresh-process replays of frozen BoTorch benchmark failures."""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import faulthandler
import hashlib
import inspect
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "results/RHPcomposition/botorch_v14_400"
OUTPUT = ROOT / "results/RHPcomposition/botorch_v14_failures"


def clean(value):
    import numpy as np
    if isinstance(value, dict):
        return {str(k): clean(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [clean(v) for v in value]
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return repr(value)


def snapshot(hp):
    fields = ("feasible", "converged", "solve_error", "constraint_errors", "numerical_retries",
              "COP", "SEC", "tol", "P_H_est", "Tref_evap_outlet", "gc_minimum_approach",
              "dT_gc_out_min_offset", "pressure_interpolation_fallbacks")
    data = {name: getattr(hp, name, None) for name in fields}
    data["pressure"] = getattr(getattr(hp, "comp_discharge", None), "pressure", None)
    data["total_heat_duty"] = sum(gc.duty for gc in getattr(hp, "gas_coolers", []))
    data["compressor_suction_volume_flows"] = [getattr(c.flow_in, "vol_flow", None)
                                               for c in getattr(hp, "compressors", [])]
    data["gas_coolers"] = [dict(duty=gc.duty, sink_mass_flow=gc.sink_in.mass_flow,
                               source_in_T=gc.source_in.temperature, source_out_T=gc.source_out.temperature,
                               sink_in_T=gc.sink_in.temperature, sink_out_T=gc.sink_out.temperature,
                               sink_in_h=gc.sink_in.enthalpy, sink_out_h=gc.sink_out.enthalpy)
                           for gc in getattr(hp, "gas_coolers", [])]
    data["unconverged_operations"] = []
    for operation in getattr(hp, "operations", []):
        if getattr(operation, "converged", False):
            continue
        details = dict(type=type(operation).__name__, name=getattr(operation, "name", None),
                       duty=getattr(operation, "duty", None), tol=getattr(operation, "tol", None))
        for port in ("source_in", "source_out", "sink_in", "sink_out"):
            flow = getattr(operation, port, None)
            if flow is not None:
                details[port] = {key: getattr(flow, key, None) for key in ("temperature", "enthalpy", "mass_flow", "pressure")}
        data["unconverged_operations"].append(details)
    return clean(data)


def worker(request_path):
    request = json.loads(Path(request_path).read_text())
    prefix = Path(request["prefix"])
    faulthandler.enable()
    faulthandler.dump_traceback_later(60, repeat=True)
    from main_RHPoptimiser_benchmark import serialise
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
    from HeatPumpLab.cycles.RecuperativeHeatPump.solver import cycle_solver as cycle
    from HeatPumpLab.cycles.RecuperativeHeatPump.solver import _discharge_pressure as pressure
    from HeatPumpLab.utils.solve_profile import profile_solve
    stream = prefix.with_suffix(".trace.jsonl").open("w")
    started = time.perf_counter()

    def emit(event, **data):
        stream.write(json.dumps(clean(dict(event=event, seconds=time.perf_counter()-started, **data))) + "\n")
        stream.flush()

    if request.get("trace_flash"):
        from HeatPumpLab.utils.properties_package import PropertyPackage
        original_flash = PropertyPackage._flash

        def flash(self, pair, first, second, comp):
            emit("flash_start", fluid=self.fld_name, pair=pair, first=first, second=second, composition=comp)
            if self.fld_name == "Air" and pair == "PH" and abs(second) > 1e6:
                frame = inspect.currentframe().f_back
                while frame is not None:
                    if "hp" in frame.f_locals:
                        emit("extreme_air_enthalpy", state=snapshot(frame.f_locals["hp"]))
                        break
                    frame = frame.f_back
                del frame
            try:
                result = original_flash(self, pair, first, second, comp)
            except Exception as exc:
                emit("flash_exception", fluid=self.fld_name, pair=pair, message=str(exc))
                raise
            emit("flash_return", fluid=self.fld_name, pair=pair, ierr=result.ierr)
            return result

        PropertyPackage._flash = flash
        original_recovery = PropertyPackage._recover_two_phase_ph

        def recovery(self, pressure, enthalpy, comp):
            def trace(frame, event, arg):
                if frame.f_code is not original_recovery.__code__:
                    return None
                if event == "return":
                    values = {k: v for k, v in frame.f_locals.items()
                              if k in ("bubble", "dew", "root", "state", "errors", "liquid", "vapor", "fl", "fv")}
                    emit("recovery_return", recovered=arg is not None, line=frame.f_lineno, values=values)
                if event == "exception":
                    emit("recovery_exception", message=str(arg[1]), line=frame.f_lineno)
                return trace
            old_trace = sys.gettrace()
            sys.settrace(trace)
            try:
                return original_recovery(self, pressure, enthalpy, comp)
            finally:
                sys.settrace(old_trace)

        PropertyPackage._recover_two_phase_ph = recovery

    original_solve = cycle._solve_case

    def solve(hp):
        emit("solve_start", state=snapshot(hp))
        try:
            failure = original_solve(hp)
        except Exception:
            emit("solve_exception", traceback=traceback.format_exc(), state=snapshot(hp))
            raise
        emit("solve_return", failure=repr(failure), state=snapshot(hp))
        return failure

    cycle._solve_case = solve
    original_probe = pressure.find_gc_pinch_residual

    def probe(value, hp, *args, **kwargs):
        emit("pressure_start", pressure=value)
        try:
            residual = original_probe(value, hp, *args, **kwargs)
        except Exception:
            emit("pressure_exception", pressure=value, traceback=traceback.format_exc(), state=snapshot(hp))
            raise
        emit("pressure_return", pressure=value, residual=residual)
        return residual

    pressure.find_gc_pinch_residual = probe
    if request.get("guard_invalid_duty"):
        from HeatPumpLab.cycles.RecuperativeHeatPump.solver import _cycle_analysis as analysis
        original_sink = analysis.set_gas_cooler_sink_temperatures

        def guarded_sink(hp):
            duties = [gc.duty for gc in hp.gas_coolers]
            if sum(duties) <= 0 or any(gc.duty < -gc.tol for gc in hp.gas_coolers):
                emit("invalid_duty_guard", state=snapshot(hp))
                raise ValueError("Diagnostic guard: invalid gas-cooler duty before sink reconstruction")
            return original_sink(hp)

        analysis.set_gas_cooler_sink_temperatures = guarded_sink
    original_record = _OuterDesignObjective._record_from_hp
    final_state = {}

    def record(self, **kwargs):
        final_state.update(snapshot(kwargs["hp"]))
        return original_record(self, **kwargs)

    _OuterDesignObjective._record_from_hp = record
    objective = _OuterDesignObjective(max_evals=1, refrigerant=request["ref"],
                                     base_case=request["base_case"] | request.get("overrides", {}),
                                     use_fast_objective=request.get("fast", True))
    with profile_solve() as profile:
        result = serialise(objective.evaluate(request["x"]))
    data = dict(request=request, record=result, state=final_state, profile=profile.as_dict(),
                wall_seconds=time.perf_counter()-started)
    prefix.with_suffix(".result.json").write_text(json.dumps(clean(data), indent=2))
    stream.close()
    faulthandler.cancel_dump_traceback_later()


def run_one(request):
    prefix = Path(request["prefix"])
    request_path = prefix.with_suffix(".request.json")
    if prefix.with_suffix(".status.json").exists():
        return json.loads(prefix.with_suffix(".status.json").read_text())
    request_path.write_text(json.dumps(request, indent=2))
    env = dict(os.environ, MPLCONFIGDIR="/tmp/heatpumplab-botorch-mpl",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    started = time.perf_counter()
    with prefix.with_suffix(".log").open("w") as log:
        try:
            result = subprocess.run([sys.executable, str(Path(__file__).resolve()), "--worker", str(request_path)],
                                    stdout=log, stderr=log, timeout=request.get("timeout", 200), env=env)
            status = dict(id=prefix.name, returncode=result.returncode, timed_out=False)
        except subprocess.TimeoutExpired:
            status = dict(id=prefix.name, returncode=None, timed_out=True)
    status["wall_seconds"] = time.perf_counter()-started
    prefix.with_suffix(".status.json").write_text(json.dumps(status, indent=2))
    return status


def run_batch(requests, workers=8):
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(run_one, request) for request in requests]):
            print(json.dumps(future.result()), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker")
    parser.add_argument("--requests")
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
        return
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.requests:
        requests = json.loads(Path(args.requests).read_text())
    else:
        requests = []
        model_hash = hashlib.sha256(b"".join(p.read_bytes() for p in sorted((ROOT / "HeatPumpLab").rglob("*.py")))).hexdigest()
        for path in sorted(SOURCE.glob("*_*.json")):
            data = json.loads(path.read_text())
            if "records" not in data:
                continue
            assert data["metadata"]["model_source_sha256"] == model_hash
            for index, row in enumerate(data["records"], 1):
                if row["success"]:
                    continue
                requests.append(dict(prefix=str(OUTPUT / f"{path.stem}_{index:03d}_baseline"),
                                     ref=data["metadata"]["refrigerant"], x=row["x"],
                                     base_case=data["metadata"]["base_case"], source=str(path),
                                     evaluation=index, original=row, model_source_sha256=model_hash))
        assert len(requests) == 44
        (OUTPUT / "baseline_requests.json").write_text(json.dumps(requests, indent=2))
    run_batch(requests, args.workers)


if __name__ == "__main__":
    main()
