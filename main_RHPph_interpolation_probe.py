"""Diagnostic PH interpolation experiment; does not change production recovery."""

import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "results/RHPcomposition/ph_interpolation_probe"


def interpolate(a, b, weight=.5):
    fields = {}
    for name in ("T", "P", "D", "Dl", "Dv", "q", "h", "s", "e"):
        fields[name] = (1-weight)*getattr(a, name) + weight*getattr(b, name)
    for name in ("x", "y"):
        fields[name] = [(1-weight)*x + weight*y for x, y in zip(getattr(a, name), getattr(b, name))]
    # Interpolate specific volume, since phase volumes add linearly.
    fields["D"] = 1 / ((1-fields["q"])/fields["Dl"] + fields["q"]/fields["Dv"])
    return a._replace(**fields)


def inspect_state(pp, state, p, h, comp):
    import math
    if not 0 < state.q < 1:
        return {"valid": False, "quality": state.q}
    try:
        liquid = pp.RP.THERMdll(state.T, state.Dl, state.x)
        vapor = pp.RP.THERMdll(state.T, state.Dv, state.y)
        fl = pp.RP.FUGCOFdll(state.T, state.Dl, state.x)
        fv = pp.RP.FUGCOFdll(state.T, state.Dv, state.y)
        return dict(valid=pp._valid_two_phase_ph(state, p, h, comp),
                    temperature_C=state.T-273.15, quality=state.q,
                    enthalpy_error_kJ_kg=((1-state.q)*liquid.h+state.q*vapor.h-h)/pp.molar_mass,
                    pressure_relative_error=max(abs(liquid.P/p-1), abs(vapor.P/p-1)),
                    fugacity_log_error=max(abs(math.log(state.x[i]*fl.phi[i]/(state.y[i]*fv.phi[i])))
                                            for i in range(len(comp))),
                    material_error=max(abs((1-state.q)*state.x[i]+state.q*state.y[i]-comp[i])
                                       for i in range(len(comp))))
    except (ValueError, ZeroDivisionError, OverflowError) as exc:
        return {"valid": False, "error": str(exc)}


def probe(request):
    from HeatPumpLab.utils.properties_package import PropertyPackage
    pp = PropertyPackage(request["ref"])
    comp = request["context"]["composition"]
    pp._ensure_refprop_context(comp)
    p, h = pp._convert_props_from_SI_with_C("PH", *request["context"]["inputs"])
    direct = pp.RP.ABFLSHdll("PH", p, h, comp, 0)
    results = []
    for axis, steps in (("h", (1e-5, 1e-4, 1e-3, .01, .1, 1., 5.)),
                        ("p", (1e-6, 1e-5, 1e-4, .001, .01))):
        valid_sides = {-1: [], 1: []}
        for step in steps:
            delta = step*pp.molar_mass if axis == "h" else step*p
            states, endpoints = [], []
            for sign in (-1, 1):
                pn = p + sign*delta if axis == "p" else p
                hn = h + sign*delta if axis == "h" else h
                state = pp.RP.ABFLSHdll("PH", pn, hn, comp, 0)
                endpoint = dict(ierr=state.ierr)
                if state.ierr == 0 and 0 < state.q < 1:
                    refined = pp._refine_binary_ph(state, pn, hn, comp)
                    state = refined if refined is not None else state
                    endpoint.update(inspect_state(pp, state, pn, hn, comp))
                else:
                    endpoint["valid"] = False
                endpoints.append(endpoint)
                states.append(state)
                if endpoint["valid"]:
                    valid_sides[sign].append((step, state))
            row = dict(axis=axis, step=step, endpoints=endpoints)
            if all(e["valid"] for e in endpoints):
                mid = interpolate(*states)
                row["interpolation"] = inspect_state(pp, mid, p, h, comp)
                solved = pp._refine_binary_ph(mid, p, h, comp)
                row["refined"] = inspect_state(pp, solved, p, h, comp) if solved is not None else None
            results.append(row)
        for low_step, low in valid_sides[-1]:
            for high_step, high in valid_sides[1]:
                if low_step == high_step:
                    continue
                mid = interpolate(low, high, low_step/(low_step+high_step))
                solved = pp._refine_binary_ph(mid, p, h, comp)
                results.append(dict(axis=axis, steps=[low_step, high_step],
                                    interpolation=inspect_state(pp, mid, p, h, comp),
                                    refined=inspect_state(pp, solved, p, h, comp) if solved is not None else None))
    return dict(request=request, direct_ierr=direct.ierr, trials=results)


def run(request):
    path = OUT / (request["id"] + ".json")
    req = path.with_suffix(".request.json")
    req.write_text(json.dumps(request, indent=2))
    env = dict(os.environ, MPLCONFIGDIR="/tmp/heatpumplab-botorch-mpl",
               OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
    try:
        with path.with_suffix(".log").open("w") as log:
            result = subprocess.run([sys.executable, __file__, "--worker", str(req)],
                                    stdout=log, stderr=log, env=env, timeout=90)
        print(request["id"], result.returncode, flush=True)
    except subprocess.TimeoutExpired:
        print(request["id"], "timeout", flush=True)


def interpolation_recovery(pp, p, h, comp):
    """Experimental fallback, installed only in the diagnostic child process."""
    for axis, steps in (("h", (1e-5, 1e-4, .001, .01, .1, 1., 5.)),
                        ("p", (1e-6, 1e-5, 1e-4, .001, .01))):
        sides = {-1: [], 1: []}
        for step in steps:
            for sign in (-1, 1):
                pn = p*(1+sign*step) if axis == "p" else p
                hn = h+sign*step*pp.molar_mass if axis == "h" else h
                state = pp.RP.ABFLSHdll("PH", pn, hn, comp, 0)
                if state.ierr or not 0 < state.q < 1:
                    continue
                corrected = pp._refine_binary_ph(state, pn, hn, comp)
                state = corrected if corrected is not None else state
                if not pp._valid_two_phase_ph(state, pn, hn, comp):
                    continue
                sides[sign].append((step, state))
                for other_step, other in sides[-sign]:
                    low, high = (state, other) if sign == -1 else (other, state)
                    low_step = step if sign == -1 else other_step
                    middle = interpolate(low, high, low_step/(step+other_step))
                    if pp._valid_two_phase_ph(middle, p, h, comp):
                        return middle
                    solved = pp._refine_binary_ph(middle, p, h, comp)
                    if solved is not None:
                        return solved
    return None


def cycle_worker(path):
    from HeatPumpLab.utils.properties_package import PropertyPackage
    from main_RHPfailure_investigation import worker
    original = PropertyPackage._recover_two_phase_ph
    events = []

    def recover(self, p, h, comp):
        result = original(self, p, h, comp)
        if result is not None:
            return result
        try:
            result = interpolation_recovery(self, p, h, comp)
            events.append(dict(pressure_MPa=p/1000, enthalpy_kJ_kg=h/self.molar_mass,
                               recovered=result is not None))
            return result
        except (ValueError, RuntimeError, ZeroDivisionError, OverflowError) as exc:
            events.append(dict(error=str(exc), recovered=False))
            return None

    PropertyPackage._recover_two_phase_ph = recover
    try:
        worker(path)
    finally:
        prefix = Path(json.loads(Path(path).read_text())["prefix"])
        prefix.with_suffix(".interpolation.json").write_text(json.dumps(events, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker")
    parser.add_argument("--cycle")
    args = parser.parse_args()
    if args.cycle:
        cycle_worker(args.cycle)
        return
    if args.worker:
        path = Path(args.worker)
        result = probe(json.loads(path.read_text()))
        Path(str(path).replace(".request.json", ".json")).write_text(json.dumps(result, indent=2))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    summary = json.loads((ROOT / "results/RHPcomposition/v15_failure_fixes/verified/summary.json").read_text())
    requests = []
    for t in summary["trials"]:
        if t["role"] != "old_failure" or t["result"]["record"]["success"]:
            continue
        contexts = [r["property_context"] for r in t["result"]["state"]["numerical_retries"]
                    if r.get("property_context")]
        requests.append(dict(id=t["id"], ref=t["request"]["ref"], context=contexts[0]))
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(run, requests))


if __name__ == "__main__":
    main()
