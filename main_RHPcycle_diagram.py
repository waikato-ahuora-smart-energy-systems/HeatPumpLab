"""Re-solve one saved design and export its T-s / P-h cycle with dome, isotherms and isobars.

Usage:
  .venv/bin/python main_RHPcycle_diagram.py --run-dir results/RHPcomposition/fixed200_constrained_v22_5000/Propane_Isopentane \
      --evaluation 937 --output results/RHPcomposition/cycle_diagrams/Propane_Isopentane_v22_e937.json

Writes a JSON with the solved cycle paths (per component, sampled states), the mixture's bubble and dew
lines, isotherms (P-h) and the cycle's isobars (T-s). Units: MPa, Celsius, kJ/kg, kJ/(kg K).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from HeatPumpLab.components.stream import Stream
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import DEFAULT_BASE_CASE

ISOTHERMS_C = list(range(0, 276, 25))


def state(s):
    vals = dict(T=s.temperature, P=s.pressure, h=s.enthalpy, s=s.entropy, q=s.quality)
    return {k: (float(v) if v is not None and np.isfinite(v) and v > -1e6 else None) for k, v in vals.items()}


def flash(fluid, **kw):
    try:
        st = Stream(fluid=fluid, **kw)
        out = state(st)
        return out if out['h'] is not None and out['s'] is not None else None
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--run-dir', required=True)
    parser.add_argument('--evaluation', type=int, required=True, help='1-based evaluation number in result.json')
    parser.add_argument('--output', required=True)
    parser.add_argument('--divisions', type=int, default=40)
    args = parser.parse_args()

    R = json.loads((Path(args.run_dir)/'result.json').read_text())
    meta, rec = R['metadata'], R['records'][args.evaluation-1]
    first, partner = meta.get('first', 'Propane'), meta['partner']
    temperature = float(meta.get('temperature_C', 200))
    objective = _OuterDesignObjective(max_evals=1, refrigerant=f'{first};{partner}|0.5;0.5|mass',
                                      base_case=dict(meta.get('base_case') or DEFAULT_BASE_CASE, T_air_tar=temperature))
    case = objective.vector_to_case(rec['x'])
    hp, _ = objective._analyse_case(case)
    report = hp.get_report()
    fluid = case['refrigerant']
    cop, vhc = float(report['COP']), float(report.get('VHC', 0))/1000
    print(f'{fluid}  x={np.round(rec["x"], 4).tolist()}')
    print(f'COP {cop:.5f} (saved {rec["COP"]:.5f})  VHC {vhc:.4f} MJ/m3 (saved {rec["VHC_MJ_per_m3"]:.4f})')

    cycle_sets, sink_sets = hp.get_components_sets_to_graph()
    paths = []
    for comp, side in cycle_sets:
        pts = [state(s) for s in comp.get_graph_data(divisions=args.divisions)[side]]
        paths.append(dict(component=str(getattr(comp, 'name', 'Component')), side=side, points=pts))
    sink = []
    for comp, side in sink_sets:
        data = comp.get_graph_data(divisions=args.divisions)
        sink.append(dict(component=str(comp.name), sink=[state(s) for s in data['sink']],
                         source=[state(s) for s in data['source']]))

    probe = Stream(fluid=fluid)
    p_crit = float(probe.P_crit)
    pressures = np.geomspace(0.05, 0.995*p_crit, 90)
    bubble = [r for r in (flash(fluid, pressure=float(P), quality=0.0) for P in pressures) if r]
    dew = [r for r in (flash(fluid, pressure=float(P), quality=1.0) for P in pressures) if r]

    p_max = max(p['P'] for path in paths for p in path['points'] if p['P'])
    sweep = np.geomspace(0.05, 1.3*p_max, 160)
    isotherms = [dict(T=T, points=[r for r in (flash(fluid, pressure=float(P), temperature=float(T)) for P in sweep) if r])
                 for T in ISOTHERMS_C]
    levels = sorted({round(p['P'], 4) for path in paths for p in path['points'][:1]+path['points'][-1:] if p['P']})
    merged = []
    for P in levels:
        if not merged or P/merged[-1] > 1.03:
            merged.append(P)
    temps = np.linspace(-20, 260, 141)
    isobars = [dict(P=P, points=[r for r in (flash(fluid, pressure=float(P), temperature=float(T)) for T in temps) if r])
               for P in merged]

    out = dict(fluid=fluid, run_dir=args.run_dir, evaluation=args.evaluation, x=rec['x'], temperature_C=temperature,
               COP=cop, VHC_MJ_per_m3=vhc, saved_COP=rec['COP'], saved_VHC=rec['VHC_MJ_per_m3'],
               P_crit_MPa=p_crit, T_crit_C=float(probe.T_crit), units='MPa, C, kJ/kg, kJ/(kg K)',
               report={k: (float(v) if isinstance(v, (int, float, np.floating)) and np.isfinite(v) else v)
                       for k, v in report.items() if isinstance(v, (int, float, str, np.floating))},
               paths=paths, sink=sink, bubble=bubble, dew=dew, isotherms=isotherms, isobars=isobars)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    Path(args.output).write_text(json.dumps(out, default=str))
    print('wrote', args.output, f'({len(paths)} paths, {len(bubble)}/{len(dew)} dome points, {len(isobars)} isobars)')


if __name__ == '__main__':
    main()
