"""Split-stream heat pump vs Design A (series RHX1), Propane/Isopentane.

Grid over the two design variables shared with Design A (RHX1 cold-side lift, stage-2 power
share); the gas-cooler pressure is optimised inside each design. Results are compared with the
Design A runs of main_RHPdesignA_caps.py at the same composition, temperature and approach.

  .venv/bin/python main_SLTHP_compare.py run --workers 4
  .venv/bin/python main_SLTHP_compare.py one --temperature 150 --lift 60 --split 0.5
  .venv/bin/python main_SLTHP_compare.py summary

Output: results/SLTHP_compare_v26/<composition>/A<approach>/T<T>/L<lift>_S<split>.json, summary.md.
"""
from __future__ import annotations

import argparse
import json
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT / 'results' / 'SLTHP_compare_v26'
DESIGN_A = ROOT / 'results' / 'RHPdesignA_caps_v26'
COMPOSITIONS = ((0.225, 0.775),)
TEMPERATURES = (150, 175, 200)
LIFTS = tuple(range(0, 141, 10))
SPLITS = (0.40, 0.50)
APPROACHES = (20.0,)
RHX1_COLD_END = (25.0, 15.0)   # RHX1 hot-branch outlet above suction-gas inlet; 25 K ~ Design A's LMTD
PMAX_CAP_BAR = 27.0
TDIS_CAPS_C = (180.0, 190.0, 200.0, 210.0, None)


def fluid(w):
    return 'Propane;Isopentane|' + ';'.join(f'{v:g}' for v in w) + '|mass'


def label(w):
    return '-'.join(f'{v:.3f}' for v in w)


def path_for(w, a, T, lift, split, ce=25.0, cap=PMAX_CAP_BAR):
    # cap=None: gas-cooler pressure searched over the full P_H_bounds; stored under CE<ce>_nocap
    tag = '' if cap == PMAX_CAP_BAR else '_nocap' if cap is None else f'_P{cap:g}'
    return OUTPUT / label(w) / f'A{a:g}' / f'CE{ce:g}{tag}' / f'T{T:g}' / f'L{lift:g}_S{split:.2f}.json'


def solve_one(w, a, T, lift, split, ce=25.0, cap=PMAX_CAP_BAR):
    from HeatPumpLab.cycles.SplitHeatPump import SplitHeatPump, SplitHPConfig
    cfg = SplitHPConfig(refrigerant=fluid(w), T_air_tar=float(T), dT_GC_min=float(a),
                        dT_rhx1_cold=float(lift), comp2_power_frac=float(split), Pmax_cap=None if cap is None else cap / 10.0,
                        dT_rhx1_cold_end=float(ce), **({} if cap is not None else dict(n_PH_scan=36)))
    started = time.time()
    res = SplitHeatPump(cfg).solve()
    out = dict(composition=list(w), approach_K=a, rhx1_cold_end_K=ce, Pmax_cap_bar=cap, T_air_tar=T, lift_K=lift, split=split,
               seconds=time.time() - started, **res.to_dict())
    p = path_for(w, a, T, lift, split, ce, cap)
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(out, indent=1, default=float))
    return out


def run(workers):
    jobs = [(w, a, T, l, s, ce) for w in COMPOSITIONS for a in APPROACHES for ce in RHX1_COLD_END for T in TEMPERATURES
            for l in LIFTS for s in SPLITS if not path_for(w, a, T, l, s, ce).exists()]
    print(f'{len(jobs)} designs to solve with {workers} workers', flush=True)
    with ProcessPoolExecutor(workers) as pool:
        futs = {pool.submit(solve_one, *j): j for j in jobs}
        for k, f in enumerate(as_completed(futs), 1):
            w, a, T, l, s, ce = futs[f]
            try:
                r = f.result()
                msg = f"COP {r['COP']:.3f} P {r['P_H'] * 10:.1f} bar" if r['success'] else r['message'][:70]
            except Exception as exc:  # noqa: BLE001
                msg = f'ERROR {exc}'
            print(f"{time.strftime('%H:%M')} [{k}/{len(jobs)}] CE{ce:g} T{T} lift {l} split {s}: {msg}", flush=True)
    summary()


def design_a_best(w, a, T, pcap, tcap):
    tag = '' if a == 20.0 else f'_A{a:g}'
    f = DESIGN_A / 'Propane_Isopentane' / label(w) / f'T{T:g}{tag}' / 'result.json'
    if not f.exists():
        return None
    best = None
    for r in json.loads(f.read_text())['records']:
        if not r.get('success') or r.get('COP') is None or r.get('Pmax_bar') is None or r['Pmax_bar'] > pcap:
            continue
        if tcap is not None and (r.get('T_comp_discharge') or 1e9) > tcap:
            continue
        if best is None or r['COP'] > best['COP']:
            best = r
    return best


def summary():
    lines = ['# Split-stream heat pump vs Design A (Propane/Isopentane)', '',
             f'Pmax <= {PMAX_CAP_BAR:g} bar abs. Split cycle: RHX1 in parallel with the low-temperature gas cooler; '
             'the LT gas-cooler branch leaves at air inlet + approach, the RHX1 branch at suction-gas inlet + CE; '
             'Design A: RHX1 in series (main_RHPdesignA_caps.py; RHX1 LMTD about 19-29 K). Compressor model 5, motor 0.93.', '']
    for w in COMPOSITIONS:
        for a in APPROACHES:
            for ce in RHX1_COLD_END:
                lines += [f'## {label(w)}, gas-cooler approach {a:g} K, RHX1 cold end {ce:g} K', '',
                          '| T air | Tdis cap | split COP | VHC | Pmax | Tdis | lift / split | RHX1 share | T split | RHX1 LMTD | '
                          'Design A COP | VHC | Pmax | Tdis | change |',
                          '|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|']
                for T in TEMPERATURES:
                    recs = []
                    for f in (OUTPUT / label(w) / f'A{a:g}' / f'CE{ce:g}' / f'T{T:g}').glob('*.json'):
                        d = json.loads(f.read_text())
                        if d['success'] and d['P_H'] * 10 <= PMAX_CAP_BAR + 1e-9:
                            recs.append(d)
                    for tc in TDIS_CAPS_C:
                        ok = [d for d in recs if tc is None or d['report']['T_comp_discharge'] <= tc + 1e-9]
                        s = max(ok, key=lambda d: d['COP']) if ok else None
                        da = design_a_best(w, a, T, PMAX_CAP_BAR, tc)
                        tct = 'none' if tc is None else f'{tc:g}'
                        if s:
                            rp = s['report']
                            ts = '' if rp['T_split'] is None else f"{rp['T_split']:.0f}"
                            lm = '' if rp['rhx1_LMTD_K'] is None else f"{rp['rhx1_LMTD_K']:.1f}"
                            sc = (f"{s['COP']:.3f} | {s['VHC']:.2f} | {s['P_H']*10:.1f} | {rp['T_comp_discharge']:.0f} | "
                                  f"{s['lift_K']:g} / {s['split']:.2f} | {rp['alpha_RHX1']:.2f} | {ts} | {lm}")
                        else:
                            sc = 'none | | | | | | | '
                        dc = (f"{da['COP']:.3f} | {da['VHC_MJ_per_m3']:.2f} | {da['Pmax_bar']:.1f} | {da['T_comp_discharge']:.0f}"
                              if da else 'n/a | | |')
                        ch = f"{100 * (s['COP'] / da['COP'] - 1):+.1f} %" if (s and da) else ''
                        lines.append(f'| {T} | {tct} | {sc} | {dc} | {ch} |')
                lines.append('')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / 'summary.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='cmd', required=True)
    r = sub.add_parser('run'); r.add_argument('--workers', type=int, default=4)
    o = sub.add_parser('one')
    o.add_argument('--temperature', type=float, required=True); o.add_argument('--lift', type=float, required=True)
    o.add_argument('--split', type=float, default=0.5); o.add_argument('--approach', type=float, default=20.0)
    o.add_argument('--composition', default='0.225,0.775')
    o.add_argument('--cold-end', type=float, default=25.0, help='RHX1 hot-branch outlet above suction-gas inlet (K)')
    o.add_argument('--no-cap', action='store_true', help='no Pmax cap (search the full gas-cooler pressure range)')
    sub.add_parser('summary')
    a = p.parse_args()
    if a.cmd == 'run':
        run(a.workers)
    elif a.cmd == 'one':
        w = tuple(float(v) for v in a.composition.split(','))
        res = solve_one(w, a.approach, a.temperature, a.lift, a.split, a.cold_end, None if a.no_cap else PMAX_CAP_BAR)
        keys = ('COP', 'VHC', 'Pmax_bar', 'P_evap_bar', 'P_int_bar', 'PR_1', 'PR_2', 'T_suction', 'Comp_1_T_out',
                'T_comp_discharge', 'T_split', 'T_valve_in', 'T_air_mid', 'alpha_RHX1', 'q_rhx1_kJkg', 'q_air_kJkg',
                'rhx1_LMTD_K', 'approaches_K')
        print(res['message'], f"({res['seconds']:.0f} s)")
        for k in keys:
            print(f'  {k}: {res["report"].get(k)}')
    else:
        summary()


if __name__ == '__main__':
    main()
