"""Weekend sweep: every binary and ternary from the fluid pool at 200 C, then a second seed for the top mixtures.

Runs main_RHPfixed_mobo.py as one subprocess per mixture, --slots at a time (each uses --batch-size
REFPROP workers). Completed mixtures are skipped and interrupted ones resume, so the sweep can simply
be restarted. Touch <output>/PAUSE to stop launching new mixtures (running ones also pause cleanly).

  .venv/bin/python main_RHPfluid_sweep.py --output results/RHPfixed200_v25 --slots 3

Phase 1: all mixtures, seed 20260915, 1,000 evaluations, ordered ternaries -> binaries -> CO2 mixtures.
Phase 2: the --top best by COP at VHC >= 3 MJ/m3, rerun with seed 20260916 in <output>_seed2.
Writes <output>/sweep_summary.{json,md} after each phase.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from main_RHPfixed_mobo import PARTNERS

ROOT = Path(__file__).resolve().parent
LEVELS = (2.0, 3.0, 4.0, 5.0)


def mixtures(pool=PARTNERS):
    """All binaries and ternaries; components in pool order (Propane first). CO2 mixtures last."""
    mixes = [list(c) for c in combinations(pool, 3)]+[list(c) for c in combinations(pool, 2)]
    return sorted(mixes, key=lambda m: ('CO2' in m, -len(m)))


def folder(output, mix):
    return Path(output)/'_'.join(mix)


def complete(output, mix):
    path = folder(output, mix)/'result.json'
    try:
        return bool(json.loads(path.read_text()).get('complete'))
    except (OSError, ValueError):
        return False


def run(output, mix, seed, budget, batch_size, extra=()):
    out = folder(output, mix)
    out.mkdir(parents=True, exist_ok=True)
    if complete(output, mix):
        return mix, 'already complete', 0.
    if (Path(output)/'PAUSE').exists():
        return mix, 'not started (PAUSE)', 0.
    args = [sys.executable, str(ROOT/'main_RHPfixed_mobo.py'), '--first', mix[0], '--partner', mix[1],
            '--budget', str(budget), '--seed', str(seed), '--batch-size', str(batch_size), '--output', str(output)]
    if len(mix) == 3:
        args += ['--third', mix[2]]
    args += list(extra)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1')
    started = time.time()
    with (Path(output)/f"{'_'.join(mix)}.log").open('a') as log:
        log.write(f"\n=== {time.strftime('%Y-%m-%d %H:%M:%S')} {' '.join(args)}\n")
        log.flush()
        code = subprocess.run(args, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=ROOT).returncode
    hours = (time.time()-started)/3600
    status = 'complete' if complete(output, mix) else f'incomplete (exit {code})'
    print(f"{time.strftime('%H:%M')} {'/'.join(mix)}: {status}, {hours:.2f} h", flush=True)
    return mix, status, hours


def metrics(output, mix):
    path = folder(output, mix)/'validated_frontier.json'
    if not path.exists():
        return None
    v = json.loads(path.read_text())
    pts = [(r['COP'], r['VHC_MJ_per_m3'], r.get('Pmax_MPa')) for r in v['frontier']]
    cop_at = {f'COP_VHC_ge_{lvl:g}': max((c for c, w, _ in pts if w >= lvl), default=None) for lvl in LEVELS}
    result = json.loads((folder(output, mix)/'result.json').read_text())
    return dict(mixture='/'.join(mix), components=len(mix), evaluations=len(result['records']),
                stop=result.get('stop_reason'), failures=sum(not r['success'] for r in result['records']),
                hypervolume=v['hypervolume'], max_COP=max(c for c, _, _ in pts), max_VHC=max(w for _, w, _ in pts),
                front_points=len(pts), validation_rejections=v.get('validation_rejections'),
                min_Pmax_MPa=min((p for *_, p in pts if p), default=None), **cop_at)


def summarise(output, mixes, title):
    rows = [m for m in (metrics(output, mix) for mix in mixes) if m]
    rows.sort(key=lambda r: -(r['COP_VHC_ge_3'] or 0))
    (Path(output)/'sweep_summary.json').write_text(json.dumps(rows, indent=1)+'\n')
    fmt = lambda v, d=4: '-' if v is None else f'{v:.{d}f}'
    lines = [f'# {title}', '', f"{len(rows)} of {len(mixes)} mixtures complete; ranked by COP at VHC >= 3 MJ/m3.", '',
             '| Mixture | Max COP | VHC>=2 | VHC>=3 | VHC>=4 | VHC>=5 | Max VHC | HV | Evals | Fails | Rej |',
             '|---|---|---|---|---|---|---|---|---|---|---|']
    for r in rows:
        lines.append(f"| {r['mixture']} | {fmt(r['max_COP'])} | {fmt(r['COP_VHC_ge_2'])} | {fmt(r['COP_VHC_ge_3'])} | "
                     f"{fmt(r['COP_VHC_ge_4'])} | {fmt(r['COP_VHC_ge_5'])} | {fmt(r['max_VHC'], 2)} | "
                     f"{fmt(r['hypervolume'], 3)} | {r['evaluations']} | {r['failures']} | {r['validation_rejections']} |")
    (Path(output)/'sweep_summary.md').write_text('\n'.join(lines)+'\n')
    return rows


def sweep(output, mixes, seed, budget, batch_size, slots, extra=()):
    pending = [m for m in mixes if not complete(output, m)]
    print(f'{len(mixes)-len(pending)} complete, {len(pending)} to run, {slots} at a time', flush=True)
    with ThreadPoolExecutor(slots) as pool:
        return list(pool.map(lambda m: run(output, m, seed, budget, batch_size, extra), pending))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--output', default='results/RHPfixed200_v25')
    parser.add_argument('--slots', type=int, default=3, help='Mixtures run in parallel')
    parser.add_argument('--batch-size', type=int, default=4)
    parser.add_argument('--budget', type=int, default=1000)
    parser.add_argument('--top', type=int, default=5, help='Mixtures rerun with a second seed (0 to skip)')
    parser.add_argument('--list', action='store_true', help='Print the mixtures and their status, then exit')
    parser.add_argument('--mixtures', nargs='+', default=None,
                        help='Run only these mixtures, e.g. Propane/Butane/Isopentane (default: the full pool)')
    parser.add_argument('--comp-model', type=int, default=None, help='Passed to main_RHPfixed_mobo.py')
    parser.add_argument('--motor-eff', type=float, default=None, help='Passed to main_RHPfixed_mobo.py')
    parser.add_argument('--min-fraction', type=float, default=None,
                        help='Passed to main_RHPfixed_mobo.py (its default: 0.01; 0 allows pure fluids)')
    args = parser.parse_args()
    output = Path(args.output)
    output.mkdir(parents=True, exist_ok=True)
    mixes = [m.split('/') for m in args.mixtures] if args.mixtures else mixtures()
    for m in mixes:
        if len(m) not in (2, 3) or any(c not in PARTNERS for c in m) or len(set(m)) != len(m):
            parser.error(f"Bad mixture {'/'.join(m)}: 2-3 distinct components from {', '.join(PARTNERS)}")
    extra = ((['--comp-model', str(args.comp_model)] if args.comp_model is not None else [])
             + (['--motor-eff', str(args.motor_eff)] if args.motor_eff is not None else [])
             + (['--min-fraction', str(args.min_fraction)] if args.min_fraction is not None else []))
    if args.list:
        for m in mixes:
            print('done ' if complete(output, m) else 'to run', '/'.join(m))
        return
    sweep(output, mixes, 20260915, args.budget, args.batch_size, args.slots, extra)
    rows = summarise(output, mixes, 'Fluid sweep at 200 C, seed 20260915')
    if args.top and not (output/'PAUSE').exists():
        best = [r['mixture'].split('/') for r in rows if r['COP_VHC_ge_3'] is not None][:args.top]
        second = Path(f'{output}_seed2')
        second.mkdir(parents=True, exist_ok=True)
        print('Second seed for:', ', '.join('/'.join(m) for m in best), flush=True)
        sweep(second, best, 20260916, args.budget, args.batch_size, args.slots, extra)
        summarise(second, best, 'Top mixtures, second seed 20260916')
    print('SWEEP DONE', flush=True)


if __name__ == '__main__':
    main()
