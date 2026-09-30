"""Compare fixed-split Pareto runs (e.g. ridge-seeded vs unseeded) across seeds and fluid pairs.

Usage: main_RHPridge_benchmark_report.py LABEL=DIR [LABEL=DIR ...] --output OUT
Each DIR holds [seed_*/]<First>_<Partner>/result.json (main_RHPfixed_mobo.py output).

Per run: evaluations, stop reason, failures, validated hypervolume (and at 600 evaluations),
max COP, COP at VHC >= 3 MJ/m3, the RHX2 = 0 ("ridge B") branch on the final front and the
evaluation that first found it, and the worst COP shortfall against the pooled front of every
run of that pair (a direct measure of a missed branch).
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

BRANCH_RHX2 = .5     # K: RHX2 at most this is the RHX2 = 0 face
BRANCH_RHX1 = 30.    # K: and RHX1 above this (not the zero-lift corner)


def front(points):
    """Non-dominated (COP, VHC) maximisation front, sorted by increasing VHC."""
    pts = sorted(points, key=lambda p: (-p[0], -p[1]))
    out, best = [], -np.inf
    for p in pts:
        if p[1] > best:
            out.append(p)
            best = p[1]
    return sorted(out, key=lambda p: p[1])


def cop_at(points, vhc):
    return max((p[0] for p in points if p[1] >= vhc), default=None)


def run_metrics(path):
    R = json.loads((path/'result.json').read_text())
    recs = R['records']
    vpath = path/'validated_frontier.json'
    V = json.loads(vpath.read_text()) if vpath.exists() else None
    ok = [(r['COP'], r['VHC_MJ_per_m3'], i+1, r['x']) for i, r in enumerate(recs)
          if r['success'] and r['COP'] is not None and r['VHC_MJ_per_m3'] is not None]
    fr = front(ok)
    branch = [p for p in fr if p[3][1] <= BRANCH_RHX2 and p[3][0] >= BRANCH_RHX1]
    c600 = path/'check_00600.json'
    seeded = [r for r in recs if r.get('method') == 'ridge_seed']
    body = [r for r in recs if r.get('method') != 'ridge_seed']
    med = lambda v: float(np.median(v)) if v else None
    return dict(
        seed=R['metadata']['seed'], pair=path.name, status=R.get('status'), n=len(recs),
        stop=R.get('stop_reason'), seed_evaluations=len(seeded), failures=sum(not r['success'] for r in recs),
        hv=V['hypervolume'] if V else None, hv600=json.loads(c600.read_text())['hypervolume'] if c600.exists() else None,
        maxcop=max(p[0] for p in fr), cop3=cop_at(fr, 3.), max_vhc=max(p[1] for p in fr),
        validation_rejections=V.get('validation_rejections') if V else None,
        branch_points=len(branch), branch_first_evaluation=min((p[2] for p in branch), default=None),
        branch_best=max(((p[0], p[1], p[3]) for p in branch), default=None),
        median_proposal_seconds=med([r.get('proposal_seconds', 0) for r in body if r.get('proposal_seconds')]),
        median_evaluation_seconds=med([r.get('evaluation_wall_seconds') for r in recs if r.get('evaluation_wall_seconds')]),
        _front=[(p[0], p[1]) for p in fr])


def shortfall(fr, pooled, grid):
    """Max relative COP shortfall of a run's front against the pooled front over a VHC grid."""
    worst = (0., None)
    for v in grid:
        a, b = cop_at(fr, v), cop_at(pooled, v)
        if b is None:
            continue
        gap = 1.-(a or 0.)/b
        if gap > worst[0]:
            worst = (gap, float(v))
    return worst


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sets', nargs='+', help='LABEL=DIR')
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    runs = []
    for item in args.sets:
        label, directory = item.split('=', 1)
        for path in sorted(set(Path(directory).glob('seed_*/*/result.json'))|set(Path(directory).glob('*/result.json'))):
            runs.append(dict(run_metrics(path.parent), label=label))
    for pair in sorted({r['pair'] for r in runs}):
        group = [r for r in runs if r['pair'] == pair]
        pooled = front([p for r in group for p in r['_front']])
        grid = np.linspace(min(p[1] for p in pooled), max(p[1] for p in pooled), 201)
        for r in group:
            r['max_cop_shortfall'], r['shortfall_at_vhc'] = shortfall(r['_front'], pooled, grid)
    summary = []
    for key in sorted({(r['pair'], r['label']) for r in runs}):
        group = [r for r in runs if (r['pair'], r['label']) == key]
        stat = {}
        for m in ('hv', 'hv600', 'maxcop', 'cop3', 'max_cop_shortfall', 'n', 'failures', 'seed_evaluations',
                  'median_proposal_seconds', 'median_evaluation_seconds'):
            v = [r[m] for r in group if r[m] is not None]
            if v:
                stat[m] = dict(mean=float(np.mean(v)), sd_percent=float(100*np.std(v, ddof=1)/np.mean(v)) if len(v) > 1 and np.mean(v) else None,
                               worst=float(min(v)) if m in ('hv', 'hv600', 'maxcop', 'cop3') else float(max(v)))
        stat['branch_found'] = f"{sum(r['branch_points'] > 0 for r in group)} of {len(group)}"
        summary.append(dict(pair=key[0], label=key[1], runs=len(group), **stat))
    out = dict(runs=[{k: v for k, v in r.items() if k != '_front'} for r in runs], summary=summary)
    Path(args.output).write_text(json.dumps(out, indent=1)+'\n')
    for s in summary:
        f = lambda m, d=4: (f"{s[m]['mean']:.{d}f} ({s[m]['sd_percent']:.2f}%)" if s.get(m) and s[m]['sd_percent'] is not None
                            else (f"{s[m]['mean']:.{d}f}" if s.get(m) else '-'))
        print(f"{s['pair']:20s} {s['label']:12s} n={s['runs']} HV {f('hv')} HV600 {f('hv600')} "
              f"maxCOP {f('maxcop')} COP@3 {f('cop3')} shortfall {f('max_cop_shortfall')} "
              f"branch {s['branch_found']} evals {f('n', 0)} seedE {f('seed_evaluations', 0)} "
              f"fail {f('failures', 0)} prop {f('median_proposal_seconds', 2)} eval {f('median_evaluation_seconds', 2)}")


if __name__ == '__main__':
    main()
