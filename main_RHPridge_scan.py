"""Ridge scan: COP and VHC on a grid of total lift S = RHX1 + RHX2 and split theta = RHX1 / S.

Tests the observation that RHX1 + RHX2 = const is an approximate ridge: for each propane mass
fraction, does the best S depend on theta, and is the middle of the ridge worse than its ends?
Fixed compressor split 0.5, 200 C sink, 0 K superheat margin. Each solve runs in a fresh
isolated worker (same oracle as the Pareto runs). Resumable: results.jsonl is appended.

Run: .venv/bin/python main_RHPridge_scan.py --partner Pentane --n-cpu 6
"""
import argparse
import json
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

from main_RHPtemperature_mobo import exact_evaluation
from rhp_local_pareto_search import max_lift_for

FRACTIONS = (0.20, 0.43, 0.79, 0.97)
THETAS = (0.0, 0.25, 0.5, 0.75, 1.0)
LIFTS = tuple(float(s) for s in range(10, 121, 10))


def grid(max_lift):
    points = []
    for f in FRACTIONS:
        points.append(dict(fraction=f, S=0.0, theta=None))
        for theta in THETAS:
            for S in LIFTS:
                # Each lift is also bounded by DEFAULT_BOUNDS (0-100 K) on its own.
                if S <= max_lift and theta*S <= 100 and (1-theta)*S <= 100:
                    points.append(dict(fraction=f, S=S, theta=theta))
    return points


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--partner', default='Pentane')
    parser.add_argument('--first', default='Propane')
    parser.add_argument('--n-cpu', type=int, default=6)
    args = parser.parse_args()
    out = Path(f'results/RHPcomposition/ridge_scan/{args.first}_{args.partner}')
    out.mkdir(parents=True, exist_ok=True)
    path = out/'results.jsonl'
    done = set()
    if path.exists():
        for line in path.read_text().splitlines():
            r = json.loads(line)
            done.add((r['fraction'], r['S'], r['theta']))
    todo = [p for p in grid(max_lift_for({})) if (p['fraction'], p['S'], p['theta']) not in done]
    print(f'{len(done)} done, {len(todo)} to run', flush=True)

    def solve(p):
        theta = p['theta'] or 0.0
        x = np.array([theta*p['S'], (1-theta)*p['S'], .5, p['fraction'], 1-p['fraction']])
        item = exact_evaluation(args.partner, x, 200., out, first=args.first)
        keep = ('success', 'COP', 'VHC_MJ_per_m3', 'Pmax_MPa', 'message', 'constraint_margins',
                'gc_pinch_active', 'PR_1', 'PR_2', 'T_comp_discharge')
        return dict(p, x=x.tolist(), **{k: item.get(k) for k in keep})

    with ThreadPoolExecutor(args.n_cpu) as pool, path.open('a') as fh:
        futures = {pool.submit(solve, p): p for p in todo}
        for i, fut in enumerate(as_completed(futures), 1):
            try:
                row = fut.result()
            except Exception as exc:  # keep the scan going; record the cause
                import traceback
                traceback.print_exc()
                row = dict(futures[fut], success=False, error=f'{type(exc).__name__}: {exc}')
            fh.write(json.dumps(row)+'\n')
            fh.flush()
            if i % 20 == 0 or i == len(futures):
                print(f'{i}/{len(futures)}', flush=True)
    (out/'DONE').write_text('done')


if __name__ == '__main__':
    main()
