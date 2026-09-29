"""Pinch-balance ridge seeding for the fixed-split Pareto search.

At a fixed composition and lift share theta = RHX1 / S, COP against total lift
S = RHX1 + RHX2 peaks where two adjacent gas coolers are pinched at once (ridge
scan, Sep 2026). Below a ridge, gas cooler k sets the pinch and k + 1 has slack;
the slack closes as S grows. The ridge is the root of

    r(S) = A[k+1](S) - A[k](S)

where A[j] is the minimum approach in gas cooler j (report 'gc_approaches',
hot-to-cold order). Each ridge point is the high-COP end of a COP-VHC trade-off
line, so solving ridges seeds every branch of the front, including RHX2 = 0.
"""
from __future__ import annotations

import math
from concurrent.futures import ThreadPoolExecutor

import numpy as np

SEED_FRACTIONS = (0.0, 1/3, 2/3, 1.0)
SEED_THETAS = (0.0, 0.5, 1.0)
COARSE_POINTS = 6          # S samples per line, including S = 0
LIFT_CAP = 100.0           # each lift's own bound (DEFAULT_BOUNDS)
X_TOLERANCE = 0.5          # K in S
R_TOLERANCE = 0.02         # K in approach balance
MAX_ROOT_STEPS = 6


def line_limit(theta, max_lift, cap=LIFT_CAP):
    limit = max_lift
    if theta > 0:
        limit = min(limit, cap/theta)
    if theta < 1:
        limit = min(limit, cap/(1-theta))
    return limit


def design(fraction, theta, S, split):
    return np.array([theta*S, (1-theta)*S, split, fraction, 1-fraction])


def active_index(record):
    label = record.get('gc_pinch_active')
    if not label or not label.startswith('gc'):
        return None
    return int(label[2:].split('_')[0]) - 1


def balance(record, k):
    """A[k+1] - A[k]; None when the approaches are unavailable."""
    a = record.get('gc_approaches') or []
    if len(a) <= k+1 or not all(math.isfinite(v) for v in a[k:k+2]):
        return None
    return a[k+1] - a[k]


def solve_line(evaluate, fraction, theta, max_lift, split):
    """Coarse sweep in S, then a safeguarded secant (Illinois) solve on each pinch switch.

    Returns (records, ridges): records in evaluation order, ridges as dicts with the
    index of the record closest to balance.
    """
    records, ridges = [], []

    def run(S):
        records.append(evaluate(design(fraction, theta, S, split)))
        records[-1]['ridge_line'] = dict(fraction=fraction, theta=theta, S=float(S))
        return records[-1]

    grid = np.linspace(0.0, line_limit(theta, max_lift), COARSE_POINTS)
    coarse = [(S, run(S)) for S in grid]
    for (s0, r0), (s1, r1) in zip(coarse, coarse[1:]):
        k0, k1 = active_index(r0), active_index(r1)
        if k0 is None or k1 is None or k1 <= k0:
            continue
        for k in range(k0, k1):          # one ridge per adjacent pair crossed
            f0, f1 = balance(r0, k), balance(r1, k)
            if f0 is None or f1 is None or not (f0 >= 0 >= f1):
                continue
            # A coarse point already on the ridge needs no root solve.
            on_ridge = [(i, f, s) for i, f, s in ((coarse.index((s0, r0)), f0, s0), (coarse.index((s1, r1)), f1, s1))
                        if abs(f) < R_TOLERANCE]
            if on_ridge:
                i, f, s = on_ridge[0]
                ridges.append(dict(fraction=fraction, theta=theta, pair=[k, k+1], S=float(s),
                                   balance=float(f), record=records.index(coarse[i][1])))
                continue
            lo, hi, flo, fhi, side = s0, s1, f0, f1, 0
            best = None
            for _ in range(MAX_ROOT_STEPS):
                S = hi - fhi*(hi-lo)/(fhi-flo)
                S = min(max(S, lo+0.1*(hi-lo)), hi-0.1*(hi-lo))
                rec = run(S)
                f = balance(rec, k)
                if f is None:
                    break
                if best is None or abs(f) < abs(best[1]):
                    best = (len(records)-1, f, S)
                if abs(f) < R_TOLERANCE:
                    break
                if f > 0:
                    lo, flo = S, f
                    if side == 1:
                        fhi /= 2
                    side = 1
                else:
                    hi, fhi = S, f
                    if side == -1:
                        flo /= 2
                    side = -1
                if hi-lo < X_TOLERANCE:
                    break
            if best is not None:
                ridges.append(dict(fraction=fraction, theta=theta, pair=[k, k+1],
                                   S=float(best[2]), balance=float(best[1]), record=best[0]))
    return records, ridges


def seed(evaluate, max_lift, split, workers=4, fractions=SEED_FRACTIONS, thetas=SEED_THETAS):
    """Solve every line in parallel; return records in a fixed (line, step) order."""
    lines = [(f, t) for f in fractions for t in thetas]
    cache = {}

    def cached(x):
        key = tuple(np.round(x, 9))
        if key not in cache:
            cache[key] = evaluate(x)
        return dict(cache[key])

    # S = 0 is the same design for every theta: evaluate it once per composition first.
    for f in fractions:
        cached(design(f, 0.0, 0.0, split))
    with ThreadPoolExecutor(workers) as pool:
        results = list(pool.map(lambda line: solve_line(cached, line[0], line[1], max_lift, split), lines))
    records, ridges, seen = [], [], set()
    for line_records, line_ridges in results:
        offset = len(records)
        index_map = {}
        for i, r in enumerate(line_records):
            key = tuple(np.round(r['x'], 9))
            if key in seen:
                continue
            seen.add(key)
            index_map[i] = len(records)
            records.append(r)
        for ridge in line_ridges:
            ridge = dict(ridge, record=index_map.get(ridge['record']))
            ridges.append(ridge)
    return records, ridges
