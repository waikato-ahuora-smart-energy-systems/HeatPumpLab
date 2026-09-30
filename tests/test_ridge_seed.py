"""Pinch-balance ridge seeding (rhp_ridge.py) and its use in the fixed-split runner."""
import json

import numpy as np

import main_RHPfixed_mobo as runner
import rhp_ridge as rr


def synthetic(x):
    """Three gas coolers; ridges at S = 40 + 40 theta (gc1/gc2) and S = 90 (gc2/gc3)."""
    S = x[0]+x[1]
    theta = x[0]/S if S > 0 else 0.
    s1, s2 = 40+40*theta, 90.
    active = 0 if S < s1 else (1 if S < s2 else 2)
    A = [19.95]*3
    if active == 0:
        A[1], A[2] = 19.95+.25*(s1-S), 19.95+.25*(s2-S)
    elif active == 1:
        A[0], A[2] = 19.95+.25*(S-s1), 19.95+.25*(s2-S)
    else:
        A[0], A[1] = 19.95+.25*(S-s1), 19.95+.25*(S-s2)
    cop = 2.4+.2*x[3]-.002*min(abs(S-s1), abs(S-s2))
    return dict(x=list(map(float, x)), success=True, COP=cop, VHC_MJ_per_m3=5.-.02*S+2*x[4], SEC=1/(cop-1),
                Pmax_MPa=5., gc_pinch_active=f'gc{active+1}_hot', gc_approaches=A,
                constraint_margins=None, message='success')


def test_ridges_found_on_every_line_including_exact_grid_hits():
    records, ridges = rr.seed(synthetic, 123.333, .5, workers=2, fractions=(.5,), thetas=(0., .3, 1.))
    found = {(r['theta'], tuple(r['pair'])): r['S'] for r in ridges}
    assert found[(0., (0, 1))] == 40. and found[(1., (0, 1))] == 80.
    assert abs(found[(.3, (0, 1))]-52.) < rr.X_TOLERANCE
    assert all(abs(S-90.) < rr.X_TOLERANCE for (t, p), S in found.items() if p == (1, 2))
    keys = [tuple(np.round(r['x'], 9)) for r in records]
    assert len(keys) == len(set(keys))           # S = 0 shared across lines, evaluated once
    assert all(ridge['record'] is not None for ridge in ridges)


def test_runner_seeds_then_searches(tmp_path, monkeypatch):
    monkeypatch.setattr(runner, 'exact_evaluation', lambda partner, x, t, out, **kw: synthetic(x))
    monkeypatch.setattr(runner, 'render', lambda *a, **k: None)
    runner.run_family('Butane', 200, 25, tmp_path, fixed_comp_frac=.5, batch_size=2)  # seeding is the default
    result = json.loads((tmp_path/'Propane_Butane'/'result.json').read_text())
    n_seed = runner.seed_count(result)
    assert result['ridge_seed_done'] and n_seed > 40 and len(result['ridge_seeds']) >= 12
    assert len(result['records']) == 200 and result['complete']
    assert all(r['method'] != 'ridge_seed' for r in result['records'][n_seed:])
    xs = [tuple(np.round(r['x'], 9)) for r in result['records']]
    assert len(xs) == len(set(xs))


def test_ridge_seeding_refuses_budgets_it_would_exhaust(tmp_path):
    import pytest
    with pytest.raises(ValueError, match='Ridge seeding'):
        runner.run_family('Butane', 150, 25, tmp_path, fixed_comp_frac=.5)


def test_refinement_fills_the_widest_cop_vhc_gap():
    # Propane/Pentane ridge B (theta = 1) from the Sep 2026 benchmark: the 2/3-1 pair trades off most.
    ridge_b = [(0., 2.3505, .438), (1/3, 2.5189, 1.437), (2/3, 2.5142, 2.452), (1., 2.4307, 4.006)]
    assert abs(rr.refine_fraction(ridge_b)-5/6) < 1e-9
    # One composition dominating the next leaves no gap.
    assert rr.refine_fraction([(0., 2.4, 1.), (1., 2.5, 2.)]) is None
    assert rr.refine_fraction([(0., 2.4, 1.)]) is None


def test_seed_adds_one_refinement_line_per_theta():
    records, ridges = rr.seed(synthetic, 123.333, .5, workers=2, thetas=(0., 1.))
    extra = [r for r in ridges if r.get('refinement')]
    assert {r['theta'] for r in extra} == {0., 1.}
    assert all(min(abs(r['fraction']-m) for m in (1/6, .5, 5/6)) < 1e-9 for r in extra)
    keys = [tuple(np.round(r['x'], 9)) for r in records]
    assert len(keys) == len(set(keys))
