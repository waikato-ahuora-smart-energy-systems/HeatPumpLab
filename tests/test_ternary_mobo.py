"""Ternary mixtures in the fixed-temperature Pareto search, ridge seeding and runner."""
import json

import numpy as np
import pytest

import main_RHPfixed_mobo as runner
import rhp_ridge as rr
from rhp_mobo_search import ParetoBOSearch
from tests.test_ridge_seed import synthetic


def test_encoding_round_trips_and_simplex_is_enforced():
    search = ParetoBOSearch(123.333, components=3)
    assert search.dimension == 5 and search.comp.tolist() == [3, 4]
    x = np.array([20., 30., .5, .5, .3, .2])
    u = search.encode(x)
    np.testing.assert_allclose(search.decode(u), x, atol=1e-12)
    assert len(search._admissible([[.1, .1, .5, .7, .4]])) == 0      # z1 + z2 > 1
    assert len(search._admissible([[.1, .1, .5, .6, .4]])) == 1      # on the edge without component 3
    with pytest.raises(ValueError):
        search.encode([20., 30., .5, .5, .5])                        # binary vector for a ternary search


def test_initial_design_covers_the_composition_triangle():
    search = ParetoBOSearch(123.333, components=3)
    z = np.asarray(search.initial)[:, search.comp]
    assert (z.sum(1) <= 1+1e-12).all()
    full = np.column_stack((z, 1-z.sum(1)))
    for vertex in np.eye(3):
        assert np.isclose(full, vertex).all(1).any()
    assert len(search.initial) > 32


def test_binary_search_is_unchanged_by_the_generalisation():
    search = ParetoBOSearch(123.333)
    assert search.dimension == 4 and len(search.initial) == 32
    np.testing.assert_allclose(search.decode([.2, .3, .5, .7]), [20., 30., .5, .7, .3])


def test_candidates_include_the_third_component_free_edge():
    search = ParetoBOSearch(123.333, components=3)
    search.tell([0., 0., .5, 1/3, 1/3, 1/3], 2.5, 3., True)
    pool, _ = search._candidates()
    z = pool[:, search.comp]
    assert np.isclose(z.sum(1), 1.).any() and (z.sum(1) <= 1+1e-12).all()


def test_ternary_ridge_seeding_runs_every_composition_line():
    records, ridges = rr.seed(synthetic, 123.333, .5, workers=2, fractions=rr.SEED_COMPOSITIONS_3)
    lines = {(tuple(np.round(r['ridge_line']['fraction'], 6)), r['ridge_line']['theta']) for r in records
             if 'ridge_line' in r}
    assert len(lines) == len(rr.SEED_COMPOSITIONS_3)*len(rr.SEED_THETAS)
    assert all(len(r['x']) == 6 and np.isclose(sum(r['x'][3:]), 1.) for r in records)
    assert not any(r.get('refinement') for r in ridges)


def test_ternary_runner_seeds_searches_and_validates(tmp_path, monkeypatch):
    seen = []

    def fake(partner, x, t, out, **kw):
        seen.append(tuple(kw.get('components') or ()))
        return synthetic(x)
    monkeypatch.setattr(runner, 'exact_evaluation', fake)
    monkeypatch.setattr(runner, 'render', lambda *a, **k: None)
    runner.run_family('Isobutane', 400, 25, tmp_path, third='Isopentane', batch_size=2)
    out = tmp_path/'Propane_Isobutane_Isopentane'
    result = json.loads((out/'result.json').read_text())
    assert result['metadata']['components'] == ['Propane', 'Isobutane', 'Isopentane']
    assert result['metadata']['search_dimensions'] == 5
    assert result['complete'] and len(result['records']) == 400
    assert set(seen) == {('Propane', 'Isobutane', 'Isopentane')}
    z = np.asarray([r['x'][3:] for r in result['records']])
    assert z.shape[1] == 3 and np.allclose(z.sum(1), 1.) and (z >= -1e-12).all()
    assert (z.min(1) > .01).any()                                    # genuinely ternary designs were tried
    assert all(.4-1e-12 <= r['x'][2] <= .6+1e-12 for r in result['records'])
