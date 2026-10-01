"""Joint-temperature runner: batching, constraint-aware failures, stop rule and ridge seeding."""
import json

import numpy as np

import main_RHPtemperature_mobo as runner
from rhp_temperature_mobo import TemperatureParetoSearch, admissible, decode, encode


def synthetic(x, temperature):
    """Two gas coolers pinch-balanced at S = 40 + 40 theta; COP falls with temperature."""
    S = x[0]+x[1]
    theta = x[0]/S if S > 0 else 0.
    ridge = 40+40*theta
    A = [19.95, 19.95+.25*(ridge-S)] if S < ridge else [19.95+.25*(S-ridge), 19.95]
    cop = 3.4-.005*(temperature-100)+.2*x[3]-.002*abs(S-ridge)
    return dict(x=list(map(float, x)), T_air_tar=float(temperature), success=True, COP=cop,
                VHC_MJ_per_m3=5.-.02*S+2*x[4]+.005*(temperature-100), SEC=1/(cop-1), Pmax_MPa=5.,
                gc_pinch_active='gc1_hot' if S < ridge else 'gc2_hot', gc_approaches=A,
                constraint_margins=None, message='success')


def record(temperature, x, success=True, margins=None):
    return dict(T_air_tar=float(temperature), COP=3., VHC_MJ_per_m3=2., x=list(x), success=success,
                constraint_margins=margins)


def test_only_marginless_failures_count_as_numerical():
    margins = dict(Comp_1_discharge_superheat=-2., Comp_1_suction_superheat=5., Comp_1_pressure_ratio=1.)
    search = TemperatureParetoSearch([record(100, [0., 0., .5, .2, .8]),
                                      record(105, [0., 0., .5, .3, .7], success=False, margins=margins),
                                      record(110, [0., 0., .5, .4, .6], success=False)])
    assert search.numerical_failures().tolist() == [False, False, True]
    assert search.margins[1] is not None and search.margins[1][0] < 0


def test_batch_proposals_are_distinct_and_pending_clears_on_tell():
    search = TemperatureParetoSearch(seed=77)
    batch = search.ask_batch(4)
    assert len(search.pending) == 4 and len(search.batch_diagnostics) == 4
    u = np.array([encode(x, t) for x, t in batch])
    assert admissible(u).all() and len(np.unique(np.round(u, 9), axis=0)) == 4
    search.tell(synthetic(*batch[0]))
    assert len(search.pending) == 3


def test_initial_points_already_observed_are_replaced():
    x, temperature = decode(TemperatureParetoSearch(seed=5).initial[0])
    seeded = dict(synthetic(x, temperature), method='ridge_seed')
    search = TemperatureParetoSearch([seeded], seed=5)
    assert search.offset == 1 and search.step_index() == 0
    search.ask()
    assert search.diagnostic['method'] == 'initial_replaced_maximin'


def test_stability_check_needs_consecutive_flat_blocks(monkeypatch):
    monkeypatch.setattr(runner, 'STOP_CHECK_INTERVAL', 2)
    monkeypatch.setattr(runner, 'STOP_EARLIEST', 4)
    monkeypatch.setattr(runner, 'STOP_BLOCKS', 2)
    records = [synthetic([0., 0., .5, .5, .5], 150.), synthetic([10., 0., .5, .5, .5], 150.)]
    records += [dict(records[0], x=[0., 0., .5, .5, .5]) for _ in range(6)]
    checks = []
    assert not runner.stability_check(records[:4], checks)     # second check has a change but only one block
    assert runner.stability_check(records, checks)
    assert [c['evaluations'] for c in checks] == [2, 4, 6, 8]


def patch_runner(monkeypatch):
    import rhp_temperature_report
    monkeypatch.setattr(runner, 'exact_evaluation', lambda partner, x, t, out, **kw: synthetic(x, t))
    monkeypatch.setattr(rhp_temperature_report, 'render_front', lambda *a, **k: None)


def test_batched_runner_pauses_and_resumes(tmp_path, monkeypatch):
    patch_runner(monkeypatch)
    kw = dict(batch_size=3, ridge_seed=False, stop_rule=False)
    monkeypatch.setattr(runner, 'verify_final', lambda *a, **k: True)
    (tmp_path/'PAUSE').touch()
    runner.run_family('Butane', 12, 3, tmp_path, **kw)
    (tmp_path/'PAUSE').unlink()
    state = json.loads((tmp_path/'Propane_Butane'/'state.json').read_text())
    assert state['status'] == 'paused' and state['evaluations'] == 0
    runner.run_family('Butane', 12, 3, tmp_path, **kw)
    records = runner.read_journal(tmp_path/'Propane_Butane'/'records.jsonl')
    assert len(records) == 12 and all(r['batch_size'] == 3 for r in records)
    assert json.loads((tmp_path/'Propane_Butane'/'state.json').read_text())['complete']


def test_ridge_seeded_runner(tmp_path, monkeypatch):
    patch_runner(monkeypatch)
    monkeypatch.setattr(runner, 'SEED_TEMPERATURES', (220.,))
    monkeypatch.setattr(runner, 'RIDGE_SEED_MIN_BUDGET', 0)
    monkeypatch.setattr(runner, 'verify_final', lambda *a, **k: True)
    runner.run_family('Butane', 120, 3, tmp_path, batch_size=2, stop_rule=False)
    out = tmp_path/'Propane_Butane'
    records = runner.read_journal(out/'records.jsonl')
    state = json.loads((out/'state.json').read_text())
    n_seed = sum(r.get('method') == 'ridge_seed' for r in records)
    assert state['ridge_seed_done'] and state['ridge_seeded_temperatures'] == [220.]
    assert 60 < n_seed < 120 and all(r.get('method') == 'ridge_seed' for r in records[:n_seed])
    assert all(r['T_air_tar'] == 220. for r in records[:n_seed])
    assert len(records) == 120
    keys = [tuple(np.round(r['x'], 9))+(r['T_air_tar'],) for r in records]
    assert len(keys) == len(set(keys))
    # Resume after completion is a no-op.
    runner.run_family('Butane', 120, 3, tmp_path, batch_size=2, stop_rule=False)
    assert len(runner.read_journal(out/'records.jsonl')) == 120
