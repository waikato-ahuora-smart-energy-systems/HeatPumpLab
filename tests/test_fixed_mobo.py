import copy
import json

import numpy as np
import pytest

from rhp_frontier_stability import assess, frontier
import main_RHPfixed_mobo as runner


def record(i, **overrides):
    point = dict(x=[0., 0., .1, .5, .5], T_air_tar=200., success=True,
                 COP=3.-i*.01, VHC_MJ_per_m3=1.+i*.1, SEC=.5, Pmax_MPa=5.)
    point.update(overrides)
    return point


def block():
    return [record(i) for i in range(20)]+[record(0, success=False) for _ in range(10)]+[record(0)]*470


def test_stability_requires_three_complete_blocks_and_2000_evaluations():
    records, previous = [], None
    for n in (500, 1000, 1500, 2000):
        records.extend(block())
        previous = assess(records, previous)
        assert previous['stop_for_stability'] == (n == 2000)
        assert previous['stable_blocks'] == n//500-1
    assert len(frontier(records)) == 20  # Duplicate observations do not create apparent coverage.


def test_curve_change_resets_streak_even_for_small_hypervolume_gain():
    previous = assess(block())
    previous['stable_blocks'] = 2
    records = block()+block()
    records[-1] = record(19, COP=2.82)
    changed = assess(records, previous)
    assert changed['hypervolume_relative_gain'] < .001
    assert changed['curve_max_relative_change'] > .001
    assert changed['stable_blocks'] == 0
    assert changed['vhc_thresholds'] == previous['vhc_thresholds']


def test_vhc_extension_and_failed_or_sparse_searches_do_not_stop():
    previous = assess(block())
    records = block()+block()
    records[-1] = record(19, COP=.001, VHC_MJ_per_m3=3.5)
    check = assess(records, previous)
    assert check['vhc_extent_relative_change'] > .001 and not check['stable_block']
    check = assess(block()+[record(0, success=False)]*500, previous)
    assert not check['stable_block']
    previous = None
    for n in (500, 1000, 1500, 2000, 2500):
        previous = assess([record(0)]*n, previous)
    assert not previous['stop_for_stability']


def test_frontier_hypervolume_uses_observed_points_and_excludes_rejections():
    records = [record(0, COP=3., VHC_MJ_per_m3=1.), record(1, COP=2., VHC_MJ_per_m3=3.),
               record(2, COP=9., validation_rejected=True)]
    check = assess(records)
    assert check['hypervolume'] == pytest.approx(7.)
    assert check['cop_at_vhc_thresholds'][50] == 2.


@pytest.mark.parametrize('first', ['Propane', 'Isobutane'])
def test_runner_pause_resume_cap_and_fresh_replay(tmp_path, monkeypatch, first):
    calls = []
    pause_after = [2]

    def evaluate(partner, x, temperature, output, **kwargs):
        assert temperature == 200.
        assert kwargs['first'] == first
        calls.append(list(x))
        if len(calls) == pause_after[0]:
            (output.parent/'PAUSE').touch()
        return record(0, x=list(x), COP=3.-x[3], VHC_MJ_per_m3=1.+x[3])

    monkeypatch.setattr(runner, 'exact_evaluation', evaluate)
    monkeypatch.setattr(runner, 'render', lambda *args, **kwargs: None)
    runner.run_family('Butane', 6, 25, tmp_path, first=first, batch_size=1)
    path = tmp_path/f'{first}_Butane'/'result.json'
    paused = json.loads(path.read_text())
    assert paused['status'] == 'paused' and len(paused['records']) == 2
    (tmp_path/'PAUSE').unlink()
    runner.run_family('Butane', 6, 25, tmp_path, first=first, batch_size=1)
    result = json.loads(path.read_text())
    assert result['complete'] and result['stop_reason'] == 'evaluation_cap'
    assert len(result['records']) == 6 and len(calls) > 6
    assert all(r['matches'] for r in result['frontier_replays'].values())
    assert result['metadata']['components'] == [first, 'Butane']
    assert (tmp_path/f'{first}_Butane'/'validated_frontier.json').exists()
    before = len(calls)
    runner.run_family('Butane', 6, 25, tmp_path, first=first, batch_size=1)
    assert len(calls) == before


def test_batched_runner_pauses_on_a_batch_boundary(tmp_path, monkeypatch):
    calls = []

    def evaluate(partner, x, temperature, output, **kwargs):
        calls.append(list(x))
        if len(calls) == 2:
            (output.parent/'PAUSE').touch()
        return record(0, x=list(x), COP=3.-x[3], VHC_MJ_per_m3=1.+x[3])

    monkeypatch.setattr(runner, 'exact_evaluation', evaluate)
    monkeypatch.setattr(runner, 'render', lambda *args, **kwargs: None)
    runner.run_family('Butane', 10, 25, tmp_path, batch_size=4)
    paused = json.loads((tmp_path/'Propane_Butane'/'result.json').read_text())
    # The PAUSE request arrives mid-batch; the whole batch is still recorded.
    assert paused['status'] == 'paused' and len(paused['records']) == 4
    assert [r['batch_size'] for r in paused['records']] == [4]*4
    (tmp_path/'PAUSE').unlink()
    runner.run_family('Butane', 10, 25, tmp_path, batch_size=4)
    result = json.loads((tmp_path/'Propane_Butane'/'result.json').read_text())
    assert result['complete'] and len(result['records']) == 10
    assert [r['batch_size'] for r in result['records']] == [4]*8+[2]*2


def test_local_checkpoint_resume_reproduces_next_proposal():
    result = dict(metadata=dict(seed=29), records=[
        record(i, x=[i*.25, i*.5, .2+i*.002, .02+i*.01, .98-i*.01]) for i in range(20)],
        local128_phase=dict(start_evaluations=20, seed=29))
    search, epoch = runner.restore(result)
    x = search.ask()
    item = record(2, x=x.tolist(), COP=2.7, VHC_MJ_per_m3=1.1)
    result['records'].append(item)
    search.tell(x, item['COP'], item['VHC_MJ_per_m3'], True)
    result['local128_state'] = dict(epoch_start=epoch,
        model_states={str(k): copy.deepcopy(v) for k, v in search.states.items()})
    restored, restored_epoch = runner.restore(result)
    assert restored_epoch == epoch
    np.testing.assert_allclose(search.ask(), restored.ask(), rtol=0, atol=1e-10)
