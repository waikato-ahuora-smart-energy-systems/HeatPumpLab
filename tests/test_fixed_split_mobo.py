"""A fixed stage allocation must stay fixed through every Pareto search path."""
import copy
import json

import numpy as np
import pytest

from rhp_mobo_search import ParetoBOSearch
import main_RHPfixed_mobo as runner
from tests.test_fixed_mobo import record


def observe(search, n=32):
    for _ in range(n):
        x = search.ask()
        assert x[2] == .5
        search.tell(x, 3.-x[3]-.001*x[0], 1.+2*x[3]+.001*x[1], True)


def test_fixed_split_is_removed_from_initial_global_gp_and_maximin_search():
    search = ParetoBOSearch(123.333, 18, fixed_comp_frac=.5)
    observe(search)
    assert np.asarray(search.x).shape == (32, 3)
    assert len(np.unique(search.x, axis=0)) == 32
    assert {0., 1.}.issubset(set(np.asarray(search.x)[:, 2]))
    pool, distance = search._candidates()
    assert pool.shape[1] == 3 and np.all(distance > 1e-7)
    assert search.decode(pool[np.argmax(distance)])[2] == .5
    proposal = search.ask()  # Actual fitted GP/acquisition optimisation in three dimensions.
    assert proposal[2] == .5
    assert not any(d.get('stage') == 'acquisition fallback' for d in search.diagnostics)
    search.tell(proposal, 2.5, 2., True)
    result = dict(metadata=dict(seed=18, fixed_comp_frac=.5), records=[
        dict(x=search.decode(u).tolist(), COP=y[0], VHC_MJ_per_m3=y[1], success=ok)
        for u,y,ok in zip(search.x,search.y,search.success)],
        kernel_state=search.kernel_state, mean_state=search.mean_state)
    restored, _ = runner.restore(json.loads(json.dumps(result)))
    np.testing.assert_allclose(search.ask(), restored.ask(), atol=1e-9, rtol=0)
    with pytest.raises(ValueError, match='fixed_comp_frac'):
        search.tell([0., 0., .4, .5, .5], 2., 2., True)


def test_fixed_split_survives_local_gp_and_checkpoint_resume():
    result = dict(metadata=dict(seed=29,fixed_comp_frac=.5), records=[
        record(i,x=[i*.25,i*.5,.5,.02+i*.01,.98-i*.01]) for i in range(20)],
        local128_phase=dict(start_evaluations=20,seed=29))
    search, epoch = runner.restore(result)
    assert search.centres.shape == (6,3)
    x = search.ask()
    assert x[2] == .5
    assert 'acquisition_error' not in search.last_diagnostic
    item=record(2,x=x.tolist(),COP=2.7,VHC_MJ_per_m3=1.1)
    result['records'].append(item)
    search.tell(x,item['COP'],item['VHC_MJ_per_m3'],True)
    result['local128_state']=dict(epoch_start=epoch,
        model_states={str(k):copy.deepcopy(v) for k,v in search.states.items()})
    restored, restored_epoch=runner.restore(result)
    assert restored_epoch == epoch and restored.dimension == 3
    np.testing.assert_allclose(search.ask(),restored.ask(),rtol=0,atol=1e-9)


def test_fixed_split_runner_preserves_setting_during_pause_resume_and_replay(tmp_path,monkeypatch):
    calls=[]
    def evaluate(partner,x,temperature,output,**kwargs):
        assert x[2] == .5
        calls.append(list(x))
        if len(calls)==2: (output.parent/'PAUSE').touch()
        return record(0,x=list(x),COP=3.-x[3],VHC_MJ_per_m3=1.+x[3])
    monkeypatch.setattr(runner,'exact_evaluation',evaluate)
    monkeypatch.setattr(runner,'render',lambda *args,**kwargs:None)
    runner.run_family('Pentane',6,25,tmp_path,fixed_comp_frac=.5, ridge_seed=False)
    (tmp_path/'PAUSE').unlink()
    runner.run_family('Pentane',6,25,tmp_path,fixed_comp_frac=.5, ridge_seed=False)
    result=json.loads((tmp_path/'Propane_Pentane/result.json').read_text())
    assert result['complete'] and len(calls)>6
    assert result['metadata']['fixed_comp_frac']==.5
    assert result['metadata']['search_dimensions']==3
    assert all(r['matches'] for r in result['frontier_replays'].values())
    with pytest.raises(ValueError,match='mismatch'):
        runner.run_family('Pentane',6,25,tmp_path,fixed_comp_frac=None, ridge_seed=False)


@pytest.mark.parametrize('value',[0.,1.,np.nan,np.inf])
def test_invalid_fixed_split_is_rejected(value):
    with pytest.raises(ValueError,match='fixed_comp_frac'):
        ParetoBOSearch(123.333,fixed_comp_frac=value)


def test_native_fixed_half_split_gives_equal_stage_pressure_ratios(tmp_path):
    from main_RHPtemperature_mobo import exact_evaluation
    r=exact_evaluation('Pentane',[20.,40.,.5,.5,.5],200.,tmp_path)
    assert r['success'],r['message']
    assert r['x'][2]==.5
    assert r['PR_1']==pytest.approx(r['PR_2'],rel=1e-6)
    assert r['PR_1']<=10 and r['PR_2']<=10
