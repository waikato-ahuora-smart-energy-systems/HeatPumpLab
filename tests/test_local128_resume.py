"""Continuation epochs must survive checkpoints without changing local proposals."""
import copy

import numpy as np

from rhp_local_pareto_search import restore_search


def dataset(count):
    return dict(records=[dict(x=[i*.25, i*.5, .45+i*.002, .02+i*.01, .98-i*.01],
                              COP=3-i*.005, VHC_MJ_per_m3=1+i*.02, success=True)
                         for i in range(count)],
                local128_phase=dict(start_evaluations=20, seed=20260912))


def test_epoch_refresh_uses_updated_front_and_resets_model_cache():
    result=dataset(67)
    search, epoch=restore_search(result)
    assert epoch==20 and search.baseline_count==20 and len(search.x)==67
    assert max(search.anchor_indices)<20
    result=dataset(68)
    result['local128_state']=dict(epoch_start=20, model_states={'0': {'visits':99}})
    search, epoch=restore_search(result)
    assert epoch==68 and search.baseline_count==68
    assert 67 in search.anchor_indices and search.states=={}
    assert len(search.training_indices(search.centres[0]))==68


def test_resume_preserves_next_proposal_and_cached_parameters():
    result=dataset(20)
    search, epoch=restore_search(result)
    proposal=search.ask()
    result['records'].append(dict(x=proposal.tolist(), COP=2.95, VHC_MJ_per_m3=1.4, success=True))
    search.tell(proposal, 2.95, 1.4, True)
    result['local128_state']=dict(epoch_start=epoch,
                                  model_states={str(k):copy.deepcopy(v) for k,v in search.states.items()})
    restored, restored_epoch=restore_search(result)
    assert restored_epoch==epoch and restored.states==search.states
    assert np.array_equal(restored.centres,search.centres)
    assert np.allclose(restored.ask(),search.ask(),rtol=0,atol=1e-10)
