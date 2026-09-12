"""Direction, anchor coverage and neighbourhood selection checks."""
import numpy as np
from rhp_local_pareto_search import LocalParetoSearch, choose_anchors, front_indices, hypervolume


def test_front_and_hypervolume_keep_ties_and_reject_invalid_values():
    values=[[2.,4.],[3.,3.],[1.,1.],[2.,4.],[0.,0.],[np.nan,5.]]
    assert sorted(front_indices(values)) == [0,1,3]
    assert hypervolume(values) == 11.


def baseline():
    records=[]
    for i in range(20):
        records.append(dict(x=[i,2*i,.2+i*.02,.05+i*.045,.95-i*.045],success=True,
                            COP=3-i*.02,VHC_MJ_per_m3=1+i*.2))
    return dict(records=records)


def test_anchors_cover_both_ends_and_subset_uses_nearest_successes():
    data=baseline()
    anchors=choose_anchors(data['records'])
    assert len(set(anchors))==6 and 0 in anchors and 19 in anchors
    search=LocalParetoSearch(data,subset_size=5)
    centre=search.centres[0]
    expected=np.argsort(np.linalg.norm(np.asarray(search.x)-centre,axis=1))[:5]
    assert set(search.training_indices(centre))==set(expected)
    assert len(front_indices(search.y))==20  # Full front survives the smaller training set.


def test_pool_stays_in_small_box_and_excludes_previous_candidates():
    search=LocalParetoSearch(baseline(),subset_size=5)
    centre=search.centres[2]
    lower,upper=np.maximum(centre-.05,0),np.minimum(centre+.05,1)
    pool,distance=search.local_pool(lower,upper,42)
    assert len(pool)>0 and (distance>1e-7).all()
    assert (pool>=lower).all() and (pool<=upper).all()
    assert (100*pool[:,:2].sum(1)<=search.max_lift).all()
