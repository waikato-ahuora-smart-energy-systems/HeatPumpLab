"""Phase 4: zero-lift face coverage, batched proposals and the hypervolume stop rule."""
import numpy as np

from rhp_frontier_stability import assess
from rhp_local_pareto_search import choose_anchors
from rhp_mobo_search import ParetoBOSearch


def record(x, cop, vhc, success=True):
    return dict(x=x, COP=cop, VHC_MJ_per_m3=vhc, success=success)


def test_each_zero_lift_face_gets_a_local_region_even_when_dominated():
    records = [record([10.+i, 20.+i, .5, .2+.03*i, .8-.03*i], 3-.02*i, 1+.2*i) for i in range(12)]
    # A dominated design on the RHX2 = 0 face, and one on the RHX1 = 0 face.
    records.append(record([70., 0., .5, .8, .2], 2.5, 2.0))
    records.append(record([0., 40., .5, .6, .4], 2.6, 1.9))
    anchors = choose_anchors(records)
    assert 12 in anchors and 13 in anchors and len(set(anchors)) == 6


def test_face_restricted_steps_alternate():
    search = ParetoBOSearch(123.333, fixed_comp_frac=.5)
    faces = []
    for n in range(8):
        search.pending = [[0., 0., 0.]]*n
        faces.append(search.face())
    assert faces == [None, 0, None, 1]*2


def test_initial_design_seeds_both_faces():
    search = ParetoBOSearch(123.333, fixed_comp_frac=.5)
    u = np.asarray(search.initial)
    assert ((u[:, 0] == 0) & (u[:, 1] > 0)).sum() >= 4
    assert ((u[:, 1] == 0) & (u[:, 0] > 0)).sum() >= 4


def test_batch_proposals_are_distinct_and_pending_clears_on_tell():
    search = ParetoBOSearch(123.333, 5, fixed_comp_frac=.5)
    for _ in range(40):
        x = search.ask()
        u = search.encode(x)
        search.tell(x, 2+.7*u[2]-.1*u[0]-.05*u[1], 1+2*(1-u[2])+.1*u[1], True)
    batch = search.ask_batch(4)
    encoded = np.asarray([search.encode(x) for x in batch])
    assert len(np.unique(np.round(encoded, 9), axis=0)) == 4
    assert len(search.pending) == 4
    for x in batch:
        assert x[2] == .5 and x[0]+x[1] <= 123.333+1e-9
        search.tell(x, 2.5, 2.5, True)
    assert search.pending == []


def test_hypervolume_only_rule_ignores_curve_changes():
    def block(shift):
        return [record([0, 0, .5, .5, .5], 2.5-.001*i+shift, 1+.01*i) for i in range(100)]
    base = block(0.)
    previous = None
    for k in range(1, 5):
        records = base*k
        previous = assess(records, previous, tolerance=2e-4, minimum_evaluations=300, hv_only=True)
    assert previous['stable_blocks'] == 3 and previous['stop_for_stability']
