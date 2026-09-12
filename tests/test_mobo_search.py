"""Domain, Pareto bookkeeping and checkpoint reproducibility for multi-objective BO."""

import json

import numpy as np
import pytest

from main_RHPmobo import frontier_indices
from rhp_mobo_search import ParetoBOSearch


def observe(search, n=32):
    for _ in range(n):
        x = search.ask()
        u = search.encode(x)
        search.tell(x, 2+.7*u[3]-.1*u[0]-.05*u[1], 1+2*(1-u[3])+.1*u[2], True)


def test_initial_design_includes_pure_endpoints_and_respects_joint_domain():
    search = ParetoBOSearch(123.333)
    observe(search)
    assert len(np.unique(search.x, axis=0)) == 32
    assert {0., 1.}.issubset(set(np.asarray(search.x)[:, 3]))
    for u in search.x:
        x = search.decode(u)
        assert x[0]+x[1] <= 123.333
        assert .1 <= x[2] <= .9 and x[3]+x[4] == 1.


def test_failed_evaluations_cannot_enter_observed_pareto_front():
    def r(cop, vhc, success=True):
        return dict(COP=cop, VHC_MJ_per_m3=vhc, success=success)
    assert frontier_indices([r(2, 4), r(3, 3), r(2, 2), r(100, 100, False),
                             r(3, 3), r(None, None, False)]) == [0, 1, 4]


def test_success_requires_both_positive_objectives():
    search = ParetoBOSearch(123.333)
    search.tell(search.ask(), 2., 0., True)
    assert search.success == [False]
    with pytest.raises(ValueError, match="Duplicate"):
        search.tell(search.decode(search.x[0]), 2., 1., True)


def test_botorch_proposal_and_checkpoint_resume_are_reproducible():
    search = ParetoBOSearch(123.333, 8)
    observe(search)
    proposal = search.ask()
    assert len(search._admissible(search.encode(proposal)[None])) == 1
    assert not any(d.get("stage") == "acquisition fallback" for d in search.diagnostics)
    search.tell(proposal, 2.4, 2., True)
    restored = ParetoBOSearch(123.333, 8)
    for u, values, success in zip(search.x, search.y, search.success):
        restored.tell(restored.decode(u), *values, success)
    restored.kernel_state = json.loads(json.dumps(search.kernel_state, allow_nan=False))
    restored.mean_state = json.loads(json.dumps(search.mean_state, allow_nan=False))
    assert np.allclose(search.ask(), restored.ask(), atol=1e-9, rtol=0)
