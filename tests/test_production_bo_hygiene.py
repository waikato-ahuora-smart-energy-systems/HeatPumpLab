"""Phase 5: production-BO initial design, lift-constraint scaling, and constraint-GP warm starts."""
import numpy as np
import pytest
import torch

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._bayesian import (
    BayesianDesignSearch,
    DesignDomain,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._constraint_model import ConstraintModel
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
    OuterDesignOptimisationConfig as Config,
)
from rhp_local_pareto_search import max_lift_for


def test_initial_design_is_spread_out_and_keeps_origin_probes_and_vertices():
    config = Config()
    domain = DesignDomain([0.2, 0.3, 0.5], 123.0, config)
    search = BayesianDesignSearch(domain, config)
    u = np.asarray(search.initial)
    assert len(u) == config.initial_evaluations
    assert np.allclose(u[0], domain.encode(domain.base))
    near_origin = np.linalg.norm(u - u[0], axis=1) <= 0.05 + 1e-12
    assert near_origin.sum() <= 1 + domain.dimension
    for vertex in np.eye(3):
        assert any(np.allclose(domain.decode(p)[3:], vertex) for p in u)
    assert all(domain.contains(p, tolerance=0) for p in u)
    # Some of the fill lies on a zero-lift face.
    assert ((u[:, 0] == 0) ^ (u[:, 1] == 0)).any()


def test_lift_constraint_matches_domain_membership():
    config = Config(vary_composition=False)
    domain = DesignDomain([1.0], 90.0, config)
    (indices, coefficients, rhs), = domain.constraints()
    for lifts in ([40., 50.], [60., 30.1], [0., 90.], [89., 2.]):
        x = domain.base.copy()
        x[:2] = lifts
        u = domain.encode(x)
        inside = float((coefficients * torch.tensor(u[indices.numpy()])).sum()) >= rhs - 1e-12
        assert inside == (sum(lifts) <= 90.0)


def test_constraint_model_warm_start_skips_refit_and_reuses_hyperparameters():
    rng = np.random.default_rng(1)
    x = rng.random((20, 2))
    m = np.c_[0.5 - x[:, 0], np.full(20, 2.0), np.full(20, 2.0)]
    first = ConstraintModel(x, m)
    warm = ConstraintModel(x, m, fit=False, state=first.state())
    assert first.fitted and not warm.fitted
    np.testing.assert_allclose(warm.probability([[0.2, 0.5], [0.8, 0.5]]),
                               first.probability([[0.2, 0.5], [0.8, 0.5]]), rtol=1e-8)
    # A changed set of active margin groups cannot reuse the old state.
    m2 = np.c_[0.5 - x[:, 0], 0.5 - x[:, 1], np.full(20, 2.0)]
    assert ConstraintModel(x, m2, fit=False, state=first.state()).fitted


def test_max_lift_follows_the_base_case():
    assert max_lift_for({}) == pytest.approx((200 - 15) * 2 / 3)
    assert max_lift_for({'base_case': {'T_air_tar': 150.0, 'T_air_inlet': 15.0,
                                       'dT_rhx3_cold': 0.0}}) == pytest.approx(90.0)
