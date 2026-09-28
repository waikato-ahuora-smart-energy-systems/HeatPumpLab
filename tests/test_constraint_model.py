"""Continuous constraint margins and the constraint-aware acquisitions (Phase 2)."""

import numpy as np
import pytest
import torch

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._bayesian import (
    BayesianDesignSearch,
    DesignDomain,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._constraint_model import (
    ConstraintModel,
    aggregate_margins,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import (
    SUPERCRITICAL_MARGIN_K,
    constraint_margins,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
    OuterDesignOptimisationConfig as Config,
)
from rhp_mobo_search import ParetoBOSearch


def margins(discharge, suction=20.0, pr=3.0):
    return {"stage_1_pressure_ratio": pr, "stage_1_discharge_superheat": discharge,
            "stage_1_suction_superheat": suction, "stage_2_pressure_ratio": pr + 1,
            "stage_2_discharge_superheat": discharge + 5, "stage_2_suction_superheat": suction}


def test_aggregate_takes_worst_stage_scales_and_clips():
    assert aggregate_margins(margins(-5.0)) == pytest.approx([-0.5, 2.0, 3.0])
    assert aggregate_margins(margins(-500.0))[0] == -3.0
    assert aggregate_margins(None) is None
    assert aggregate_margins(dict(margins(1.0), stage_2_suction_superheat=None)) is None


def test_constraint_margins_from_report():
    report = {"PR_1": 7.5, "PR_2": 11.0,
              "Comp_1_dT_sh_in": 4.0, "Comp_1_dT_sh_out": 12.0, "Comp_1_phase_out": "superheated_vapour",
              "Comp_2_dT_sh_in": 3.0, "Comp_2_dT_sh_out": -2.0, "Comp_2_phase_out": "supercritical"}
    m = constraint_margins(report, 2, 10.0, min_discharge_superheat_K=5.0)
    assert m["stage_1_pressure_ratio"] == pytest.approx(2.5)
    assert m["stage_2_pressure_ratio"] == pytest.approx(-1.0)
    assert m["stage_1_discharge_superheat"] == pytest.approx(7.0)
    # The phase check accepts a supercritical discharge; its margin agrees.
    assert m["stage_2_discharge_superheat"] == pytest.approx(SUPERCRITICAL_MARGIN_K)
    assert m["stage_1_suction_superheat"] == pytest.approx(4.0)
    assert constraint_margins({"PR_1": 2.0}, 1, 10.0)["stage_1_discharge_superheat"] is None


def test_constraint_model_separates_feasible_side_with_finite_gradient():
    rng = np.random.default_rng(0)
    x = rng.random((24, 2))
    m = np.c_[0.5 - x[:, 0], np.full(24, 2.0), np.full(24, 2.0)]
    model = ConstraintModel(x, m)
    assert model.active == [0]
    p = model.probability([[0.1, 0.5], [0.9, 0.5]])
    assert p[0] > 0.9 and p[1] < 0.1
    X = torch.tensor([[[0.45, 0.5]]], dtype=torch.double, requires_grad=True)
    model.log_probability(X).sum().backward()
    assert torch.isfinite(X.grad).all()


def run_production(use_margins, steps=40):
    """Minimise SEC, which improves towards a lift that violates the margin.

    The optimum lies on the constraint (u0 = 0.6), as in the real problem, so
    some violations while refining the boundary are expected either way.
    """
    config = Config(vary_composition=False, initial_evaluations=8, candidate_pool_size=128,
                    acquisition_restarts=2, seed=11)
    domain = DesignDomain([1.0], 123.0, config)
    search = BayesianDesignSearch(domain, config)
    log = []
    for step in range(steps):
        x = search.ask()
        u0 = x[0] / 100
        feasible = u0 <= 0.6
        sec = 1.0 - 0.5 * u0 + 0.05 * (x[2] - 0.5) ** 2
        search.tell(x, sec if feasible else 1e6, feasible,
                    margins(30 * (0.6 - u0)) if use_margins else None)
        exploring = (step >= config.initial_evaluations and
                     (step - config.initial_evaluations) % config.exploration_interval
                     == config.exploration_interval - 1)
        log.append((u0, feasible, exploring))
    best = min((y, x) for y, x, ok in zip(search.y, search.x, search.success) if ok)
    return log, domain.decode(best[1])


def test_production_bo_learns_the_constraint_boundary():
    with_margins, best = run_production(True)
    without, _ = run_production(False)
    late = lambda log: [(u0, ok) for u0, ok, exploring in log[16:] if not exploring]
    violations = lambda log: sum(not ok for _, ok in late(log))
    assert violations(with_margins) < violations(without)
    # Guided violations stay close to the boundary instead of deep inside it.
    assert all(u0 <= 0.7 for u0, ok in late(with_margins) if not ok)
    assert 50.0 <= best[0] <= 60.0


def test_mobo_counts_only_marginless_failures_as_numerical():
    search = ParetoBOSearch(123.333, 3, fixed_comp_frac=0.5)
    for i in range(12):
        x = search.ask()
        violated = i % 4 == 1
        numerical = i == 6
        ok = not (violated or numerical)
        search.tell(x, 2.5 + 0.1 * x[3] if ok else None, 3.0 - x[3] if ok else None, ok,
                    None if numerical else margins(-4.0 if violated else 8.0))
    assert search.numerical_failures().tolist() == [i == 6 for i in range(12)]
    assert search.constraint_model() is not None
    x = search.ask()
    assert x[2] == 0.5 and x[0] + x[1] <= 123.333 + 1e-9
