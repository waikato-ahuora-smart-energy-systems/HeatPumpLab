"""Production BO contract: joint geometry, exact budgets and failure handling."""

from dataclasses import replace

import numpy as np
import pytest
from botorch.exceptions.errors import ModelFittingError

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._bayesian import (
    BayesianDesignSearch,
    DesignDomain,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
    OuterDesignOptimisationConfig as Config,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import outer_design
from tests.test_outer_design_optimisation import SyntheticAnalyser
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import (
    _OuterDesignObjective,
)


@pytest.mark.parametrize("n", [1, 2, 3, 5])
@pytest.mark.parametrize("fixed", [False, True])
def test_joint_domain_roundtrip_and_boundary_proposals(n, fixed):
    config = Config(
        vary_dT_rhx1_cold=not fixed,
        fixed_dT_rhx1_cold=17.0,
        vary_composition=not fixed,
        initial_evaluations=8,
        candidate_pool_size=64,
    )
    fractions = np.full(n, 1 / n)
    domain = DesignDomain(fractions, 70.0, config)
    search = BayesianDesignSearch(domain, config)
    for _ in range(20):
        x = search.ask()
        assert np.allclose(domain.decode(domain.encode(x)), x)
        assert min(x[3:]) >= -1e-12
        assert sum(x[3:]) == pytest.approx(1.0)
        assert sum(x[:2]) <= 70.0 + 1e-10
        if fixed:
            assert x[0] == 17.0
            assert x[3:] == pytest.approx(fractions)
        # All failures must still explore, without fitting invalid penalty values.
        search.tell(x, 1e6, False)
    assert len({tuple(x) for x in search.x}) == 20


def test_composition_vertices_are_available_without_material_repair():
    domain = DesignDomain([0.2, 0.3, 0.5], 123.0, Config())
    search = BayesianDesignSearch(domain, Config())
    initial = [domain.decode(u)[3:] for u in search.initial]
    for vertex in np.eye(3):
        assert any(np.array_equal(vertex, p) for p in initial)


def test_all_fixed_and_zero_lift_domain():
    config = Config(
        vary_dT_rhx1_cold=False,
        vary_dT_rhx2_cold=False,
        vary_comp2_power_frac=False,
        vary_composition=False,
    )
    domain = DesignDomain([0.3, 0.7], 0.0, config)
    search = BayesianDesignSearch(domain, config)
    x = search.ask()
    search.tell(x, 1.0, True)
    assert search.ask() is None
    constrained = DesignDomain([1.0], 0.0, Config(vary_comp2_power_frac=False))
    assert constrained.dimension == 0
    with pytest.raises(ValueError, match="combined lift"):
        DesignDomain([1.0], 10.0, replace(config, fixed_dT_rhx1_cold=11.0))


def test_joint_gp_improves_coupled_objective_and_repeats():
    config = Config(
        initial_evaluations=6, candidate_pool_size=64, acquisition_restarts=2
    )
    domain = DesignDomain([0.5, 0.5], 123.0, config)

    def run():
        search = BayesianDesignSearch(domain, config)
        for _ in range(22):
            x = search.ask()
            u = domain.encode(x)
            # The loading optimum changes with composition, requiring joint movement.
            y = (
                (u[0] - 0.3) ** 2
                + (u[1] - 0.4) ** 2
                + (u[2] - 0.2 - 0.5 * u[3]) ** 2
                + (u[3] - 0.75) ** 2
            )
            search.tell(x, y, True)
        return search

    a, b = run(), run()
    assert min(a.y) < a.y[0] * 0.12
    assert np.allclose(a.x, b.x, atol=1e-10)
    assert not a.diagnostics or all(d["stage"] == "warning" for d in a.diagnostics)


def test_gp_fitting_failure_returns_distinct_feasible_candidate(monkeypatch):
    config = Config(initial_evaluations=4, candidate_pool_size=32)
    search = BayesianDesignSearch(DesignDomain([1.0], 123.0, config), config)
    for _ in range(4):
        x = search.ask()
        search.tell(x, sum(x[:2]) + 1, True)

    def fail():
        raise ModelFittingError("synthetic fitting failure")

    monkeypatch.setattr(search, "_fit", fail)
    x = search.ask()
    assert not search.observed(x)
    assert search.domain.contains(search.domain.encode(x))
    assert search.diagnostics[-1]["stage"] == "GP fallback"


def test_acquisition_candidate_failure_falls_back_to_pool(monkeypatch):
    from botorch.exceptions.errors import CandidateGenerationError

    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _bayesian

    config = Config(initial_evaluations=4, candidate_pool_size=32)
    search = BayesianDesignSearch(DesignDomain([1.0], 123.0, config), config)
    for _ in range(4):
        x = search.ask()
        search.tell(x, sum(x[:2]) + 1, True)

    def fail(*args, **kwargs):
        raise CandidateGenerationError("all restarts infeasible")

    monkeypatch.setattr(_bayesian, "optimize_acqf", fail)
    x = search.ask()
    assert not search.observed(x)
    assert search.domain.contains(search.domain.encode(x))
    assert any(d["stage"] == "acquisition" for d in search.diagnostics)


@pytest.mark.parametrize(
    "budget,fail", [(0, False), (1, False), (12, False), (12, True)]
)
def test_production_budget_and_resource_cleanup(monkeypatch, budget, fail):
    analyser = SyntheticAnalyser(fail=fail)
    instances = []

    def factory(**kwargs):
        objective = _OuterDesignObjective(**kwargs, analyser=analyser)
        objective.closed = False

        def close():
            objective.closed = True

        objective.close = close
        instances.append(objective)
        return objective

    monkeypatch.setattr(outer_design, "_OuterDesignObjective", factory)
    result = outer_design.optimise_recuperative_heat_pump(
        config=Config(max_exact_evaluations=budget)
    )
    assert result.exact_evaluations == budget
    assert len(analyser.calls) == budget
    assert all(o.closed for o in instances)
    assert result.failed_evaluations == (budget if fail else 0)
    assert (result.best is None) == (fail or budget == 0)
    if budget:
        assert result.evaluations[0].phase == "origin"
        assert result.evaluations[0].x[:3] == (0.0, 0.0, 0.5)
    assert "line_budget" not in result.to_summary_dict()


def test_legacy_configuration_is_rejected():
    with pytest.raises(TypeError):
        Config(line_budget=8)
    assert Config().max_exact_evaluations == 200
    assert "botorch" in outer_design.OPTIMISER_VERSION


def test_single_free_coordinate_uses_bo_and_honours_fixed_inputs():
    config = Config(
        vary_dT_rhx1_cold=False,
        vary_dT_rhx2_cold=False,
        fixed_dT_rhx1_cold=12.0,
        fixed_dT_rhx2_cold=8.0,
        vary_composition=False,
        initial_evaluations=4,
        candidate_pool_size=32,
        acquisition_restarts=2,
    )
    search = BayesianDesignSearch(DesignDomain([0.3, 0.7], 123.0, config), config)
    for _ in range(12):
        x = search.ask()
        assert x[:2] == pytest.approx([12.0, 8.0])
        assert x[3:] == pytest.approx([0.3, 0.7])
        search.tell(x, (x[2] - 0.53) ** 2, True)
    assert min(search.y) < 1e-4


def test_historical_candidate_cannot_shift_fixed_composition_by_repair(monkeypatch):
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
        HistoricalIncumbent,
    )

    analyser = SyntheticAnalyser()
    monkeypatch.setattr(
        outer_design,
        "_OuterDesignObjective",
        lambda **kwargs: _OuterDesignObjective(**kwargs, analyser=analyser),
    )
    result = outer_design.optimise_recuperative_heat_pump(
        config=Config(max_exact_evaluations=2, vary_composition=False),
        incumbents=(
            HistoricalIncumbent(
                "old", "old.csv", 0.5, (10.0, 10.0, 0.4, 0.500000001, 0.499999999)
            ),
        ),
    )
    assert not result.incumbent_replays
    assert result.exact_evaluations == 2
    assert all(r.x[3:] == (0.5, 0.5) for r in result.evaluations)
