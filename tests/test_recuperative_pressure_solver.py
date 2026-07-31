from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.common import ScalarSearchResult
from HeatPumpLab.cycles.RecuperativeHeatPump.solver import (
    _discharge_pressure as pressure_solver,
)


@pytest.mark.parametrize("compute_jac", [False, True])
def test_pressure_bounds_reuse_the_available_local_domain(monkeypatch, compute_jac):
    observed = []
    hp = SimpleNamespace(
        compute_jac=compute_jac,
        tol=1e-6,
        comp_discharge=SimpleNamespace(pressure=10.0, pressure_bnds=(8.0, 12.0)),
    )

    def residual(pressure, solved_hp):
        observed.append(pressure)
        return -1.0 if pressure < 10.0 else 1.0

    monkeypatch.setattr(pressure_solver, "find_gc_pinch_residual", residual)

    bounds, initial, feasible, residuals = (
        pressure_solver.get_discharge_pressure_bounds(hp)
    )

    expected_bounds = (9.95, 10.05) if compute_jac else (8.0, 12.0)
    assert bounds == pytest.approx(expected_bounds)
    assert initial == pytest.approx(10.0)
    assert feasible is True
    assert residuals == pytest.approx((1.0, -1.0))
    assert observed == pytest.approx(list(expected_bounds))


def test_pressure_bound_expansion_stops_at_property_limit(monkeypatch):
    evaluated = []
    hp = SimpleNamespace(
        tol=1e-6,
        comp_discharge=SimpleNamespace(P_max=3.0),
    )

    def residual(pressure, solved_hp):
        evaluated.append(pressure)
        return -1.0

    monkeypatch.setattr(pressure_solver, "find_gc_pinch_residual", residual)

    upper, upper_residual, lower, lower_residual = (
        pressure_solver.expand_discharge_pressure_bounds(
            hp,
            suction_pressure=1.0,
            pressure_start=2.0,
            residual_start=-1.0,
        )
    )

    assert max(evaluated) == pytest.approx(3.0)
    assert upper == pytest.approx(3.0)
    assert lower < upper
    assert (upper_residual, lower_residual) == pytest.approx((-1.0, -1.0))


def test_minimum_pressure_accepts_endpoint_within_target_tolerance(monkeypatch):
    hp = SimpleNamespace(
        tol=1e-6,
        feasible=True,
        comp_discharge=SimpleNamespace(pressure=1.5, pressure_bnds=None),
    )
    solved = []
    monkeypatch.setattr(
        pressure_solver,
        "get_discharge_pressure_bounds",
        lambda solved_hp: ((1.0, 2.0), 1.5, True, (1.0, 0.049)),
    )
    monkeypatch.setattr(
        pressure_solver,
        "solve_at_discharge_pressure",
        lambda pressure, solved_hp: solved.append(pressure),
    )

    failure = pressure_solver.solve_at_minimum_discharge_pressure(hp)

    assert failure is None
    assert hp.feasible is True
    assert solved == pytest.approx([1.0])


def test_infeasible_pressure_returns_structured_shortfall(monkeypatch):
    hp = SimpleNamespace(
        tol=1e-6,
        refrigerant="Butane",
        comp_discharge=SimpleNamespace(pressure=1.5, pressure_bnds=None),
    )
    monkeypatch.setattr(
        pressure_solver,
        "get_discharge_pressure_bounds",
        lambda solved_hp: ((1.0, 2.0), 1.5, False, (-0.75, -0.75)),
    )

    failure = pressure_solver.solve_at_minimum_discharge_pressure(hp)

    assert failure is not None
    assert hp.feasible is False
    assert failure.constraint_errors == {
        "gas_cooler_pinch_shortfall": pytest.approx(0.70)
    }


def test_multicomponent_pressure_uses_sampled_fallback_after_bounds_fail(monkeypatch):
    hp = SimpleNamespace(
        tol=1e-6,
        feasible=True,
        refrigerant="Propane;Isopentane|0.5;0.5|mass",
        comp1=SimpleNamespace(flow_in=SimpleNamespace(pressure=1.0)),
        comp_discharge=SimpleNamespace(
            pressure=4.0,
            pressure_bnds=None,
            P_max=10.0,
        ),
    )
    solved = []
    monkeypatch.setattr(
        pressure_solver,
        "get_discharge_pressure_bounds",
        lambda solved_hp: ((4.0, 10.0), 4.0, False, (-0.8, -0.8)),
    )
    monkeypatch.setattr(
        pressure_solver,
        "find_sampled_log_target",
        lambda *args, **kwargs: ScalarSearchResult(3.0, 0.02, True),
    )
    monkeypatch.setattr(
        pressure_solver,
        "solve_at_discharge_pressure",
        lambda pressure, solved_hp: solved.append(pressure),
    )

    failure = pressure_solver.solve_at_minimum_discharge_pressure(hp)

    assert failure is None
    assert solved == pytest.approx([3.0])


def test_multicomponent_bracket_selects_lowest_accepted_pressure(monkeypatch):
    hp = SimpleNamespace(
        tol=1e-6,
        refrigerant="Propane;Butane|0.5;0.5|mass",
    )
    monkeypatch.setattr(
        pressure_solver,
        "find_gc_pinch_residual",
        lambda pressure, solved_hp: -0.2 if pressure < 1.2 else 0.01,
    )

    pressure, converged = pressure_solver.select_gc_target_pressure(
        (1.0, 2.0),
        1.5,
        (0.01, -0.2),
        hp,
    )

    assert converged is True
    assert pressure == pytest.approx(1.2, rel=1e-3)
