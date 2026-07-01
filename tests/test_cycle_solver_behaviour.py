from types import SimpleNamespace

import numpy as np
import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.solver import cycle_solver


def test_solver_failure_marks_case_infeasible_and_records_report(monkeypatch, capsys):
    recorded = []

    class FailedHeatPump:
        debug = True

        def __init__(self):
            self.plots = []
            self.output = {}

        def get_cycle_plot(self, *args, **kwargs):
            self.plots.append((args, kwargs))

        def get_report(self):
            self.output = {"case": "failed"}
            return self.output

    hp = FailedHeatPump()

    monkeypatch.setattr(
        cycle_solver,
        "_solve_case",
        lambda solved_hp: (_ for _ in ()).throw(RuntimeError("no root")),
    )
    monkeypatch.setattr(
        cycle_solver,
        "print_result_to_csv",
        lambda output, file_name: recorded.append((output, file_name)),
    )

    result = cycle_solver.solve_heat_pump_cycle(hp)

    assert result is hp
    assert hp.converged is False
    assert hp.feasible is False
    assert hp.solve_error == "no root"
    assert hp.plots == [((), {})]
    assert recorded == [({"case": "failed"}, "Cases-Unsolved")]
    assert "Heat pump solver failed with an exception" in capsys.readouterr().out


def test_solver_debug_mode_renders_cycle_plot_before_metrics(monkeypatch):
    calls = []
    hp = SimpleNamespace(debug=True, get_cycle_plot=lambda: calls.append("plot"))

    monkeypatch.setattr(
        cycle_solver,
        "_find_hp_with_min_discharge_pressure",
        lambda solved_hp: calls.append("pressure") or solved_hp,
    )
    monkeypatch.setattr(
        cycle_solver,
        "_get_evaporator_inlet_temperature_offset",
        lambda solved_hp: 0.0,
    )
    monkeypatch.setattr(
        cycle_solver,
        "_get_metrics",
        lambda solved_hp: calls.append("metrics") or solved_hp,
    )

    assert cycle_solver._solve_case(hp) is hp
    assert calls == ["pressure", "plot", "metrics"]


def test_fixed_pressure_objective_solves_at_seed_pressure_without_pressure_search(
    monkeypatch,
):
    calls = []
    hp = SimpleNamespace(
        debug=True,
        P_H_est=4.2,
        fixed_discharge_pressure_objective=True,
        get_cycle_plot=lambda: calls.append("plot"),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_find_hp_with_min_discharge_pressure",
        lambda solved_hp: pytest.fail("fixed-pressure objective should not search"),
    )
    monkeypatch.setattr(
        cycle_solver,
        "_solve_cycle_at_discharge_pressure",
        lambda pressure, solved_hp: calls.append(("fixed", pressure)) or solved_hp,
    )
    monkeypatch.setattr(
        cycle_solver,
        "_get_evaporator_inlet_temperature_offset",
        lambda solved_hp: pytest.fail("fixed-pressure objective should not align"),
    )
    monkeypatch.setattr(
        cycle_solver,
        "_get_objective_metrics",
        lambda solved_hp: calls.append("objective_metrics") or solved_hp,
    )

    assert cycle_solver._solve_case(hp) is hp
    assert calls == [("fixed", 4.2), "plot", "objective_metrics"]


def test_fixed_pressure_objective_metrics_only_require_sec_inputs():
    hp = SimpleNamespace(
        converged=True,
        feasible=True,
        gas_coolers=[SimpleNamespace(duty=10.0)],
        compressors=[SimpleNamespace(power=2.0)],
    )

    cycle_solver._get_objective_metrics(hp)

    assert hp.COP == pytest.approx(5.0)
    assert hp.SEC == pytest.approx(0.25)


@pytest.mark.parametrize("compute_jac", [False, True])
def test_discharge_pressure_bounds_reuse_cached_bounds(monkeypatch, compute_jac):
    observed_pressures = []
    hp = SimpleNamespace(
        compute_jac=compute_jac,
        comp_discharge=SimpleNamespace(pressure=10.0, pressure_bnds=(8.0, 12.0)),
    )

    def offset(pressure, solved_hp):
        observed_pressures.append(pressure)
        return -1.0 if pressure < 10.0 else 1.0

    monkeypatch.setattr(cycle_solver, "_find_current_gc_dtmin_offset", offset)

    bounds, x0, feasible, dt_values = cycle_solver._get_discharge_pressure_bounds(hp)

    expected_bounds = (9.95, 10.049999999999999) if compute_jac else (8.0, 12.0)
    assert bounds == pytest.approx(expected_bounds)
    assert x0 == pytest.approx(10.0)
    assert feasible is True
    assert dt_values == pytest.approx((1.0, -1.0))
    assert observed_pressures == pytest.approx(list(expected_bounds))


def test_discharge_pressure_bounds_contract_until_sign_changes(monkeypatch):
    hp = SimpleNamespace(
        compute_jac=False,
        comp1=SimpleNamespace(flow_in=SimpleNamespace(pressure=1.0)),
        comp_discharge=SimpleNamespace(pressure=10.0, pressure_bnds=None),
    )

    def offset(pressure, solved_hp):
        return -1.0 if pressure < 9.5 else 1.0

    monkeypatch.setattr(cycle_solver, "_find_current_gc_dtmin_offset", offset)

    bounds, x0, feasible, dt_values = cycle_solver._get_discharge_pressure_bounds(hp)

    assert bounds[0] < 9.5
    assert bounds[1] == pytest.approx(10.0)
    assert x0 == pytest.approx(10.0)
    assert feasible is True
    assert dt_values == pytest.approx((1.0, -1.0))


def test_minimum_discharge_pressure_reports_infeasible_bounds(monkeypatch):
    recorded = []
    hp = SimpleNamespace(
        comp_discharge=SimpleNamespace(pressure_bnds=None),
        get_report=lambda: setattr(hp, "output", {"case": "infeasible"}),
        output={},
    )

    monkeypatch.setattr(
        cycle_solver,
        "_get_discharge_pressure_bounds",
        lambda solved_hp: ((1.0, 2.0), 1.5, False, (1.0, 1.0)),
    )
    monkeypatch.setattr(
        cycle_solver,
        "print_result_to_csv",
        lambda output, file_name: recorded.append((output, file_name)),
    )

    with pytest.raises(Exception, match="Infeasible pressure bounds"):
        cycle_solver._find_hp_with_min_discharge_pressure(hp)

    assert hp.comp_discharge.pressure_bnds == pytest.approx((1.0, 2.0))
    assert recorded == [({"case": "infeasible"}, "Cases-Unsolved")]


@pytest.mark.parametrize(
    ("dt_bounds", "expected_pressure"),
    [
        ((1.0, 0.049), 1.0),
        ((0.049, 1.0), 2.0),
    ],
)
def test_minimum_discharge_pressure_accepts_endpoint_within_target_tolerance(
    monkeypatch, dt_bounds, expected_pressure
):
    solved_pressures = []
    hp = SimpleNamespace(
        tol=1e-6,
        feasible=True,
        comp_discharge=SimpleNamespace(pressure_bnds=None),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_get_discharge_pressure_bounds",
        lambda solved_hp: ((1.0, 2.0), 1.5, True, dt_bounds),
    )
    monkeypatch.setattr(
        cycle_solver,
        "root_scalar",
        lambda **kwargs: pytest.fail("endpoint already satisfies the pressure target"),
    )
    monkeypatch.setattr(
        cycle_solver,
        "_do_sequential_iterative_cycle_analysis",
        lambda pressure, solved_hp: solved_pressures.append(pressure) or solved_hp,
    )
    monkeypatch.setattr(
        cycle_solver,
        "_update_pinch_props_for_gcs",
        lambda P_comp_out, hp: hp,
    )

    result = cycle_solver._find_hp_with_min_discharge_pressure(hp)

    assert result is hp
    assert hp.converged is True
    assert hp.comp_discharge.pressure_bnds == pytest.approx((1.0, 2.0))
    assert solved_pressures == pytest.approx([expected_pressure])


def test_minimum_discharge_pressure_prefers_bracketed_root_to_endpoint_tolerance(
    monkeypatch,
):
    solved_pressures = []
    hp = SimpleNamespace(
        tol=1e-6,
        feasible=True,
        comp_discharge=SimpleNamespace(pressure_bnds=None),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_get_discharge_pressure_bounds",
        lambda solved_hp: ((1.0, 2.0), 1.5, True, (0.049, -1.0)),
    )
    monkeypatch.setattr(
        cycle_solver,
        "root_scalar",
        lambda **kwargs: SimpleNamespace(converged=True, root=1.25),
    )
    monkeypatch.setattr(
        cycle_solver,
        "_do_sequential_iterative_cycle_analysis",
        lambda pressure, solved_hp: solved_pressures.append(pressure) or solved_hp,
    )
    monkeypatch.setattr(
        cycle_solver,
        "_update_pinch_props_for_gcs",
        lambda P_comp_out, hp: hp,
    )

    result = cycle_solver._find_hp_with_min_discharge_pressure(hp)

    assert result is hp
    assert hp.converged is True
    assert solved_pressures == pytest.approx([1.25])


def test_multicomponent_pressure_search_uses_lowest_pressure_within_target(
    monkeypatch,
):
    hp = SimpleNamespace(
        tol=1e-6,
        refrigerant="Propane;Butane|0.5;0.5|mass",
    )
    sampled_pressures = []

    monkeypatch.setattr(
        cycle_solver,
        "root_scalar",
        lambda **kwargs: pytest.fail("mixture pressure search should stop by target"),
    )

    def offset(pressure, solved_hp):
        sampled_pressures.append(pressure)
        return -0.2 if pressure < 1.2 else 0.01

    monkeypatch.setattr(cycle_solver, "_find_current_gc_dtmin_offset", offset)

    pressure, converged = cycle_solver._find_pressure_for_gc_target(
        (1.0, 2.0), 1.5, (0.01, -1.0), hp
    )

    assert converged is True
    assert 1.0 < pressure < 2.0
    assert pressure == pytest.approx(sampled_pressures[-1])


def test_minimum_discharge_pressure_falls_back_to_near_target_endpoint(
    monkeypatch,
):
    solved_pressures = []
    hp = SimpleNamespace(
        tol=1e-6,
        feasible=True,
        comp_discharge=SimpleNamespace(pressure_bnds=None),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_get_discharge_pressure_bounds",
        lambda solved_hp: ((1.0, 2.0), 1.5, True, (1.0, -0.049)),
    )
    monkeypatch.setattr(
        cycle_solver,
        "root_scalar",
        lambda **kwargs: (_ for _ in ()).throw(ValueError("stale bracket")),
    )
    monkeypatch.setattr(
        cycle_solver,
        "_do_sequential_iterative_cycle_analysis",
        lambda pressure, solved_hp: solved_pressures.append(pressure) or solved_hp,
    )
    monkeypatch.setattr(
        cycle_solver,
        "_update_pinch_props_for_gcs",
        lambda P_comp_out, hp: hp,
    )

    result = cycle_solver._find_hp_with_min_discharge_pressure(hp)

    assert result is hp
    assert hp.converged is True
    assert solved_pressures == pytest.approx([1.0])


def test_evaporator_alignment_accepts_domain_temperature_tolerance(monkeypatch):
    hp = SimpleNamespace(
        tol=1e-6,
        Tref_evap_inlet=10.0,
        Tref_evap_outlet=0.0,
        evaporator=SimpleNamespace(
            sink_in=SimpleNamespace(temperature=10.00005),
            sink_out=SimpleNamespace(temperature=21.0),
        ),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_find_hp_with_min_discharge_pressure",
        lambda solved_hp: pytest.fail("evaporator is already aligned"),
    )

    result = cycle_solver._align_evaporator_inlet_temperature(hp)

    assert result is hp
    assert hp.Tref_evap_outlet == pytest.approx(21.0)


def test_discharge_pressure_bounds_accept_degenerate_near_target_endpoint(
    monkeypatch,
):
    hp = SimpleNamespace(
        compute_jac=False,
        tol=1e-6,
        comp_discharge=SimpleNamespace(pressure=10.0, pressure_bnds=(10.0, 10.0)),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_find_current_gc_dtmin_offset",
        lambda pressure, solved_hp: 0.02,
    )

    bounds, x0, feasible, dt_values = cycle_solver._get_discharge_pressure_bounds(hp)

    assert bounds == pytest.approx((10.0, 10.0))
    assert x0 == pytest.approx(10.0)
    assert feasible is True
    assert dt_values == pytest.approx((0.02, 0.02))


def test_evaporator_alignment_reports_non_convergence(monkeypatch):
    class SinkOutlet:
        temperature = 20.0

        def props_in_ab(self, *args):
            self.temperature = args[1]

    hp = SimpleNamespace(
        tol=1e-6,
        Tref_evap_inlet=10.0,
        comp_discharge=SimpleNamespace(pressure_bnds=(1.0, 2.0)),
        evaporator=SimpleNamespace(
            sink_in=SimpleNamespace(temperature=12.0),
            sink_out=SinkOutlet(),
        ),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_find_hp_with_min_discharge_pressure",
        lambda solved_hp: solved_hp,
    )

    with pytest.raises(Exception, match="failed to converge"):
        cycle_solver._align_evaporator_inlet_temperature(hp)

    assert hp.comp_discharge.pressure_bnds is None


def test_gas_cooler_dtmin_offset_ignores_inlet_air_pinch_candidates(monkeypatch):
    hp = SimpleNamespace(dT_GC_min_tar=15.0)
    monkeypatch.setattr(
        cycle_solver,
        "_get_all_possible_gc_dtmin_pos",
        lambda pressure, solved_hp: [
            [25.0, 1.0, True],
            [12.0, 0.0, True],
            [20.0, 0.0, True],
        ],
    )

    assert cycle_solver._find_current_gc_dtmin_offset(pressure := 10.0, hp) == (
        pytest.approx(10.0)
    )
    assert cycle_solver._find_current_gc_dtmin_offset(
        pressure, hp, return_abs=True
    ) == pytest.approx(10.0)


def test_all_possible_gas_cooler_pinch_positions_include_debug_plot(monkeypatch):
    plotted = []
    gas_cooler = SimpleNamespace(
        source_in=SimpleNamespace(temperature=120.0),
        sink_out=SimpleNamespace(temperature=90.0),
        source_out=SimpleNamespace(temperature=70.0),
        sink_in=SimpleNamespace(temperature=40.0),
    )
    hp = SimpleNamespace(debug=True, dT_GC_min_tar=15.0, gas_coolers=[gas_cooler])

    monkeypatch.setattr(
        cycle_solver,
        "_do_sequential_iterative_cycle_analysis",
        lambda pressure, solved_hp: solved_hp,
    )
    monkeypatch.setattr(
        cycle_solver,
        "_find_gc_dt_approaches",
        lambda gc: (20.0, 0.4, True),
    )
    monkeypatch.setattr(
        cycle_solver,
        "plot_driving_force",
        lambda gas_coolers, dtmin: plotted.append((gas_coolers, dtmin)),
    )

    positions = cycle_solver._get_all_possible_gc_dtmin_pos(10.0, hp)

    assert positions == [(30.0, 1.0, True), (20.0, 0.4, True), (30.0, 0.0, True)]
    assert plotted == [([gas_cooler], 15.0)]


def test_sequential_cycle_analysis_runs_debug_cycle_plots(monkeypatch):
    plots = []
    hp = SimpleNamespace(
        tol=1e-3,
        debug=True,
        comp_discharge=SimpleNamespace(pressure=0.0),
        exp_valve=SimpleNamespace(flow_in=SimpleNamespace(enthalpy=100.0)),
        gas_cooler1=SimpleNamespace(source_out=SimpleNamespace(pressure=8.0)),
        compressors=[SimpleNamespace(power=5.0)],
        gas_coolers=[SimpleNamespace(duty=20.0)],
        operations=[SimpleNamespace(solve=lambda: None)],
        get_cycle_plot=lambda **kwargs: plots.append(kwargs),
    )

    monkeypatch.setattr(
        cycle_solver,
        "_set_gas_coolers_sink_temperatures",
        lambda solved_hp: solved_hp,
    )

    assert cycle_solver._do_sequential_iterative_cycle_analysis(12.0, hp) is hp
    assert plots == [
        {"x_axis": "s", "y_axis": "T"},
        {"x_axis": "h", "y_axis": "P"},
        {"x_axis": "h", "y_axis": "T"},
    ]


def test_gas_cooler_midpoint_search_is_skipped_when_no_duty_is_available():
    gc = SimpleNamespace(
        duty=0.0,
        source_out=SimpleNamespace(pressure=5.0),
        sink_in=SimpleNamespace(pressure=1.0),
        source_pinch=SimpleNamespace(),
        sink_pinch=SimpleNamespace(),
    )

    dtmin, fraction, success = cycle_solver._find_gc_dt_approaches(gc)

    assert dtmin == np.inf
    assert fraction == pytest.approx(0.0)
    assert success is False


def test_gas_cooler_sink_temperatures_are_back_calculated_from_total_duty():
    class SinkStream:
        def __init__(self, enthalpy):
            self.enthalpy = enthalpy
            self.mass_flow = 0.0
            self.heat_flow = 0.0

        def props_in_mPT(self):
            self.enthalpy = self.temperature

        def props_in_mPH(self):
            pass

    first = SimpleNamespace(
        duty=30.0,
        sink_out=SinkStream(enthalpy=0.0),
        sink_in=SinkStream(enthalpy=0.0),
        sink_pinch=SimpleNamespace(mass_flow=0.0),
    )
    second = SimpleNamespace(
        duty=10.0,
        sink_out=SinkStream(enthalpy=0.0),
        sink_in=SinkStream(enthalpy=20.0),
        sink_pinch=SimpleNamespace(mass_flow=0.0),
    )
    hp = SimpleNamespace(T_air_tar=100.0, gas_coolers=[first, second])

    cycle_solver._set_gas_coolers_sink_temperatures(hp)

    assert first.sink_out.mass_flow == pytest.approx(0.5)
    assert first.sink_out.heat_flow == pytest.approx(50.0)
    assert first.sink_in.mass_flow == pytest.approx(0.5)
    assert first.sink_in.enthalpy == pytest.approx(40.0)
    assert second.sink_in.mass_flow == pytest.approx(0.5)
    assert second.sink_in.enthalpy == pytest.approx(-20.0)
