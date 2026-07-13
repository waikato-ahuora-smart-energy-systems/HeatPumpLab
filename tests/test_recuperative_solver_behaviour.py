from types import SimpleNamespace

import numpy as np
import pytest

from HeatPumpLab.cycles.common import CycleFailure
from HeatPumpLab.cycles.RecuperativeHeatPump.solver import cycle_solver
from HeatPumpLab.cycles.RecuperativeHeatPump.solver import (
    _cycle_analysis as cycle_analysis,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.solver import (
    _discharge_pressure as discharge_pressure,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.solver import (
    _evaporator_alignment as evaporator_alignment,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.solver import (
    _gas_cooler_pinch as gas_cooler_pinch,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.solver import _metrics as metrics


def test_expected_failure_marks_case_infeasible_and_records_report(monkeypatch, capsys):
    recorded = []

    class FailedHeatPump:
        debug = True
        record_failed_solve = True

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
        lambda solved_hp: CycleFailure("no root", {"pinch": 2.0}),
    )
    monkeypatch.setattr(
        cycle_solver,
        "print_result_to_csv",
        lambda output, file_name: recorded.append((output, file_name)),
    )

    result = cycle_solver.solve_heat_pump_cycle(hp)

    assert result is hp
    assert hp.feasible is False
    assert hp.converged is False
    assert hp.solve_error == "no root"
    assert hp.constraint_errors == {"pinch": 2.0}
    assert hp.plots == [((), {})]
    assert recorded == [({"case": "failed"}, "Cases-Unsolved")]
    assert "Heat pump solver failed." in capsys.readouterr().out


def test_unexpected_exception_is_caught_without_physical_diagnostics(
    monkeypatch,
    capsys,
):
    hp = SimpleNamespace(debug=False, record_failed_solve=False)
    monkeypatch.setattr(
        cycle_solver,
        "_solve_case",
        lambda solved_hp: (_ for _ in ()).throw(RuntimeError("property failure")),
    )

    result = cycle_solver.solve_heat_pump_cycle(hp)

    assert result.constraint_errors == {}
    assert result.solve_error == "property failure"
    assert "unexpected exception" in capsys.readouterr().out


def test_fixed_pressure_objective_uses_seed_and_objective_metrics(monkeypatch):
    calls = []
    hp = SimpleNamespace(
        debug=True,
        P_H_est=4.2,
        fixed_discharge_pressure_objective=True,
        get_cycle_plot=lambda: calls.append("plot"),
    )
    monkeypatch.setattr(
        cycle_solver,
        "solve_at_discharge_pressure",
        lambda pressure, solved_hp: calls.append(("fixed", pressure)),
    )
    monkeypatch.setattr(
        cycle_solver,
        "solve_at_minimum_discharge_pressure",
        lambda solved_hp: pytest.fail("fixed objective must not search pressure"),
    )
    monkeypatch.setattr(
        cycle_solver,
        "populate_objective_metrics",
        lambda solved_hp: calls.append("objective_metrics"),
    )

    assert cycle_solver.solve_heat_pump_cycle(hp) is hp
    assert calls == [("fixed", 4.2), "plot", "objective_metrics"]


def test_initial_evaporator_outlet_reports_superheater_reversal():
    hp = SimpleNamespace(
        tol=1e-6,
        Tref_evap_outlet=36.0,
        evap_sup=SimpleNamespace(sink_out=SimpleNamespace(temperature=35.0)),
    )

    failure = evaporator_alignment.validate_initial_evaporator_outlet(hp)

    assert failure is not None
    assert failure.constraint_errors == {
        "superheater_temperature_reversal": pytest.approx(1.0)
    }


def test_objective_metrics_require_only_cycle_energy_inputs():
    hp = SimpleNamespace(
        feasible=True,
        iteration_converged=True,
        operations=[SimpleNamespace(converged=True)],
        gas_coolers=[SimpleNamespace(duty=10.0)],
        compressors=[SimpleNamespace(power=2.0)],
    )

    metrics.populate_objective_metrics(hp)

    assert hp.COP == pytest.approx(5.0)
    assert hp.SEC == pytest.approx(0.25)


def test_sequential_cycle_analysis_reports_iteration_limit(monkeypatch):
    compressor = SimpleNamespace(power=1.0)
    operation = SimpleNamespace(converged=True)

    def solve_operation():
        compressor.power += 1.0

    operation.solve = solve_operation
    hp = SimpleNamespace(
        tol=1e-6,
        debug=False,
        comp_discharge=SimpleNamespace(pressure=3.0),
        exp_valve=SimpleNamespace(flow_in=SimpleNamespace(enthalpy=100.0)),
        gas_cooler1=SimpleNamespace(source_out=SimpleNamespace(pressure=2.0)),
        compressors=[compressor],
        gas_coolers=[SimpleNamespace(duty=10.0)],
        operations=[operation],
    )
    monkeypatch.setattr(
        cycle_analysis,
        "set_gas_cooler_sink_temperatures",
        lambda solved_hp: solved_hp,
    )

    result = cycle_analysis.run_sequential_cycle_analysis(3.0, hp)

    assert result.iteration_converged is False
    assert result.converged is False


def test_final_pressure_solution_returns_component_convergence_failure(monkeypatch):
    hp = SimpleNamespace(
        feasible=True,
        iteration_converged=True,
        compressors=[SimpleNamespace(power=1.0)],
        operations=[SimpleNamespace(converged=False)],
    )
    monkeypatch.setattr(
        discharge_pressure,
        "run_sequential_cycle_analysis",
        lambda pressure, solved_hp: solved_hp,
    )
    monkeypatch.setattr(
        discharge_pressure,
        "update_gc_pinch_properties",
        lambda pressure, solved_hp: solved_hp,
    )

    failure = discharge_pressure.solve_at_discharge_pressure(3.0, hp)

    assert failure is not None
    assert failure.constraint_errors == {"component_convergence": 1.0}


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
        evaporator_alignment,
        "solve_at_minimum_discharge_pressure",
        lambda solved_hp: pytest.fail("already aligned"),
    )

    failure = evaporator_alignment.align_evaporator_inlet_temperature(hp)

    assert failure is None
    assert hp.Tref_evap_outlet == pytest.approx(21.0)


def test_evaporator_alignment_stops_at_superheater_limit(monkeypatch):
    class SinkOutlet:
        temperature = 20.0

        def props_in_ab(self, input_pair, temperature, quality):
            self.temperature = temperature

    hp = SimpleNamespace(
        tol=1e-6,
        Tref_evap_outlet=20.0,
        Tref_evap_inlet=10.0,
        comp_discharge=SimpleNamespace(pressure_bnds=(1.0, 2.0)),
        evaporator=SimpleNamespace(
            sink_in=SimpleNamespace(temperature=0.0),
            sink_out=SinkOutlet(),
        ),
        evap_sup=SimpleNamespace(sink_out=SimpleNamespace(temperature=25.0)),
    )
    searches = []
    monkeypatch.setattr(
        evaporator_alignment,
        "solve_at_minimum_discharge_pressure",
        lambda solved_hp: searches.append(solved_hp.Tref_evap_outlet),
    )

    failure = evaporator_alignment.align_evaporator_inlet_temperature(hp)

    assert failure is None
    assert hp.Tref_evap_outlet == pytest.approx(24.9999)
    assert searches == pytest.approx([24.9999])


def test_evaporator_alignment_exhaustion_returns_failure(monkeypatch):
    class SinkOutlet:
        temperature = 20.0

        def props_in_ab(self, input_pair, temperature, quality):
            self.temperature = temperature

    hp = SimpleNamespace(
        tol=1e-6,
        Tref_evap_outlet=20.0,
        Tref_evap_inlet=10.0,
        comp_discharge=SimpleNamespace(pressure_bnds=(1.0, 2.0)),
        evaporator=SimpleNamespace(
            sink_in=SimpleNamespace(temperature=12.0),
            sink_out=SinkOutlet(),
        ),
    )
    monkeypatch.setattr(
        evaporator_alignment,
        "solve_at_minimum_discharge_pressure",
        lambda solved_hp: None,
    )

    failure = evaporator_alignment.align_evaporator_inlet_temperature(hp)

    assert failure is not None
    assert failure.message == (
        "Evaporator inlet temperature correction failed to converge."
    )
    assert failure.constraint_errors == {}


def test_gas_cooler_midpoint_search_skips_nonpositive_duty():
    gas_cooler = SimpleNamespace(
        duty=0.0,
        source_pinch=SimpleNamespace(),
        source_out=SimpleNamespace(pressure=2.0),
        sink_pinch=SimpleNamespace(),
        sink_in=SimpleNamespace(pressure=1.0),
    )

    approach = gas_cooler_pinch.find_internal_gc_approach(gas_cooler)

    assert approach == (np.inf, 0.0, False)


def test_gas_cooler_sink_temperatures_follow_total_cycle_duty():
    class Stream:
        def __init__(self, enthalpy):
            self.enthalpy = enthalpy
            self.temperature = 0.0
            self.mass_flow = 0.0
            self.heat_flow = 0.0

        def props_in_mPT(self):
            self.enthalpy = 220.0

        def props_in_mPH(self):
            pass

    first = SimpleNamespace(
        duty=20.0,
        sink_out=Stream(220.0),
        sink_in=Stream(200.0),
        sink_pinch=Stream(0.0),
    )
    second = SimpleNamespace(
        duty=30.0,
        sink_out=first.sink_in,
        sink_in=Stream(170.0),
        sink_pinch=Stream(0.0),
    )
    hp = SimpleNamespace(T_air_tar=200.0, gas_coolers=[first, second])

    cycle_analysis.set_gas_cooler_sink_temperatures(hp)

    assert first.sink_out.mass_flow == pytest.approx(1.0)
    assert first.sink_in.enthalpy == pytest.approx(200.0)
    assert second.sink_in.enthalpy == pytest.approx(170.0)
