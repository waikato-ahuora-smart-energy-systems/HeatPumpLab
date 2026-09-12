"""Regression cases for invalid pressure probes and bounded numerical recovery."""

from types import SimpleNamespace

import pytest

from HeatPumpLab.components.stream import Stream
from HeatPumpLab.components.common.base_heat_exchanger import BaseHeatExchanger
from HeatPumpLab.cycles.RecuperativeHeatPump.solver._cycle_analysis import (
    InvalidSinkStateError, set_gas_cooler_sink_temperatures,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.utils import analysis
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective


@pytest.mark.parametrize("duties", [(-10., 10.0001), (0., 0.), (-1., -1.), (float("nan"), 2.), (float("inf"), 2.)])
def test_invalid_sink_duty_never_reaches_property_flash(duties):
    def native():
        pytest.fail("Invalid duties must be rejected before a native call")
    coolers = [SimpleNamespace(duty=d, tol=1e-6, sink_out=SimpleNamespace(props_in_mPT=native)) for d in duties]
    with pytest.raises(InvalidSinkStateError) as failure:
        set_gas_cooler_sink_temperatures(SimpleNamespace(gas_coolers=coolers))
    assert failure.value.constraint_errors


@pytest.mark.parametrize("rise", [0., -10., float("nan")])
def test_invalid_sink_enthalpy_rise_never_reaches_ph(rise):
    def ph():
        pytest.fail("Invalid enthalpy rise must not reach PH")
    cooler = SimpleNamespace(duty=10., sink_out=SimpleNamespace(enthalpy=100.+rise, props_in_mPT=lambda: None),
                             sink_in=SimpleNamespace(enthalpy=100., props_in_mPH=ph))
    with pytest.raises(InvalidSinkStateError, match="enthalpy rise"):
        set_gas_cooler_sink_temperatures(SimpleNamespace(gas_coolers=[cooler], T_air_tar=200.))


def test_tiny_tolerated_negative_duty_cannot_amplify_sink_enthalpy():
    def ph():
        pytest.fail("Amplified air enthalpy must be rejected before PH")
    first = SimpleNamespace(duty=-5e-7, tol=1e-6,
                            sink_out=SimpleNamespace(enthalpy=220., props_in_mPT=lambda: None),
                            sink_in=SimpleNamespace(enthalpy=200., props_in_mPH=ph), sink_pinch=SimpleNamespace())
    last = SimpleNamespace(duty=5.00000001e-7, tol=1e-6, sink_in=SimpleNamespace(enthalpy=170.))
    with pytest.raises(InvalidSinkStateError, match="outside"):
        set_gas_cooler_sink_temperatures(SimpleNamespace(gas_coolers=[first, last], T_air_tar=200.))


def make_flow(temperature, enthalpy, pressure):
    flow = Stream()
    flow.temperature, flow.enthalpy, flow.pressure = temperature, enthalpy, pressure
    flow.entropy = flow.heat_capacity = flow.density = flow.specific_volume = 1.
    flow.quality = 998.
    flow.set_flows(1.)
    return flow


@pytest.mark.parametrize("duty", [0., 2e-7, 1e-6])
def test_negligible_terminal_duty_conserves_both_streams_without_ph(duty):
    exchanger = BaseHeatExchanger.__new__(BaseHeatExchanger)
    exchanger.source_in = make_flow(100., 1000., 2.)
    exchanger.sink_in = make_flow(90., 100., 1.)
    exchanger.source_out, exchanger.sink_out = Stream(), Stream()
    exchanger.source_out_max = SimpleNamespace(enthalpy=1000.-duty, props_in_ab=lambda *a: None)
    exchanger.sink_out_max = SimpleNamespace(enthalpy=100.+duty, props_in_ab=lambda *a: None)
    exchanger.dP_h = exchanger.dP_c = 0.
    exchanger.dTmin = 5.
    exchanger.tol = 1e-6
    exchanger.get_dTmin_terminal_model()
    assert exchanger.converged and exchanger.duty == 0.
    for before, after in ((exchanger.source_in, exchanger.source_out), (exchanger.sink_in, exchanger.sink_out)):
        assert (after.temperature, after.pressure, after.enthalpy, after.mass_flow, after.heat_flow) == (
            before.temperature, before.pressure, before.enthalpy, before.mass_flow, before.heat_flow)


@pytest.mark.parametrize("duty,pressure_drop", [(-2e-7, 0.), (2e-7, .01), (2e-6, 0.)])
def test_zero_duty_treatment_does_not_hide_reversal_pressure_drop_or_finite_heat(duty, pressure_drop):
    exchanger = BaseHeatExchanger.__new__(BaseHeatExchanger)
    exchanger.name = "test"
    exchanger.source_in = make_flow(100., 1000., 2.)
    exchanger.sink_in = make_flow(90., 100., 1.)
    calls = []
    def flash():
        calls.append(1)
        raise RuntimeError("reached normal property path")
    exchanger.source_out = SimpleNamespace(props_in_mPH=flash)
    exchanger.sink_out = Stream()
    exchanger.source_out_max = SimpleNamespace(enthalpy=1000.-duty, props_in_ab=lambda *a: None)
    exchanger.sink_out_max = SimpleNamespace(enthalpy=100.+duty, props_in_ab=lambda *a: None)
    exchanger.dP_h, exchanger.dP_c = pressure_drop, 0.
    exchanger.dTmin, exchanger.tol = 5., 1e-6
    with pytest.raises(ValueError):
        exchanger.get_dTmin_terminal_model()
    assert calls == [1]


def test_negligible_fixed_temperature_rise_conserves_energy():
    exchanger = BaseHeatExchanger.__new__(BaseHeatExchanger)
    exchanger.source_in = make_flow(100., 1000., 2.)
    exchanger.sink_in = make_flow(90., 100., 1.)
    exchanger.source_out, exchanger.sink_out = Stream(), make_flow(90., 100., 1.)
    def pt():
        exchanger.sink_out.enthalpy = 100.+exchanger.sink_out.temperature-90.
        exchanger.sink_out.set_flows(1.)
    exchanger.sink_out.props_in_mPT = pt
    exchanger.dP_h = exchanger.dP_c = 0.
    exchanger.dT_c, exchanger.tol = 2e-7, 1e-6
    exchanger.get_fixed_cold_delta_T_model()
    assert exchanger.converged and exchanger.duty == 0.
    assert exchanger.source_out.heat_flow == exchanger.source_in.heat_flow
    assert exchanger.sink_out.heat_flow == exchanger.sink_in.heat_flow
    assert abs(exchanger.sink_out.temperature-exchanger.sink_in.temperature-exchanger.dT_c) <= exchanger.tol


def failed_model(kind="native_property", message="PH [226]"):
    return SimpleNamespace(feasible=False, converged=False, solve_failure_kind=kind,
                           constraint_errors={}, numerical_retries=[], solve_error=message)


@pytest.mark.parametrize("seed,expected", [(0., [0., 3., 7.]), (3., [3., 7.]), (7., [7., 3.])])
def test_pressure_retries_are_bounded_and_preserve_original_cause(monkeypatch, seed, expected):
    seen = []
    def once(cls, case, *args):
        seen.append(case.get("P_H_est", 0.))
        model = failed_model(message=f"PH failed at {seen[-1]}")
        model.numerical_retries = case.get("numerical_retry_history", [])
        return model, 1.
    monkeypatch.setattr(analysis, "_analyse_RHP_case_once", once)
    hp, elapsed = analysis.analyse_RHP_case(None, {"refrigerant": "Propane;Pentane", "P_H_est": seed}, record_result=False)
    assert seen == expected and elapsed == len(expected)
    assert hp.solve_error == f"PH failed at {seed}"
    assert [r["P_H_est"] for r in hp.numerical_retries] == expected[1:]


def test_retry_stops_at_first_valid_result_and_retains_physical_settings(monkeypatch):
    seen = []
    def once(cls, case, *args):
        seen.append(case)
        if len(seen) == 1:
            return failed_model(), 1.
        return SimpleNamespace(feasible=True, converged=True), 2.
    monkeypatch.setattr(analysis, "_analyse_RHP_case_once", once)
    case = dict(refrigerant="Propane;Pentane", dT_GC_min_tar=20., dT_rhx1_cold=17., tol=1e-6)
    hp, dt = analysis.analyse_RHP_case(None, case, record_result=False)
    assert hp.converged and dt == 3. and len(seen) == 2
    assert all(seen[1][k] == v for k, v in case.items())
    assert "P_H_est" not in case


@pytest.mark.parametrize("success,kind", [(True, None), (False, "physical")])
def test_valid_or_physical_results_do_not_trigger_pressure_retries(monkeypatch, success, kind):
    calls = []
    model = failed_model(kind)
    model.feasible = model.converged = success
    def once(*args):
        calls.append(1)
        return model, 1.
    monkeypatch.setattr(analysis, "_analyse_RHP_case_once", once)
    analysis.analyse_RHP_case(None, {"refrigerant": "Propane;Pentane"}, record_result=False)
    assert len(calls) == 1


def test_objective_retains_native_failure_context():
    hp = failed_model()
    hp.output = {"SEC": None}
    hp.property_failure_context = {"fluid": "Propane;Pentane", "input_pair": "PH", "inputs": [4., 350.]}
    hp.solve_error_stage = "discharge_pressure_search"
    hp.constraint_errors = {"native_property_failure": 1.}
    objective = _OuterDesignObjective(max_evals=1, refrigerant="Propane;Pentane|0.3;0.7|mass",
                                     molar_masses=(44., 72.), analyser=lambda *a, **k: (hp, .1))
    record = objective.evaluate([10., 20., .5, .3, .7])
    assert not record.success and record.message == "PH [226]"
    assert record.report["property_failure_context"] == hp.property_failure_context
    assert record.report["solve_failure_kind"] == "native_property"
    assert record.report["solve_error_stage"] == "discharge_pressure_search"
    from main_RHPoptimiser_benchmark import serialise
    saved = serialise(record)
    assert saved["property_failure_context"] == hp.property_failure_context
    assert saved["solve_error_stage"] == hp.solve_error_stage


@pytest.mark.parametrize("pressure,temperature", [(2.7002492339609954, 133.4524442395218),
                                                 (3.0685294327334156, 139.0634020121554)])
def test_binary_ph_refinement_preserves_requested_enthalpy_and_heat_direction(hpl, pressure, temperature):
    pp = hpl.PropertyPackage("Propane;Pentane|0.3;0.7|mass")
    composition = pp.fld_composition
    before = pp.get_std_properties("PT", pressure, temperature, composition)
    target = before["h"] - 2e-6
    native_p, native_h = pp._convert_props_from_SI_with_C("PH", pressure, target)
    raw = pp._flash("PH", native_p, native_h, composition)
    assert raw.ierr == 0 and 0 < raw.q < 1
    assert abs(raw.h-native_h) > 1e-6
    refined = pp._refine_binary_ph(raw, native_p, native_h, composition)
    assert refined is not None and pp._valid_two_phase_ph(refined, native_p, native_h, composition)
    after = pp.get_std_properties("PH", pressure, target, composition)
    assert after["h"] == pytest.approx(target, abs=1e-8, rel=0)
    assert after["T"] < before["T"]
    assert 0 < after["q"] < 1


def test_binary_ph_refinement_rejects_failed_independent_checks(hpl, monkeypatch):
    pp = hpl.PropertyPackage("Propane;Pentane|0.3;0.7|mass")
    p, h = pp._convert_props_from_SI_with_C("PH", 3.0685294327334156, 423.9590880044125)
    raw = pp._flash("PH", p, h, pp.fld_composition)
    monkeypatch.setattr(pp, "_valid_two_phase_ph", lambda *a: False)
    assert pp._refine_binary_ph(raw, p, h, pp.fld_composition) is None
