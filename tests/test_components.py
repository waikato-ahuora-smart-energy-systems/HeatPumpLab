import math
from types import SimpleNamespace

import pytest

from HeatPumpLab.components.common.base_heat_exchanger import (
    BaseHeatExchanger,
    _get_LMTD,
)


def set_state(stream, **state):
    for name, value in state.items():
        setattr(stream, name, value)

    if stream.pressure is not None and stream.temperature is not None:
        stream.props_in_mPT()
    elif stream.pressure is not None and stream.enthalpy is not None:
        stream.props_in_mPH()
    elif stream.pressure is not None and stream.quality is not None:
        stream.props_in_mPQ()
    elif stream.temperature is not None and stream.quality is not None:
        stream.props_in_mTQ()
    return stream


def test_compressor_preserves_flow_and_raises_refrigerant_pressure(hpl):
    compressor = hpl.Compressor(fluid="Butane", model_index=1, motor_eff=1.0)
    compressor.flow_in = hpl.Stream(
        fluid="Butane", temperature=10.0, quality=1.0, mass_flow=2.0
    )
    compressor.flow_out = hpl.Stream(
        fluid="Butane",
        pressure=compressor.flow_in.pressure * 3.0,
        mass_flow=2.0,
    )

    compressor.solve()

    assert compressor.converged is True
    assert compressor.flow_out.mass_flow == pytest.approx(compressor.flow_in.mass_flow)
    assert compressor.pressure_ratio == pytest.approx(3.0)
    assert 0.4 <= compressor.isentropic_eff <= 1.0
    assert compressor.flow_out.pressure == pytest.approx(
        compressor.flow_in.pressure * 3.0
    )
    assert compressor.flow_out.temperature > compressor.flow_in.temperature
    assert compressor.flow_out.enthalpy > compressor.flow_in.enthalpy
    assert compressor.power > 0.0
    assert compressor.delta_entropy > 0.0


@pytest.mark.parametrize("model_index", [0, 1, 2, 3])
def test_compressor_efficiency_models_remain_physical(hpl, model_index):
    compressor = hpl.Compressor(fluid="Butane", model_index=model_index, motor_eff=1.0)
    compressor.flow_in = hpl.Stream(
        fluid="Butane", temperature=10.0, quality=1.0, mass_flow=1.0
    )
    compressor.flow_out = hpl.Stream(
        fluid="Butane", pressure=compressor.flow_in.pressure * 2.5
    )

    compressor.solve()

    assert compressor.converged is True
    assert 0.4 <= compressor.isentropic_eff <= 1.0
    assert 0.4 <= compressor.overall_eff <= 1.0
    assert math.isfinite(compressor.power)
    assert compressor.power > 0.0


def test_brendel_compressor_efficiency_handles_no_pressure_lift(hpl):
    compressor = hpl.Compressor(fluid="Butane", model_index=3, motor_eff=1.0)
    compressor.flow_in = hpl.Stream(
        fluid="Butane", temperature=10.0, quality=1.0, mass_flow=1.0
    )
    compressor.flow_out = hpl.Stream(
        fluid="Butane",
        pressure=compressor.flow_in.pressure,
        mass_flow=compressor.flow_in.mass_flow,
    )

    assert compressor.get_isentropic_efficiency() == pytest.approx(0.4)


def test_compressor_graph_path_spans_inlet_to_outlet_state(hpl):
    compressor = hpl.Compressor(fluid="Butane", model_index=1)
    compressor.flow_in = hpl.Stream(
        fluid="Butane", temperature=10.0, quality=1.0, mass_flow=1.0
    )
    compressor.flow_out = hpl.Stream(
        fluid="Butane", pressure=compressor.flow_in.pressure * 3.0
    )
    compressor.solve()

    path = compressor.get_graph_data(divisions=5)["siso"]

    assert len(path) == 6
    assert path[0].pressure == pytest.approx(compressor.flow_in.pressure)
    assert path[-1].pressure == pytest.approx(compressor.flow_out.pressure)
    assert path[0].enthalpy == pytest.approx(compressor.flow_in.enthalpy)
    assert path[-1].enthalpy == pytest.approx(compressor.flow_out.enthalpy)
    assert [state.pressure for state in path] == sorted(state.pressure for state in path)


def test_compressor_motor_efficiency_only_changes_electrical_power(hpl):
    def solve_compressor(motor_eff):
        compressor = hpl.Compressor(
            fluid="Butane",
            model_index=1,
            motor_eff=motor_eff,
        )
        compressor.flow_in = hpl.Stream(
            fluid="Butane", temperature=10.0, quality=1.0, mass_flow=1.0
        )
        compressor.flow_out = hpl.Stream(
            fluid="Butane", pressure=compressor.flow_in.pressure * 3.0
        )
        compressor.solve()
        return compressor

    reference = solve_compressor(motor_eff=1.0)
    lossy_motor = solve_compressor(motor_eff=0.8)

    assert lossy_motor.isentropic_eff == pytest.approx(reference.isentropic_eff)
    assert lossy_motor.flow_out.enthalpy == pytest.approx(
        reference.flow_out.enthalpy
    )
    assert lossy_motor.shaft_work == pytest.approx(reference.shaft_work)
    assert lossy_motor.power == pytest.approx(reference.power / 0.8)
    assert lossy_motor.overall_eff == pytest.approx(
        lossy_motor.isentropic_eff * lossy_motor.motor_eff
    )


def test_expansion_valve_isenthalpic_drop_preserves_mass_flow(hpl):
    valve = hpl.ExpansionValve(fluid="Butane", specific_heat_loss=2.0)
    valve.flow_in = hpl.Stream(
        fluid="Butane", temperature=60.0, quality=0.0, mass_flow=1.5
    )
    valve.flow_out = hpl.Stream(
        fluid="Butane", pressure=valve.flow_in.pressure / 3.0, mass_flow=1.5
    )

    valve.solve()

    assert valve.converged is True
    assert valve.flow_out.mass_flow == pytest.approx(valve.flow_in.mass_flow)
    assert valve.flow_out.pressure == pytest.approx(valve.flow_in.pressure / 3.0)
    assert valve.flow_out.enthalpy == pytest.approx(valve.flow_in.enthalpy - 2.0)
    assert valve.power == pytest.approx(0.0)


def test_expansion_valve_graph_path_tracks_pressure_drop(hpl):
    valve = hpl.ExpansionValve(fluid="Butane")
    valve.flow_in = hpl.Stream(
        fluid="Butane", temperature=60.0, quality=0.0, mass_flow=1.0
    )
    valve.flow_out = hpl.Stream(
        fluid="Butane", pressure=valve.flow_in.pressure / 2.0, mass_flow=1.0
    )
    valve.solve()

    path = valve.get_graph_data(divisions=4)["siso"]

    assert len(path) == 5
    assert path[0].pressure == pytest.approx(valve.flow_in.pressure)
    assert path[-1].pressure == pytest.approx(valve.flow_out.pressure)
    assert [state.pressure for state in path] == sorted(
        (state.pressure for state in path), reverse=True
    )


def test_mixer_balances_mass_pressure_and_weighted_enthalpy(hpl):
    mixer = hpl.Mixer(fluid="Butane", specific_heat_loss=1.0)
    mixer.flow_in_1 = hpl.Stream(
        fluid="Butane", temperature=20.0, quality=1.0, mass_flow=1.0
    )
    mixer.flow_in_2 = hpl.Stream(
        fluid="Butane", temperature=40.0, quality=1.0, mass_flow=3.0
    )
    mixer.flow_out = hpl.Stream(fluid="Butane")

    mixer.solve()

    expected_enthalpy = (
        mixer.flow_in_1.enthalpy * 1.0 + mixer.flow_in_2.enthalpy * 3.0
    ) / 4.0 - 1.0
    assert mixer.converged is True
    assert mixer.flow_out.mass_flow == pytest.approx(4.0)
    assert mixer.flow_out.pressure == pytest.approx(
        min(mixer.flow_in_1.pressure, mixer.flow_in_2.pressure)
    )
    assert mixer.flow_out.enthalpy == pytest.approx(expected_enthalpy)
    assert mixer.flow_in_1.temperature < mixer.flow_out.temperature
    assert mixer.flow_out.temperature < mixer.flow_in_2.temperature


def test_mixer_zero_flow_returns_first_inlet_state_and_reports_it(hpl, capsys):
    mixer = hpl.Mixer(fluid="Butane")
    mixer.flow_in_1 = hpl.Stream(
        fluid="Butane", temperature=20.0, quality=1.0, mass_flow=0.0
    )
    mixer.flow_in_2 = hpl.Stream(
        fluid="Butane", temperature=40.0, quality=1.0, mass_flow=0.0
    )
    mixer.flow_out = hpl.Stream(fluid="Butane")

    mixer.solve()

    assert mixer.flow_out.mass_flow == pytest.approx(0.0)
    assert mixer.flow_out.enthalpy == pytest.approx(mixer.flow_in_1.enthalpy)
    assert "No mass flow through the mixer" in capsys.readouterr().out


def test_splitter_conserves_mass_and_copies_state_to_each_branch(hpl):
    splitter = hpl.Splitter(fluid="Butane", r_ratio=3.0)
    splitter.flow_in = hpl.Stream(
        fluid="Butane", temperature=30.0, quality=1.0, mass_flow=8.0
    )

    splitter.solve()

    assert splitter.converged is True
    assert splitter.flow_out_1.mass_flow == pytest.approx(2.0)
    assert splitter.flow_out_2.mass_flow == pytest.approx(6.0)
    assert splitter.flow_out_1.mass_flow + splitter.flow_out_2.mass_flow == (
        pytest.approx(splitter.flow_in.mass_flow)
    )
    for outlet in (splitter.flow_out_1, splitter.flow_out_2):
        assert outlet.pressure == pytest.approx(splitter.flow_in.pressure)
        assert outlet.temperature == pytest.approx(splitter.flow_in.temperature)
        assert outlet.enthalpy == pytest.approx(splitter.flow_in.enthalpy)


def test_splitter_zero_ratio_sends_all_flow_to_primary_branch(hpl):
    splitter = hpl.Splitter(fluid="Butane", r_ratio=0.0)
    splitter.flow_in = hpl.Stream(
        fluid="Butane", temperature=30.0, quality=1.0, mass_flow=5.0
    )

    splitter.solve()

    assert splitter.flow_out_1.mass_flow == pytest.approx(5.0)
    assert splitter.flow_out_2.mass_flow == pytest.approx(0.0)
    assert splitter.flow_out_1.enthalpy == pytest.approx(splitter.flow_in.enthalpy)
    assert splitter.flow_out_2.enthalpy == pytest.approx(splitter.flow_in.enthalpy)


def test_turbine_generates_power_from_pressure_drop(hpl):
    turbine = hpl.Turbine(fluid="Butane", eff_isen=0.85, eff_gen=0.9)
    turbine.flow_in = hpl.Stream(
        fluid="Butane", pressure=1.0, temperature=80.0, mass_flow=2.0
    )
    turbine.flow_out = hpl.Stream(fluid="Butane", pressure=0.25, mass_flow=2.0)

    turbine.solve()

    assert turbine.converged is True
    assert turbine.flow_out.enthalpy < turbine.flow_in.enthalpy
    assert turbine.flow_out.mass_flow == pytest.approx(turbine.flow_in.mass_flow)
    assert turbine.power > 0.0
    assert turbine.power == pytest.approx(
        turbine.dh * turbine.flow_out.mass_flow * turbine.eff_gen
    )


def test_turbine_conserves_mass_and_lower_efficiency_extracts_less_power(hpl):
    def solve_turbine(eff_isen):
        turbine = hpl.Turbine(fluid="Butane", eff_isen=eff_isen, eff_gen=1.0)
        turbine.flow_in = hpl.Stream(
            fluid="Butane", pressure=1.0, temperature=80.0, mass_flow=2.0
        )
        turbine.flow_out = hpl.Stream(fluid="Butane", pressure=0.25)
        turbine.solve()
        return turbine

    reversible = solve_turbine(eff_isen=1.0)
    lossy = solve_turbine(eff_isen=0.7)

    assert lossy.flow_out.mass_flow == pytest.approx(lossy.flow_in.mass_flow)
    assert reversible.flow_out.mass_flow == pytest.approx(
        reversible.flow_in.mass_flow
    )
    assert lossy.power == pytest.approx(reversible.power * 0.7)
    assert lossy.flow_out.enthalpy > reversible.flow_out.enthalpy


def test_stream_state_pairs_recreate_the_same_thermodynamic_state(hpl):
    reference = hpl.Stream(
        fluid="Butane", pressure=0.6, temperature=40.0, mass_flow=2.5
    )

    recreated_from_enthalpy = hpl.Stream(
        fluid="Butane",
        pressure=reference.pressure,
        enthalpy=reference.enthalpy,
        mass_flow=reference.mass_flow,
    )
    recreated_from_entropy = hpl.Stream(
        fluid="Butane",
        pressure=reference.pressure,
        entropy=reference.entropy,
        mass_flow=reference.mass_flow,
    )

    for recreated in (recreated_from_enthalpy, recreated_from_entropy):
        assert recreated.pressure == pytest.approx(reference.pressure)
        assert recreated.temperature == pytest.approx(reference.temperature, abs=1e-8)
        assert recreated.enthalpy == pytest.approx(reference.enthalpy, abs=1e-8)
        assert recreated.entropy == pytest.approx(reference.entropy, abs=1e-8)
        assert recreated.heat_flow == pytest.approx(
            reference.enthalpy * reference.mass_flow
        )


def test_stream_density_pair_recreates_state_and_transport_properties(hpl):
    reference = hpl.Stream(fluid="Butane", pressure=0.6, temperature=40.0)
    recreated = hpl.Stream(fluid="Butane", pressure=reference.pressure)
    recreated.density = reference.density

    recreated.props_in_mPD()
    recreated.props_extend()

    assert recreated.temperature == pytest.approx(reference.temperature, abs=1e-8)
    assert recreated.enthalpy == pytest.approx(reference.enthalpy, abs=1e-8)
    assert recreated.prandtl > 0.0
    assert recreated.viscosity > 0.0
    assert recreated.thermal_conductivity > 0.0


def test_stream_generic_property_update_can_rename_and_reflow_stream(hpl):
    stream = hpl.Stream(fluid="Butane")

    stream.props_in_ab("PT", 0.7, 35.0, mass_flow=2.4, name="updated stream")

    assert stream.name == "updated stream"
    assert stream.pressure == pytest.approx(0.7)
    assert stream.temperature == pytest.approx(35.0)
    assert stream.mass_flow == pytest.approx(2.4)
    assert stream.heat_flow == pytest.approx(stream.enthalpy * stream.mass_flow)
    assert stream.vol_flow == pytest.approx(stream.specific_volume * stream.mass_flow)


def test_stream_temperature_entropy_pair_recreates_state(hpl):
    reference = hpl.Stream(fluid="Butane", pressure=0.6, temperature=35.0)
    recreated = hpl.Stream(
        fluid="Butane",
        temperature=reference.temperature,
        entropy=reference.entropy,
        mass_flow=1.2,
    )

    assert recreated.pressure == pytest.approx(reference.pressure, abs=1e-6)
    assert recreated.enthalpy == pytest.approx(reference.enthalpy, abs=1e-5)
    assert recreated.heat_flow == pytest.approx(recreated.enthalpy * 1.2)


def test_stream_saturation_helpers_return_critical_state_above_critical_pressure(hpl):
    stream = hpl.Stream(fluid="Butane", pressure=0.6, temperature=40.0)
    stream.pressure = stream.P_crit * 1.05

    assert math.isfinite(stream.get_Tsat(1.0))
    assert math.isfinite(stream.get_hsat(1.0))


def test_stream_copy_and_flow_updates_keep_extensive_properties_consistent(hpl):
    source = hpl.Stream(fluid="water", pressure=1.2, temperature=45.0, mass_flow=3.0)
    target = hpl.Stream(fluid="water")

    target.copy_fld_props(source)
    target.set_flows(1.25)

    assert target.pressure == pytest.approx(source.pressure)
    assert target.temperature == pytest.approx(source.temperature)
    assert target.enthalpy == pytest.approx(source.enthalpy)
    assert target.mass_flow == pytest.approx(1.25)
    assert target.heat_flow == pytest.approx(target.enthalpy * target.mass_flow)
    assert target.vol_flow == pytest.approx(target.specific_volume * target.mass_flow)


def test_binary_stream_saturation_helpers_preserve_temperature_glide(hpl):
    fluid = "Propane;Butane|0.5;0.5"
    dew = hpl.Stream(fluid=fluid, temperature=10.0, quality=1.0)

    bubble_temperature_at_dew_pressure = dew.get_Tsat(0.0)

    assert dew.get_Tsat(1.0) == pytest.approx(dew.temperature, abs=1e-8)
    assert dew.get_Psat(1.0) == pytest.approx(dew.pressure, abs=1e-8)
    assert dew.get_hsat(1.0) == pytest.approx(dew.enthalpy, abs=1e-8)
    assert dew.temperature - bubble_temperature_at_dew_pressure > 8.0


def test_binary_compressor_preserves_mixture_identity_while_raising_pressure(hpl):
    fluid = "Propane;Butane|0.5;0.5"
    compressor = hpl.Compressor(fluid=fluid, model_index=1, motor_eff=1.0)
    compressor.flow_in = hpl.Stream(
        fluid=fluid, temperature=10.0, quality=1.0, mass_flow=1.4
    )
    compressor.flow_out = hpl.Stream(
        fluid=fluid,
        pressure=compressor.flow_in.pressure * 2.5,
        mass_flow=compressor.flow_in.mass_flow,
    )

    compressor.solve()

    assert compressor.converged is True
    assert compressor.flow_out.fld_composition == pytest.approx(
        compressor.flow_in.fld_composition
    )
    assert compressor.flow_out.mass_flow == pytest.approx(
        compressor.flow_in.mass_flow
    )
    assert compressor.flow_out.pressure == pytest.approx(
        compressor.flow_in.pressure * 2.5
    )
    assert compressor.flow_out.temperature > compressor.flow_in.temperature
    assert compressor.power > 0.0


def test_binary_expansion_valve_preserves_enthalpy_and_composition(hpl):
    fluid = "Propane;Butane|0.5;0.5"
    valve = hpl.ExpansionValve(fluid=fluid)
    valve.flow_in = hpl.Stream(fluid=fluid, temperature=30.0, quality=0.0, mass_flow=2.0)
    valve.flow_out = hpl.Stream(
        fluid=fluid,
        pressure=valve.flow_in.pressure / 2.0,
        mass_flow=valve.flow_in.mass_flow,
    )

    valve.solve()

    assert valve.converged is True
    assert valve.flow_out.fld_composition == pytest.approx(valve.flow_in.fld_composition)
    assert valve.flow_out.mass_flow == pytest.approx(valve.flow_in.mass_flow)
    assert valve.flow_out.enthalpy == pytest.approx(valve.flow_in.enthalpy)
    assert valve.flow_out.temperature < valve.flow_in.temperature


def test_heat_exchanger_fixed_cold_outlet_delivers_requested_temperature(hpl):
    exchanger = hpl.HeatExchanger(
        "water",
        "water",
        model="fixed cold temperature outlet",
        dP_h=0.01,
        dP_c=0.02,
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)
    exchanger.source_out.mass_flow = exchanger.source_in.mass_flow
    exchanger.sink_out.mass_flow = exchanger.sink_in.mass_flow
    exchanger.sink_out.temperature = 50.0

    exchanger.solve()

    assert exchanger.converged is True
    assert exchanger.duty > 0.0
    assert exchanger.sink_out.temperature == pytest.approx(50.0)
    assert exchanger.source_out.temperature < exchanger.source_in.temperature
    assert exchanger.source_out.temperature > exchanger.sink_in.temperature
    assert exchanger.source_out.pressure == pytest.approx(1.99)
    assert exchanger.sink_out.pressure == pytest.approx(1.98)


def test_heat_exchanger_effectiveness_model_moves_both_streams_toward_each_other(hpl):
    exchanger = hpl.HeatExchanger(
        "water",
        "water",
        model="effectiveness",
        HX_Eff=0.5,
        dP_h=0.01,
        dP_c=0.02,
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)
    exchanger.source_out.mass_flow = 2.0
    exchanger.sink_out.mass_flow = 3.0

    exchanger.solve()

    assert exchanger.converged is True
    assert exchanger.duty > 0.0
    assert exchanger.source_out.temperature < exchanger.source_in.temperature
    assert exchanger.sink_out.temperature > exchanger.sink_in.temperature
    assert exchanger.source_out.pressure == pytest.approx(1.99)
    assert exchanger.sink_out.pressure == pytest.approx(1.98)


def test_heat_exchanger_dtmin_terminal_limits_hot_outlet_approach(hpl):
    exchanger = hpl.HeatExchanger(
        "Butane",
        "Air",
        model="dTmin_terminal",
        dTmin=10.0,
        dP_h=0.01,
        dP_c=0.02,
    )
    set_state(exchanger.source_in, pressure=1.0, temperature=100.0, mass_flow=1.0)
    set_state(exchanger.sink_in, pressure=1.0, temperature=20.0, mass_flow=2.0)
    exchanger.source_out.mass_flow = 1.0
    exchanger.sink_out.mass_flow = 2.0

    exchanger.solve()

    assert exchanger.converged is True
    assert exchanger.duty > 0.0
    assert exchanger.source_out.temperature == pytest.approx(
        exchanger.sink_in.temperature + exchanger.dTmin, abs=1e-5
    )
    assert exchanger.sink_out.temperature > exchanger.sink_in.temperature
    assert exchanger.source_out.pressure == pytest.approx(0.99)
    assert exchanger.sink_out.pressure == pytest.approx(0.98)


def test_heat_exchanger_fixed_duty_models_balance_source_and_sink_energy(hpl):
    cold_duty = hpl.HeatExchanger(
        "Butane",
        "Air",
        model="fixed cold duty",
        dP_h=0.01,
        dP_c=0.02,
    )
    set_state(cold_duty.source_in, pressure=1.0, temperature=100.0, mass_flow=1.0)
    set_state(cold_duty.sink_in, pressure=1.0, temperature=20.0, mass_flow=2.0)
    cold_duty.source_out.mass_flow = 1.0
    cold_duty.sink_out.pressure = cold_duty.sink_in.pressure
    cold_duty.sink_out.enthalpy = cold_duty.sink_in.enthalpy + 10.0
    cold_duty.sink_out.props_in_mPH()
    cold_duty.sink_out.set_flows(2.0)

    cold_duty.solve()

    hot_duty = hpl.HeatExchanger(
        "Butane",
        "Air",
        model="fixed hot duty",
        dP_h=0.01,
        dP_c=0.02,
    )
    set_state(hot_duty.source_in, pressure=1.0, temperature=100.0, mass_flow=1.0)
    set_state(hot_duty.sink_in, pressure=1.0, temperature=20.0, mass_flow=2.0)
    hot_duty.source_out.pressure = hot_duty.source_in.pressure
    hot_duty.source_out.enthalpy = hot_duty.source_in.enthalpy - 20.0
    hot_duty.source_out.props_in_mPH()
    hot_duty.source_out.set_flows(1.0)
    hot_duty.sink_out.mass_flow = 2.0

    hot_duty.solve()

    for exchanger in (cold_duty, hot_duty):
        assert exchanger.converged is True
        assert exchanger.duty == pytest.approx(20.0)
        assert exchanger.source_out.temperature < exchanger.source_in.temperature
        assert exchanger.sink_out.temperature > exchanger.sink_in.temperature


def test_heat_exchanger_hot_outlet_constraints_set_temperature_and_quality(hpl):
    delta_t = hpl.HeatExchanger(
        "water",
        "water",
        model="fixed hot delta temperature",
        dP_h=0.01,
        dP_c=0.02,
    )
    set_state(delta_t.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(delta_t.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)
    delta_t.source_out.mass_flow = 2.0
    delta_t.sink_out.mass_flow = 3.0
    delta_t.dT_h = 25.0

    delta_t.solve()

    quality = hpl.HeatExchanger("Butane", "water", model="fixed hot quality outlet")
    set_state(quality.source_in, pressure=0.7, quality=1.0, mass_flow=1.0)
    set_state(quality.sink_in, pressure=2.0, temperature=20.0, mass_flow=4.0)
    quality.source_out.mass_flow = 1.0
    quality.sink_out.mass_flow = 4.0
    quality.source_out.quality = 0.0

    quality.solve()

    assert delta_t.source_out.temperature == pytest.approx(65.0)
    assert delta_t.sink_out.temperature > delta_t.sink_in.temperature
    assert quality.source_out.quality == pytest.approx(0.0)
    assert quality.sink_out.temperature > quality.sink_in.temperature


def test_heat_exchanger_bypass_and_zero_flow_do_not_exchange_heat(hpl):
    bypassed = hpl.HeatExchanger(
        "water",
        "water",
        model="fixed cold delta temperature",
        bypassed=True,
    )
    set_state(bypassed.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(bypassed.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)

    bypassed.solve()

    zero_flow = hpl.HeatExchanger(
        "water",
        "water",
        model="fixed cold delta temperature",
    )
    set_state(zero_flow.source_in, pressure=2.0, temperature=90.0, mass_flow=0.0)
    set_state(zero_flow.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)

    zero_flow.solve()

    assert bypassed.source_out is bypassed.source_in
    assert bypassed.sink_out is bypassed.sink_in
    assert bypassed.duty == pytest.approx(0.0)
    assert zero_flow.source_out.mass_flow == pytest.approx(0.0)
    assert zero_flow.source_out.temperature == pytest.approx(
        zero_flow.source_in.temperature
    )
    assert zero_flow.sink_out.temperature == pytest.approx(zero_flow.sink_in.temperature)
    assert zero_flow.duty == pytest.approx(0.0)


def test_heat_exchanger_temperature_cross_is_reported_as_not_converged(hpl):
    exchanger = hpl.HeatExchanger(
        "water",
        "water",
        model="fixed cold temperature outlet",
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=70.0, mass_flow=1.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=60.0, mass_flow=1.0)
    exchanger.source_out.mass_flow = 1.0
    exchanger.sink_out.mass_flow = 1.0
    exchanger.sink_out.temperature = 80.0

    exchanger.solve()

    assert exchanger.duty > 0.0
    assert exchanger.source_in.temperature < exchanger.sink_out.temperature
    assert exchanger.converged is False


def test_heat_exchanger_rejects_unknown_model(hpl):
    exchanger = hpl.HeatExchanger("water", "water", model="unknown")
    set_state(exchanger.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)

    with pytest.raises(ValueError):
        exchanger.solve()


def test_base_heat_exchanger_specialised_model_dispatches_to_subclass():
    class SpecialisedHeatExchanger(BaseHeatExchanger):
        model = "specialised"

        def __init__(self):
            self.called = False

        def get_detailed_performance(self):
            self.called = True

    exchanger = SpecialisedHeatExchanger()

    exchanger.get_performance()

    assert exchanger.called is True


def test_base_heat_exchanger_base_constructor_is_available():
    assert BaseHeatExchanger() is not None


def test_heat_exchanger_wraps_terminal_dt_property_failures(hpl, monkeypatch):
    exchanger = hpl.HeatExchanger(
        "water",
        "water",
        name="terminal",
        model="dTmin_terminal",
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)

    def fail_props(*args, **kwargs):
        raise RuntimeError("property failure")

    monkeypatch.setattr(exchanger.source_out_max, "props_in_ab", fail_props)

    with pytest.raises(ValueError, match="terminal") as exc_info:
        exchanger.solve()

    assert isinstance(exc_info.value.__cause__, RuntimeError)


def test_heat_exchanger_wraps_effectiveness_property_failures(hpl, monkeypatch):
    exchanger = hpl.HeatExchanger(
        "water",
        "water",
        name="effectiveness",
        model="effectiveness",
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)

    def fail_props():
        raise RuntimeError("property failure")

    monkeypatch.setattr(exchanger.source_out_max, "props_in_mPT", fail_props)

    with pytest.raises(ValueError, match="effectiveness") as exc_info:
        exchanger.solve()

    assert isinstance(exc_info.value.__cause__, RuntimeError)


def test_heat_exchanger_lmtd_average_is_positive_for_counterflow_profile(hpl):
    exchanger = hpl.HeatExchanger(
        "water",
        "water",
        model="fixed cold temperature outlet",
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=20.0, mass_flow=3.0)
    exchanger.source_out.mass_flow = 2.0
    exchanger.sink_out.mass_flow = 3.0
    exchanger.sink_out.temperature = 50.0
    exchanger.solve()

    assert exchanger.get_LMTD_ave() > 0.0


def test_lmtd_returns_common_temperature_drop_when_terminal_drops_match():
    assert _get_LMTD(Thin=100.0, Thout=80.0, Tcin=30.0, Tcout=50.0) == pytest.approx(
        50.0
    )


def test_lmtd_reports_temperature_cross_as_infeasible():
    with pytest.raises(ValueError, match="temperature cross"):
        _get_LMTD(Thin=100.0, Thout=40.0, Tcin=50.0, Tcout=80.0)


def test_plate_heat_exchanger_fixed_cold_delta_temperature_heats_sink(hpl):
    exchanger = hpl.PlateHeatExchanger(
        "water",
        "water",
        model="fixed cold delta temperature",
        dT_c=25.0,
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=95.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=25.0, mass_flow=2.0)
    exchanger.source_out.mass_flow = 2.0
    exchanger.sink_out.mass_flow = 2.0

    exchanger.solve()

    assert exchanger.converged is True
    assert exchanger.duty > 0.0
    assert exchanger.sink_out.temperature == pytest.approx(50.0)
    assert exchanger.source_out.temperature < exchanger.source_in.temperature
    assert exchanger.source_out.temperature > exchanger.sink_in.temperature


def test_heat_exchanger_graph_data_spans_solved_endpoints(hpl):
    exchanger = hpl.PlateHeatExchanger(
        "water",
        "water",
        model="fixed cold delta temperature",
        dT_c=20.0,
    )
    set_state(exchanger.source_in, pressure=2.0, temperature=90.0, mass_flow=2.0)
    set_state(exchanger.sink_in, pressure=2.0, temperature=25.0, mass_flow=2.0)
    exchanger.source_out.mass_flow = 2.0
    exchanger.sink_out.mass_flow = 2.0
    exchanger.solve()

    data = exchanger.get_graph_data(divisions=4)

    assert len(data["source"]) == 5
    assert len(data["sink"]) == 5
    assert data["source"][0].temperature == pytest.approx(
        exchanger.source_in.temperature
    )
    assert data["source"][-1].temperature == pytest.approx(
        exchanger.source_out.temperature
    )
    assert data["sink"][0].temperature == pytest.approx(exchanger.sink_in.temperature)
    assert data["sink"][-1].temperature == pytest.approx(exchanger.sink_out.temperature)


def test_gas_cooler_fixed_hot_outlet_heats_air_sink(hpl):
    gas_cooler = hpl.GasCooler(
        "Butane",
        "Air",
        model="fixed hot temperature outlet",
        dTmin=10.0,
    )
    set_state(gas_cooler.source_in, pressure=2.0, temperature=120.0, mass_flow=1.0)
    set_state(gas_cooler.sink_in, pressure=1.0, temperature=25.0, mass_flow=3.0)
    gas_cooler.source_out.mass_flow = 1.0
    gas_cooler.sink_out.mass_flow = 3.0
    gas_cooler.source_out.temperature = 55.0

    gas_cooler.solve()

    assert gas_cooler.converged is True
    assert gas_cooler.duty > 0.0
    assert gas_cooler.source_out.temperature == pytest.approx(55.0)
    assert gas_cooler.sink_out.temperature > gas_cooler.sink_in.temperature
    assert gas_cooler.source_in.temperature > gas_cooler.sink_out.temperature


def test_evaporator_saturated_gas_model_reaches_vapour_outlet(hpl):
    evaporator = hpl.Evaporator("water", "Butane")
    set_state(evaporator.source_in, pressure=2.0, temperature=40.0, mass_flow=10.0)
    evaporator.sink_out.temperature = 5.0
    evaporator.sink_out.quality = 1.0
    evaporator.sink_out.mass_flow = 1.0
    evaporator.sink_out.props_in_mTQ()
    evaporator.sink_in.mass_flow = 1.0
    evaporator.sink_in.pressure = evaporator.sink_out.pressure
    evaporator.sink_in.enthalpy = evaporator.sink_out.get_hsat(0.0)
    evaporator.sink_in.props_in_mPH()
    evaporator.source_out.mass_flow = evaporator.source_in.mass_flow

    evaporator.solve()

    assert evaporator.converged is True
    assert evaporator.duty > 0.0
    assert evaporator.sink_out.temperature == pytest.approx(5.0)
    assert evaporator.sink_out.quality == pytest.approx(1.0)
    assert evaporator.source_out.temperature < evaporator.source_in.temperature
    assert evaporator.source_out.temperature > evaporator.sink_out.temperature


def test_binary_evaporator_heats_gliding_mixture_to_dew_point(hpl):
    fluid = "Propane;Butane|0.5;0.5"
    evaporator = hpl.Evaporator("water", fluid)
    set_state(evaporator.source_in, pressure=2.0, temperature=40.0, mass_flow=20.0)
    evaporator.sink_out.temperature = 10.0
    evaporator.sink_out.quality = 1.0
    evaporator.sink_out.mass_flow = 0.4
    evaporator.sink_out.props_in_mTQ()
    evaporator.sink_in.mass_flow = evaporator.sink_out.mass_flow
    evaporator.sink_in.pressure = evaporator.sink_out.pressure
    evaporator.sink_in.enthalpy = evaporator.sink_out.get_hsat(0.0)
    evaporator.sink_in.props_in_mPH()
    evaporator.source_out.mass_flow = evaporator.source_in.mass_flow

    evaporator.solve()

    assert evaporator.converged is True
    assert evaporator.duty > 0.0
    assert evaporator.sink_out.quality == pytest.approx(1.0)
    assert evaporator.sink_out.temperature > evaporator.sink_in.temperature + 8.0
    assert evaporator.source_out.temperature < evaporator.source_in.temperature
    assert evaporator.source_out.temperature > evaporator.sink_out.temperature


def test_economiser_adds_heat_until_target_approach_is_met(hpl):
    economiser = hpl.Economiser("Butane", "Butane", dTmin=5.0)
    set_state(economiser.source_out, pressure=0.8, temperature=40.0, mass_flow=1.2)
    set_state(economiser.sink_in, pressure=0.3, quality=0.05, mass_flow=0.3)
    comp1_out = hpl.Stream(fluid="Butane", pressure=0.3, temperature=20.0)
    comp2_out = hpl.Stream(fluid="Butane", pressure=0.8, temperature=80.0)
    economiser.set_add_obj(comp1_out, comp2_out)

    economiser.solve()

    assert economiser.converged is True
    assert economiser.duty > 0.0
    assert economiser.source_in.temperature - economiser.sink_out.temperature == (
        pytest.approx(5.0, abs=1e-5)
    )
    assert economiser.sink_out.enthalpy > economiser.sink_in.enthalpy
    assert economiser.source_in.enthalpy > economiser.source_out.enthalpy


def test_economiser_graph_data_spans_solved_hot_and_cold_profiles(hpl):
    economiser = hpl.Economiser("Butane", "Butane", dTmin=5.0)
    set_state(economiser.source_out, pressure=0.8, temperature=40.0, mass_flow=1.2)
    set_state(economiser.sink_in, pressure=0.3, quality=0.05, mass_flow=0.3)
    economiser.set_add_obj(
        hpl.Stream(fluid="Butane", pressure=0.3, temperature=20.0),
        hpl.Stream(fluid="Butane", pressure=0.8, temperature=80.0),
    )
    economiser.solve()

    data = economiser.get_graph_data(divisions=4)

    assert len(data["source"]) == 5
    assert len(data["sink"]) == 5
    assert data["source"][0].temperature == pytest.approx(
        economiser.source_in.temperature
    )
    assert data["source"][-1].temperature == pytest.approx(
        economiser.source_out.temperature
    )
    assert data["sink"][0].temperature == pytest.approx(economiser.sink_in.temperature)
    assert data["sink"][-1].temperature == pytest.approx(economiser.sink_out.temperature)
    assert data["source"][0].enthalpy > data["source"][-1].enthalpy
    assert data["sink"][0].enthalpy < data["sink"][-1].enthalpy


def test_economiser_reports_zero_duty_when_sink_flow_is_zero(hpl):
    economiser = hpl.Economiser("Butane", "Butane", dTmin=5.0)
    set_state(economiser.source_out, pressure=0.8, temperature=40.0, mass_flow=1.2)
    set_state(economiser.sink_in, pressure=0.3, quality=0.05, mass_flow=0.0)
    economiser.set_add_obj(
        hpl.Stream(fluid="Butane", pressure=0.3, temperature=20.0),
        hpl.Stream(fluid="Butane", pressure=0.8, temperature=80.0),
    )

    economiser.solve()

    assert economiser.converged is True
    assert economiser.duty == pytest.approx(0.0)
    assert economiser.source_in.enthalpy == pytest.approx(
        economiser.source_out.enthalpy
    )
    assert economiser.sink_out.enthalpy == pytest.approx(economiser.sink_in.enthalpy)


def test_economiser_reports_zero_duty_for_supercritical_sink_inlet(hpl):
    economiser = hpl.Economiser("Butane", "Butane", dTmin=5.0)
    set_state(economiser.source_out, pressure=0.8, temperature=40.0, mass_flow=1.2)
    economiser.sink_in.pressure = economiser.sink_in.P_crit * 1.05
    economiser.sink_in.temperature = economiser.sink_in.T_crit + 10.0
    economiser.sink_in.mass_flow = 0.3
    economiser.sink_in.props_in_mPT()
    economiser.set_add_obj(
        hpl.Stream(fluid="Butane", pressure=0.3, temperature=20.0),
        hpl.Stream(fluid="Butane", pressure=0.8, temperature=80.0),
    )

    economiser.solve()

    assert economiser.converged is True
    assert economiser.duty == pytest.approx(0.0)
    assert economiser.source_in.enthalpy == pytest.approx(
        economiser.source_out.enthalpy
    )
    assert economiser.sink_out.enthalpy == pytest.approx(economiser.sink_in.enthalpy)


def test_set_relationship_updates_temperature_from_source_stream(hpl):
    source = hpl.Stream(fluid="water", pressure=1.5, temperature=35.0)
    target = hpl.Stream(fluid="water", pressure=1.5)

    relationship = hpl.SetRelationship(
        source_obj=source,
        target_obj=[target],
        m=1.0,
        c=10.0,
        prop="temperature",
        prop2="pressure",
    )

    assert relationship.converged is True
    assert target.temperature == pytest.approx(45.0)
    assert target.pressure == pytest.approx(source.pressure)
    assert math.isfinite(target.enthalpy)


def test_set_relationship_updates_quality_at_fixed_pressure(hpl):
    target = hpl.Stream(fluid="Butane", pressure=0.3)

    hpl.SetRelationship(
        source_obj=[1.0],
        target_obj=[target],
        prop="quality",
        prop2="pressure",
    )

    assert target.quality == pytest.approx(1.0)
    assert target.pressure == pytest.approx(0.3)
    assert math.isfinite(target.temperature)


def test_set_relationship_updates_temperature_at_fixed_quality(hpl):
    source = hpl.Stream(fluid="Butane", temperature=12.0, quality=1.0)
    target = hpl.Stream(fluid="Butane", quality=1.0)

    hpl.SetRelationship(
        source_obj=source,
        target_obj=[target],
        m=1.0,
        c=3.0,
        prop="temperature",
        prop2="quality",
    )

    assert target.temperature == pytest.approx(15.0)
    assert target.quality == pytest.approx(1.0)
    assert math.isfinite(target.pressure)


def test_set_relationship_updates_quality_at_fixed_temperature(hpl):
    target = hpl.Stream(fluid="Butane", temperature=10.0)

    hpl.SetRelationship(
        source_obj=[0.0],
        target_obj=[target],
        prop="quality",
        prop2="temperature",
    )

    assert target.temperature == pytest.approx(10.0)
    assert target.quality == pytest.approx(0.0)
    assert math.isfinite(target.pressure)


def test_set_relationship_updates_mass_flow_and_heat_flow(hpl):
    source = hpl.Stream(fluid="water", pressure=1.0, temperature=20.0, mass_flow=2.0)
    target = hpl.Stream(fluid="water", pressure=1.0, temperature=20.0, mass_flow=0.0)

    hpl.SetRelationship(
        source_obj=source,
        target_obj=[target],
        m=3.0,
        c=1.0,
        prop="Mass Flow",
    )

    assert target.mass_flow == pytest.approx(7.0)
    assert target.heat_flow == pytest.approx(target.mass_flow * target.enthalpy)


def test_set_relationship_updates_each_target_stream(hpl):
    source = hpl.Stream(fluid="water", pressure=1.0, temperature=20.0, mass_flow=1.5)
    targets = [
        hpl.Stream(fluid="water", pressure=1.0, temperature=20.0, mass_flow=0.0),
        hpl.Stream(fluid="water", pressure=1.0, temperature=30.0, mass_flow=0.0),
    ]

    relationship = hpl.SetRelationship(
        source_obj=source,
        target_obj=targets,
        m=2.0,
        c=0.5,
        prop="Mass Flow",
    )

    assert relationship.converged is True
    for target in targets:
        assert target.mass_flow == pytest.approx(3.5)
        assert target.heat_flow == pytest.approx(target.mass_flow * target.enthalpy)


def test_set_relationship_rejects_unknown_state_pair(hpl):
    source = hpl.Stream(fluid="water", pressure=1.5, temperature=35.0)
    target = hpl.Stream(fluid="water", pressure=1.5)

    with pytest.raises(ValueError):
        hpl.SetRelationship(
            source_obj=source,
            target_obj=[target],
            prop="temperature",
            prop2="enthalpy",
        )


def test_set_relationship_rejects_unknown_quality_state_pair(hpl):
    target = hpl.Stream(fluid="Butane", pressure=0.3)

    with pytest.raises(ValueError):
        hpl.SetRelationship(
            source_obj=[1.0],
            target_obj=[target],
            prop="quality",
            prop2="enthalpy",
        )


def test_set_relationship_rejects_unknown_relationship_type(hpl):
    source = hpl.Stream(fluid="water", pressure=1.5, temperature=35.0)
    target = hpl.Stream(fluid="water", pressure=1.5)

    with pytest.raises(ValueError):
        hpl.SetRelationship(
            name="unsupported",
            source_obj=source,
            target_obj=[target],
            prop="pressure",
            prop2="temperature",
        )


def test_custom_relationship_distributes_intermediate_compressor_pressures(hpl):
    suction = hpl.Stream(fluid="Butane", pressure=0.2, temperature=10.0)
    discharge = hpl.Stream(fluid="Butane", pressure=1.6, temperature=80.0)
    intermediate_1 = hpl.Stream(fluid="Butane")
    intermediate_2 = hpl.Stream(fluid="Butane")
    relationship = hpl.CustomRelationship(
        source_obj=[suction, discharge, 0.2, 0.3],
        target_obj=[intermediate_1, intermediate_2],
        relation=1,
    )

    assert relationship.solve() is True
    assert relationship.converged is True
    assert suction.pressure < intermediate_1.pressure
    assert intermediate_1.pressure < intermediate_2.pressure
    assert intermediate_2.pressure < discharge.pressure


def test_custom_relationship_temperature_target_uses_configured_offset(hpl):
    target = SimpleNamespace(temperature=None)
    relationship = hpl.CustomRelationship(
        source_obj=[
            SimpleNamespace(temperature=20.0),
            SimpleNamespace(temperature=40.0),
        ],
        target_obj=target,
        relation=0,
        var=35.0,
    )

    assert relationship.solve() is True
    assert relationship.converged is True
    assert target.temperature == pytest.approx(55.0)


def test_custom_relationship_rejects_unknown_relation_id(hpl):
    relationship = hpl.CustomRelationship(relation=99)

    with pytest.raises(ValueError):
        relationship.solve()
