from types import SimpleNamespace

import pytest

from HeatPumpLab.components.custom_relationship import CustomRelationship
from HeatPumpLab.components import custom_relationship as custom_relationship_module


class FakeStream:
    def __init__(
        self,
        fluid=None,
        pp_instance=None,
        pressure=None,
        temperature=None,
        enthalpy=None,
        quality=None,
        mass_flow=1.0,
        P_crit=100.0,
    ):
        self.fluid = fluid
        self.pp_instance = pp_instance
        self.pressure = pressure
        self.temperature = temperature
        self.enthalpy = enthalpy
        self.quality = quality
        self.mass_flow = mass_flow
        self.heat_capacity = 100.0
        self.P_crit = P_crit
        self.heat_flow = (self.enthalpy or 0.0) * self.mass_flow
        if self.enthalpy is None and self.temperature is not None:
            self.enthalpy = self.temperature
            self.heat_flow = self.enthalpy * self.mass_flow
        if self.temperature is None and self.enthalpy is not None:
            self.temperature = self.enthalpy
        if self.enthalpy is None and self.quality is not None:
            self.enthalpy = 80.0
            self.temperature = 80.0
            self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPT(self):
        if self.enthalpy is None:
            self.enthalpy = self.temperature
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPH(self):
        self.temperature = self.enthalpy
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPQ(self):
        self.enthalpy = 80.0
        self.temperature = 80.0
        self.heat_flow = self.enthalpy * self.mass_flow


def stream(**kwargs):
    return FakeStream(**kwargs)


class LinearHeatCapacityStream(FakeStream):
    def __init__(
        self,
        fluid=None,
        pp_instance=None,
        pressure=None,
        temperature=None,
        enthalpy=None,
        quality=None,
        mass_flow=1.0,
        P_crit=100.0,
    ):
        if enthalpy is None and temperature is not None:
            enthalpy = temperature
        if temperature is None and enthalpy is not None:
            temperature = enthalpy
        if enthalpy is None and quality is not None:
            enthalpy = 90.0
            temperature = 90.0
        super().__init__(
            fluid=fluid,
            pp_instance=pp_instance,
            pressure=pressure,
            temperature=temperature,
            enthalpy=enthalpy,
            quality=quality,
            mass_flow=mass_flow,
            P_crit=P_crit,
        )
        self.heat_capacity = self.temperature

    def props_in_mPT(self):
        self.enthalpy = self.temperature
        self.heat_capacity = self.temperature
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPH(self):
        self.temperature = self.enthalpy
        self.heat_capacity = self.enthalpy
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPQ(self):
        self.enthalpy = 90.0
        self.temperature = 90.0
        self.heat_capacity = 90.0
        self.heat_flow = self.enthalpy * self.mass_flow


class ExpansionValveInletStream(FakeStream):
    def props_in_mPH(self):
        self.temperature = self.enthalpy
        self.heat_flow = self.enthalpy * self.mass_flow


class ExpansionValveOutletStream(FakeStream):
    def props_in_mPQ(self):
        self.enthalpy = 5.0
        self.temperature = 2.5
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPH(self):
        self.temperature = self.enthalpy / 2.0
        self.heat_flow = self.enthalpy * self.mass_flow


def test_air_mass_flow_sets_sink_flow_from_gas_cooler_duty(monkeypatch):
    monkeypatch.setattr(custom_relationship_module, "Stream", FakeStream)
    gas_cooler = SimpleNamespace(
        duty=160.0,
        dTmin=20.0,
        source_in=stream(
            fluid="Butane",
            pressure=120.0,
            temperature=120.0,
            enthalpy=120.0,
            mass_flow=1.0,
        ),
        source_out=stream(
            fluid="Butane",
            pressure=100.0,
            temperature=80.0,
            enthalpy=80.0,
            mass_flow=1.0,
        ),
        sink_in=stream(
            fluid="Air",
            pressure=1.0,
            temperature=20.0,
            enthalpy=20.0,
            mass_flow=0.0,
        ),
        sink_out=stream(fluid="Air", pressure=1.0, enthalpy=0.0, mass_flow=0.0),
    )
    gas_cooler.source_in.heat_flow = 280.0
    gas_cooler.source_out.heat_flow = 120.0
    target = SimpleNamespace(
        sink_in=stream(fluid="Air", pressure=1.0, temperature=20.0, enthalpy=20.0),
        sink_out=stream(fluid="Air", pressure=1.0, enthalpy=0.0),
        source_in=gas_cooler.source_in,
        source_out=gas_cooler.source_out,
    )

    mass_flow_candidates = CustomRelationship(
        source_obj=[gas_cooler],
        target_obj=[target],
        relation=2,
    ).solve()

    assert mass_flow_candidates[0] == pytest.approx(2.0)
    assert target.sink_in.mass_flow == pytest.approx(2.0)
    assert target.sink_out.mass_flow == pytest.approx(2.0)
    assert target.sink_out.enthalpy == pytest.approx(100.0)


def test_air_mass_flow_considers_internal_and_saturation_pinch_candidates(monkeypatch):
    monkeypatch.setattr(custom_relationship_module, "Stream", LinearHeatCapacityStream)
    gas_cooler = SimpleNamespace(
        duty=8000.0,
        dTmin=20.0,
        source_in=LinearHeatCapacityStream(
            fluid="Butane",
            pressure=50.0,
            temperature=120.0,
            enthalpy=120.0,
            mass_flow=1.0,
            P_crit=100.0,
        ),
        source_out=LinearHeatCapacityStream(
            fluid="Butane",
            pressure=40.0,
            temperature=80.0,
            enthalpy=80.0,
            mass_flow=1.0,
            P_crit=100.0,
        ),
        sink_in=LinearHeatCapacityStream(
            fluid="Air",
            pressure=1.0,
            temperature=20.0,
            enthalpy=20.0,
            mass_flow=0.0,
        ),
        sink_out=LinearHeatCapacityStream(
            fluid="Air",
            pressure=1.0,
            enthalpy=0.0,
            mass_flow=0.0,
        ),
    )
    gas_cooler.source_in.heat_capacity = 50.0
    gas_cooler.source_out.heat_capacity = 150.0
    target = SimpleNamespace(
        sink_in=LinearHeatCapacityStream(
            fluid="Air", pressure=1.0, temperature=20.0, enthalpy=20.0
        ),
        sink_out=LinearHeatCapacityStream(fluid="Air", pressure=1.0, enthalpy=0.0),
        source_in=gas_cooler.source_in,
        source_out=gas_cooler.source_out,
    )

    mass_flow_candidates = CustomRelationship(
        source_obj=[gas_cooler],
        target_obj=[target],
        relation=2,
    ).solve()

    assert mass_flow_candidates == pytest.approx([100.0, 133.0, 159.4])
    assert target.sink_in.mass_flow == pytest.approx(159.4)
    assert target.sink_out.mass_flow == pytest.approx(159.4)


def test_air_mass_flow_uses_endpoint_bound_search_and_bisect_fallback(monkeypatch):
    monkeypatch.setattr(custom_relationship_module, "Stream", LinearHeatCapacityStream)
    root_methods = []

    def fake_root_scalar(**kwargs):
        root_methods.append(kwargs["method"])
        return {"root": 110.0 if kwargs["method"] == "brentq" else 100.0}

    monkeypatch.setattr(custom_relationship_module, "root_scalar", fake_root_scalar)

    gas_cooler = SimpleNamespace(
        duty=8000.0,
        dTmin=20.0,
        source_in=LinearHeatCapacityStream(
            fluid="Butane",
            pressure=50.0,
            temperature=120.0,
            enthalpy=120.0,
            mass_flow=1.0,
            P_crit=100.0,
        ),
        source_out=LinearHeatCapacityStream(
            fluid="Butane",
            pressure=40.0,
            temperature=80.0,
            enthalpy=80.0,
            mass_flow=1.0,
            P_crit=100.0,
        ),
        sink_in=LinearHeatCapacityStream(
            fluid="Air",
            pressure=1.0,
            temperature=20.0,
            enthalpy=20.0,
            mass_flow=0.0,
        ),
        sink_out=LinearHeatCapacityStream(fluid="Air", pressure=1.0, enthalpy=0.0),
    )
    gas_cooler.source_in.heat_capacity = 50.0
    gas_cooler.source_out.heat_capacity = 80.0
    target = SimpleNamespace(
        sink_in=LinearHeatCapacityStream(
            fluid="Air", pressure=1.0, temperature=20.0, enthalpy=20.0
        ),
        sink_out=LinearHeatCapacityStream(fluid="Air", pressure=1.0, enthalpy=0.0),
        source_in=gas_cooler.source_in,
        source_out=gas_cooler.source_out,
    )

    mass_flow_candidates = CustomRelationship(
        source_obj=[gas_cooler],
        target_obj=[target],
        relation=2,
    ).solve()

    assert root_methods == ["brentq", "bisect"]
    assert mass_flow_candidates[1] == pytest.approx(133.0)
    assert target.sink_in.mass_flow == pytest.approx(max(mass_flow_candidates))


def test_air_source_mass_flow_sets_source_flow_from_evaporator_duty(monkeypatch):
    monkeypatch.setattr(custom_relationship_module, "Stream", FakeStream)
    evaporator = SimpleNamespace(
        sink_in=stream(fluid="Butane", pressure=1.0, enthalpy=100.0, mass_flow=1.0),
        sink_out=stream(
            fluid="Butane",
            pressure=1.0,
            temperature=75.0,
            enthalpy=200.0,
            mass_flow=1.0,
        ),
        source_out=stream(
            fluid="Air",
            pp_instance="src",
            pressure=1.0,
            temperature=20.0,
            enthalpy=100.0,
        ),
    )
    evaporator.sink_in.heat_flow = 100.0
    evaporator.sink_out.heat_flow = 200.0
    target = SimpleNamespace(
        source_in=stream(fluid="Air", pressure=1.0, temperature=20.0, enthalpy=0.0),
        source_out=stream(fluid="Air", pressure=1.0, temperature=20.0, enthalpy=0.0),
    )

    result = CustomRelationship(
        source_obj=[evaporator],
        target_obj=[target],
        relation=3,
    ).solve()

    assert result is True
    assert target.source_in.mass_flow == pytest.approx(2.0)
    assert target.source_out.mass_flow == pytest.approx(2.0)


def test_economiser_expansion_valve_relationship_uses_discharge_state_when_supercritical():
    comp1 = SimpleNamespace(flow_out=stream(pressure=3.0, mass_flow=4.0))
    economiser = SimpleNamespace(dP_h=0.2, dP_c=0.3, source_out=stream(mass_flow=0.0))
    comp_discharge = stream(pressure=10.0, enthalpy=500.0)
    expansion_valve = SimpleNamespace(
        dTmin=5.0,
        flow_in=stream(P_crit=0.0, pressure=0.0, enthalpy=0.0, mass_flow=1.0),
        flow_out=stream(P_crit=1.0, pressure=0.0, enthalpy=0.0, mass_flow=1.0),
    )
    splitter = SimpleNamespace(r_ratio=0.25, flow_in=stream(pressure=0.0, enthalpy=0.0))

    result = CustomRelationship(
        source_obj=[comp1, economiser, comp_discharge, 0.5],
        target_obj=[expansion_valve, splitter],
        relation=4,
    ).solve()

    assert result is True
    assert expansion_valve.flow_in.mass_flow == pytest.approx(0.0)
    assert economiser.source_out.mass_flow == pytest.approx(comp1.flow_out.mass_flow)
    assert expansion_valve.flow_in.pressure == pytest.approx(9.3)
    assert expansion_valve.flow_out.pressure == pytest.approx(3.3)
    assert expansion_valve.flow_in.enthalpy == pytest.approx(500.0)
    assert expansion_valve.flow_out.enthalpy == pytest.approx(500.0)
    assert splitter.flow_in.pressure == pytest.approx(9.5)
    assert splitter.flow_in.enthalpy == pytest.approx(500.0)


def test_economiser_expansion_valve_relationship_solves_subcritical_temperature_drop():
    comp1 = SimpleNamespace(flow_out=stream(pressure=3.0, mass_flow=4.0))
    economiser = SimpleNamespace(dP_h=0.2, dP_c=0.3, source_out=stream(mass_flow=0.0))
    comp_discharge = stream(pressure=10.0, enthalpy=50.0)
    expansion_valve = SimpleNamespace(
        dTmin=10.0,
        flow_in=ExpansionValveInletStream(
            P_crit=100.0,
            pressure=0.0,
            enthalpy=0.0,
            mass_flow=1.0,
        ),
        flow_out=ExpansionValveOutletStream(
            P_crit=100.0,
            pressure=0.0,
            enthalpy=0.0,
            mass_flow=1.0,
        ),
    )
    splitter = SimpleNamespace(r_ratio=0.25, flow_in=stream(pressure=0.0, enthalpy=0.0))

    result = CustomRelationship(
        source_obj=[comp1, economiser, comp_discharge, 0.5],
        target_obj=[expansion_valve, splitter],
        relation=4,
    ).solve()

    assert result is True
    assert expansion_valve.flow_in.mass_flow == pytest.approx(1.0)
    assert economiser.source_out.mass_flow == pytest.approx(5.0)
    assert expansion_valve.flow_in.enthalpy == pytest.approx(20.0)
    assert expansion_valve.flow_out.enthalpy == pytest.approx(20.0)
    assert expansion_valve.flow_in.temperature - expansion_valve.flow_out.temperature == (
        pytest.approx(10.0)
    )
