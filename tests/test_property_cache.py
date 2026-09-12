from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.solver import cycle_solver
from HeatPumpLab.utils.properties_package import PropertyPackage
from HeatPumpLab.utils.property_cache import property_call_cache
from HeatPumpLab.utils.solve_profile import profile_solve


def _fake_property_package(refprop):
    package = object.__new__(PropertyPackage)
    package.RP = refprop
    package.fld_name = "Butane"
    package.fld_composition = [1.0]
    package.UNITS = 0
    package.iMass = 0
    package.iFlag = 0
    package.molar_mass = 2.0
    package.TC_to_K = 273.15
    package.DEFAULT_VAL = -9999990.0
    return package


class CountingRefprop:
    def __init__(self):
        self.std_calls = 0
        self.specified_calls = 0

    def SETFLUIDSdll(self, name):
        return 0

    def ABFLSHdll(self, *args):
        self.std_calls += 1
        return SimpleNamespace(P=1000.0, h=20.0, q=-1.0, T=300.15, s=4.0, Cp=6.0, D=2.0, ierr=0, herr="")

    def REFPROPdll(self, *args):
        self.specified_calls += 1
        return SimpleNamespace(Output=[1.0, 2.0, 3.0], ierr=0, herr="")


def test_property_cache_reuses_identical_property_requests_within_context():
    refprop = CountingRefprop()
    package = _fake_property_package(refprop)

    with property_call_cache():
        first = package.get_std_properties("PT", 1.0, 20.0, [1.0])
        second = package.get_std_properties("PT", 1.0, 20.0, [1.0])
        third = package.get_specified_properties("PRANDTL;VIS;TCX", "PH", 1.0, 2.0)
        fourth = package.get_specified_properties("PRANDTL;VIS;TCX", "PH", 1.0, 2.0)

    assert refprop.std_calls == 1
    assert refprop.specified_calls == 1
    assert second == first
    assert fourth == third


def test_property_cache_returns_independent_result_objects():
    refprop = CountingRefprop()
    package = _fake_property_package(refprop)

    with property_call_cache():
        first = package.get_std_properties("PT", 1.0, 20.0, [1.0])
        first["T"] = -100.0
        second = package.get_std_properties("PT", 1.0, 20.0, [1.0])

    assert second["T"] == pytest.approx(27.0)


def test_property_cache_does_not_share_results_between_contexts():
    refprop = CountingRefprop()
    package = _fake_property_package(refprop)

    with property_call_cache():
        package.get_std_properties("PT", 1.0, 20.0, [1.0])

    with property_call_cache():
        package.get_std_properties("PT", 1.0, 20.0, [1.0])

    assert refprop.std_calls == 2


def test_profile_counts_actual_refprop_calls_when_cache_is_active():
    refprop = CountingRefprop()
    package = _fake_property_package(refprop)

    with profile_solve() as profile:
        with property_call_cache():
            package.get_std_properties("PT", 1.0, 20.0, [1.0])
            package.get_std_properties("PT", 1.0, 20.0, [1.0])

    report = profile.as_dict()
    assert report["counts"]["refprop_std"] == 1
    assert refprop.std_calls == 1


def test_recuperative_solver_uses_one_property_cache_per_solve(monkeypatch):
    refprop = CountingRefprop()
    package = _fake_property_package(refprop)
    hp = SimpleNamespace(debug=False)

    def pressure_search(solved_hp):
        package.get_std_properties("PT", 1.0, 20.0, [1.0])
        package.get_std_properties("PT", 1.0, 20.0, [1.0])
        return None

    monkeypatch.setattr(
        cycle_solver,
        "solve_at_minimum_discharge_pressure",
        pressure_search,
    )
    monkeypatch.setattr(
        cycle_solver,
        "evaporator_inlet_temperature_offset",
        lambda solved_hp: 0.0,
    )
    monkeypatch.setattr(
        cycle_solver,
        "populate_full_metrics",
        lambda solved_hp: solved_hp,
    )

    assert cycle_solver.solve_heat_pump_cycle(hp) is hp
    assert refprop.std_calls == 1
