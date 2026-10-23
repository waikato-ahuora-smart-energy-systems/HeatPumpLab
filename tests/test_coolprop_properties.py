"""CoolProp property backend (HeatPumpLab/utils/coolprop_properties.py).

The plumbing tests (units, input pairs, phase labels, errors, backend switch) run everywhere with a
fake CoolProp module backed by the toy fluid of test_conventional_heat_pump. The CoolProp tests need
CoolProp >= 8; the comparison tests also need REFPROP.
"""
import math
import types

import pytest

from HeatPumpLab.utils import coolprop_properties as cpp
from HeatPumpLab.utils.properties_package import RefpropError
from test_conventional_heat_pump import ToyFluid


# ---------------------------------------------------------------- fake CoolProp (SI in, SI out)


def fake_coolprop():
    toy = ToyFluid()
    m = types.SimpleNamespace()
    names = ["PT_INPUTS", "HmassP_INPUTS", "PSmass_INPUTS", "PQ_INPUTS", "QT_INPUTS",
             "DmassT_INPUTS", "DmassP_INPUTS", "SmassT_INPUTS"]
    for i, n in enumerate(names, start=1):
        setattr(m, n, i)
    phases = dict(iphase_liquid=0, iphase_supercritical=1, iphase_supercritical_gas=2,
                  iphase_supercritical_liquid=3, iphase_gas=5, iphase_twophase=6)
    for k, v in phases.items():
        setattr(m, k, v)
    # like CoolProp 8: constants only, no get_input_pair_index / get_phase_index
    m.reference_calls = []
    m.set_reference_state = lambda fluid, ref: m.reference_calls.append((fluid, ref))
    m.get_global_param_string = lambda k: "8.0.0-fake"

    class AbstractState:
        def __init__(self, backend, fluids):
            if fluids != "n-Butane":
                raise ValueError(f"unknown fluid {fluids}")
            self.st = None

        def update(self, pair, a, b):
            if pair == m.PT_INPUTS:
                st = toy.ref_PT(a / 1e6, b - 273.15)
            elif pair == m.HmassP_INPUTS:
                st = toy.ref_PH(b / 1e6, a / 1e3)
            elif pair == m.PSmass_INPUTS:
                st = toy.ref_PS(a / 1e6, b / 1e3)
            elif pair == m.PQ_INPUTS:
                st = toy.ref_PQ(a / 1e6, b)
            elif pair == m.QT_INPUTS:
                st = toy.ref_TQ(b - 273.15, a)
            else:
                raise ValueError("pair not supported by the fake")
            self.st = st

        def T(self): return self.st["T"] + 273.15
        def p(self): return self.st["P"] * 1e6
        def hmass(self): return self.st["h"] * 1e3
        def smass(self): return self.st["s"] * 1e3
        def rhomass(self): return self.st["D"]
        def cpmass(self): return 1700.0
        def Q(self): return self.st["q"] if 0 <= self.st["q"] <= 1 else -1.0
        def phase(self): return {"liquid": 0, "superheated_vapour": 5}.get(self.st["phase"], 6)
        def molar_mass(self): return 0.05812
        def p_critical(self): return toy.P_crit * 1e6
        def T_critical(self): return toy.T_crit + 273.15
        def rhomass_critical(self): return 228.0
        def pmax(self): return 200e6

    m.AbstractState = AbstractState
    return m, toy


@pytest.fixture
def fake(monkeypatch):
    m, toy = fake_coolprop()
    monkeypatch.setattr(cpp, "_CP", m)
    return m, toy


def test_fluid_strings_and_names():
    assert cpp.parse_fluid("Butane") == (["n-Butane"], [1.0], "mol")
    names, x, basis = cpp.parse_fluid("Propane;Isopentane|0.2;0.8|mass")
    assert names == ["Propane", "Isopentane"] and x == pytest.approx([0.2, 0.8]) and basis == "mass"
    assert cpp.coolprop_name("ISOBUTANE.FLD") == "IsoButane"
    assert cpp.coolprop_name("R1233zd(E)") == "R1233zd(E)"
    with pytest.raises(ValueError):
        cpp.parse_fluid("Propane;Isopentane|0.2|mass")


def test_backend_switch(monkeypatch):
    monkeypatch.delenv(cpp.ENV_VAR, raising=False)
    assert cpp.property_backend() == "coolprop"                   # the default
    monkeypatch.setenv(cpp.ENV_VAR, "REFPROP")
    assert cpp.property_backend() == "refprop"
    assert cpp.property_backend("coolprop") == "coolprop"        # an explicit choice wins
    with pytest.raises(ValueError):
        cpp.property_backend("nist")


def test_missing_coolprop_raises_refprop_error(monkeypatch):
    monkeypatch.setattr(cpp, "_CP", None)
    with pytest.raises(RefpropError, match="CoolProp is not installed"):
        cpp.CoolPropPackage("Butane")


def test_package_units_pairs_and_phases(fake):
    m, toy = fake
    pkg = cpp.CoolPropPackage("Butane")
    assert ("n-Butane", "IIR") in m.reference_calls
    assert pkg.P_crit == pytest.approx(toy.P_crit) and pkg.T_crit == pytest.approx(toy.T_crit, abs=1e-9)
    assert pkg.molar_mass == pytest.approx(58.12)
    P = toy.Psat(50.0)
    sat_l = pkg.get_std_properties("TQ", 50.0, 0.0)
    assert sat_l["P"] == pytest.approx(P) and sat_l["phase"] == "saturated_liquid" and sat_l["q"] == 0.0
    sat_v = pkg.get_std_properties("PQ", P, 1.0)
    assert sat_v["T"] == pytest.approx(50.0) and sat_v["phase"] == "saturated_vapour"
    gas = pkg.get_std_properties("PT", P, 70.0)
    assert gas["phase"] == "superheated_vapour" and gas["q"] == cpp.DEFAULT_VAL and gas["q_refprop"] == 998.0
    assert gas["h"] == pytest.approx(toy.ref_PT(P, 70.0)["h"])
    liq = pkg.get_std_properties("PT", P, 30.0)
    assert liq["phase"] == "liquid" and liq["q_refprop"] == -1.0
    mid = pkg.get_std_properties("PH", P, 0.5 * (sat_l["h"] + sat_v["h"]))
    assert mid["phase"] == "two_phase" and mid["q"] == pytest.approx(0.5)
    back = pkg.get_std_properties("PS", P, gas["s"])
    assert back["T"] == pytest.approx(70.0, abs=1e-6)
    hp = pkg.get_std_properties("HP", gas["h"], P)                # reversed pair order
    assert hp["T"] == pytest.approx(70.0, abs=1e-6)


def test_flash_failures_raise_refprop_error(fake):
    pkg = cpp.CoolPropPackage("Butane")
    with pytest.raises(RefpropError, match="Unsupported"):
        pkg.get_std_properties("XY", 1.0, 2.0)
    with pytest.raises(RefpropError, match="flash failed"):
        pkg.get_std_properties("TD", 20.0, 500.0)                 # the fake has no TD flash


def test_conventional_cycle_runs_on_the_coolprop_backend(fake):
    from HeatPumpLab.cycles.ConventionalHeatPump import (
        ConventionalHeatPump, ConventionalHPConfig, EfficiencyCurve,
    )
    cfg = ConventionalHPConfig(refrigerant="Butane", T_evap=10.0, T_cond=80.0, dT_subcool=5.0,
                               dT_superheat=5.0, loss_UA_kWK=None, property_backend="coolprop",
                               compressor=EfficiencyCurve(eta_is=0.7, motor_eff=0.9))
    via_coolprop = ConventionalHeatPump(cfg).solve()
    direct = ConventionalHeatPump(cfg, props=ToyFluid()).solve()
    assert via_coolprop.success, via_coolprop.message
    assert via_coolprop.report["property_backend"] == "coolprop"
    assert via_coolprop.COP == pytest.approx(direct.COP, rel=1e-9)


# ---------------------------------------------------------------- real CoolProp


@pytest.fixture
def coolprop():
    CP = pytest.importorskip("CoolProp.CoolProp")
    version = cpp.coolprop_version() or "0"
    if int(version.split(".")[0]) < 8:
        pytest.skip(f"CoolProp {version} < 8")
    return CP


def test_coolprop_butane_saturation(coolprop):
    pkg = cpp.CoolPropPackage("Butane")
    assert pkg.T_crit == pytest.approx(151.98, abs=0.05)
    assert pkg.P_crit == pytest.approx(3.796, abs=0.005)
    sat = pkg.get_std_properties("TQ", 0.0, 0.0)
    assert sat["h"] == pytest.approx(200.0, abs=1e-6) and sat["s"] == pytest.approx(1.0, abs=1e-9)  # IIR
    assert pkg.get_std_properties("TQ", -0.5, 1.0)["P"] == pytest.approx(0.1013, rel=0.01)      # NBP


def test_coolprop_mixture_dew_and_bubble(coolprop):
    pkg = cpp.CoolPropPackage("Propane;Isopentane|0.2;0.8|mass")
    dew = pkg.get_std_properties("PQ", 0.5, 1.0)
    bub = pkg.get_std_properties("PQ", 0.5, 0.0)
    assert dew["phase"] == "saturated_vapour" and bub["phase"] == "saturated_liquid"
    assert dew["T"] > bub["T"] + 5.0                               # zeotropic glide


def test_conventional_cycle_coolprop_matches_refprop(coolprop, hpl):
    from HeatPumpLab.cycles.ConventionalHeatPump import ConventionalHeatPump, ConventionalHPConfig
    res = {b: ConventionalHeatPump(ConventionalHPConfig(property_backend=b)).solve() for b in ("refprop", "coolprop")}
    assert all(r.success for r in res.values()), {b: r.message for b, r in res.items()}
    rp, cp = res["refprop"].report, res["coolprop"].report
    for k in ("COP", "Q_heat_kW", "W_electric_kW", "m_dot_kgs", "P_cond_bar", "P_evap_bar"):
        assert cp[k] == pytest.approx(rp[k], rel=2e-3), k
    assert cp["T_discharge"] == pytest.approx(rp["T_discharge"], abs=0.05)
    assert math.isclose(cp["UA_loss_kWK"], rp["UA_loss_kWK"], rel_tol=0.02)
