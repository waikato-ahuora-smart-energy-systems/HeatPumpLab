"""Screening estimate of the refrigerant charge (HeatPumpLab/cycles/common/charge.py)."""
import pytest

from HeatPumpLab.cycles.common.charge import (
    GEOMETRY,
    ChargeSection,
    estimate_charge,
    void_fraction_density,
)
from HeatPumpLab.cycles.DualPressureHeatPump import (
    DualPressureConfig,
    DualPressureHeatPump,
)
from HeatPumpLab.cycles.EconomiserHeatPump import EconomiserConfig, EconomiserHeatPump
from HeatPumpLab.cycles.MultiPressureHeatPump import (
    MultiPressureConfig,
    MultiPressureHeatPump,
)

pytest.importorskip("CoolProp")


def _constant_section(service="air", rho=100.0, duty=10.0, dT=20.0, n=4, sides=1):
    return ChargeSection("hx", service, duty, [dT] * n, [[(rho, rho)] * n] * sides)


def test_void_fraction_density_limits():
    assert void_fraction_density(600.0, 600.0, 20.0) == pytest.approx(600.0)
    assert void_fraction_density(20.0, 600.0, 20.0) == pytest.approx(20.0)
    # slip holds back liquid, so two-phase density exceeds the homogeneous value
    assert void_fraction_density(60.0, 600.0, 20.0) > 60.0


def test_single_phase_section_matches_hand_calculation():
    e = estimate_charge([_constant_section()], Q_heat=10.0)
    U, v_per_a = GEOMETRY["typical"]["air"]
    # m = Q / (U dT) (V/A) rho, per kW of heating, in g
    assert e.exchanger_central == pytest.approx(1000.0 / (U * 20.0) * v_per_a * 100.0 * 1000.0)
    assert e.by_exchanger["hx"]["area_m2_per_kW"] == pytest.approx(1000.0 / (U * 20.0))
    assert e.exchanger_low < e.exchanger_central < e.exchanger_high
    assert e.system_low > e.exchanger_low and e.system_high > e.exchanger_high


def test_both_recuperator_sides_hold_charge():
    one = estimate_charge([_constant_section("recuperator")], Q_heat=10.0)
    two = estimate_charge([_constant_section("recuperator", sides=2)], Q_heat=10.0)
    assert two.exchanger_central == pytest.approx(2.0 * one.exchanger_central)
    assert two.by_exchanger["hx"]["volume_L_per_kW"] == pytest.approx(2.0 * one.by_exchanger["hx"]["volume_L_per_kW"])


def test_specific_charge_does_not_depend_on_capacity_and_kg_does():
    small = estimate_charge([_constant_section(duty=10.0)], Q_heat=10.0)
    large = estimate_charge([_constant_section(duty=1000.0)], Q_heat=1000.0)
    assert large.exchanger_central == pytest.approx(small.exchanger_central)
    assert small.system_kg(1000.0)[1] == pytest.approx(small.system_central)


def test_geometry_can_be_replaced():
    doubled = {case: {s: (U, 2.0 * v) for s, (U, v) in services.items()} for case, services in GEOMETRY.items()}
    base = estimate_charge([_constant_section()], Q_heat=10.0)
    assert estimate_charge([_constant_section()], Q_heat=10.0, geometry=doubled).exchanger_central == pytest.approx(
        2.0 * base.exchanger_central)


@pytest.mark.parametrize("dT", [0.0, -1.0])
def test_rejects_no_driving_force(dT):
    with pytest.raises(ValueError, match="driving force"):
        estimate_charge([_constant_section(dT=dT)], Q_heat=10.0)


def test_rejects_no_heating():
    with pytest.raises(ValueError, match="heating"):
        estimate_charge([_constant_section()], Q_heat=0.0)


def test_range_brackets_the_propane_reference_of_dawoud_et_al():
    """Sanity check, not validation: the 12 kW propane reference of that paper (A-8/W55, 14.3 g/kW for the system)
    should fall inside the range. It uses a plate condenser and a fin-and-tube evaporator."""
    from HeatPumpLab.cycles.common.charge import _pair
    from HeatPumpLab.cycles.common.dryer import Fluid

    fl = Fluid("Propane", "coolprop")
    n = 24

    def rho(P, h0, h1):
        def sat(q):
            return fl.PQ(P, q)["D"]
        return [_pair(fl.PH(P, h0 + (h1 - h0) * (i + 0.5) / n), sat) for i in range(n)]
    Pc, Pe = fl.TQ(60.0, 1.0)["P"], fl.TQ(-15.0, 1.0)["P"]
    h_liq, h_dew, h_suc = fl.TQ(58.0, 0.0)["h"], fl.TQ(60.0, 1.0)["h"], fl.TQ(-15.0, 1.0)["h"]
    Q_cond = 12.0
    Q_evap = Q_cond * (1.0 - 1.0 / 1.7)
    sections = [ChargeSection("condenser", "liquid", Q_cond, [8.0] * n, [rho(Pc, h_liq, h_dew)]),
                ChargeSection("evaporator", "air", Q_evap, [7.0] * n, [rho(Pe, h_liq, h_suc)])]
    e = estimate_charge(sections, Q_heat=Q_cond)
    assert e.system_low < 14.3 < e.system_high


@pytest.fixture(scope="module")
def multi():
    return MultiPressureHeatPump(MultiPressureConfig(property_backend="coolprop"))


@pytest.mark.parametrize("n", [2, 3])
def test_multi_pressure_has_every_exchanger(multi, n):
    e = multi.estimate_charge(multi.solve(n_pressures=n))
    assert len(e.by_exchanger) == n + 1
    assert 0.0 < e.exchanger_low < e.exchanger_central < e.exchanger_high
    assert e.exchanger_central == pytest.approx(sum(row["central"] for row in e.by_exchanger.values()))


def test_dual_pressure_matches_two_pressures(multi):
    dual = DualPressureHeatPump(DualPressureConfig(property_backend="coolprop", T_cond_low=120.0))
    a = dual.estimate_charge(dual.solve())
    b = multi.estimate_charge(multi.solve(n_pressures=2, T_cond=[120.0]))
    assert a.exchanger_central == pytest.approx(b.exchanger_central, rel=1e-9)


def test_subcooling_against_air_costs_charge():
    """Flash subcools to near the air inlet, while the single stage leaves the condenser hot and light."""
    hp = EconomiserHeatPump(EconomiserConfig(property_backend="coolprop"))
    flash = hp.estimate_charge(hp.solve(arrangement="flash"))
    single = hp.estimate_charge(hp.solve(arrangement="single"))
    assert flash.exchanger_central > single.exchanger_central


def test_supercritical_design_finds_its_gas_cooler_path():
    hp = EconomiserHeatPump(EconomiserConfig(property_backend="coolprop", fluid="n-Butane"))
    r = hp.solve(arrangement="flash", T_air_inlet=90.0, T_air_tar=140.0)
    assert r.success and r.report.get("supercritical")
    assert hp.estimate_charge(r).exchanger_central > 0.0


def test_failed_design_is_rejected(multi):
    with pytest.raises(ValueError, match="successful"):
        multi.estimate_charge(multi.solve(T_cond=[48.0, 110.0]))


def _rhp(**kwargs):
    from HeatPumpLab.cycles.RecuperativeHeatPump import RecuperativeHP

    args = {"refrigerant": "Butane", "T_source_in": 40.0, "Tref_evap_inlet": 10.0, "T_air_inlet": 15.0, "T_air_tar": 150.0,
                "dT_GC_min_tar": 20.0, "dTmin_rhx": 5.0, "n_comp": 2, "dT_rhx1_cold": 60.0, "dT_rhx2_cold": 0.0, "comp2_power_frac": 0.5,
                "sink": "Air", "source": "water"}
    return RecuperativeHP(**{**args, **kwargs}).solve()


def test_recuperative_counts_only_active_exchangers(monkeypatch):
    monkeypatch.setenv("HEATPUMPLAB_PROPERTY_BACKEND", "coolprop")
    hp = _rhp()
    assert hp.converged
    e = hp.estimate_charge()
    assert {"gas_cooler1", "gas_cooler2", "evaporator", "rhx1"} <= set(e.by_exchanger)
    assert "rhx2" not in e.by_exchanger
    assert e.exchanger_low < e.exchanger_central < e.exchanger_high
    assert not any("zeotrope" in w for w in e.warnings)


def test_recuperative_zeotrope_warns(monkeypatch):
    monkeypatch.setenv("HEATPUMPLAB_PROPERTY_BACKEND", "coolprop")
    hp = _rhp(refrigerant="Propane;Isopentane|0.2;0.8|mass", T_air_tar=200.0)
    assert hp.converged
    e = hp.estimate_charge()
    assert any("zeotrope" in w for w in e.warnings)
    assert e.exchanger_central > 0.0


def test_unconverged_recuperative_is_rejected(monkeypatch):
    monkeypatch.setenv("HEATPUMPLAB_PROPERTY_BACKEND", "coolprop")
    hp = _rhp()
    hp.converged = False
    with pytest.raises(ValueError, match="converged"):
        hp.estimate_charge()
