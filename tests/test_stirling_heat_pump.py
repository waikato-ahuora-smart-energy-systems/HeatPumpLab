"""Reverse-Stirling heat pump block (HeatPumpLab/cycles/StirlingHeatPump)."""
import pytest

from HeatPumpLab.cycles.StirlingHeatPump import StirlingConfig, StirlingHeatPump

pytest.importorskip("CoolProp")


@pytest.fixture(scope="module")
def hp():
    return StirlingHeatPump()


def test_calibration_point_rfc_maasdam():
    """The default Carnot fraction reproduces the supplier-reported point: COP about 2 for 155 C steam from about
    75 C waste heat (hot and cold sides at those temperatures)."""
    m = StirlingHeatPump(StirlingConfig(dT_hot=0, dT_cold=0, T_source_in=80, T_source_out=75, T_air_inlet=150,
                                        T_air_tar=155, air_RH=0))
    assert m.solve().COP == pytest.approx(2.0, abs=0.05)


def test_one_unit_is_the_carnot_fraction(hp):
    r = hp.solve()
    u = r.report["units"][0]
    Th, Tc = 170 + 20 + 273.15, 50 - 5 + 273.15
    assert u["T_hot"] == pytest.approx(190.0) and u["T_cold"] == pytest.approx(45.0)
    assert r.COP == pytest.approx(0.37 * Th / (Th - Tc), rel=1e-12)
    assert r.q_air == pytest.approx(r.w_shaft + r.q_source, rel=1e-12)
    assert r.q_air == pytest.approx(hp.air.h(170) - hp.air.h(15), rel=1e-12)


def test_units_in_series_follow_the_air_glide(hp):
    cops = [hp.solve(n_units=n).COP for n in (1, 2, 3, 4)]
    assert all(b > a for a, b in zip(cops, cops[1:]))
    r = hp.solve(n_units=3)
    units = r.report["units"]
    assert [u["T_air_out"] for u in units][-1] == pytest.approx(170.0)
    assert r.q_air == pytest.approx(sum(u["q_air"] for u in units), rel=1e-12)
    assert r.w_shaft == pytest.approx(sum(u["w_el"] for u in units), rel=1e-12)


def test_optimised_boundaries_are_no_worse(hp):
    start = hp.solve(n_units=3)
    best = hp.optimise(n_units=3)
    b = best.report["optimiser"]["boundaries"]
    assert best.COP >= start.COP - 1e-12 and 15 < b[0] < b[1] < 170


def test_ratio_dependent_efficiency_and_warnings():
    m = StirlingHeatPump(StirlingConfig(eta_of_ratio=lambda tr: 0.5 / tr))
    u = m.solve().report["units"][0]
    assert u["eta_carnot"] == pytest.approx(0.5 / ((190 + 273.15) / (45 + 273.15)))
    hot = StirlingHeatPump(StirlingConfig(T_air_tar=200.0)).solve()
    assert any("field data" in w for w in hot.warnings)
    assert not StirlingHeatPump(StirlingConfig(n_units=5)).solve().success
    assert not StirlingHeatPump(StirlingConfig(n_units=2, boundaries=[200.0])).solve().success
