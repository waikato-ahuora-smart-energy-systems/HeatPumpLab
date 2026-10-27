"""Ammonia-water properties (HeatPumpLab/utils/ammonia_water.py): the IAPWS G4-01 verification values, the flashes
and a comparison with REFPROP where it is installed."""
import pytest

pytest.importorskip("teqp")

from HeatPumpLab.utils.ammonia_water import (  # noqa: E402
    AmmoniaWater,
    AmmoniaWaterStream,
    mass_fraction,
    mole_fraction,
)


@pytest.fixture(scope="module")
def model():
    return AmmoniaWater()


# IAPWS Guideline G4-01 (2001), section 8: x is the ammonia mole fraction, rho in mol/dm3, p in MPa
TABLE6 = [(0.1, 600, 35, 32.1221333), (0.1, 600, 4, 12.7721090), (0.5, 500, 32, 21.3208159),
          (0.5, 500, 1, 3.6423080), (0.9, 400, 30, 22.2830797), (0.9, 400, 0.5, 1.5499708)]
TABLE7 = [(0.2, 300, 0.040710, 0.9360, 51.941, 0.01640), (0.4, 400, 2.5545, 0.9363, 43.318, 0.8608),
          (0.6, 500, 16.698, 0.7844, 25.459, 8.86)]
TABLE8 = [(0.2, 300, 0.00437062, 0.010672, 55.16434, 0.00175506), (0.4, 400, 0.394694, 0.051541, 50.83187, 0.122658),
          (0.6, 500, 6.52607, 0.22135, 39.93714, 2.00730)]


@pytest.mark.parametrize("x,T,rho,p", TABLE6)
def test_iapws_single_phase_pressures(model, x, T, rho, p):
    assert model.pressure(T, rho * 1000, x) / 1e6 == pytest.approx(p, abs=1e-7)


@pytest.mark.parametrize("x,T,p,y,rhoL,rhoV", TABLE7)
def test_iapws_bubble_points(model, x, T, p, y, rhoL, rhoV):
    pb, yb = model.bubble_p(T, x)
    assert pb / 1e6 == pytest.approx(p, rel=5e-5) and yb == pytest.approx(y, abs=5e-5)
    assert model.density(T, pb, x, "liquid") / 1000 == pytest.approx(rhoL, rel=5e-5)
    assert model.density(T, pb, yb, "vapour") / 1000 == pytest.approx(rhoV, rel=1e-3)


@pytest.mark.parametrize("y,T,p,x,rhoL,rhoV", TABLE8)
def test_iapws_dew_points(model, y, T, p, x, rhoL, rhoV):
    pd, xd = model.dew_p(T, y)
    assert pd / 1e6 == pytest.approx(p, rel=5e-6) and xd == pytest.approx(x, abs=5e-6)
    assert model.density(T, pd, xd, "liquid") / 1000 == pytest.approx(rhoL, rel=5e-6)
    assert model.density(T, pd, y, "vapour") / 1000 == pytest.approx(rhoV, rel=5e-6)


def test_composition_conversions():
    for w in (0.01, 0.3, 0.5, 0.97):
        assert mass_fraction(mole_fraction(w)) == pytest.approx(w, rel=1e-14)
    assert mole_fraction(0.5) == pytest.approx(0.5140532624, rel=1e-9)


def test_saturation_temperatures_invert(model):
    for x in (0.3, 0.5, 0.7):
        Tb, _ = model.bubble_T(2.5e6, x)
        Td, _ = model.dew_T(2.5e6, x)
        assert Tb < Td
        assert model.bubble_p(Tb, x)[0] == pytest.approx(2.5e6, rel=1e-8)
        assert model.dew_p(Td, x)[0] == pytest.approx(2.5e6, rel=1e-8)


def test_stream_flashes_round_trip(model):
    s = AmmoniaWaterStream(0.5, model)
    two = s.PT(2.5, 150.0)
    assert two["phase"] == "two_phase" and 0 < two["q"] < 1 and two["w_liq"] < 0.5 < two["w_vap"]
    # the lever rule holds on a mass basis
    assert two["q"] * two["w_vap"] + (1 - two["q"]) * two["w_liq"] == pytest.approx(0.5, rel=1e-9)
    assert s.PH(2.5, two["h"])["T"] == pytest.approx(150.0, abs=1e-6)
    assert s.PS(2.5, two["s"])["T"] == pytest.approx(150.0, abs=1e-6)
    for T in (60.0, 230.0):                              # liquid, then superheated vapour
        st = s.PT(2.5, T)
        assert s.PH(2.5, st["h"])["T"] == pytest.approx(T, abs=1e-6)
    q = s.PQ(2.5, 0.3)
    assert q["q"] == pytest.approx(0.3, abs=1e-8)
    tq = s.TQ(120.0, 0.0)
    assert tq["T"] == pytest.approx(120.0, abs=1e-6) and tq["q"] == 0.0
    # the fast isobar table follows the flash closely
    lo, hi = s.PT(2.5, 40.0)["h"], s.PT(2.5, 240.0)["h"]
    for k in range(1, 20):
        h = lo + (hi - lo) * k / 20
        assert s.T_ph_fast(2.5, h) == pytest.approx(s.PH(2.5, h)["T"], abs=0.02)


def test_close_to_refprop_where_installed():
    """REFPROP 10 uses the newer Bell et al. (2018) ammonia-water model: dew points agree within about 0.4 K, bubble
    points within 0.5 K at 30 % ammonia but about 4 K lower here at 50 %."""
    try:
        from HeatPumpLab.utils.coolprop_properties import make_property_package
        pp = make_property_package("Ammonia;Water|0.3;0.7|mass", backend="refprop")
        comp = pp.fld_composition
        rp_b = pp.get_std_properties("PQ", 2.5, 0.0, comp)["T"]
        rp_d = pp.get_std_properties("PQ", 2.5, 1.0, comp)["T"]
    except Exception as exc:
        pytest.skip(f"REFPROP not available ({exc})")
    s = AmmoniaWaterStream(0.3)
    Tb, Td, _, _ = s.bubble_dew(2.5)
    assert Tb - 273.15 == pytest.approx(rp_b, abs=0.6)
    assert Td - 273.15 == pytest.approx(rp_d, abs=0.5)
