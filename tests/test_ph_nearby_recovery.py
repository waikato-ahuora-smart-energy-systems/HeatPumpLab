"""Internal PH fallback routing, bounded search, validation and real regressions."""

from collections import namedtuple
from types import SimpleNamespace

import pytest

from HeatPumpLab.utils.properties_package import PropertyPackage, RefpropError
from HeatPumpLab.utils.property_cache import property_call_cache


Flash = namedtuple("Flash", "P T D Dl Dv h e s q x y Cp ierr herr")


def state(p=1000., h=600.):
    return Flash(p, 300.+h/100, 1.6, 4., 1., h, h-10, 4., .5,
                 [.6, .4], [.4, .6], -9999990., 0, "")


def package(native):
    pp = object.__new__(PropertyPackage)
    pp.RP = native
    pp.fld_name = "A.FLD*B.FLD"
    pp.fld_composition = [.5, .5]
    pp.molar_mass = 1.
    pp.polish_ph = False
    return pp


@pytest.mark.parametrize("pair,inputs", [("PH", (1., 600.)), ("HP", (600., 1.))])
def test_internal_call_recovers_before_conversion_and_caches_only_success(monkeypatch, pair, inputs):
    calls = []
    failed = SimpleNamespace(ierr=226, herr="original PH failure")
    pp = package(SimpleNamespace(ABFLSHdll=lambda *args: calls.append(args) or failed))
    monkeypatch.setattr(pp, "_ensure_refprop_context", lambda comp: None)
    attempts = []
    monkeypatch.setattr(pp, "_recover_two_phase_ph", lambda *args: attempts.append("PT") or None)

    def recover(p, h, comp):
        assert (p, h, comp) == (1000., 600., [.5, .5])
        attempts.append("nearby")
        return state(p, h)

    monkeypatch.setattr(pp, "_recover_nearby_binary_ph", recover)
    with property_call_cache() as cache:
        first = pp.get_std_properties(pair, *inputs, pp.fld_composition)
        second = pp.get_std_properties(pair, *inputs, pp.fld_composition)
        assert len(cache.values) == 1
    assert first == second and first["P"] == 1. and first["h"] == 600.
    assert attempts == ["PT", "nearby"] and len(calls) == 1


def test_nearby_search_uses_asymmetric_checked_endpoints_without_recursing(monkeypatch):
    calls = []

    def native(pair, p, h, comp, flag):
        calls.append((pair, p, h, flag))
        if abs(h-(600.-1e-5)) < 1e-10 or abs(h-(600.+1e-4)) < 1e-10:
            return state(p, h)
        # Missing every property: failed native output must never be inspected.
        return SimpleNamespace(ierr=226)

    pp = package(SimpleNamespace(ABFLSHdll=native))
    monkeypatch.setattr(pp, "_flash", lambda *a: pytest.fail("Recovery must not recurse"))
    monkeypatch.setattr(pp, "_refine_binary_ph", lambda *a: None)
    monkeypatch.setattr(pp, "_valid_two_phase_ph", lambda s, p, h, c: abs(s.h-h) < 1e-10)
    result = pp._recover_nearby_binary_ph(1000., 600., [.5, .5])
    assert result.h == pytest.approx(600., abs=1e-10, rel=0)
    assert result.T == pytest.approx(306., abs=1e-10, rel=0)
    assert result.P == pytest.approx(1000.)
    assert len(calls) == 4
    assert all(pair == "PH" and flag == 0 for pair, _, _, flag in calls)


@pytest.mark.parametrize("native_success", [False, True])
def test_nearby_search_is_bounded_and_rejects_unvalidated_endpoints(monkeypatch, native_success):
    calls = []

    def native(pair, p, h, comp, flag):
        calls.append(1)
        return state(p, h) if native_success else SimpleNamespace(ierr=226)

    pp = package(SimpleNamespace(ABFLSHdll=native))
    monkeypatch.setattr(pp, "_refine_binary_ph", lambda *a: None)
    monkeypatch.setattr(pp, "_valid_two_phase_ph", lambda *a: False)
    assert pp._recover_nearby_binary_ph(1000., 600., [.5, .5]) is None
    assert len(calls) == 24


def test_all_valid_endpoints_allow_at_most_74_bracket_corrections(monkeypatch):
    corrections = []
    pp = package(SimpleNamespace(ABFLSHdll=lambda pair, p, h, comp, flag: state(p, h)))

    def refine(s, p, h, comp):
        if p == 1000. and h == 600.:
            corrections.append(1)
        return None

    monkeypatch.setattr(pp, "_refine_binary_ph", refine)
    monkeypatch.setattr(pp, "_valid_two_phase_ph", lambda s, p, h, c: p != 1000. or h != 600.)
    assert pp._recover_nearby_binary_ph(1000., 600., [.5, .5]) is None
    assert len(corrections) == 74


@pytest.mark.parametrize("pair,error,comp", [("PT", 226, [.5, .5]), ("PH", 248, [.5, .5]),
                                          ("PH", 226, [1.])])
def test_other_failures_do_not_enter_binary_nearby_fallback(monkeypatch, pair, error, comp):
    pp = package(SimpleNamespace(ABFLSHdll=lambda *a: SimpleNamespace(ierr=error, herr="original error")))
    monkeypatch.setattr(pp, "_recover_nearby_binary_ph", lambda *a: pytest.fail("Unexpected fallback"))
    with pytest.raises(RefpropError, match="original error"):
        pp._flash(pair, 1000., 600., comp)


@pytest.mark.parametrize("comp", [[1.], [.3, .3, .4], [0., 1.], [float("nan"), .5]])
def test_nearby_search_requires_an_actual_binary_mixture(comp):
    pp = package(SimpleNamespace(ABFLSHdll=lambda *a: pytest.fail("Unexpected native call")))
    assert pp._recover_nearby_binary_ph(1000., 600., comp) is None


@pytest.mark.parametrize("fluid,pressure,enthalpy", [
    ("Propane;Pentane|0.6;0.4|mass", 4.8271138480928695, 547.5455623779903),
    ("Propane;Pentane|0.6;0.4|mass", 4.843474231959467, 546.8070821833591),
    ("Propane;Isobutane|0.3;0.7|mass", 3.8968903691775676, 650.5677987511655),
    ("Propane;Isobutane|0.3;0.7|mass", 3.9210507549126223, 648.2953120090376),
])
def test_failed_native_ph_recovers_validated_target_through_internal_call(hpl, fluid, pressure, enthalpy):
    pp = hpl.PropertyPackage(fluid)
    comp = pp.fld_composition
    p, h = pp._convert_props_from_SI_with_C("PH", pressure, enthalpy)
    result = pp._flash("PH", p, h, comp)
    assert pp._valid_two_phase_ph(result, p, h, comp)
    assert any(key.startswith("PH nearby recovery succeeded") for key in pp.refprop_diagnostics)
    converted = pp.get_std_properties("HP", enthalpy, pressure, comp)
    assert converted["h"] == pytest.approx(enthalpy, abs=1e-7, rel=0)
    assert converted["P"] == pytest.approx(pressure, abs=1e-9, rel=0)
