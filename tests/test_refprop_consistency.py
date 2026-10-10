"""Phase-envelope lifecycle, flash failure handling, and near-zero precision."""

from types import SimpleNamespace

import pytest

from HeatPumpLab.utils.properties_package import PropertyPackage, RefpropError
from HeatPumpLab.utils.property_cache import property_call_cache


def state(**overrides):
    return SimpleNamespace(**dict(P=1000., T=300., D=2., h=600., s=4.,
                                  Cp=2., q=-1., ierr=0, herr="") | overrides)


class NativeStub:
    def __init__(self):
        self.setups = []
        self.envelopes = []
        self.calls = []
        self.spline_error = 0
        self.result = state()

    def SETFLUIDSdll(self, name):
        self.setups.append(name)
        return 0

    def SATSPLNdll(self, comp):
        self.envelopes.append(tuple(comp))
        return SimpleNamespace(ierr=self.spline_error, herr="spline diagnostic")

    def ABFLSHdll(self, pair, a, b, comp, flag):
        self.calls.append((pair, a, b, tuple(comp), flag))
        return self.result

    def REFPROPdll(self, *args):
        return SimpleNamespace(ierr=248, herr="failed flash", Output=[float("nan")])


def package(native, composition=(.3, .7)):
    pp = object.__new__(PropertyPackage)
    pp.RP = native
    pp.fld_name = "A.FLD*B.FLD"
    pp.fld_composition = list(composition)
    pp.molar_mass = 1.
    pp.polish_ph = False
    return pp


def test_envelope_is_shared_and_rebuilt_after_composition_or_fluid_change():
    native = NativeStub()
    first = package(native)
    same = package(native)
    other_composition = package(native, (.6, .4))
    pure = package(native, (1.,))
    pure.fld_name = "C"
    for pp in (first, same, other_composition, first, pure, first):
        pp.get_std_properties("PT", 1., 26.85, pp.fld_composition)

    assert native.setups == [first.fld_name, "C", first.fld_name]
    assert native.envelopes == [(.3, .7), (.6, .4), (.3, .7), (.3, .7)]
    assert all(call[-1] == 0 for call in native.calls)  # blind molar flashes


def test_failed_envelope_setup_is_not_marked_ready_or_cached():
    native = NativeStub()
    native.spline_error = 355
    pp = package(native)
    with property_call_cache() as cache:
        for _ in range(2):
            with pytest.raises(RefpropError, match="SATSPLNdll.*355"):
                pp.get_std_properties("PT", 1., 26.85, pp.fld_composition)
        assert not cache.values
    assert len(native.envelopes) == 2
    assert not native.calls


def test_failed_composition_switch_invalidates_the_previous_envelope():
    native = NativeStub()
    first = package(native)
    second = package(native, (.6, .4))
    first.get_std_properties("PT", 1., 26.85, first.fld_composition)
    native.spline_error = 355
    with pytest.raises(RefpropError):
        second.get_std_properties("PT", 1., 26.85, second.fld_composition)
    native.spline_error = 0
    first.get_std_properties("PT", 1., 26.85, first.fld_composition)
    assert native.envelopes == [(.3, .7), (.6, .4), (.3, .7)]


@pytest.mark.parametrize("change", [dict(ierr=248, herr="iteration failed"),
                                   dict(ierr=-235, herr="phase unknown"),
                                   dict(h=float("nan")), dict(D=0), dict(P=float("inf"))])
def test_invalid_flash_never_enters_the_property_cache(change):
    native = NativeStub()
    native.result = state(**change)
    pp = package(native)
    with property_call_cache() as cache:
        for _ in range(2):
            with pytest.raises(RefpropError):
                pp.get_std_properties("PT", 1., 26.85, pp.fld_composition)
        assert not cache.values
    assert len(native.calls) == 2


def test_nonfatal_flash_warning_is_preserved():
    native = NativeStub()
    native.result = state(ierr=-1, herr="temperature outside model bounds")
    pp = package(native)
    pp.get_std_properties("PT", 1., 26.85, pp.fld_composition)
    assert pp.refprop_diagnostics["ABFLSHdll(PT):-1"] == {
        "count": 1, "message": "temperature outside model bounds",
    }


def test_high_level_flash_failure_is_rejected_before_caching():
    pp = package(NativeStub())
    with property_call_cache() as cache:
        with pytest.raises(RefpropError, match="248"):
            pp.get_specified_properties("H", "PT", 1., 26.85)
        assert not cache.values


@pytest.mark.parametrize("pair,inputs", [("PH", (1., 600.)), ("HP", (600., 1.))])
def test_single_phase_ph_roundoff_is_polished_with_blind_pt(pair, inputs):
    native = NativeStub()

    def flash(pair, a, b, comp, flag):
        native.calls.append((pair, flag))
        if pair in ("PH", "HP"):
            return state(T=300.000001)
        return state(T=b, h=2*b)

    native.ABFLSHdll = flash
    pp = package(native)
    pp.polish_ph = True
    result = pp.get_std_properties(pair, *inputs, pp.fld_composition)
    assert result["T"] == pytest.approx(26.85, abs=1e-10)
    assert result["h"] == pytest.approx(600., abs=1e-10)
    assert native.calls == [(pair, 0), ("PT", 0), ("PT", 0)]


def test_two_phase_ph_is_not_polished_as_single_phase():
    native = NativeStub()
    native.result = state(q=.5, Cp=-9999990.)
    pp = package(native)
    pp.polish_ph = True
    assert pp.get_std_properties("PH", 1., 600., pp.fld_composition)["q"] == .5
    assert len(native.calls) == 1


def test_failed_polish_retains_valid_original_and_records_diagnostic():
    native = NativeStub()
    native.ABFLSHdll = lambda pair, *args: (
        state() if pair == "PH" else state(ierr=248, herr="PT failed")
    )
    pp = package(native)
    pp.polish_ph = True
    result = pp.get_std_properties("PH", 1., 600., pp.fld_composition)
    assert result["h"] == 600.
    assert pp.refprop_diagnostics["ABFLSHdll(PT):248"]["count"] == 1
    assert pp.refprop_diagnostics["PH polish declined: PT flash failed"]["count"] == 1


def test_precision_mode_has_a_separate_cache_key():
    native = NativeStub()
    pp = package(native)
    with property_call_cache():
        pp.get_std_properties("PT", 1., 26.85, pp.fld_composition)
        pp.polish_ph = True
        pp.get_std_properties("PT", 1., 26.85, pp.fld_composition)
    assert len(native.calls) == 2


def test_failed_two_phase_ph_does_not_force_another_phase():
    native = NativeStub()

    def flash(pair, pressure, value, comp, flag):
        native.calls.append((pair, flag))
        if pair == "PH":
            return state(ierr=226, herr="two-phase iteration failed")
        assert pair == "PQ"
        return state(h=500+100*value, q=value, Cp=-9999990.)

    native.ABFLSHdll = flash
    pp = package(native)
    with pytest.raises(RefpropError, match="226"):
        pp.get_std_properties("PH", 1., 540., pp.fld_composition)
    assert pp.refprop_diagnostics["ABFLSHdll(PH):226"]["count"] == 1
    # Saturation endpoints may be checked, but an invalid bracket must not
    # lead to a guessed quality or a forced single-phase result.
    assert native.calls == [("PH", 0), ("PQ", 0), ("PQ", 0)] + [("PH", 0)]*24


@pytest.mark.parametrize("fluid,pressure,enthalpy", [
    ("Isopentane;Neopentane|0.5;0.5|mass", 3208.337199446526, 36439.25748946073),
    ("Propane;Pentane|0.6;0.4|mass", 4850.4647826692335, 28448.223314311726),
])
def test_checked_two_phase_recovery_reproduces_requested_enthalpy(hpl, fluid, pressure, enthalpy):
    pp = PropertyPackage(fluid, "ref")
    result = pp._recover_two_phase_ph(pressure, enthalpy, pp.fld_composition)
    assert result is not None
    assert abs(result.h-enthalpy) < 1e-6
    assert 0 < result.q < 1 and result.Dl > result.Dv


def test_two_phase_recovery_rejects_inconsistent_phase_pressure(hpl, monkeypatch):
    pp = PropertyPackage("Isopentane;Neopentane|0.5;0.5|mass", "ref")
    original = pp.RP.THERMdll
    def inconsistent(*args):
        result = original(*args)
        return result._replace(P=result.P*1.01)
    monkeypatch.setattr(pp.RP, "THERMdll", inconsistent)
    result = pp._recover_two_phase_ph(3208.337199446526, 36439.25748946073, pp.fld_composition)
    assert result is None


def test_unrecoverable_ph_error_is_not_cached():
    native = NativeStub()
    native.result = state(ierr=226, herr="failed PH and PQ")
    pp = package(native)
    with property_call_cache() as cache:
        with pytest.raises(RefpropError, match="226"):
            pp.get_std_properties("PH", 1., 540., pp.fld_composition)
        assert not cache.values


def test_liquid_water_ph_avoids_the_spurious_solid_root(hpl):
    pp = hpl.PropertyPackage("water")
    expected = pp.get_std_properties("PT", 2., 20., pp.fld_composition)
    recreated = pp.get_std_properties("PH", 2., expected["h"], pp.fld_composition)
    assert recreated["T"] == pytest.approx(20., abs=1e-7)
    assert recreated["h"] == pytest.approx(expected["h"], abs=1e-7)
    with pytest.raises(RefpropError):
        pp.get_std_properties("PH", 2., -10000., pp.fld_composition)
    with pytest.raises(RefpropError):
        pp.get_std_properties("PH", 2., -20., pp.fld_composition)


def test_near_critical_ph_never_uses_failed_output(hpl):
    pp = hpl.PropertyPackage("Propane;Isobutane|0.3;0.7|mass")
    pressure = 3.946390399004728
    enthalpy = 34257.983530203266 / pp.molar_mass
    try:
        recovered = pp.get_std_properties("PH", pressure, enthalpy, pp.fld_composition)
    except RefpropError as error:
        # REFPROP 10 currently fails this two-phase flash. A newer library
        # may solve it; either way, the old erroneous enthalpy must not pass.
        assert "226" in str(error)
        return
    assert recovered["h"] == pytest.approx(enthalpy, abs=1e-7)
    assert recovered["T"] > 120.
    assert recovered["P"] == pytest.approx(pressure, abs=1e-9)


def test_real_mixture_roundtrip_survives_library_rebinding(hpl):
    pp = hpl.PropertyPackage("Propane;Isobutane|0.3;0.7|mass")
    pressure, enthalpy = 4.4062458384733825, 865.3571439882853
    for instance in ("default", "source", "sink", "default"):
        pp._set_RP_lib(instance)
        inlet = pp.get_std_properties("PH", pressure, enthalpy, pp.fld_composition)
        equal = pp.get_std_properties("PT", pressure, inlet["T"], pp.fld_composition)
        hotter = pp.get_std_properties("PT", pressure, inlet["T"]+1e-8, pp.fld_composition)
        assert equal["h"] == pytest.approx(enthalpy, abs=1e-9)
        assert hotter["h"] - equal["h"] == pytest.approx(inlet["Cp"]*1e-8, rel=1e-4)


@pytest.mark.parametrize("quality", [0., .5, 1.])
def test_real_mixture_saturation_keeps_two_phase_properties(hpl, quality):
    pp = hpl.PropertyPackage("Propane;Pentane|0.3;0.7|mass")
    saturated = pp.get_std_properties("TQ", 10., quality, pp.fld_composition)
    recreated = pp.get_std_properties("PH", saturated["P"], saturated["h"], pp.fld_composition)
    assert recreated["T"] == pytest.approx(10., abs=1e-4)
    assert recreated["h"] == pytest.approx(saturated["h"], rel=1e-6)


@pytest.mark.parametrize("lift", [0., 1e-8, 1e-6, .00060080442336426, .02, .05])
def test_near_zero_rhx2_gap_is_feasible(solve_recuperative_case, assert_successful_cycle, lift):
    hp = solve_recuperative_case(
        "Propane;Isobutane|0.3;0.7|mass",
        dT_rhx1_cold=36.16489790906914, dT_rhx2_cold=lift,
        comp2_power_frac=.19148358339296268,
        record_failed_solve=False, report_failed_solve=False,
    )
    assert_successful_cycle(hp)
    assert hp.gc_minimum_approach >= 19.95
    if lift > 0:
        assert hp.rhx2.duty >= -1e-6


def test_mixture_without_critical_line_falls_back_to_true_critical_point(monkeypatch):
    class CriticalStub:
        def REFPROPdll(self, fld, hin, hout, *args):
            if hout == "M":
                return SimpleNamespace(ierr=0, herr="", Output=[60.0])
            return SimpleNamespace(ierr=313, herr="[CRTHMX error 313] no coefficients",
                                   Output=[-9999990.0] * 4)

        def INFOdll(self, icomp):
            return SimpleNamespace(Tc=(370.0, 511.7)[icomp - 1], Pc=(4250.0, 4570.0)[icomp - 1])

        def CRTPNTdll(self, z, t_guess, p_guess):
            assert 370.0 < t_guess < 511.7 and 4250.0 < p_guess < 4570.0
            return SimpleNamespace(Tc=450.0, Pc=5000.0, Dc=3.0, ierr=0, herr="")

        def WMOLdll(self, z):
            return 60.0

    pp = package(CriticalStub(), composition=(.4, .6))
    pp.n_components = 2
    monkeypatch.setattr(pp, "_ensure_refprop_context", lambda comp: None)

    M, Pc, Tc, Dc = pp._get_critical_properties()

    assert M == pytest.approx(60.0)
    assert Pc == pytest.approx(5.0)
    assert Tc == pytest.approx(450.0 - 273.15)
    assert Dc == pytest.approx(180.0)


def test_pure_fluid_critical_failure_is_not_masked(monkeypatch):
    class FailingStub:
        def REFPROPdll(self, *args):
            return SimpleNamespace(ierr=313, herr="failed", Output=[float("nan")] * 4)

    pp = package(FailingStub(), composition=(1.0,))
    pp.n_components = 1
    monkeypatch.setattr(pp, "_ensure_refprop_context", lambda comp: None)

    with pytest.raises(RefpropError):
        pp._get_critical_properties()
