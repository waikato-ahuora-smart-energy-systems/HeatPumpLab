"""TT CO2 cascade (Kong TT-4) against the dryer web app's JS model, plus behaviour checks."""
import math

import pytest

from dryer_cascade_reference import close, compare_result, compare_scalars, coolprop_or_skip, reference_cases

coolprop_or_skip()

from HeatPumpLab.cycles.TTCascadeHeatPump import TTCascadeConfig, TTCascadeHeatPump  # noqa: E402

CASES = reference_cases("tt")
_MODELS = {}


def model(**over) -> TTCascadeHeatPump:
    """One model per air specification (the property caches are shared, as in the app)."""
    cfg = TTCascadeConfig(property_backend="coolprop", **over)
    key = (cfg.air_RH, cfg.T_air_inlet, cfg.P_air_kPa)
    if key not in _MODELS:
        _MODELS[key] = TTCascadeHeatPump(TTCascadeConfig(property_backend="coolprop", air_RH=cfg.air_RH,
                                                         T_air_inlet=cfg.T_air_inlet, P_air_kPa=cfg.P_air_kPa))
    base = _MODELS[key]
    return TTCascadeHeatPump(cfg, bottom=base.bottom, top=base.top, air=base.air)


@pytest.mark.parametrize("case", [c for c in CASES if not c["optimise"]], ids=lambda c: c["name"])
def test_matches_js(case):
    res = model(**case["config"]).solve()
    exp = case["expected"]
    compare_result(res, exp)
    if exp["success"]:
        for loop in ("top", "bot"):     # CO2 keeps CoolProp's default reference state, so h and s compare too
            for name, st in exp["states"][loop].items():
                compare_scalars(res.states[loop][name], st, f"states.{loop}.{name}")


def test_app_default_values():
    """The app's default design: 150 bar both cycles, r 0.43, 170 C air, no secondary air."""
    res = model().solve()
    assert res.success, res.message
    r = res.report
    assert res.COP == pytest.approx(2.57999, abs=1e-4)
    assert r["P_evap_top_bar"] == pytest.approx(57.2905, abs=1e-3)       # CO2 dew at 20 C
    assert r["P_evap_bot_bar"] == pytest.approx(45.0218, abs=1e-3)       # and at 10 C
    assert r["T_dis_top"] == pytest.approx(235.28, abs=0.01) and r["T_dis_bot"] == pytest.approx(162.92, abs=0.01)
    assert r["gap_MT"] == pytest.approx(20.0, abs=0.01) and r["gap_HT"] > 20.0
    assert r["bottom_supercritical"] is True and r["share_sec"] == 0.0
    # energy balance per kg top flow: heat to air = compressor work + source heat (to the 0.01 K closure of the
    # gas-cooler outlet iteration, as in the app)
    assert res.q_air == pytest.approx(res.w_shaft + res.q_source, rel=1e-4)
    assert res.q_air == pytest.approx(r["q_HT"] + r["q_MT"], rel=1e-9)
    assert set(res.paths) == {"Top cycle", "Bottom cycle"}
    assert all(len(p["pts"]) > 1 for loop in res.paths.values() for p in loop)


def test_expander_credit_and_two_phase_warning():
    valve = model(P_b=7.0, r=2.5, gc_off=5.0).solve()
    exp = model(P_b=7.0, r=2.5, gc_off=5.0, expander_eta=0.75).solve()
    assert valve.success and exp.success
    assert valve.report["w_exp_kJkg"] == 0.0 and exp.report["w_exp_kJkg"] > 0.0
    assert exp.COP > valve.COP
    assert exp.report["w_el_kJkg"] == pytest.approx(exp.w_shaft / 0.93 - exp.report["w_exp_kJkg"] * 0.93)
    assert any("two-phase expansion" in w for w in exp.warnings)
    assert exp.report["exp_out_quality"][1] is not None and 0 < exp.report["exp_out_quality"][1] < 1
    assert exp.paths["Top cycle"][0]["name"] == "Expander"


def test_secondary_air_share():
    res = model(air_ratio=0.5, T_mt=102.0).solve()
    assert res.success
    r = res.report
    assert 0 < r["share_sec"] < 0.5
    assert r["q_primary"] + r["q_secondary"] == pytest.approx(res.q_air, rel=1e-6)


@pytest.mark.parametrize("over, msg", [
    (dict(r=0.0), "the bottom/top mass flow must be positive"),
    (dict(P_t=5.0), "discharge pressure below the evaporating pressure"),
    (dict(air_ratio=0.5, T_mt=180.0), "the secondary air target must lie between"),
    (dict(P_t=8.0, r=0.1), "cannot heat the primary air"),
    (dict(r=5.0), "the bottom cycle brings more heat"),
])
def test_failures(over, msg):
    res = model(**over).solve()
    assert not res.success and msg in res.message and res.COP is None


@pytest.mark.parametrize("case", [c for c in CASES if c["optimise"]], ids=lambda c: c["name"])
def test_optimiser_matches_js(case):
    res = model(**case["config"]).optimise(max_eval=case["optimise"]["max_eval"])
    exp = case["expected"]
    assert res.success, res.message
    close(res.COP, exp["COP"], "COP")
    opt = res.report["optimiser"]
    for k in ("P_b", "P_t", "r", "gc_off", "start"):
        close(opt[k], exp["optimiser"][k], k)
    assert opt["nEval"] == exp["optimiser"]["nEval"]
    assert res.COP >= opt["start"] - 1e-9


def test_refprop_backend_agrees(hpl):
    """REFPROP (when installed) gives nearly the same design: the same equations of state for CO2."""
    res = TTCascadeHeatPump(TTCascadeConfig(property_backend="refprop", air_RH=0.0)).solve()
    ref = TTCascadeHeatPump(TTCascadeConfig(property_backend="coolprop", air_RH=0.0)).solve()
    assert res.success and ref.success
    assert res.COP == pytest.approx(ref.COP, rel=5e-3)
    assert math.isfinite(res.VHC)
