"""sCO2 reverse Brayton heat pump (HeatPumpLab/cycles/BraytonHeatPump) against the "Advanced dryer heat pumps" web app.

The reference values are the page's JavaScript model (DRYCASC.brayton and braytonOptimiseSteps, artifact version
1792578914-d233) run in Node over CoolProp 8.0.0 (WebAssembly), with the page's inputs: humid dryer air (60 % RH at the
air inlet, 101.325 kPa), eta_c = eta_t = 0.85, eta_m 0.99, IHX and preheater effectiveness 0.85, pressure drops
50 kPa (high side, per exchanger) and 10 kPa (low side), 5 K source approach, 20 K air approach, Pmax 300 bar.

- App defaults: air 15 -> 170 C, source water 40 -> 30 C, motor 0.93, P_L 80 bar, P_H 165 bar, x 0.37. That design
  has no solution (compressor 2 stays below the target); the page's optimiser finds one at a subcritical suction.
- Paper basis (Euro sCO2 2027, Dalal, Klöppel, Rose and Walmsley): air 35 -> 170 C, sensible source 80 -> 60 C,
  COP on shaft power with eta_m (motor 1.0), P_L 80 bar. The JavaScript model gives architecture 1 a COP of 3.83 and
  architecture 2 3.99 at a 20 K approach (P_H and x optimised), and 3.52 / 3.64 at a 30 K approach; at the paper's
  EBSILON pressures (188.85 / 147 bar) and 30 K, 3.47 / 3.63. The paper's own EBSILON values (3.64 / 3.72) are not
  reproduced exactly by this model (see the co-author notes of the Euro sCO2 study).

Architecture 1 reproduces the JavaScript to about 1e-9. Architecture 2 to about 1e-4 in the COP: its bisection on
the air temperature between the HTHXs (0.02 K tolerance) reacts to the ~1e-10 differences between CoolProp's native
and WebAssembly builds. The full optimiser runs (a minute or more each) run with HEATPUMPLAB_SLOW_TESTS=1.
"""
from __future__ import annotations

import math
import os

import pytest

from HeatPumpLab.cycles.BraytonHeatPump import BraytonHeatPump, BraytonHPConfig
from HeatPumpLab.utils import coolprop_properties as cpp

pytest.importorskip("CoolProp.CoolProp")
if int((cpp.coolprop_version() or "0").split(".")[0]) < 8:
    pytest.skip("CoolProp < 8", allow_module_level=True)

slow = pytest.mark.skipif(not os.environ.get("HEATPUMPLAB_SLOW_TESTS"), reason="set HEATPUMPLAB_SLOW_TESTS=1")

PAPER = dict(T_air_inlet=35.0, T_air_tar=170.0, T_source_in=80.0, T_source_out=60.0, motor_eff=1.0)
KEYS = ("COP", "P_H_bar", "T_suction", "T_c2_out", "T_hthx_out", "T_turb_out", "T_lthx_out", "T_air_mid", "q_air_1",
        "q_air_2", "q_lthx")

# (case, config overrides (MPa), reference report values from the JavaScript model, then VHC)
CASES = [
    ("app default, arch 1 optimum", dict(P_L=6.625, P_H=19.1975, x=0.36),
     (3.2427171345031827, 191.975, 44.09509489273, 233.88235308712444, 39.50485484033322, 21.937056612998845,
      34.63541485363942, None, 235.9140023125396, 0, 168.25972777294396), 28.14916384191099),
    ("app default, arch 2 optimum", dict(arch=2, P_L=6.5, P_H=17.1268, x=0.3225),
     (3.2484341272884785, 171.268, 58.81480752220773, 227.2624174574363, 64.4789572698686, 23.502972205173137,
      35.00000011913107, 35.214331662372274, 200.0791178737934, 29.833531823597536, 164.960548913776),
     25.149869782136772),
    ("paper, arch 1, 20 K, optimum at 80 bar", dict(PAPER, P_H=16.453125, x=0.36875),
     (3.8297076262881546, 164.53125, 76.15623442966375, 204.3413620735028, 72.02821118891859, 36.057328748977795,
      74.99999999999994, None, 151.71535278458558, 0, 112.81020279042787), 20.636241540424198),
    ("paper, arch 1, 30 K, 188.85 bar", dict(PAPER, dT_GC_min=30.0, P_H=18.885, x=0.37),
     (3.4657498866966354, 188.85, 77.1008514463494, 234.1116299505208, 71.94815782473324, 35.20261983549341,
      74.99999999999994, None, 180.6200533433044, 0, 129.37772164739295), 23.669433802729888),
    ("paper, arch 2, 30 K, 147 bar", dict(PAPER, dT_GC_min=30.0, arch=2, P_H=14.7, x=0.41),
     (3.6281381427957027, 147.0, 101.2071302455821, 214.29264426510025, 107.86597277788428, 37.755324747658335,
      74.99999999999994, 76.70840328146501, 90.66868758154034, 40.28525193282872, 95.51054255400959),
     16.146010519843106),
    ("paper, arch 2, 20 K, 136 bar", dict(PAPER, arch=2, P_H=13.6, x=0.41),
     (3.990959949738941, 136.0, 99.86127128721614, 196.19676890167358, 107.5081924983956, 38.600260578708,
      74.99999999999994, 85.1487982199335, 74.86443671188938, 43.96530350734585, 89.65488102081332), 15.00050581354671),
    ("paper, arch 1, lowest discharge pressure", dict(PAPER, P_H_min=True),
     (3.625034654464999, 146.970703125, 84.5055873462693, 190.0007021065549, 83.9424136621036, 40.19227344354948,
      74.99999999999994, None, 110.74294421617213, 0, 80.81924083798953), 14.858047651430475),
    ("paper, arch 1, LTHX duty capped", dict(PAPER, P_H=16.453125, x=0.36875, q_lthx_max=100.0),
     (3.363436514365356, 164.53125, 70.95140300478306, 198.5802020146786, 74.88393537166564, 34.337728272721336,
      46.580916124246585, None, 141.37791868270617, 0, 100.0), 19.784538037524282),
]


@pytest.fixture(scope="module")
def shared():
    """Property objects shared by the cases (their memo and isobar tables make the suite fast), one air model per
    air inlet temperature."""
    m = BraytonHeatPump(BraytonHPConfig(backend="coolprop"))
    return {"fluid": m.pr, "air": {15.0: m.air}}


def model(shared, **kw):
    cfg = BraytonHPConfig(backend="coolprop", **kw)
    air = shared["air"].get(cfg.T_air_inlet)
    m = BraytonHeatPump(cfg, fluid=shared["fluid"], air=air)
    shared["air"][cfg.T_air_inlet] = m.air
    return m


def check(res, ref, vhc, arch):
    # arch 2's bisection on the air temperature between the two HTHXs amplifies last-digit property differences
    # (native vs WebAssembly CoolProp, input values returned exactly): its temperatures move by up to ~0.1 K and its
    # COP by up to ~0.1 %, as much as REFPROP differs from the JavaScript on the same case
    rel, dT = (1e-6, 1e-4) if arch == 1 else (1e-3, 0.1)
    assert res.success, res.message
    for key, value in zip(KEYS, ref):
        if value is None:
            assert res.report[key] is None
        elif key == "COP" or key.startswith("q_"):
            assert res.report[key] == pytest.approx(value, rel=rel * 10, abs=1e-6), key
        else:
            assert res.report[key] == pytest.approx(value, abs=dT), key
    assert res.COP == pytest.approx(ref[0], rel=rel)
    assert res.VHC == pytest.approx(vhc, rel=rel * 10)


@pytest.mark.parametrize("name, kw, ref, vhc", CASES, ids=[c[0] for c in CASES])
def test_reproduces_the_web_app(shared, name, kw, ref, vhc):
    check(model(shared, **kw).solve(), ref, vhc, kw.get("arch", 1))


def test_cycle_closes(shared):
    res = model(shared, **PAPER, P_H=16.453125, x=0.36875).solve()
    r = res.report
    assert r["converged"] and abs(r["balance_kJkg"]) < 0.005 * res.q_air + 0.2
    assert res.q_air == pytest.approx(res.q_source + r["w_c1"] + r["w_c2"] - r["w_turb"], abs=0.01)
    assert r["COP_shaft"] == pytest.approx(res.COP)                      # motor 1.0
    assert r["gap_hthx"] >= 20.0 - 0.1 and r["gap_lthx"] >= 5.0 - 0.1
    st = res.states
    assert st["suction"]["P"] == pytest.approx(8.0) and st["c2_in"]["P"] == pytest.approx(7.99)
    assert st["turb_in"]["P"] == pytest.approx(16.453125 - 0.1) and st["turb_out"]["P"] == pytest.approx(8.02)
    assert st["lthx_out"]["T"] == pytest.approx(75.0, abs=1e-6)          # source in 80 C less the 5 K approach
    assert len(res.paths) == 11 and len(res.profiles) == 4


def test_app_defaults_fail_as_in_the_web_app(shared):
    res = model(shared).solve()
    assert not res.success and res.P_H == pytest.approx(16.5)
    assert res.message == ("compressor 2 discharges at 69.3 °C, below the 170 °C target plus the 20 K approach: "
                           "raise the discharge pressure")
    low = model(shared, P_H_min=True).solve()
    assert low.message == ("compressor 2 reaches only 125 °C at 300 bar, short of the 170 °C target plus the 20 K "
                           "approach: raise the max pressure or the preheater/IHX effectiveness")


def test_failures(shared):
    low = model(shared, **PAPER, P_H=12.0, x=0.37).solve()
    assert low.message == ("compressor 2 discharges at 140.6 °C, below the 170 °C target plus the 20 K approach: "
                           "raise the discharge pressure to at least 147.0 bar")
    assert model(shared, P_H=8.2).solve().message == "the discharge pressure must be above the suction pressure"
    src = model(shared, T_source_in=25.0, T_source_out=20.0, P_L=12.0, P_H=25.0).solve()
    assert not src.success
    assert src.message.startswith("LTHX: the CO₂ leaves the turbine at 16.7 °C, above the source water outlet "
                                  "(20.0 °C) less the 5 K source approach")
    sub = model(shared, P_L=6.625, P_H=19.1975, x=0.36).solve()
    assert sub.warnings[0].startswith("suction 66 bar is below the critical pressure (73.8 bar)")


# ---------------------------------------------------------------- optimiser


def test_optimiser_at_fixed_suction_pressure(shared):
    for kw, ref in ((dict(PAPER), (3.8297076262881546, 16.453125, 0.36875, 38)),
                    (dict(PAPER, dT_GC_min=30.0), (3.5177910087946715, 16.46875, 0.37, 32))):
        opt = model(shared, **kw).optimise(opt_P_L=False)
        assert opt.success and opt.P_L == 8.0
        assert opt.COP == pytest.approx(ref[0], rel=1e-7)
        assert (opt.P_H, opt.x, opt.n_eval) == pytest.approx(ref[1:])
        assert opt.result.COP == pytest.approx(opt.COP)


def test_optimiser_with_a_solve_budget(shared):
    # P_L free from the paper basis: the search climbs the suction pressure until the budget stops it
    opt = model(shared, **PAPER).optimise(budget=80)
    assert opt.COP == pytest.approx(4.164504944978695, rel=1e-7)
    assert (opt.P_L, opt.P_H, opt.x, opt.n_eval) == pytest.approx((10.0, 20.789184027777782, 0.40475, 108))
    assert opt.budget_hit
    # app defaults: no design at 80 bar, so the suction pressure steps down (to 70 bar) before the budget runs out
    opt = model(shared).optimise(budget=60)
    assert opt.COP == pytest.approx(3.0977428852176967, rel=1e-7)
    assert (opt.P_L, opt.P_H, opt.x, opt.n_eval) == pytest.approx((7.0, 24.41375, 0.37, 69))
    assert opt.start == -math.inf and not opt.start_ok


@slow
def test_optimiser_at_the_app_defaults(shared):
    opt = model(shared).optimise()
    assert opt.COP == pytest.approx(3.2427169571859022, rel=1e-6)
    assert (opt.P_L, opt.P_H, opt.x, opt.n_eval) == pytest.approx((6.625, 19.19750123626374, 0.36, 499), rel=1e-6)
    assert not opt.budget_hit and opt.at_bound == {"P_L": False, "P_H": False, "x": False}


@slow
def test_optimiser_at_the_app_defaults_arch2(shared):
    opt = model(shared, arch=2).optimise()
    assert opt.COP == pytest.approx(3.248393364875174, rel=1e-3)
    assert opt.P_L == pytest.approx(6.5, abs=0.125) and opt.P_H == pytest.approx(17.126777380952383, abs=0.2)
    assert opt.x == pytest.approx(0.3225, abs=0.01)


@slow
def test_arch2_paper_optimum_and_lowest_pressure(shared):
    for kw, ref in ((dict(PAPER), (3.9947618704086554, 13.625, 0.4075)),
                    (dict(PAPER, dT_GC_min=30.0), (3.6422070381067986, 14.125, 0.415))):
        opt = model(shared, arch=2, **kw).optimise(opt_P_L=False)
        assert opt.COP == pytest.approx(ref[0], rel=5e-4)
        # the optimum is flat in P_H: the compass search can stop a step away from the JavaScript's (142.8 bar
        # against 141.25 bar at 30 K) for a COP within 1e-4
        assert opt.P_H == pytest.approx(ref[1], abs=0.25) and opt.x == pytest.approx(ref[2], abs=0.02)
    low = model(shared, **PAPER, arch=2, P_H_min=True).solve()
    assert low.COP == pytest.approx(3.899403323317905, rel=1e-3)
    assert low.report["P_H_bar"] == pytest.approx(131.1396484375, abs=0.05)


def test_refprop_backend_agrees():
    try:
        res = BraytonHeatPump(BraytonHPConfig(backend="refprop", **PAPER, P_H=16.453125, x=0.36875)).solve()
    except Exception as exc:                                  # REFPROP not installed
        pytest.skip(f"REFPROP unavailable: {exc}")
    assert res.success, res.message
    assert res.COP == pytest.approx(3.8297076262881546, rel=1e-6)
