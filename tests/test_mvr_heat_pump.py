"""Heat pump + MVR steam (HeatPumpLab/cycles/MVRHeatPump) against the "Advanced dryer heat pumps" web app.

The reference values are the page's JavaScript model (DRYCASC.booster and boosterOptimise, artifact version
1792578914-d233) run in Node over CoolProp 8.0.0 (WebAssembly), with the page's inputs: humid dryer air
(15 C, 60 % RH, 101.325 kPa), source water 40 -> 30 C, 5 K source approach, 20 K air approach, ammonia bottom cycle,
Pareto front compressor, motor 0.93, MVR fans on the Pareto front at PR 1.5 with 5 K exit superheat and a 0.93 motor. The
Python port runs on CoolProp 8 and reproduces them to about 1e-9.
"""
from __future__ import annotations

import pytest

from HeatPumpLab.cycles.MVRHeatPump import MVRHeatPump, MVRHPConfig
from HeatPumpLab.utils import coolprop_properties as cpp

pytest.importorskip("CoolProp.CoolProp")
if int((cpp.coolprop_version() or "0").split(".")[0]) < 8:
    pytest.skip("CoolProp < 8", allow_module_level=True)

KEYS = ("COP", "P_evap_bar", "P_cond_bar", "T_dis_bottom", "P_steam_bar", "P_final_bar", "steam_per_kg", "gap_bottom",
        "gap_steam", "sub_raise", "q_ihx", "w_exp_kJkg")

# (case, config overrides, reference report values from the JavaScript model)
CASES = [
    ("app default", {},
     (2.302367110634495, 10.026949730931909, 51.16415040525371, 184.31291335036076, 0.5786697165209151,
      22.246007871642917, 0.49438164766096365, 19.99999999990007, 56.86383770482678, 0, 0, 0)),
    ("air 200 C, 12 fans, condensing 100 C", dict(T_air_tar=200.0, n_mvr=12, T_cond=100.0),
     (1.6573606839756128, 10.026949730931909, 62.55124503032476, 208.45020139327158, 0.8460846631954946,
      109.77638659503832, 0.8017231249266, 19.99999999912677, 136.2544376470561, 0, 0, 0)),
    ("split fill", dict(split_fill=True),
     (2.994450708819638, 10.026949730931909, 51.16415040525371, 184.31291335036076, 0.5786697165209151,
      22.246007871642917, 0.43662867320933046, 19.99999999990007, 20.517752737600432, 0, 0, 0)),
    ("bottom expander 0.75", dict(expander_eta=0.75),
     (2.31723402436787, 10.026949730931909, 51.16415040525371, 184.31291335036076, 0.5786697165209151,
      22.246007871642917, 0.49438164766096365, 19.99999999990007, 56.86383770482678, 0, 0, 5.689306619970637)),
    ("bottom IHX, no subcooler", dict(ihx=True, sub_air=False),
     (2.042243360359039, 10.026949730931909, 51.16415040525371, 263.1491547258098, 0.5786697165209151,
      22.246007871642917, 0.47221101292302603, 20.000000000001407, 56.86383770482681, 0, 147.61986019004462, 0)),
    ("propane, subcooler pinch inside", dict(refrigerant="Propane", T_cond=85.0),
     (2.4268643375360477, 9.520745233424922, 34.36097694710335, 99.34016727155694, 0.47414474030233716,
      18.227716647208794, 0.13774053555553548, 19.992706037176106, 47.22347477776856, 1.1568124312283516, 0, 0)),
    ("constant-efficiency fans", dict(fan_model=1, fan_eta=0.75, PR_mvr=1.4, sh_inject=10.0, fan_motor=0.95),
     (2.449982165063315, 10.026949730931909, 51.16415040525371, 184.31291335036076, 0.5786697165209151,
      11.95592208552264, 0.45630736379526204, 19.99999999990007, 26.830305768900644, 0, 0, 0)),
    ("process splits", dict(split=(0, 0, 0, 0, 0.2, 0.3, 0.3, 0.5)),
     (2.617974641982939, 10.026949730931909, 51.16415040525371, 184.31291335036076, 0.5786697165209151,
      22.246007871642917, 0.45976854543781637, 19.99999999990007, 49.150247329353675, 0, 0, 0)),
    ("split fill with a common margin", dict(dT_GC_min=30.0, split_fill=True, T_cond=80.0),
     (2.778384645251874, 10.026949730931909, 41.41294727322482, 160.2662858092413, 0.38595395847525743,
      14.837366727868147, 0.5062883923203969, 29.99999999999295, 29.950065515159878, 0, 0, 0)),
]


@pytest.fixture(scope="module")
def shared():
    """One set of property objects for all cases (the memo caches make the suite fast)."""
    base = MVRHeatPump(MVRHPConfig(backend="coolprop"))
    return base


def model(shared, **kw):
    cfg = MVRHPConfig(backend="coolprop", **kw)
    fluid = shared.fb if cfg.refrigerant == shared.c.refrigerant else None
    return MVRHeatPump(cfg, fluid=fluid, water=shared.water, air=shared.air)


@pytest.mark.parametrize("name, kw, ref", CASES, ids=[c[0] for c in CASES])
def test_reproduces_the_web_app(shared, name, kw, ref):
    res = model(shared, **kw).solve()
    assert res.success, res.message
    for key, value in zip(KEYS, ref):
        assert res.report[key] == pytest.approx(value, rel=1e-6, abs=1e-6), key
    assert res.COP == pytest.approx(ref[0], rel=1e-7)


def test_default_case_details(shared):
    res = model(shared).solve()
    r = res.report
    assert res.warnings == ["steam raised at 0.58 bar (vacuum)"]
    # steam at 85 C, nine fans at PR 1.5, air after the bottom cycle at 90 - 20 = 70 C
    assert r["T_steam"] == 85.0 and len(r["stages"]) == 9 and r["T_air_mid"] == 70.0
    assert r["P_final_bar"] == pytest.approx(r["P_steam_bar"] * 1.5 ** 9)
    assert r["splits"] == [0.0] * 8
    assert res.states["suction"]["T"] == pytest.approx(30.0)          # source out 30 - 5 K, + 5 K superheat
    # energy balance per kg of bottom refrigerant: air = source + compressor + fans (- expander)
    assert res.q_air == pytest.approx(res.q_source + r["w_bottom"] + r["w_mvr"] - r["w_exp_kJkg"], abs=1e-6)
    assert res.q_air == pytest.approx(r["q_air_bottom"] + r["q_air_steam"])
    assert res.COP == pytest.approx(res.q_air / r["w_el_kJkg"])
    assert {p["name"] for p in res.paths} == {"Bottom cycle", "Steam"}


def test_split_fill_margin(shared):
    res = model(shared, dT_GC_min=30.0, split_fill=True, T_cond=80.0).solve()
    assert res.report["fill_margin"] == pytest.approx(3.4140625, abs=1e-9)
    assert res.report["gap_steam"] >= 30.0 - 0.05


def test_expander_warns_of_a_two_phase_outlet(shared):
    res = model(shared, expander_eta=0.75).solve()
    assert res.report["exp_out_quality"] == pytest.approx(0.04, abs=0.01)
    assert res.warnings[0].startswith("two-phase expansion: the expander outlet is 4 % vapour (25.0 C)")


def test_failures(shared):
    low = model(shared, T_cond=30.0).solve()
    assert not low.success and low.message == "condensing temperature too low to heat the air"
    cross = model(shared, T_air_tar=200.0, n_mvr=2).solve()
    assert not cross.success
    assert cross.message.startswith("the process steam cannot heat the air to 200 C (temperature cross)")


def test_optimiser_reproduces_the_web_app(shared):
    opt = model(shared).optimise()
    assert opt.success
    assert opt.COP == pytest.approx(2.643304452233928, rel=1e-7)
    assert (opt.T_cond, opt.n_mvr, opt.n_eval) == (102.4, 6, 440)
    assert opt.converged and opt.at_bound == {"T_cond": False, "n_mvr": False}
    assert opt.start == pytest.approx(2.302367110634495, rel=1e-7) and opt.start_ok
    assert opt.result.COP == opt.COP and opt.splits == [0.0] * 5


def test_refprop_backend_agrees():
    try:
        res = MVRHeatPump(MVRHPConfig(backend="refprop")).solve()
    except Exception as exc:                                  # REFPROP not installed
        pytest.skip(f"REFPROP unavailable: {exc}")
    assert res.success, res.message
    assert res.COP == pytest.approx(2.302367110634495, rel=1e-6)
