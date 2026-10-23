"""n-level cascade (Hoang et al.) against the dryer web app's JS model, plus behaviour checks."""
import pytest

from dryer_cascade_reference import close, compare_result, compare_scalars, coolprop_or_skip, reference_cases

coolprop_or_skip()

from HeatPumpLab.cycles.CascadeHeatPump import CascadeConfig, CascadeHeatPump  # noqa: E402
from HeatPumpLab.cycles.common.dryer import Fluid  # noqa: E402

CASES = reference_cases("cascade")
_FLUIDS: dict = {}


def model(**over) -> CascadeHeatPump:
    cfg = CascadeConfig(property_backend="coolprop", **over)
    fluids = []
    for name in cfg.fluids[:cfg.n]:
        if name not in _FLUIDS:
            _FLUIDS[name] = Fluid(name, "coolprop")
        fluids.append(_FLUIDS[name])
    return CascadeHeatPump(cfg, fluids=fluids)


def compare_levels(res, exp):
    assert len(res.report["levels"]) == len(exp["report"]["levels"])
    for k, (got, want) in enumerate(zip(res.report["levels"], exp["report"]["levels"])):
        compare_scalars(got, want, f"levels[{k}]")


@pytest.mark.parametrize("case", [c for c in CASES if not c["optimise"]], ids=lambda c: c["name"])
def test_matches_js(case):
    res = model(**case["config"]).solve()
    exp = case["expected"]
    compare_result(res, exp)
    if exp["success"]:
        compare_levels(res, exp)


def test_app_default_values():
    """The app's default: isobutane condensing at 100 C, cyclopentane on top, 170 C air."""
    res = model().solve()
    assert res.success, res.message
    l1, l2 = res.report["levels"]
    assert res.COP == pytest.approx(2.97819, abs=1e-4)
    assert l1["Te"] == pytest.approx(25.0) and l1["Tc"] == pytest.approx(100.0, abs=1e-9)
    assert l2["Te"] == pytest.approx(95.0) and l2["Tc"] == pytest.approx(185.008, abs=1e-3)
    assert l1["Pc_bar"] == pytest.approx(19.8652, abs=1e-3) and l2["Pc_bar"] == pytest.approx(21.1014, abs=1e-3)
    # level 1 sends ~9 % of its condensing heat to the air, the air leaving it at 100 - 20 = 80 C
    assert l1["f_air"] == pytest.approx(0.0909, abs=1e-3) and l2["air_in"] == pytest.approx(80.0, abs=0.01)
    assert res.q_air == pytest.approx(res.w_shaft + res.q_source, rel=1e-9)
    assert res.q_air == pytest.approx(l1["q_air"] + l2["q_air"], rel=1e-9)
    assert l2["gap_air"] == pytest.approx(20.0, abs=0.05)
    assert set(res.paths) == {"Level 1", "Level 2"}


def test_ihx_raises_the_suction():
    off = model(n=3, T_cond=[80.0, 130.0], T_air_tar=200.0).solve()
    on = model(n=3, T_cond=[80.0, 130.0], T_air_tar=200.0, ihx=[False, True, True, False]).solve()
    assert off.success and on.success
    assert on.report["levels"][1]["superheat"] > off.report["levels"][1]["superheat"]
    assert any(e["name"] == "Level 2 IHX" and e["gap"] >= 5.0 - 0.05 for e in on.exchangers)


def test_expander_credit():
    valve = model(n=4, fluids=["Isobutane", "Pentane", "Cyclopentane", "Cyclopentane"], T_cond=[70.0, 110.0, 150.0])
    exp = model(n=4, fluids=["Isobutane", "Pentane", "Cyclopentane", "Cyclopentane"], T_cond=[70.0, 110.0, 150.0],
                expander_eta=0.7)
    rv, re_ = valve.solve(), exp.solve()
    assert rv.success and re_.success
    assert re_.report["w_exp_kJkg"] > 0 and re_.COP > rv.COP
    assert re_.report["exp_out_quality"] is not None and len(re_.report["exp_out_quality"]) == 4
    assert any("two-phase expansion" in w for w in re_.warnings)


@pytest.mark.parametrize("over, msg", [
    (dict(T_air_tar=250.0), "no top condensing temperature up to 238 C"),
    (dict(n=5, fluids=["Isobutane"] * 5, T_cond=[60.0, 80.0, 100.0, 120.0]), "a cascade has 2 to 4 levels"),
    (dict(T_cond=[20.0]), "no top condensing temperature"),
])
def test_failures(over, msg):
    res = model(**over).solve()
    assert not res.success and msg in res.message and res.COP is None


def test_missing_fluid():
    m = CascadeHeatPump(CascadeConfig(property_backend="coolprop", n=3), fluids=[Fluid("Isobutane", "coolprop")])
    res = m.solve()
    assert not res.success and res.message == "a refrigerant is needed for every level"


@pytest.mark.parametrize("case", [c for c in CASES if c["optimise"]], ids=lambda c: c["name"])
def test_optimiser_matches_js(case):
    res = model(**case["config"]).optimise(max_eval=case["optimise"]["max_eval"])
    exp = case["expected"]
    assert res.success, res.message
    close(res.COP, exp["COP"], "COP")
    opt = res.report["optimiser"]
    assert opt["T_cond"] == pytest.approx(exp["optimiser"]["T_cond"], abs=1e-9)
    close(opt["T_cond_top"], exp["optimiser"]["T_cond_top"], "T_cond_top")
    assert opt["nEval"] == exp["optimiser"]["nEval"]
    assert res.COP >= opt["start"] - 1e-9


def test_refprop_backend_agrees(hpl):
    res = CascadeHeatPump(CascadeConfig(property_backend="refprop", air_RH=0.0)).solve()
    ref = CascadeHeatPump(CascadeConfig(property_backend="coolprop", air_RH=0.0)).solve()
    assert res.success and ref.success
    assert res.COP == pytest.approx(ref.COP, rel=1e-2)
