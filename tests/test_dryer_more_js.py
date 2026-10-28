"""The dryer app's JavaScript ports of the economiser, multi-pressure, Stirling and ammonia-water hybrid cycles
(apps/src/dryer/dryer_more.js) against the Python models, run in Node over CoolProp 8 (WebAssembly) and on the Python
CoolProp backend: the same equations, property library and iterations, so they agree to round-off. The hybrid runs on
the page's IAPWS G4-01 port (apps/src/dryer/ammonia_water.js) and the Python model on teqp. The air is the app's default humid air (the
moisture of 15 C, 60 % RH air), which the Python models build once from their configured inlet."""
from __future__ import annotations

import importlib.util
import json

import pytest

from test_dryer_live import LIVE, _run_node

pytest.importorskip("CoolProp")

from HeatPumpLab.cycles.EconomiserHeatPump import EconomiserConfig, EconomiserHeatPump  # noqa: E402
from HeatPumpLab.cycles.MultiPressureHeatPump import MultiPressureConfig, MultiPressureHeatPump  # noqa: E402
from HeatPumpLab.cycles.StirlingHeatPump import StirlingConfig, StirlingHeatPump  # noqa: E402

REL = 2e-6

# (name, model, overrides); the JS keys are the Python config fields
CASES = [
    ("flash", "E", dict(arrangement="flash")),
    ("economiser", "E", dict(arrangement="economiser")),
    ("parallel", "E", dict(arrangement="parallel")),
    ("single", "E", dict(arrangement="single")),
    ("flash hot inlet", "E", dict(arrangement="flash", T_air_inlet=110.0)),
    ("economiser hot inlet", "E", dict(arrangement="economiser", T_air_inlet=110.0, T_int=95.0)),
    ("parallel hot inlet", "E", dict(arrangement="parallel", T_air_inlet=110.0)),
    ("flash butane", "E", dict(arrangement="flash", T_air_inlet=90.0, T_air_tar=140.0, fluid="n-Butane")),
    ("too hot", "E", dict(arrangement="flash", T_air_tar=230.0, P_max_bar=30.0)),
    ("dual", "D", dict(n_pressures=2, T_cond=[110.0])),
    ("dual 80", "D", dict(n_pressures=2, T_cond=[80.0])),
    ("three", "D", dict()),
    ("four", "D", dict(n_pressures=4)),
    ("five", "D", dict(n_pressures=5)),
    ("three butane", "D", dict(T_cond=[80.0, 110.0], T_air_tar=140.0, fluid="n-Butane")),
    ("three falling", "D", dict(T_cond=[130.0, 110.0])),
    ("three too low", "D", dict(T_cond=[48.0, 110.0])),
    ("stirling 1", "R", dict()),
    ("stirling 3", "R", dict(n_units=3)),
    ("stirling bounds", "R", dict(n_units=2, boundaries=[70.0], eta_carnot=0.4)),
    ("stirling bad", "R", dict(n_units=2, boundaries=[200.0])),
    ("hybrid", "Y", dict()),
    ("hybrid dry", "Y", dict(wet_compression=False)),
    ("hybrid hot", "Y", dict(T_air_tar=150.0)),
    ("hybrid capped", "Y", dict(T_air_tar=150.0, P_max_bar=28.0)),
    ("hybrid lean", "Y", dict(w_rich=0.35, q_des=0.25)),
]
OPT = [("flash", "E", dict(arrangement="flash", T_air_inlet=110.0)), ("dual", "D", dict(n_pressures=2, T_cond=[110.0])),
       ("three", "D", dict(max_eval=40)), ("stirling", "R", dict(n_units=3)),
       ("hybrid", "Y", dict(max_eval=12))]
TEQP = importlib.util.find_spec("teqp") is not None

FLUID = {"Cyclopentane": "Cyclopentane", "n-Butane": "n-Butane"}
SCALARS = ("COP", "VHC", "Pmax_bar", "T_dis_max", "q_air", "w_shaft", "q_source")


def _python(model: str, over: dict, optimise: bool = False):
    over = dict(over)
    if model == "Y":
        if not TEQP:
            pytest.skip("teqp not installed")
        from HeatPumpLab.cycles.HybridHeatPump import HybridHeatPump
        return HybridHeatPump().optimise(**over) if optimise else HybridHeatPump().solve(**over)
    if model == "R":
        hp = StirlingHeatPump(StirlingConfig(**over))
    elif model == "E":
        hp = EconomiserHeatPump(EconomiserConfig(property_backend="coolprop", fluid=over.pop("fluid", "Cyclopentane")))
    else:
        hp = MultiPressureHeatPump(MultiPressureConfig(property_backend="coolprop", fluid=over.pop("fluid", "Cyclopentane")))
    if model == "R":
        return hp.optimise() if optimise else hp.solve()
    if not optimise:
        over.pop("max_eval", None)
    return hp.optimise(**over) if optimise else hp.solve(**over)


@pytest.fixture(scope="module")
def js(tmp_path_factory):
    """Every case and optimiser in one Node run."""
    cases = [dict(name=n, model=m, over=o) for n, m, o in CASES]
    opts = [dict(name=n, model=m, over=o) for n, m, o in OPT]
    return _run_node(tmp_path_factory.mktemp("more"), f"""
const M = req({json.dumps(str(LIVE / "dryer_more.js"))});
const W = CP.HAPropsSI('W', 'T', 288.15, 'P', 101325, 'R', 0.6);
const humid = W0 => {{ const a = DRY.airModel(T => CP.HAPropsSI('H', 'T', T + 273.15, 'P', 101325, 'W', W0) / 1000); a.W = W0; return a; }};
const fls = {{}}, fl = n => fls[n] || (fls[n] = new Fluid(CP, [[n, 100]], 'mass'));
DRY.setCompModel(5);
const run = (m, o) => {{
  const {{ fluid = 'Cyclopentane', ...c }} = o, a = humid(W);
  return m === 'E' ? M.economiser(fl(fluid), a, c) : m === 'D' ? M.multiPressure(fl(fluid), a, c) : m === 'Y' ? M.hybrid(a, c) : M.stirling(a, c);
}};
const gen = (m, o) => {{
  const {{ fluid = 'Cyclopentane', ...c }} = o, a = humid(W);
  const {{ max_eval, ...cc }} = c, opt = max_eval ? {{ maxEval: max_eval }} : {{}};
  const g = m === 'E' ? M.economiserOptimiseSteps(fl(fluid), a, cc, opt) : m === 'D' ? M.multiPressureOptimiseSteps(fl(fluid), a, cc, opt)
    : m === 'Y' ? M.hybridOptimiseSteps(a, cc, opt) : M.stirlingOptimiseSteps(a, cc, opt);
  let s; do {{ s = g.next(); }} while (!s.done); return s.value;
}};
const pick = r => ({{ success: r.success, message: r.message, warnings: r.warnings, COP: r.COP, VHC: r.VHC, Pmax_bar: r.Pmax_bar, T_dis_max: r.T_dis_max,
  q_air: r.q_air, w_shaft: r.w_shaft, q_source: r.q_source, report: r.report }});
const out = {{ cases: {{}}, opts: {{}} }};
for (const c of {json.dumps(cases)}) out.cases[c.name] = pick(run(c.model, c.over));
for (const c of {json.dumps(opts)}) out.opts[c.name] = gen(c.model, c.over);
console.log(JSON.stringify(out));
""")


def _close(a, b, what):
    if b is None or isinstance(b, (bool, str)):
        assert a == b, what
    else:
        assert a == pytest.approx(b, rel=REL, abs=1e-6), what


@pytest.mark.parametrize("name,model,over", CASES, ids=[c[0] for c in CASES])
def test_matches_python(js, name, model, over):
    got, want = js["cases"][name], _python(model, over)
    assert got["success"] == want.success, (got["message"], want.message)
    if not want.success:
        assert got["message"].split(" (")[0].replace(".0 ", " ") == want.message.split(" (")[0].replace(".0 ", " ")
        return
    for k in SCALARS:
        if getattr(want, k) is not None:
            _close(got[k], getattr(want, k), k)
    assert got["warnings"] == want.warnings
    for k, v in want.report.items():
        if k not in got["report"]:
            continue
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            _close(got["report"][k], v, f"report.{k}")
        elif isinstance(v, list) and v and all(isinstance(x, (int, float)) for x in v):
            assert len(got["report"][k]) == len(v), k
            for i, x in enumerate(v):
                _close(got["report"][k][i], x, f"report.{k}[{i}]")


@pytest.mark.parametrize("name,model,over", OPT, ids=[c[0] for c in OPT])
def test_optimisers_match_python(js, name, model, over):
    got, want = js["opts"][name], _python(model, over, optimise=True)
    assert got["success"] and want.success
    o = want.report["optimiser"]
    _close(got["COP"], want.COP, "COP")
    assert got["nEval"] == o["nEval"]
    for key in {"E": ["T_int"], "D": ["T_cond"], "R": ["boundaries"], "Y": ["w_rich", "q_des"]}[model]:
        assert got[key] == pytest.approx(o[key], abs=1e-9)
    assert got["COP"] >= got["start"] - 1e-12
