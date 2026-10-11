"""The three manufacturer R600a compressor maps (Recip-S, Recip-M, Recip-L): swept volumes, the model list,
the browser copy of the fits, and the cycle on each compressor (Python and JavaScript agree)."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from HeatPumpLab.cycles.ConventionalHeatPump import (
    DISPLACEMENT_50HZ_M3H, R600A_COEFFICIENTS, recip_models, recip_r600a,
)
from test_steam_hot_water import JS_TOY, WATER, steam_model

ROOT = Path(__file__).resolve().parents[1]
JS = ROOT / "local_scripts" / "sghp_live" / "sghp_model.js"
MODELS = ("Recip-S", "Recip-M", "Recip-L")


def node_json(script):
    node = shutil.which("node")
    if node is None or not JS.exists():
        pytest.skip("node or sghp_model.js not available")
    out = subprocess.run([node, "-e", f"const S = require({json.dumps(str(JS))});" + script],
                         capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


def test_every_model_has_a_swept_volume():
    assert recip_models() == list(MODELS)                                   # smallest first
    assert set(DISPLACEMENT_50HZ_M3H) == set(R600A_COEFFICIENTS) == set(MODELS)
    assert DISPLACEMENT_50HZ_M3H == {"Recip-S": 13.42, "Recip-M": 56.25, "Recip-L": 151.6}
    for m in MODELS:
        comp = recip_r600a(m)
        assert comp.name == m and comp.displacement_m3h == DISPLACEMENT_50HZ_M3H[m]
    assert recip_r600a("Recip-S", displacement_m3h=None).displacement_m3h is None   # size-less still possible
    with pytest.raises(KeyError, match="available"):
        recip_r600a("4NES-20")


def test_browser_copy_of_the_fits():
    js = node_json("console.log(JSON.stringify({c: S.COEFFS, d: S.DISPLACEMENT_50HZ, "
                   "e: Object.keys(S.COEFFS).map(m => [m, S.compressorMap(m).displacement_m3h, "
                   "[[-5, 40], [20, 90], [45, 115]].map(([a, b]) => [S.compressorMap(m).eta_is(a, b), S.compressorMap(m).eta_vol(a, b)])])}));")
    assert list(js["c"]) == list(MODELS) and js["d"] == DISPLACEMENT_50HZ_M3H
    for m in MODELS:
        assert js["c"][m]["vol"] == list(R600A_COEFFICIENTS[m]["volumetric"])
        assert js["c"][m]["is"] == list(R600A_COEFFICIENTS[m]["isentropic"])
    for m, disp, pts in js["e"]:
        py = recip_r600a(m)
        assert disp == DISPLACEMENT_50HZ_M3H[m]
        for (a, b), (eis, evol) in zip([(-5, 40), (20, 90), (45, 115)], pts):
            assert eis == pytest.approx(py.eta_is(a, b), rel=1e-12) and evol == pytest.approx(py.eta_vol(a, b), rel=1e-12)


def test_the_models_differ():
    eta = {m: recip_r600a(m).eta_is(20.0, 90.0) for m in MODELS}
    assert len({round(v, 4) for v in eta.values()}) == 3


def test_cycle_on_each_compressor_matches_the_browser():
    for model in MODELS:
        check_cycle(model)


def check_cycle(model):
    comp = recip_r600a(model, basis="temperature")
    py = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp).solve()
    assert py.success, py.message
    r = py.report
    # mass flow = swept volume x volumetric efficiency x suction density
    rho = py.states["suction"]["D"]
    assert r["m_dot_kgs"] == pytest.approx(comp.displacement_m3h / 3600 * r["eta_vol"] * rho, rel=1e-9)
    js = node_json(JS_TOY + f"const base = {json.dumps(WATER | {'dT_superheat': 5.0})};"
                   f"const r = S.solve(toy, base, S.compressorMap({json.dumps(model)}, {{map_basis: 'temperature'}}), 0.0);"
                   "console.log(JSON.stringify(r.success ? r.report : {msg: r.message}));")
    for k in ("COP", "Q_heat_kW", "m_dot_kgs", "eta_is_overall", "eta_vol", "W_electric_kW"):
        assert js[k] == pytest.approx(r[k], rel=1e-9), k


def test_heat_output_scales_with_size():
    q = {m: steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=recip_r600a(m, basis="temperature")).solve().report["Q_heat_kW"]
         for m in MODELS}
    assert q["Recip-S"] < q["Recip-M"] < q["Recip-L"]
