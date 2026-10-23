"""Cases on the SGHP live page (local_scripts/sghp_live/sghp_cases.js): the link code, input
cleaning, and the page wiring. The shared store and the browser behaviour are checked by hand."""
import importlib.util
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LIVE = ROOT / "local_scripts" / "sghp_live"


def run_cases_js(script):
    node = shutil.which("node")
    if node is None or not (LIVE / "sghp_cases.js").exists():
        pytest.skip("node or sghp_cases.js not available")
    live = (LIVE / "sghp_live.js").read_text()
    fluids = live[live.index("const FLUIDS"):live.index("class CoolPropFluid")]
    src = (f"const SGHP = require({json.dumps(str(LIVE / 'sghp_model.js'))});\n" + fluids
           + (LIVE / "sghp_cases.js").read_text() + "\n" + script)
    out = subprocess.run([node, "-e", src], capture_output=True, text=True, timeout=60)
    assert out.returncode == 0, out.stderr[-800:]
    return json.loads(out.stdout)


INPUTS = dict(comps=[["Propane", 20], ["Isopentane", 80]], basis="mass", compressor="Recip-L", T_source=50, T_steam=115,
              T_water_in=80, T_water_out=95, compressor_Hz=45, dTmin_source=5, dTmin_steam=5, dTmin_water=5,
              dT_superheat=15, dTmin_ihx=5, T_ambient=15, UA_W_K=15, UA_cond_W_K=20, UA_subcool_W_K=10, UA_evap_W_K=5, preheat_subcool_pct=50,
              dPK_cond=0.05, dPK_subcool=0.5, dPK_evap=0.02)


def test_link_round_trip():
    o = run_cases_js(f"const inp = {json.dumps(INPUTS)};"
                     "const link = caseLink('Propane/isopentane · 115 °C — test', inp);"
                     "const code = link.split('#case=')[1];"
                     "const back = JSON.parse(unb64url(code));"
                     "console.log(JSON.stringify({link, back, clean: cleanInputs(back.i), label: fluidLabel(inp), urlsafe: /^[A-Za-z0-9_-]+$/.test(code)}));")
    assert o["link"].startswith("https://waikato-ahuora-smart-energy-systems.github.io/HeatPumpLab/steam/#case=") and o["urlsafe"]
    assert o["back"]["n"] == "Propane/isopentane · 115 °C — test"            # non-ASCII names survive
    assert o["clean"] == INPUTS and o["back"]["v"] == 1
    assert o["label"] == "Propane/Isopentane 20/80 (mass)"
    old = {k: v for k, v in INPUTS.items() if not k.startswith(("UA_cond", "UA_subcool", "UA_evap", "preheat_subcool_pct", "dPK_"))}   # saved before the exchanger losses
    o2 = run_cases_js(f"console.log(JSON.stringify(cleanInputs({json.dumps(old)})));")
    assert (o2["UA_cond_W_K"], o2["UA_subcool_W_K"], o2["UA_evap_W_K"], o2["preheat_subcool_pct"]) == (0, 0, 0, 0)
    assert (o2["dPK_cond"], o2["dPK_subcool"], o2["dPK_evap"]) == (0, 0, 0)


def test_untrusted_inputs_are_cleaned():
    o = run_cases_js(
        "const bad = Object.assign({}, %s, {comps: [['Propane', 150], ['NotAFluid', 10], ['Isopentane', -5]], basis: 'x',"
        " compressor: 'nope', T_steam: 500, compressor_Hz: 5, extra: '<script>'});"
        "const c = cleanInputs(bad); let e1 = null, e2 = null;"
        "try { cleanInputs({comps: [['NotAFluid', 100]]}); } catch (e) { e1 = e.message; }"
        "try { cleanInputs(Object.assign({}, bad, {T_source: 'abc'})); } catch (e) { e2 = e.message; }"
        "console.log(JSON.stringify({c, e1, e2}));" % json.dumps(INPUTS))
    c = o["c"]
    assert c["comps"] == [["Propane", 100], ["Isopentane", 0]] and c["basis"] == "mass" and c["compressor"] == "Recip-M"
    assert c["T_steam"] == 140 and c["compressor_Hz"] == 30 and "extra" not in c
    assert "no known refrigerant" in o["e1"] and "T_source is missing" in o["e2"]


def test_summary_and_rows():
    o = run_cases_js("const r = {COP: 3.7, Q_heat_kW: 50, Q_cond_kW: 32, Q_subcool_kW: 18, W_electric_kW: 13.5, P_cond_bar: 24,"
                     " P_evap_bar: 5, PR: 4.8, T_discharge: 130, eta_is_overall: 0.69, eta_vol: 0.79, m_dot_kgs: 0.14,"
                     " m_water_kgs: 0.28, VHC: 4.05, rho_suction_kgm3: 18.9, V_suction_m3h: 26.7, V_discharge_m3h: 5.1, in_envelope: true};"
                     "console.log(JSON.stringify({s: summarize(r), keys: RESULT_ROWS.map(x => x[0]), max: MAX_COMPARE,"
                     " colors: CASE_COLORS, sweeps: Object.keys(SWEEPS)}));")
    s = o["s"]
    assert s["ok"] and s["m_dot_kgh"] == pytest.approx(504) and s["m_water_kgh"] == pytest.approx(1008) and s["in_envelope"]
    assert set(o["keys"]) <= set(s) and o["max"] == 4 and len(o["colors"]) == 4
    assert o["sweeps"] == ["T_steam", "T_source", "compressor_Hz"]


@pytest.mark.skipif(not (LIVE / "build_sghp_live.py").exists(), reason="local_scripts/sghp_live not present (gitignored)")
def test_page_wiring(tmp_path):
    spec = importlib.util.spec_from_file_location("build_sghp_live", LIVE / "build_sghp_live.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    html = mod.build(tmp_path / "out", None, mock=True).read_text()
    for needle in ('id="cases-panel"', 'id="compare" hidden', 'id="case-save"', 'id="case-copy"', 'id="case-open"',
                   'id="sweep-x"', 'id="sweep-y"', 'id="cmp-inputs"', 'id="cmp-ts"', 'id="cmp-ph"',
                   "boot();\ncasesBoot();", 'claude.use("db")', 'collection("cases")', "function setCompressor(name)",
                   'if (typeof casesOnRender === "function") casesOnRender();',
                   'id="cmp-side"', "function sideBySide(cs)", "dialSVG(r.P_discharge_bar ?? r.P_cond_bar, pTop,", 'id="dPK_cond"', "dP_K_cond: inp.dPK_cond",
                   ".side-kv dd.win,.win-chip{", "function startDrag(ev, col)", "function moveCompare(id, to)",
                   "const caseColor = id =>", ".grip{cursor:grab;touch-action:none",
                   "const base = okRs[0];", 'r => r.rho_suction_kgm3', '["Discharge flow, m³/h", r => r.V_discharge_m3h, 1, -1]', '["Suction flow, m³/h", r => r.V_suction_m3h, 1, -1]', 'id="case-all"', 'id="case-none"', "function selectAll()", 'id="case-del-sel"', "async function deleteSelected()"):
        assert needle in html, needle
    assert "set baseline" not in html and "CS.baseline" not in html          # the left-most case is the baseline
    assert 'id="cmp-results"' not in html                                   # the side-by-side view replaces the results table
    for needle in ():
        assert needle in html, needle
    assert '"lose"' not in html and "--lose" not in html and ".lose" not in html   # only the best value is highlighted
    for needle in ():
        assert needle in html, needle
    assert "<select" not in html                                          # buttons only (embedded viewers)
    assert 'id="cmp-order"' in html and 'id="cmp-order-list"' in html and "function applyOrder(cs)" in html
    assert "localStorage.setItem(LOCAL_KEY" in html                       # works without the shared store
    # an automatic case name follows the inputs (refrigerant, temperatures...) until someone types one
    assert "if (CS.nameAuto && document.activeElement !== $(\"case-name\")) $(\"case-name\").value = \"\";" in html
    assert "CS.nameAuto = !name || name === suggestName(inp);" in html


def test_order_by():
    """The comparison can be kept sorted: by name, or by a result best first (lower first where lower is better,
    higher first otherwise); cases without a result go last in their current order."""
    o = run_cases_js(
        "const cs = [{id: 'a', name: 'Case 10'}, {id: 'b', name: 'Case 2'}, {id: 'c', name: 'Base'}, {id: 'd', name: 'Failed'}];"
        "const val = {COP: {a: 3.5, b: 3.9, c: 3.7, d: null}, P_discharge_bar: {a: 24, b: 30, c: 20, d: null}, m_dot_kgh: {a: 500, b: 400, c: 600, d: null}};"
        "const by = k => sortedCompare(cs, k, c => val[k][c.id]);"
        "console.log(JSON.stringify({none: sortedCompare(cs, '', null), name: sortedCompare(cs, 'name', null), cop: by('COP'),"
        " p: by('P_discharge_bar'), m: by('m_dot_kgh'), opts: orderOptions().map(o => o[0])}));")
    assert o["none"] == ["a", "b", "c", "d"]
    assert o["name"] == ["c", "b", "a", "d"]                 # numeric-aware: Case 2 before Case 10
    assert o["cop"] == ["b", "c", "a", "d"]                  # highest COP first
    assert o["p"] == ["c", "a", "b", "d"]                    # lowest max pressure first
    assert o["m"] == ["c", "a", "b", "d"]                    # neither better: highest first
    assert o["opts"][:3] == ["", "name", "COP"] and "P_discharge_bar" in o["opts"]
