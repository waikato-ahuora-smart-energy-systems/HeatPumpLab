"""Build checks for the SGHP live page (local_scripts/sghp_live/build_sghp_live.py).

These catch a page that loads but can never solve: the CoolProp loader missing or not handed to the
page, or coolprop.wasm not copied next to index.html. They need no CoolProp and no browser.
"""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "local_scripts" / "sghp_live" / "build_sghp_live.py"

FAKE_ES_LOADER = """var Module = (() => { return async function(opts) { return {}; }; })();
var s = "</script>";
export default Module;
"""
FAKE_CLASSIC_LOADER = "var Module = typeof Module != 'undefined' ? Module : {};\n"


@pytest.fixture(scope="module")
def build_mod():
    if not BUILD.exists():
        pytest.skip("local_scripts/sghp_live not present")
    spec = importlib.util.spec_from_file_location("build_sghp_live", BUILD)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fake_coolprop(tmp_path, loader):
    d = tmp_path / "coolprop_js"
    d.mkdir()
    (d / "coolprop.js").write_text(loader)
    (d / "coolprop.wasm").write_bytes(b"\0asm\1\0\0\0")
    return d / "coolprop.js"


def test_missing_loader_fails_the_build(build_mod, tmp_path):
    # a page without the loader shows "Loading CoolProp…" forever, so the build must refuse
    with pytest.raises(SystemExit, match="coolprop.js not found"):
        build_mod.build(tmp_path / "out", tmp_path / "nowhere" / "coolprop.js")
    assert not (tmp_path / "out" / "index.html").exists()


def test_es_module_loader_is_inlined_and_handed_to_the_page(build_mod, tmp_path):
    js = fake_coolprop(tmp_path, FAKE_ES_LOADER)
    page = build_mod.build(tmp_path / "out", js)
    html = page.read_text()
    assert '<script type="module">' in html
    assert "window.__cpFactoryReady(Module);" in html
    assert "export default" not in html                       # not valid in an inline module that hands off
    assert 'var s = "<\\/script>";' in html                   # loader text cannot close its own tag
    assert "window.__cpFactoryPromise" in html and "boot();" in html
    assert (page.parent / "coolprop.wasm").read_bytes() == (js.parent / "coolprop.wasm").read_bytes()


def test_classic_loader_gets_locate_file(build_mod, tmp_path):
    js = fake_coolprop(tmp_path, FAKE_CLASSIC_LOADER)
    html = build_mod.build(tmp_path / "out", js).read_text()
    assert "locateFile: p => (p.endsWith('.wasm') ? 'coolprop.wasm' : p)" in html
    assert "__cpFactoryReady(Module)" not in html


def test_mock_build_needs_no_coolprop(build_mod, tmp_path):
    page = build_mod.build(tmp_path / "out", None, mock=True)
    assert page.exists() and not (page.parent / "coolprop.wasm").exists()


def test_page_controls_and_headline(build_mod, tmp_path):
    html = build_mod.build(tmp_path / "out", None, mock=True).read_text()
    # compressor speed bar runs the full width of its panel, readout underneath
    assert '<div class="pickrow wide temp"><div class="slhead"><label class="lab" for="hz">Compressor speed</label>' in html
    assert ".pickrow.wide{grid-template-columns:minmax(0,1fr)" in html
    # five evenly spaced headline results; dials 1.25x the original 150 px; COP number sized to the stars
    assert ".headline{display:grid;grid-template-columns:repeat(5,minmax(0,1fr))" in html     # COP + four gauges
    assert 'dialSVG(pmax, dialTop(pmax, rm.P_discharge_bar ?? rm.P_cond_bar), "--s6", "bar maximum pressure")' in html and "Max pressure, bar" in html
    # refrigerant pressure-loss coefficients sit in Approaches and losses, before the ambient (last)
    i = [html.index(f'id="{k}"') for k in ("UA_evap", "dPK_cond", "dPK_subcool", "dPK_evap", "T_ambient")]
    assert i == sorted(i) and "dP_K_cond: Math.max(0, v(\"dPK_cond\") || 0)" in html
    assert ".headline > .dial{width:100%;max-width:188px}" in html
    assert "function fitCop()" in html
    assert ".hl-cop{padding-right:28px;border-right" not in html                # no vertical divider
    assert 'id="hz" min="30" max="70"' in html                                   # 30 Hz minimum speed
    # source, steam and hot-water temperatures on sliders; no condensing or subcooling inputs
    for sid, lo, hi, val in (("T_source", 20, 80, 55), ("T_steam", 90, 140, 120), ("T_water_in", 20, 120, 80), ("T_water_out", 30, 135, 95)):
        assert f'<input type="range" id="{sid}" min="{lo}" max="{hi}" step="1" value="{val}">' in html, sid
    assert 'id="T_cond"' not in html and 'id="dT_subcool"' not in html and 'id="T_evap"' not in html
    assert 'id="dTmin_source"' in html and 'id="dTmin_steam"' in html and 'id="dTmin_water"' in html
    assert 'id="ihx"' not in html                          # IHX after the subcooler, on whenever superheat > 0
    # no preset fluids; refrigerant and exchanger settings stacked beside the operating point
    assert 'id="presets"' not in html and "PRESETS" not in html and "chipbtn" not in html
    # no native drop-downs (some embedded viewers do not open them): fluid picker and mass/mole toggle are buttons
    assert "<select" not in html and "createElement(\"select\")" not in html
    assert 'class="tbtn" data-basis="mass" aria-pressed="true"' in html and 'id="frac-basis" value="mass"' in html
    assert 'list.className = "fluid-list"' in html
    # left column: compressor and operating point; right column: refrigerant and approaches and losses
    left = html.index('class="ctl-col"'); right = html.index('class="ctl-col"', left + 1)
    assert left < html.index('id="h-comp"') < html.index('id="h-point"') < right < html.index('id="h-fluid"') < html.index('id="h-set"')
    assert html.index('id="h-comp"') < html.index('id="compressor-pick"') < html.index('id="h-set"')   # the compressor has its own box
    assert '<div class="temps">' in html and ".temps{display:grid;grid-template-columns:repeat(2,minmax(0,1fr))" in html
    # compressor map on the density basis: R600a vapour-density table and per-fluid maps in the page
    assert "R600A_SAT_VAPOUR" in html and "function fluidMaps(fl, sh)" in html and "R600A_SUCTION_15K" in html and "map_basis: \"density\"" in html
    assert "W/K (fixed)" not in html and "Efficiency maps are 50 Hz data" not in html      # chips removed
    assert html.index('id="chips"') < html.index('class="controls"')                       # warnings above the inputs
    # fixed UA of 15 W/K; no UA calibration anywhere in the page
    assert 'id="preheat_sub" min="0" max="100" step="5" value="100"' in html and 'id="btn-design"' not in html
    pair = html.index('<div class="temps pair">')                             # preheat and suction superheat side by side
    assert pair < html.index('id="preheat_sub"') < html.index('id="dT_superheat"') < html.index('id="hz"')
    assert '<input type="range" id="dT_superheat" min="0" max="40" step="1" value="15">' in html
    assert 'id="dT_subcool_cond"' not in html                                  # a share of the preheat, not K
    assert 'for="dT_superheat" title=' in html and ">Superheat, K<" not in html and 'type="number" id="dT_superheat"' not in html
    assert '>Compressor heat loss UA, W/K<input type="number" id="UA_manual" value="15"' in html and "Shell UA" not in html
    for k in ("UA_cond", "UA_subcool", "UA_evap"):                      # exchanger heat losses, 0 W/K by default
        assert f'<input type="number" id="{k}" value="0" step="1" min="0">' in html, k
    assert '["Total heat loss", fmt(totalLoss(r), 2), "kW"' in html and '"Shell heat loss"' not in html
    assert "calibrat" not in html.lower() and 'id="ua-mode"' not in html
    assert ">Steam-Generating Heat Pump | Ahuora Centre</h1>" in html
    # headline: COP stars and three dials, no split bar
    assert "function copStars(cop) { return Math.max(0, Math.min(5, cop)); }" in html   # 1 star = 1 COP from 0, partial, no rounding
    assert "/5</small>" not in html and "% of Carnot · " not in html          # stars only, no score or subtext
    assert "hot water from subcooling" not in html
    # one shared dial scale: largest value rounded up to the next 10 kW
    assert "const dialTop = (...v) => Math.max(10, Math.ceil(Math.max(...v.filter(isFinite)) / 10) * 10);" in html
    assert "mx = { steam: top, hw: top, el: top }" in html and "DIAL_MAX" not in html
    assert "function maxHzCase()" in html and "LIVE.resultMaxHz =" in html            # scale set at the top speed
    assert "cooling COP" not in html
    # heat exchange profiles: IHX drawn after the evaporator, shaded areas between the streams
    assert "const shift = pts => pts.map(p => [p[0] + xE, p[1]]);" in html and 'fill-opacity": 0.16' in html
    assert "`Dmolar|${ph}`" in html and "DmolarT_INPUTS" not in html   # phase-imposed (rho, T) calls: no stability analysis
    # compressor picker: the three manufacturer maps as buttons; header follows the choice
    assert 'id="compressor-pick"' in html and "function compressorPick()" in html
    assert '<div class="eyebrow" id="eyebrow">Single-stage heat pump · solved in your browser with CoolProp</div>' in html
    assert "internal heat exchanger · solved" not in html
    assert "(COMPRESSOR.displacement_m3h * cur / 50).toFixed(1)} m³/h swept" in html and "56.25 * cur" not in html
    assert "const COMPRESSORS = Object.keys(SGHP.COEFFS);" in html and "LIVE.fluidMeta.mapSH = undefined;" in html
    for m in ("Recip-S", "Recip-M", "Recip-L"):
        assert f'"{m}": {{' in html, m
    assert '["Isentropic efficiency", ' in html and '["Volumetric efficiency", ' in html and '"η isentropic"' not in html
    assert '["Mass flow", fmt(r.m_dot_kgs * 3600, 1), "kg/h",' in html and "hot water ${fmt(r.m_water_kgs * 3600, 0)} kg/h" in html   # flows in kg/h
    for name in ("function copStars", "function dialSVG", "function headline(c)", 'id="headline"'):
        assert name in html, name


def test_coolprop_stack_helpers_are_exposed(build_mod, tmp_path):
    # a C++ exception leaving CoolProp's WebAssembly does not unwind its stack pointer; the page
    # restores it (cpCall) using two helpers the build exports from the loader
    loader = ('async function Module(m={}){var Module=m;var stackRestore=val=>0;var stackSave=()=>0;'
              'Module["decrementExceptionRefcount"]=decrementExceptionRefcount;return Module}\nexport default Module;\n')
    js = fake_coolprop(tmp_path, loader)
    html = build_mod.build(tmp_path / "out", js).read_text()
    assert 'Module["sghpStackSave"]=stackSave;Module["sghpStackRestore"]=stackRestore;' in html
    assert "function cpCall(fn)" in html and "CP.sghpStackRestore(sp)" in html
    assert 'guarded(cpCall(() => M.factory("HEOS", fluids)))' in html


def test_real_coolprop_survives_repeated_errors(build_mod, tmp_path):
    """With the stack restored after each CoolProp exception, hundreds of failing flashes leave
    CoolProp working (without it, about 60 abort the module)."""
    import shutil, subprocess
    js, node = build_mod.CP_DIR / "coolprop.js", shutil.which("node")
    if not js.exists() or node is None:
        pytest.skip("coolprop.js or node not available")
    anchor = 'Module["decrementExceptionRefcount"]=decrementExceptionRefcount;'
    (tmp_path / "coolprop.mjs").write_text(js.read_text().replace(
        anchor, anchor + 'Module["sghpStackSave"]=stackSave;Module["sghpStackRestore"]=stackRestore;', 1))
    shutil.copy(js.parent / "coolprop.wasm", tmp_path / "coolprop.wasm")
    (tmp_path / "t.mjs").write_text(
        "import F from './coolprop.mjs';\n"
        "const M = await F({ locateFile: p => new URL(p, import.meta.url).pathname });\n"
        "const a = M.factory('HEOS', 'Propane');\n"
        "for (let i = 0; i < 300; i++) { const sp = M.sghpStackSave();\n"
        "  try { a.update(M.input_pairs.QT_INPUTS, 0, 398.15); } catch (e) { M.sghpStackRestore(sp); } }\n"
        "console.log(M.PropsSI('T', 'P', 101325, 'Q', 0, 'Water'));\n")
    out = subprocess.run([node, str(tmp_path / "t.mjs")], capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr[-500:]
    assert float(out.stdout.strip().splitlines()[-1]) == pytest.approx(373.124, abs=0.01)


def test_real_coolprop_build(build_mod, tmp_path):
    js = build_mod.CP_DIR / "coolprop.js"
    if not js.exists():
        pytest.skip("local_scripts/coolprop_js/coolprop.js not present")
    page = build_mod.build(tmp_path / "out", js)
    html = page.read_text()
    assert len(html) > len(js.read_text())                    # loader inlined, not just referenced
    assert 'Module["sghpStackSave"]=stackSave' in html          # the real loader can be patched
    assert "__cpFactoryReady(Module)" in html or "onRuntimeInitialized" in html
    assert (page.parent / "coolprop.wasm").stat().st_size > 1_000_000


def test_mixture_near_critical_fallbacks(build_mod, tmp_path):
    """n-Butane/Propane 50/50 (mass): CoolProp's bubble point by temperature fails above about 105 C
    (critical point near 125.6 C). The page re-solves it with pressure-quality flashes, estimates the
    critical point from the phase envelope, and can interpolate saturation states from the envelope and
    evaluate single-phase states from (density, temperature)."""
    import json, shutil, subprocess
    js, node = build_mod.CP_DIR / "coolprop.js", shutil.which("node")
    if not js.exists() or node is None:
        pytest.skip("coolprop.js or node not available")
    anchor = 'Module["decrementExceptionRefcount"]=decrementExceptionRefcount;'
    (tmp_path / "coolprop.mjs").write_text(js.read_text().replace(
        anchor, anchor + 'Module["sghpStackSave"]=stackSave;Module["sghpStackRestore"]=stackRestore;', 1))
    shutil.copy(js.parent / "coolprop.wasm", tmp_path / "coolprop.wasm")
    live = ROOT / "local_scripts" / "sghp_live"
    (tmp_path / "t.mjs").write_text(f"""
import F from './coolprop.mjs'; import fs from 'fs'; import {{ createRequire }} from 'module';
const SGHP = createRequire(import.meta.url)({json.dumps(str(live / "sghp_model.js"))});
const src = fs.readFileSync({json.dumps(str(live / "sghp_live.js"))}, 'utf8');
const cut = (a, b) => src.slice(src.indexOf(a), src.indexOf(b));
const code = cut('function cpCall', '/* ---------------------------------------------------------------- property adapter')
  + cut('const LABEL', '/* ---------------------------------------------------------------- state */');
const CP = await F({{ locateFile: p => new URL(p, import.meta.url).pathname }});
const Fluid = new Function('CP', 'SGHP', code + '\\nreturn CoolPropFluid;')(CP, SGHP);
const fl = new Fluid(CP, [['n-Butane', 50], ['Propane', 50]], 'mass');
const out = {{ Tc: fl.T_crit, Pc: fl.P_crit, notes: fl.notes }};
const b = fl.TQ(115, 0); out.bub = b; out.check = fl._flash(fl.pairs.PQ_INPUTS, b.P * 1e6, 0, 0);
out.envBub = fl._envSat(0, 'P', b.P);
const P = 3.0; out.gasFast = fl._singleFast('gas', P, 'T', 120); out.gasDT = fl._dtSingle('gas', P, 'T', 120);
out.gasDTh = fl._dtSingle('gas', P, 'h', out.gasFast.h);
out.liqFast = fl._singleFast('liquid', P, 'T', 60); out.liqDT = fl._dtSingle('liquid', P, 'T', 60);
const cfg = {{ T_source: 55, dTmin_source: 5, T_steam: 110, dTmin_steam: 5, T_water_in: 80, T_water_out: 95, dTmin_water: 5, dT_superheat: 5, Q_heat_kW: null }};
const r = SGHP.solve(fl, cfg, SGHP.compressorMap('Recip-M'), 0.015); out.ok = r.success; out.msg = r.message;
// n-Butane/Isobutane 50/50: dew points above about 128 C need the envelope; each must stay quick
// (a phase-stability analysis on every density-temperature call once made these take 6-12 s)
const fl2 = new Fluid(CP, [['n-Butane', 50], ['IsoButane', 50]], 'mass');
const t0 = Date.now(); out.iso = [130, 140].map(T => fl2.TQ(T, 1)); out.isoMs = Date.now() - t0;
console.log(JSON.stringify(out));
""")
    res = subprocess.run([node, str(tmp_path / "t.mjs")], capture_output=True, text=True, timeout=600)
    assert res.returncode == 0, res.stderr[-800:]
    o = json.loads(res.stdout.strip().splitlines()[-1])
    assert 124.5 < o["Tc"] < 126.5 and 4.2 < o["Pc"] < 4.5                       # envelope critical point
    assert any("phase envelope" in n for n in o["notes"])
    assert o["bub"]["fallback"] == "PQ" and o["bub"]["T"] == 115                 # exact, via PQ
    assert o["check"]["T"] == pytest.approx(115.0, abs=1e-5) and o["bub"]["D"] == pytest.approx(o["check"]["D"], rel=1e-6)
    assert o["envBub"]["T"] == pytest.approx(115.0, abs=0.3)                     # envelope interpolation
    assert o["envBub"]["D"] == pytest.approx(o["check"]["D"], rel=0.03)
    for k in ("h", "s", "D"):
        assert o["gasDT"][k] == pytest.approx(o["gasFast"][k], rel=1e-6), k
        assert o["liqDT"][k] == pytest.approx(o["liqFast"][k], rel=1e-6), k
    assert o["gasDTh"]["T"] == pytest.approx(120.0, abs=1e-5)
    assert o["ok"], o["msg"]
    assert all(st["fallback"] is True for st in o["iso"]) and o["iso"][0]["P"] < o["iso"][1]["P"]
    assert o["isoMs"] < 3000, o["isoMs"]
