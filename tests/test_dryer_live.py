"""Advanced dryer heat pumps page (local_scripts/dryer_live): the build, and the JavaScript cycle models over
CoolProp against the REFPROP study.

The regression cases are study designs (results/dryer_compare_v26/dashboard_data.json) with their REFPROP COP and
gas-cooler pressure written in here, so the test needs node and the CoolProp loader but not the study files:
RHP-A optima for every refrigerant family and approach, and SLTHP optima including supercritical gas coolers
(above the mixture's 41.7 bar critical pressure), discharge gas past CoolProp's 270 C inversion bound, a study
optimum on the evaporator-dew bound and one on the RHX1-branch bound (its rounded flow share is raised).
"""
from __future__ import annotations

import importlib.util
import json
import shutil
import subprocess
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
LIVE = ROOT / "local_scripts" / "dryer_live"
SGHP = ROOT / "local_scripts" / "sghp_live"
CP_JS = ROOT / "local_scripts" / "coolprop_js" / "coolprop.js"

# (fluid key, approach K, air C, lift K, split, GC outlet offset K, REFPROP COP, REFPROP gas-cooler pressure bar)
RHPA = [
    ("0.200-0.800", 20, 150, 80, 0.5, 0.0, 3.198, 16.45),
    ("0.200-0.800", 20, 200, 100, 0.5, 0.0, 2.761, 23.455),
    ("0.200-0.800", 30, 170, 0, 0.5, 0.0, 2.733, 64.496),          # no RHX1 lift, IHX on, supercritical
    ("0.200-0.800", 10, 180, 90, 0.5, 0.0, 3.173, 19.275),         # evaporator dew capped at 35 C
    ("Propane_Isobutane_Isopentane:0.325-0.200-0.475", 20, 170, 80, 0.5, 0.0, 2.974, 27.423),
    ("Propane_Butane_Isopentane:0.325-0.200-0.475", 20, 170, 85, 0.5, 0.0, 2.983, 25.747),
    ("Propane_Butane:0.250-0.750", 20, 170, 75, 0.5, 0.0, 2.827, 29.038),
]
# Propane/cyclopentane: CoolProp's mixing parameters differ from REFPROP's; about 0.15 % in COP
RHPA_C5 = ("Propane_Cyclopentane:0.400-0.600", 20, 170, 80, 0.5, 0.0, 2.685, 18.776)
# (fluid key, approach, air, lift, gas-cooler pressure bar, RHX1 flow share, REFPROP COP), power share 0.50
SLTHP = [
    ("0.200-0.800", 20, 170, 30, 42.714, 0.11, 2.761),
    ("0.200-0.800", 20, 200, 50, 55.23, 0.165, 2.46),
    ("0.200-0.800", 30, 190, 40, 59.959, 0.133, 2.335),            # flow share raised from 13.30 to 13.34 %
    ("0.200-0.800", 10, 150, 30, 26.69, 0.179, 3.25),               # evaporator dew on the 35 C bound
    ("0.200-0.800", 20, 150, 80, 15.134, 0.375, 2.976),
]


@pytest.fixture(scope="module")
def build_mod():
    spec = importlib.util.spec_from_file_location("build_dryer_live", LIVE / "build_dryer_live.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_build_page(build_mod, tmp_path):
    if not CP_JS.exists():
        pytest.skip("local_scripts/coolprop_js/coolprop.js not present")
    data = tmp_path / "data"
    data.mkdir()
    for name in build_mod.DATA_FILES:
        (data / name).write_text("{}")
    page = build_mod.build(tmp_path / "out", CP_JS, data)
    html = page.read_text()
    assert "<title>Advanced dryer heat pumps</title>" in html
    assert "Advanced dryer heat pumps | Ahuora Centre" in html
    assert "Barlow Semi Condensed" in html and "--accent:#1f5f8b" in html          # the steam page's style
    assert 'Module["sghpStackSave"]=stackSave' in html
    assert "class CoolPropFluid" in html and "root.DRYER = api" in html
    assert "bindLive();\nstudyBoot();" in html
    assert "const D = {" not in html                                                # data fetched, not inlined
    for tag in ('id="case-save"', 'id="cmp-results"', 'id="cmp-side"', 'id="cmp-order"', 'id="sweep-run"', 'casesBoot();', 'function solveCfg(cfg)', 'function cleanInputs(raw)',
                'id="lvTmt"', 'id="lvPhi"', 'id="lvPrMax"', 'id="lvStages"', 'id="lvIhx"', 'id="h-lvcycle"', "['H', 'Conventional']",
                'id="lvOpt"', 'function optimise(c)', 'cascadeOptimiseSteps', 'boosterOptimiseSteps',
                'id="lvSrc"', 'id="lvSrcT"', 'id="lvSrcM"', 'id="lvSrcW"', 'function withSource(cfg, solveOne)', 'id="lvFanSec"', 'id="lvFanPR"', 'function splitFields(p)'):
        assert tag in html, tag                                                     # secondary air, PR limit, conventional, Cycle box
    for name in build_mod.DATA_FILES + ("coolprop.wasm",):
        assert (page.parent / name).exists(), name
    page2 = build_mod.build(tmp_path / "out", CP_JS, data)                          # rebuild over read-only copies
    assert page2 == page


def test_build_needs_coolprop(build_mod, tmp_path):
    with pytest.raises(SystemExit, match="coolprop.js not found"):
        build_mod.build(tmp_path / "out", tmp_path / "nowhere" / "coolprop.js", None)


def _run_node(tmp_path, body: str, timeout=900):
    node = shutil.which("node")
    if not CP_JS.exists() or node is None:
        pytest.skip("coolprop.js or node not available")
    anchor = 'Module["decrementExceptionRefcount"]=decrementExceptionRefcount;'
    (tmp_path / "coolprop.mjs").write_text(CP_JS.read_text().replace(
        anchor, anchor + 'Module["sghpStackSave"]=stackSave;Module["sghpStackRestore"]=stackRestore;', 1))
    shutil.copy(CP_JS.parent / "coolprop.wasm", tmp_path / "coolprop.wasm")
    (tmp_path / "t.mjs").write_text(f"""
import F from './coolprop.mjs'; import fs from 'fs'; import {{ createRequire }} from 'module';
const req = createRequire(import.meta.url);
const SGHP = req({json.dumps(str(SGHP / "sghp_model.js"))}), DRY = req({json.dumps(str(LIVE / "dryer_model.js"))});
const src = fs.readFileSync({json.dumps(str(SGHP / "sghp_live.js"))}, 'utf8');
const code = src.slice(0, src.indexOf('/* ---------------------------------------------------------------- state */'));
const CP = await F({{ locateFile: p => new URL(p, import.meta.url).pathname }});
globalThis.document = {{ getElementById: () => ({{}}) }};
const Fluid = new Function('CP_', 'SGHP', code.replace('let CP = null', 'let CP = CP_') + '\\nreturn CoolPropFluid;')(CP, SGHP);
const NAME = {{ Propane: 'Propane', Isobutane: 'IsoButane', Butane: 'n-Butane', Isopentane: 'Isopentane', Cyclopentane: 'Cyclopentane' }};
const comps = w => {{ const [mix, fr] = w.includes(':') ? w.split(':') : ['Propane_Isopentane', w]; const f = fr.split('-').map(Number);
  return mix.split('_').map((n, i) => [NAME[n], f[i] * 100]); }};
const air = DRY.airModel(T => CP.PropsSI('H', 'P', 101000, 'T', T + 273.15, 'Air') / 1000);
const fluids = {{}}; const fluid = w => fluids[w] || (fluids[w] = new Fluid(CP, comps(w), 'mass'));
{body}
""")
    res = subprocess.run([node, str(tmp_path / "t.mjs")], capture_output=True, text=True, timeout=timeout)
    assert res.returncode == 0, res.stderr[-1500:]
    return json.loads(res.stdout.strip().splitlines()[-1])


def test_rhpa_matches_refprop(tmp_path):
    cases = RHPA + [RHPA_C5]
    out = _run_node(tmp_path, f"""
const out = [];
for (const [w, a, T, lift, split, off, cop, P] of {json.dumps(cases)}) {{
  const r = DRY.rhpa(fluid(w), air, {{ T_air_tar: T, dT_GC_min: a, lift, split, gc_out_offset: off, P_H_guess: P / 10 }});
  out.push(r.success ? {{ ok: true, COP: r.COP, P: r.P_H * 10, paths: DRY.paths(fluid(w), r).length }} : {{ ok: false, msg: r.message }});
}}
console.log(JSON.stringify(out));
""")
    for case, r in zip(cases, out):
        assert r["ok"], (case, r)
        tol = 0.002 if case is RHPA_C5 else 3e-4
        assert r["COP"] == pytest.approx(case[6], rel=tol), case       # study COP is stored to 3 dp
        assert r["P"] == pytest.approx(case[7], rel=2e-3), case
        assert r["paths"] >= 6


def test_slthp_matches_refprop(tmp_path):
    out = _run_node(tmp_path, f"""
const out = [];
for (const [w, a, T, lift, P, alpha, cop] of {json.dumps(SLTHP)}) {{
  const r = DRY.slthp(fluid(w), air, {{ T_air_tar: T, dT_GC_min: a, lift, split: 0.5, P_H: P / 10, alpha }});
  out.push(r.success ? {{ ok: true, COP: r.COP, warnings: r.warnings, paths: DRY.paths(fluid(w), r).length }} : {{ ok: false, msg: r.message }});
}}
console.log(JSON.stringify(out));
""")
    for case, r in zip(SLTHP, out):
        assert r["ok"], (case, r)
        assert r["COP"] == pytest.approx(case[6], rel=5e-4), case
        assert r["paths"] >= 7
    assert any("raised" in w for w in out[2]["warnings"])


def test_supercritical_roots(tmp_path):
    """Above the critical pressure a phase-imposed solve can return the other density root; the adapter takes
    the liquid-like root below the critical temperature and the gas-like root above it, and inverts hot gas
    past CoolProp's mole-averaged maximum temperature itself."""
    o = _run_node(tmp_path, """
const fl = fluid('0.200-0.800'), P = 5.523;
const liq = fl.PT(P, 35), gas = fl.PT(P, 250), hot = fl.PT(P, 300);
const back = fl.PH(P, hot.h), mid = fl.PH(P, liq.h);
console.log(JSON.stringify({ Tc: fl.T_crit, Pc: fl.P_crit, liqD: liq.D, gasD: gas.D, hotT: back.T, midT: mid.T }));
""", timeout=300)
    assert 170 < o["Tc"] < 171.5 and 4.1 < o["Pc"] < 4.25
    assert o["liqD"] > 450 and o["gasD"] < 150
    assert o["hotT"] == pytest.approx(300.0, abs=1e-4)
    assert o["midT"] == pytest.approx(35.0, abs=1e-4)


def test_source_water_sets_the_evaporator(tmp_path):
    """With a finite source-water flow the suction pressure is the highest at which the evaporator and superheater
    stay dTmin_rhx (5 K) below the water line; water 40 -> 15 C pinches inside the glide, so the suction pressure
    falls below the study's (evaporator inlet 10 C against an isothermal 40 C source) and COP with it."""
    o = _run_node(tmp_path, """
const fl = fluid('0.200-0.800'), cfg = { T_air_tar: 170, dT_GC_min: 20, lift: 90, split: 0.4, gc_out_offset: 0, P_H_guess: 1.93 };
const study = DRY.rhpa(fl, air, cfg), w15 = DRY.rhpa(fl, air, { ...cfg, T_source_out: 15 }), w30 = DRY.rhpa(fl, air, { ...cfg, T_source_out: 30 });
const s = DRY.slthp(fl, air, { T_air_tar: 170, dT_GC_min: 20, lift: 30, split: 0.5, P_H: 4.2714, alpha: 0.11, T_source_out: 15 });
const pick = r => ({ ok: r.success, msg: r.message, COP: r.COP, PL: r.P_L, pinch: r.report && r.report.evap_pinch_K, frac: r.report && r.report.power_frac_2 });
console.log(JSON.stringify({ study: pick(study), w15: pick(w15), w30: pick(w30), s: pick(s) }));
""")
    assert all(o[k]["ok"] for k in o), o
    assert o["study"]["pinch"] is None
    for k in ("w15", "w30", "s"):
        assert o[k]["pinch"] == pytest.approx(5.0, abs=0.1), k
    assert o["w15"]["PL"] < o["w30"]["PL"] and o["w15"]["PL"] < o["study"]["PL"]
    assert o["w15"]["COP"] < o["study"]["COP"]
    assert 0.3 < o["study"]["frac"] < 0.5                    # shaft-power share from the pressure split


def test_multi_loop_cycles(tmp_path):
    """The TT cascade (CO2/CO2), the 2-level cascade (isobutane/cyclopentane) and the booster (ammonia + MVR steam)
    solve at their defaults, keep their approaches and balance energy (heat to air = electricity + source heat)."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const co2 = new Fluid(CP, [['CarbonDioxide', 100]], 'mass'), ib = new Fluid(CP, [['IsoButane', 100]], 'mass');
const cp5 = new Fluid(CP, [['Cyclopentane', 100]], 'mass'), nh3 = new Fluid(CP, [['Ammonia', 100]], 'mass'), w = new Fluid(CP, [['Water', 100]], 'mass');
const base = {{ T_air_tar: 200, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30 }};
const k = C.kong(co2, co2, air, {{ ...base, T_source_out: 20 }}), n = C.cascade([ib, cp5], air, base), b = C.booster(nh3, w, air, base);
const pick = r => ({{ ok: r.success, msg: r.message, COP: r.COP, bal: r.success ? (r.q_air - r.w_shaft - r.q_source) / r.q_air : null, rep: r.report }});
console.log(JSON.stringify({{ k: pick(k), n: pick(n), b: pick(b) }}));
""", timeout=600)
    for key, lo, hi in (("k", 1.8, 3.2), ("n", 2.0, 3.5), ("b", 1.8, 3.2)):
        r = o[key]
        assert r["ok"], (key, r["msg"])
        assert lo < r["COP"] < hi, (key, r["COP"])
        assert abs(r["bal"]) < 0.02, (key, r["bal"])          # energy balance on the air, work and source heat
    assert o["k"]["rep"]["gap_MT"] == pytest.approx(20, abs=0.2)
    assert o["b"]["rep"]["P_steam_bar"] < 1.0                # vacuum steam at 85 C
    assert o["n"]["rep"]["levels"][1]["superheat"] > 5         # cyclopentane needs IHX superheat for a dry discharge


def test_secondary_air_composite(tmp_path):
    """Primary + secondary air as one stream per kg of primary air: the enthalpy rises (1 + phi) times as fast up to
    the secondary target and at the primary rate above it; phi = 0 (or a target outside the heating) leaves the air
    unchanged, and the inverse holds."""
    o = _run_node(tmp_path, """
const A = DRY.airSecondary(air, 15, 102, 0.5, 200), dh = (a, T1, T2) => a.h(T2) - a.h(T1);
console.log(JSON.stringify({
  below: dh(A, 15, 102) / dh(air, 15, 102), above: dh(A, 102, 200) / dh(air, 102, 200),
  inv: A.T(A.h(150)), inv2: A.T(A.h(60)),
  same0: DRY.airSecondary(air, 15, 102, 0, 200) === air, sameHi: DRY.airSecondary(air, 15, 210, 0.5, 200) === air,
  share: DRY.secondaryShare(air, 15, 102, 200, 0.5), shareCheck: 0.5 * dh(air, 15, 102) / (0.5 * dh(air, 15, 102) + dh(air, 15, 200)),
}));
""", timeout=300)
    assert o["below"] == pytest.approx(1.5, rel=1e-3)
    assert o["above"] == pytest.approx(1.0, rel=1e-3)
    assert o["inv"] == pytest.approx(150, abs=0.05) and o["inv2"] == pytest.approx(60, abs=0.05)
    assert o["same0"] and o["sameHi"]
    assert o["share"] == pytest.approx(o["shareCheck"], rel=1e-9) and 0.18 < o["share"] < 0.2


def test_secondary_air_every_cycle(tmp_path):
    """Every design takes the secondary air: phi = 0 reproduces the primary-only result exactly; with phi = 0.5 the
    heat per kg of primary air follows the composite curve and the COP rises (more of the heat is low-temperature);
    the TT cascade gives the same answer whether it splits the air itself or is handed the composite air."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const fl = fluid('0.200-0.800'), A = DRY.airSecondary(air, 15, 102, 0.5, 170);
const cfg = {{ T_air_tar: 170, dT_GC_min: 20, lift: 80, split: 0.5, gc_out_offset: 0, P_H_guess: 2.0 }};
const r0 = DRY.rhpa(fl, air, cfg), rA = DRY.rhpa(fl, DRY.airSecondary(air, 15, 102, 0, 170), cfg), r5 = DRY.rhpa(fl, A, cfg);
const s0 = DRY.slthp(fl, air, {{ T_air_tar: 170, dT_GC_min: 20, lift: 30, split: 0.5, P_H: 4.2714, alpha: 0.11 }});

const s5 = DRY.slthp(fl, A, {{ T_air_tar: 170, dT_GC_min: 20, lift: 30, split: 0.5, P_H: 4.2714, alpha: 0.11 }});
const co2 = new Fluid(CP, [['CarbonDioxide', 100]], 'mass');
DRY.setCompModel(1, 0.7);
const kb = {{ T_air_tar: 200, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 20, motor_eff: 0.96, P_b: 15, P_t: 15, r: 0.43, eta: null }};
const kOwn = C.kong(co2, co2, air, {{ ...kb, air_ratio: 0.5, T_mt: 102 }});
const kComp = C.kong(co2, co2, DRY.airSecondary(air, 15, 102, 0.5, 200), {{ ...kb, air_ratio: 0, T_mt: 108 }});
const kOne = C.kong(co2, co2, air, {{ ...kb, air_ratio: 0, T_mt: 108 }});
DRY.setCompModel(5);
const n0 = C.cascade([new Fluid(CP, [['IsoButane', 100]], 'mass'), new Fluid(CP, [['Cyclopentane', 100]], 'mass')], air, {{ T_air_tar: 200, T_source_out: 30 }});
const n5 = C.cascade([new Fluid(CP, [['IsoButane', 100]], 'mass'), new Fluid(CP, [['Cyclopentane', 100]], 'mass')], DRY.airSecondary(air, 15, 102, 0.5, 200), {{ T_air_tar: 200, T_source_out: 30 }});
const pick = r => ({{ ok: r.success, msg: r.message, COP: r.COP, q: r.success ? (r.report.q_air_kJkg ?? r.q_air) : null, m: r.success && r.report.air_per_ref_kgkg, marg: r.success ? r.report.min_approach_margin_K : null }});
console.log(JSON.stringify({{ r0: pick(r0), rA: pick(rA), r5: pick(r5), s0: pick(s0), s5: pick(s5), kOwn: pick(kOwn), kComp: pick(kComp), kOne: pick(kOne), n0: pick(n0), n5: pick(n5), hP: air.h(170) - air.h(15),
  hA: A.h(170) - A.h(15) }}));
""", timeout=900)
    for k, r in o.items():
        if isinstance(r, dict):
            assert r["ok"], (k, r["msg"])
    assert o["rA"]["COP"] == o["r0"]["COP"]                                      # phi = 0 is the plain air
    assert o["r5"]["COP"] > o["r0"]["COP"] and o["n5"]["COP"] > o["n0"]["COP"]
    # the SLTHP takes its gas-cooler pressure (and so its refrigerant cycle) as an input: the same duty now heats
    # less primary air, by the ratio of the air enthalpy rises
    assert o["s5"]["COP"] == pytest.approx(o["s0"]["COP"], rel=1e-9)
    assert o["s5"]["m"] / o["s0"]["m"] == pytest.approx(o["hP"] / o["hA"], rel=2e-3)
    assert o["r5"]["q"] / o["r5"]["m"] == pytest.approx(o["hA"], rel=1e-3)        # heat per kg primary air
    assert o["kComp"]["COP"] == pytest.approx(o["kOwn"]["COP"], rel=3e-3)         # two representations agree
    assert o["kOne"]["COP"] == pytest.approx(2.22, rel=0.01)                      # Kong thesis Table 12, TT-4 at 200 C
    assert o["kOwn"]["COP"] == pytest.approx(2.466, rel=0.005)                    # two gas coolers, a third of the air at 102 C


def test_conventional_cycle(tmp_path):
    """The conventional heat pump is the RHP model with no recuperator: one stage (split 0) or two, IHX optional.
    n-butane to 200 C air from water 40 -> 15 C needs a single-stage pressure ratio above 20 (1.48 bar suction,
    supercritical gas cooler), so it stops at the default limit and solves when the limit is raised; the same model
    as RHP-B (RHX2 65 K) on pure n-butane solves with the superheater outlet on the dew point."""
    o = _run_node(tmp_path, """
const nb = new Fluid(CP, [['n-Butane', 100]], 'mass');
DRY.setCompModel(5);
const base = { T_air_tar: 200, dT_GC_min: 20, T_source_in: 40, T_air_inlet: 15, T_source_out: 15, motor_eff: 0.93, lift: 0, lift2: 0, split: 0, ihx_on: true, gc_out_offset: 0 };
const capped = DRY.rhpa(nb, air, base), open = DRY.rhpa(nb, air, { ...base, max_stage_pr: 35 });
const fl = fluid('0.200-0.800'), conv = { T_air_tar: 170, dT_GC_min: 20, T_source_out: 30, lift: 0, lift2: 0, ihx_on: true };
const one = DRY.rhpa(fl, air, { ...conv, split: 0, gc_out_offset: 40 }), two = DRY.rhpa(fl, air, { ...conv, split: 0.45, gc_out_offset: 15 });
const dryNo = DRY.rhpa(fl, air, { ...conv, split: 0.45, gc_out_offset: 15, ihx_on: false });
DRY.setCompModel(1);
const rhpb = DRY.rhpa(nb, air, { T_air_tar: 170, dT_GC_min: 20, T_source_out: 30, motor_eff: 0.93, lift: 0, lift2: 65, split: 0.45, gc_out_offset: 0 });
DRY.setCompModel(5);
const pick = r => ({ ok: r.success, msg: r.message, COP: r.COP, rep: r.success ? r.report : null, pr: r.stage_pr });
console.log(JSON.stringify({ capped: pick(capped), open: pick(open), one: pick(one), two: pick(two), dryNo: pick(dryNo), rhpb: pick(rhpb) }));
""", timeout=900)
    assert not o["capped"]["ok"] and "stage pressure ratio above 20" in o["capped"]["msg"]
    assert o["open"]["ok"] and 30 < o["open"]["rep"]["PR_1"] < 35 and o["open"]["rep"]["PR_2"] == pytest.approx(1, abs=1e-6)
    assert o["open"]["COP"] == pytest.approx(1.94, abs=0.02)
    assert o["one"]["ok"] and o["one"]["rep"]["PR_2"] == pytest.approx(1, abs=1e-6) and o["one"]["rep"]["q_rhx1_kJkg"] == 0
    assert o["two"]["ok"] and o["two"]["COP"] > o["one"]["COP"]
    assert not o["dryNo"]["ok"] and "wet compression" in o["dryNo"]["msg"]       # isopentane-rich blend needs the IHX
    assert o["rhpb"]["ok"] and o["rhpb"]["COP"] == pytest.approx(2.845, abs=0.005)


def test_tt_optimiser(tmp_path):
    """The TT optimiser never returns a worse design than it started from, returns one that re-solves, and on the
    thesis case settles near Kong's optimum (both pressures at the 150 bar limit, bottom/top flow about 0.43)."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const co2 = new Fluid(CP, [['CarbonDioxide', 100]], 'mass');
DRY.setCompModel(1, 0.7);
const kb = {{ T_air_tar: 200, T_mt: 102, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 20, motor_eff: 0.96, P_b: 13, P_t: 14, r: 0.6, air_ratio: 0.5, eta: null }};
const opt = C.kongOptimise(co2, co2, air, kb, {{ Pmax: 15 }});
const again = C.kong(co2, co2, air, {{ ...kb, P_b: opt.P_b, P_t: opt.P_t, r: opt.r, gc_off: opt.gc_off }});
console.log(JSON.stringify({{ opt, again: {{ ok: again.success, COP: again.COP }} }}));
""", timeout=900)
    opt = o["opt"]
    assert opt["success"] and opt["COP"] >= opt["start"] - 1e-9
    assert o["again"]["ok"] and o["again"]["COP"] == pytest.approx(opt["COP"], rel=1e-9)
    assert opt["P_t"] == pytest.approx(15, abs=1e-6) and opt["P_b"] > 14.0
    assert 0.38 < opt["r"] < 0.48
    assert opt["COP"] == pytest.approx(2.467, abs=0.01)


def test_cascade_and_booster_optimisers(tmp_path):
    """The n-level cascade and booster optimisers never return a worse design than a feasible start, return one that
    re-solves to the same COP with every air approach kept, and beat a coarse scan of the same variables."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const ib = new Fluid(CP, [['IsoButane', 100]], 'mass'), cp5 = new Fluid(CP, [['Cyclopentane', 100]], 'mass');
const nh3 = new Fluid(CP, [['Ammonia', 100]], 'mass'), w = new Fluid(CP, [['Water', 100]], 'mass');
const base = {{ T_air_tar: 200, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30 }};
const n = C.cascadeOptimise([ib, cp5], air, {{ ...base, T_cond: [90] }});
const nA = C.cascade([ib, cp5], air, {{ ...base, T_cond: n.T_cond }});
const nScan = [];
for (let T = 60; T <= 130; T += 5) {{ const r = C.cascade([ib, cp5], air, {{ ...base, T_cond: [T] }}); if (r.success && !r.warnings.some(x => /approach|IHX/.test(x))) nScan.push(r.COP); }}
const n3 = C.cascadeOptimise([ib, cp5, cp5], air, {{ ...base, n: 3, T_cond: [80, 140] }});
const b = C.boosterOptimise(nh3, w, air, base);
const bA = C.booster(nh3, w, air, {{ ...base, T_cond: b.T_cond, n_mvr: b.n_mvr }});
const bScan = [];
for (let m = 1; m <= 12; m++) for (let T = 60; T <= 125; T += 5) {{ const r = C.booster(nh3, w, air, {{ ...base, T_cond: T, n_mvr: m }});
  if (r.success && r.report.gap_steam >= 19.95 && r.report.gap_bottom >= 19.95) bScan.push(r.COP); }}
console.log(JSON.stringify({{ n, nA: {{ ok: nA.success, COP: nA.COP, w: nA.warnings }}, nScan: Math.max(...nScan), n3,
  b, bA: {{ ok: bA.success, COP: bA.COP, w: bA.warnings, g1: bA.report.gap_bottom, g2: bA.report.gap_steam }}, bScan: Math.max(...bScan) }}));
""", timeout=900)
    n, b = o["n"], o["b"]
    assert n["success"] and n["startOk"] and n["COP"] >= n["start"] - 1e-9
    assert o["nA"]["ok"] and o["nA"]["COP"] == pytest.approx(n["COP"], rel=1e-9) and not o["nA"]["w"]
    assert n["COP"] >= o["nScan"] - 1e-3                    # within the 0.1 K rounding of the result
    assert o["n3"]["success"] and o["n3"]["COP"] >= o["n3"]["start"] - 1e-9
    assert o["n3"]["T_cond"][0] < o["n3"]["T_cond"][1]
    assert b["success"] and b["COP"] >= b["start"] - 1e-9 and isinstance(b["n_mvr"], int)
    assert o["bA"]["ok"] and o["bA"]["COP"] == pytest.approx(b["COP"], rel=1e-9)
    assert o["bA"]["g1"] >= 19.95 and o["bA"]["g2"] >= 19.95
    assert b["COP"] >= o["bScan"] - 1e-6


def test_mvr_fan_exit_superheat(tmp_path):
    """HP + MVR: water injected after every fan but the last brings the steam to the set superheat above saturation;
    a larger superheat injects less water, and turning injection off leaves the steam to superheat stage on stage."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const nh3 = new Fluid(CP, [['Ammonia', 100]], 'mass'), w = new Fluid(CP, [['Water', 100]], 'mass');
const base = {{ T_air_tar: 200, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30 }};
const run = x => {{ const r = C.booster(nh3, w, air, {{ ...base, ...x }}); return {{ ok: r.success, COP: r.COP,
  inj: r.report.stages.map(s => s.inj), sat: r.report.stages.map(s => w.PQ(s.P_out, 1).T), Tout: r.report.stages.map(s => s.T_out) }}; }};
// the steam entering stage i+1 (after injection) = the next stage's suction; recover its temperature from the path
const suc = x => C.booster(nh3, w, air, {{ ...base, ...x }}).loops[1].paths.filter(p => p.name.startsWith('MVR')).map(p => p.pts[0][0]);
console.log(JSON.stringify({{ s5: run({{ sh_inject: 5 }}), s20: run({{ sh_inject: 20 }}), off: run({{ inject: false }}), suc20: suc({{ sh_inject: 20 }}) }}));
""", timeout=600)
    s5, s20, off = o["s5"], o["s20"], o["off"]
    assert s5["ok"] and s20["ok"] and off["ok"]
    assert all(i > 0 for i in s5["inj"][:-1]) and s5["inj"][-1] == 0      # every fan but the last
    assert s20["inj"][0] < s5["inj"][0]                                  # same fan-1 discharge: more superheat, less water
    for k in range(1, len(o["suc20"])):                                  # suction of fan k+1 = 20 K above saturation
        assert o["suc20"][k] == pytest.approx(s20["sat"][k - 1] + 20, abs=0.05)
    assert max(off["Tout"]) > max(s5["Tout"]) + 100 and off["COP"] < s5["COP"]


def test_mvr_fans_and_process_splits(tmp_path):
    """HP + MVR: the fans have their own efficiency model and motor (the compressor model leaves the fan work alone);
    with no process split all the steam goes through every fan as before; Auto sizes each fan's process split to the
    air approach, balances energy and beats no split; the same splits entered by hand give the same COP; and the
    optimiser (Auto) returns a design that re-solves with every approach kept."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const nh3 = new Fluid(CP, [['Ammonia', 100]], 'mass'), w = new Fluid(CP, [['Water', 100]], 'mass');
const base = {{ T_air_tar: 200, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, T_cond: 80, n_mvr: 11 }};
const sm = r => r.success ? {{ ok: true, COP: r.COP, w: r.warnings, splits: r.report.splits, g1: r.report.gap_bottom, g2: r.report.gap_steam,
  w_mvr: r.report.w_mvr / r.report.steam_per_kg, bal: (r.q_air - r.w_shaft - r.q_source) / r.q_air, tp: r.report.stages.map(s => s.to_process) }} : {{ ok: false, msg: r.message }};
const B = x => sm(C.booster(nh3, w, air, {{ ...base, ...x }}));
const none = B({{}}), fill = B({{ split_fill: true }});
const man = fill.ok ? B({{ split: fill.splits }}) : null;
DRY.setCompModel(1, 0.6); const fanA = B({{ fan_model: 1, fan_eta: 0.8 }});
DRY.setCompModel(1, 0.9); const fanB = B({{ fan_model: 1, fan_eta: 0.8 }});
const fan9 = B({{ fan_model: 1, fan_eta: 0.9 }}), mot1 = B({{ fan_model: 1, fan_eta: 0.8, fan_motor: 1 }});
DRY.setCompModel(5);
const opt = C.boosterOptimise(nh3, w, air, {{ ...base, split_fill: true }});
const again = sm(C.booster(nh3, w, air, {{ ...base, split_fill: true, T_cond: opt.T_cond, n_mvr: opt.n_mvr }}));
const optNone = C.boosterOptimise(nh3, w, air, base);
console.log(JSON.stringify({{ none, fill, man, fanA, fanB, fan9, mot1, opt, again, optNone }}));
""", timeout=900)
    none, fill, man = o["none"], o["fill"], o["man"]
    assert none["ok"] and fill["ok"] and man["ok"]
    assert all(t == 0 for t in none["tp"][:-1]) and none["tp"][-1] > 0          # no split: all to the last fan
    assert all(0 < x < 1 for x in fill["splits"]) and all(t > 0 for t in fill["tp"])
    assert abs(fill["bal"]) < 1e-6 and fill["g2"] >= 19.95 and fill["COP"] > none["COP"] + 0.2
    assert man["COP"] == pytest.approx(fill["COP"], rel=1e-9)                    # Auto = the same splits by hand
    assert o["fanA"]["w_mvr"] == pytest.approx(o["fanB"]["w_mvr"], rel=1e-9)    # fan work per kg steam: compressor model leaves it alone
    assert o["fan9"]["w_mvr"] < o["fanB"]["w_mvr"]
    assert o["mot1"]["COP"] > o["fanB"]["COP"]
    opt, again = o["opt"], o["again"]
    assert opt["success"] and again["ok"] and again["COP"] == pytest.approx(opt["COP"], rel=1e-9)
    assert again["g1"] >= 19.95 and again["g2"] >= 19.95
    assert opt["COP"] > o["optNone"]["COP"] + 0.2


def test_source_curve_every_model(tmp_path):
    """An air source reaches every model as a source curve (temperature against the fraction of the evaporator duty
    from its cold end): a straight curve gives exactly the water results; a curve that sags below the straight line
    (humid air condensing) tightens the RHP evaporator pinch; and the cascade, booster and TT evaporator lines follow it."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const co2 = new Fluid(CP, [['CarbonDioxide', 100]], 'mass'), ib = new Fluid(CP, [['IsoButane', 100]], 'mass'), cp5 = new Fluid(CP, [['Cyclopentane', 100]], 'mass');
const nh3 = new Fluid(CP, [['Ammonia', 100]], 'mass'), w = new Fluid(CP, [['Water', 100]], 'mass'), mix = fluid('0.200-0.800');
const base = {{ T_air_tar: 170, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30 }};
const line = f => 30 + 10 * f, sag = f => 30 + 10 * f * f;          // same ends; sag sits below the line
const A = x => DRY.rhpa(mix, air, {{ ...base, lift: 80, split: 0.45, gc_out_offset: 5, ...x }});
const a0 = A({{}}), a1 = A({{ source_curve: line }}), a2 = A({{ source_curve: sag }});
const n2 = C.cascade([ib, cp5], air, {{ ...base, source_curve: sag, source_name: 'source air' }});
const b2 = C.booster(nh3, w, air, {{ ...base, source_curve: sag, source_name: 'source air' }});
const kb = {{ ...base, T_air_tar: 200, T_mt: 102, T_source_out: 20, P_b: 15, P_t: 15, r: 0.43, air_ratio: 0.5, eta: null }};
DRY.setCompModel(1, 0.7);
const k2 = C.kong(co2, co2, air, {{ ...kb, source_curve: f => 20 + 20 * f * f }});
DRY.setCompModel(5);
const ev = r => r.hx.filter(h => /vaporator/.test(h.name)).map(h => ({{ name: h.name, other: h.otherName, hot: h.hot }}));
console.log(JSON.stringify({{ a0: [a0.success, a0.COP, a0.P_L], a1: [a1.success, a1.COP, a1.P_L], a2: [a2.success, a2.COP, a2.P_L],
  n2: n2.success ? ev(n2) : n2.message, b2: b2.success ? ev(b2) : b2.message,
  k2: k2.success ? {{ mid: k2.report.T_water_mid, fMid: k2.report.q_evap_bot / (k2.report.q_evap_bot + k2.report.q_evap_top), ev: ev(k2) }} : k2.message }}));
""", timeout=900)
    a0, a1, a2 = o["a0"], o["a1"], o["a2"]
    assert a0[0] and a1[0] and a2[0]
    assert a1[1] == pytest.approx(a0[1], rel=1e-9) and a1[2] == pytest.approx(a0[2], rel=1e-9)
    assert a2[2] < a0[2] and a2[1] < a0[1]                          # the sagging source forces a lower evaporator pressure
    for key in ("n2", "b2"):
        (e,) = o[key]
        assert e["other"] == "source air" and len(e["hot"]) == 16
        f = [p[0] / e["hot"][-1][0] for p in e["hot"]]
        assert all(abs(p[1] - (30 + 10 * x * x)) < 1e-9 for p, x in zip(e["hot"], f))
    k = o["k2"]
    assert k["mid"] == pytest.approx(20 + 20 * k["fMid"] ** 2, abs=1e-9)
    assert all(len(e["hot"]) == 16 for e in k["ev"])


def test_air_source_page(build_mod, tmp_path):
    """In the page, an air source (temperature, flow, moisture) is solved on its outlet: the evaporator duty at the
    swept volume equals the source air's enthalpy drop, moisture condenses below the dew point, switching back
    to water gives the water result again, and sizing by a target sink duty sets the swept volume to meet it."""
    sync_api = pytest.importorskip("playwright.sync_api")
    if not CP_JS.exists():
        pytest.skip("coolprop.js not present")
    import os, threading, http.server, functools
    data = Path(os.environ.get("DRYER_DATA", str(build_mod.DATA_DIR)))      # the page needs the real study data to start
    if not all((data / n).exists() for n in build_mod.DATA_FILES):
        pytest.skip(f"study data not found in {data} (set DRYER_DATA)")
    page = build_mod.build(tmp_path / "out", CP_JS, data)
    handler = functools.partial(http.server.SimpleHTTPRequestHandler, directory=str(page.parent))
    handler.log_message = lambda *a, **k: None
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    url = f"http://127.0.0.1:{srv.server_address[1]}/{page.name}"
    done = "(document.getElementById('solve-time').textContent||'').includes(' in ') && !LIVE.busy"
    try:
        with sync_api.sync_playwright() as p:
            b = p.chromium.launch(); pg = b.new_page(); errs = []
            pg.on("pageerror", lambda e: errs.append(str(e)))
            pg.goto(url); pg.wait_for_function(done, timeout=600000)
            pg.evaluate("LV.cyc = 'N'; LV.src = 'water'; renderLiveControls(); changed()")
            pg.wait_for_timeout(300); pg.wait_for_function(done, timeout=600000)
            water = pg.evaluate("LIVE.res.COP")
            pg.evaluate("LV.src = 'air'; LV.srcT = 20; LV.srcM = 4; LV.srcW = 13; renderLiveControls(); changed()")
            pg.wait_for_timeout(300); pg.wait_for_function(done, timeout=600000)
            o = pg.evaluate("""(() => { const r = LIVE.res, a = r.srcAir, A = LIVE.srcAir;
              return {ok: r.success, a, drop: a.mDry * (A.h(a.Tin) - A.h(a.Tout)), hidden: [$('lvTsRow').hidden, $('lvSrcT').closest('.pickrow').hidden],
                      line: document.querySelector('#lvHeadline .hl-s').textContent}; })()""")
            pg.evaluate("LV.src = 'water'; renderLiveControls(); changed()")
            pg.wait_for_timeout(300); pg.wait_for_function(done, timeout=600000)
            back = pg.evaluate("[LIVE.res.COP, LIVE.res.srcAir || null]")
            # sizing by a target sink duty: the swept volume follows so the heat to air equals the target
            pg.evaluate("state.sizeBy = 'duty'; state.duty = 300; sizeChanged()")
            pg.wait_for_timeout(300); pg.wait_for_function(done, timeout=600000)
            duty = pg.evaluate("[kwAt(LIVE.res, LIVE.cfg).air, state.vsw, document.getElementById('vsw').disabled]")
            b.close()
    finally:
        srv.shutdown()
    assert not errs, errs
    a = o["a"]
    assert o["ok"] and abs(a["resid"]) < 0.05
    assert a["kw"] == pytest.approx(o["drop"], rel=2e-3)            # duty = the source air's enthalpy drop
    assert a["Tout"] < a["Td"] and a["cond"] > 0 and a["Wout"] < a["W"]   # 13 g/kg at 20 C: condenses
    assert o["hidden"] == [True, False] and "kg/s of source air" in o["line"]
    assert back[1] is None and back[0] == pytest.approx(water, rel=1e-12)
    assert duty[0] == pytest.approx(300, rel=1e-9) and duty[1] > 0 and duty[2]


def test_mvr_bottom_subcooler_keeps_approach(tmp_path):
    """HP + MVR: a bottom fluid whose subcooled liquid bends (propane) would pinch inside the subcooler. Either the
    liquid outlet is raised or the air leaving the bottom cycle is lowered until the approach holds; the higher-COP
    fix is kept (for propane, lowering the air), with Auto splits too."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const c3 = new Fluid(CP, [['Propane', 100]], 'mass'), w = new Fluid(CP, [['Water', 100]], 'mass');
const out = [];
for (const split_fill of [false, true]) for (const T_cond of [66, 75, 85]) {{
  const r = C.booster(c3, w, air, {{ T_air_tar: 170, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, split_fill, T_cond, n_mvr: 10 }});
  const raiseOnly = C.booster(c3, w, air, {{ T_air_tar: 170, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, split_fill, T_cond, n_mvr: 10, t1_drop: 0 }});
  out.push({{ ok: r.success, g1: r.report.gap_bottom, raise: r.report.sub_raise, drop: r.report.t1_drop, COP: r.COP, COPraise: raiseOnly.COP, g1raise: raiseOnly.report.gap_bottom,
    w: r.warnings.filter(x => /bottom-cycle/.test(x)) }});
}}
console.log(JSON.stringify(out));
""", timeout=600)
    for r in o:
        assert r["ok"] and r["g1"] >= 19.95 and not r["w"]
        assert r["g1raise"] >= 19.95 and r["COP"] >= r["COPraise"] - 1e-12   # the better of the two fixes
    assert max(r["drop"] for r in o) > 5                                     # lowering the bottom air wins for propane


def test_cascade_lower_level_subcooling_capped(tmp_path):
    """n-level cascade: a lower level with a large flow could heat the air past its condensing temperature less the
    approach with its subcooled liquid alone; the subcooling is capped so every level keeps the air approach."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const ib = new Fluid(CP, [['IsoButane', 100]], 'mass'), cp5 = new Fluid(CP, [['Cyclopentane', 100]], 'mass');
const out = [];
for (const [n, tc] of [[2, [100]], [2, [120]], [3, [120, 170]]]) {{
  const r = C.cascade([ib, cp5, cp5], air, {{ T_air_tar: 200, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, n, T_cond: tc }});
  out.push({{ ok: r.success, msg: r.message, w: r.success ? r.warnings.filter(x => /approach/.test(x)) : [], lv: r.success ? r.report.levels.map(l => [l.Tc, l.air_in]) : [],
    bal: r.success ? (r.q_air - r.w_shaft - r.q_source) / r.q_air : null }});
}}
console.log(JSON.stringify(out));
""", timeout=600)
    for r in o:
        assert r["ok"], r["msg"]
        assert not r["w"], r["w"]
        assert abs(r["bal"]) < 0.02
        for k in range(len(r["lv"]) - 1):                      # air into the next level <= this level's Tc - approach
            assert r["lv"][k + 1][1] <= r["lv"][k][0] - 20 + 0.1


def test_cascade_level_ihx(tmp_path):
    """n-level cascade: an IHX on a level heats its suction to the liquid entering the IHX less the IHX approach; the
    exchanger is reported (key ihx) with its approach kept, energy still balances, and on the 3-level case the top
    level's IHX raises the COP."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const ib = new Fluid(CP, [['IsoButane', 100]], 'mass'), cp5 = new Fluid(CP, [['Cyclopentane', 100]], 'mass');
const base = {{ T_air_tar: 170, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, n: 3, T_cond: [90, 140] }};
const pick = r => ({{ ok: r.success, COP: r.COP, bal: (r.q_air - r.w_shaft - r.q_source) / r.q_air, w: r.warnings,
  ihx: r.hx.filter(h => h.key === 'ihx').map(h => ({{ name: h.name, gap: h.gap }})), sh: r.report.levels.map(l => l.superheat) }});
console.log(JSON.stringify({{ off: pick(C.cascade([ib, cp5, cp5], air, base)), on: pick(C.cascade([ib, cp5, cp5], air, {{ ...base, ihx: [false, false, true], dT_ihx: 5 }})) }}));
""", timeout=600)
    off, on = o["off"], o["on"]
    assert off["ok"] and on["ok"]
    assert abs(on["bal"]) < 1e-6
    assert on["COP"] > off["COP"] + 0.05 and on["sh"][2] > off["sh"][2] + 5
    top = [h for h in on["ihx"] if h["name"] == "Level 3 IHX"]
    assert top and top[0]["gap"] == pytest.approx(5, abs=0.1)
    assert not [w for w in on["w"] if "IHX" in w]


def test_mvr_bottom_ihx(tmp_path):
    """HP + MVR: the bottom IHX heats the suction to the liquid entering it less its approach; with the liquid not
    subcooled to the air it raises the COP for isobutane, the energy balance holds, and it is reported as an IHX."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const ib = new Fluid(CP, [['IsoButane', 100]], 'mass'), w = new Fluid(CP, [['Water', 100]], 'mass');
const base = {{ T_air_tar: 170, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, T_cond: 80, n_mvr: 10, split_fill: true, sub_air: false }};
const pick = r => ({{ ok: r.success, COP: r.COP, bal: (r.q_air - r.w_shaft - r.q_source) / r.q_air, rep: r.report, keys: r.hx.map(h => h.key) }});
console.log(JSON.stringify({{ off: pick(C.booster(ib, w, air, base)), on: pick(C.booster(ib, w, air, {{ ...base, ihx: true, dT_ihx: 5 }})) }}));
""", timeout=600)
    off, on = o["off"], o["on"]
    assert off["ok"] and on["ok"] and abs(on["bal"]) < 1e-6
    assert on["COP"] > off["COP"] + 0.1
    assert on["rep"]["T_suction"] == pytest.approx(80 - 5, abs=0.2)         # bubble-point liquid at 80 C less 5 K
    assert on["rep"]["gap_ihx"] >= 4.95 and "ihx" in on["keys"] and "ihx" not in off["keys"]


def test_tt_evaporators_in_series(tmp_path):
    """The source water passes the TT cascade's top evaporator, then its bottom one: one water line through both,
    the water between them set by the energy balance, and each evaporator's approach taken to the water leaving it."""
    CASC = LIVE / "dryer_cascades.js"
    o = _run_node(tmp_path, f"""
const C = req({json.dumps(str(CASC))});
const co2 = new Fluid(CP, [['CarbonDioxide', 100]], 'mass');
DRY.setCompModel(1, 0.7);
const kb = {{ T_air_tar: 200, T_mt: 102, T_air_inlet: 15, dT_GC_min: 20, T_source_in: 40, T_source_out: 20, motor_eff: 0.96, P_b: 15, P_t: 15, r: 0.43, air_ratio: 0.5, eta: null }};
const k = C.kong(co2, co2, air, kb), cold = C.kong(co2, co2, air, {{ ...kb, T_source_out: 11 }});
const ev = n => k.hx.find(h => h.name === n);
console.log(JSON.stringify({{ ok: k.success, rep: k.report, top: ev('Top evaporator').hot, bot: ev('Bottom evaporator').hot, topCold: ev('Top evaporator').cold,
  coldWarn: cold.success && cold.warnings.some(w => /bottom evaporator/.test(w)) }}));
""", timeout=600)
    rep = o["rep"]
    assert o["ok"]
    qt, qb, tm = rep["q_evap_top"], rep["q_evap_bot"], rep["T_water_mid"]
    assert (40 - tm) / (tm - 20) == pytest.approx(qt / qb, rel=1e-9)            # one water flow through both
    assert o["bot"][0] == [0, 20] and o["bot"][-1][0] == pytest.approx(qb) and o["bot"][-1][1] == pytest.approx(tm)
    assert o["top"][0][0] == pytest.approx(qb) and o["top"][0][1] == pytest.approx(tm) and o["top"][-1][1] == 40
    assert o["topCold"][0][0] == pytest.approx(qb)                               # top evaporator stacked after the bottom one
    assert rep["gap_evap_top"] == pytest.approx(tm - 20) and rep["gap_evap_bot"] == pytest.approx(10)
    assert o["coldWarn"]                                                         # water out 11 C against 10 C evaporating
