"""Supercritical high side of the SGHP model (ConventionalHeatPump and the live page's sghp_model.js).

At or above the critical temperature there is no bubble point to pin: the gas-cooler exit sits at the steam
temperature plus dTmin_steam on the critical isochore (rho = rho_crit), which sets the high pressure; the
compressor map reads the critical-isochore state at the discharge pressure. Checked with real CoolProp
(the page's WebAssembly build, driven from Node): the Python and browser models agree, the pressure and the
duties are continuous through the critical point, and the high-side states are where they should be.
"""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from HeatPumpLab.cycles.ConventionalHeatPump import ConventionalHeatPump, ConventionalHPConfig, recip_r600a

ROOT = Path(__file__).resolve().parents[1]
LIVE = ROOT / "local_scripts" / "sghp_live"
CP_JS = ROOT / "local_scripts" / "coolprop_js" / "coolprop.js"

SERVER = """
import F from './coolprop.mjs'; import fs from 'fs'; import readline from 'readline'; import { createRequire } from 'module';
const SGHP = createRequire(import.meta.url)(%(model)s);
const src = fs.readFileSync(%(live)s, 'utf8');
const cut = (a, b) => src.slice(src.indexOf(a), src.indexOf(b));
const code = cut('function cpCall', '/* ---------------------------------------------------------------- property adapter')
  + cut('const LABEL', '/* ---------------------------------------------------------------- state */');
const CP = await F({ locateFile: p => new URL(p, import.meta.url).pathname });
const Fluid = new Function('CP', 'SGHP', code + '\\nreturn CoolPropFluid;')(CP, SGHP);
const fl = new Fluid(CP, JSON.parse(process.argv[2]), 'mass');
const out = o => process.stdout.write(JSON.stringify(o) + '\\n');
out({ T_crit: fl.T_crit, P_crit: fl.P_crit, rho_crit: fl.rho_crit });
readline.createInterface({ input: process.stdin }).on('line', line => {
  const m = JSON.parse(line);
  try {
    if (m.op === 'solve') { const r = SGHP.solve(fl, m.cfg, SGHP.compressorMap('Recip-M'), m.UA); out({ ok: r.success, msg: r.message, report: r.report }); return; }
    const s = fl[m.op](m.a, m.b); out({ ok: true, st: { T: s.T, P: s.P, h: s.h, s: s.s, D: s.D, phase: s.phase, q: s.q } });
  } catch (e) { out({ ok: false, err: String(e.message || e) }); }
});
"""


class NodeFluid:
    """HeatPumpLab property interface over the page's CoolProp fluid (one Node process)."""

    def __init__(self, tmp, comps):
        node = shutil.which("node")
        if node is None or not CP_JS.exists():
            pytest.skip("node or local_scripts/coolprop_js/coolprop.js not available")
        anchor = 'Module["decrementExceptionRefcount"]=decrementExceptionRefcount;'
        (tmp / "coolprop.mjs").write_text(CP_JS.read_text().replace(
            anchor, anchor + 'Module["sghpStackSave"]=stackSave;Module["sghpStackRestore"]=stackRestore;', 1))
        shutil.copy(CP_JS.parent / "coolprop.wasm", tmp / "coolprop.wasm")
        (tmp / "server.mjs").write_text(SERVER % dict(model=json.dumps(str(LIVE / "sghp_model.js")),
                                                      live=json.dumps(str(LIVE / "sghp_live.js"))))
        self.p = subprocess.Popen([node, str(tmp / "server.mjs"), json.dumps(comps)], stdin=subprocess.PIPE,
                                  stdout=subprocess.PIPE, text=True, cwd=tmp)
        h = json.loads(self.p.stdout.readline())
        self.T_crit, self.P_crit, self.rho_crit = h["T_crit"], h["P_crit"], h["rho_crit"]

    def call(self, m):
        self.p.stdin.write(json.dumps(m) + "\n")
        self.p.stdin.flush()
        return json.loads(self.p.stdout.readline())

    def _st(self, op, a, b):
        r = self.call(dict(op=op, a=a, b=b))
        if not r["ok"]:
            raise ValueError(r["err"])
        st = r["st"]
        st["q"] = -1.0 if st["q"] is None else st["q"]
        return st

    def ref_PT(self, P, T): return self._st("PT", P, T)
    def ref_PH(self, P, h): return self._st("PH", P, h)
    def ref_PS(self, P, s): return self._st("PS", P, s)
    def ref_PQ(self, P, q): return self._st("PQ", P, q)
    def ref_TQ(self, T, q): return self._st("TQ", T, q)
    def ref_TD(self, T, D): return self._st("TD", T, D)

    def close(self):
        self.p.stdin.close()
        self.p.wait(timeout=30)


@pytest.fixture(scope="module")
def propane(tmp_path_factory):
    fl = NodeFluid(tmp_path_factory.mktemp("cp"), [["Propane", 100]])
    yield fl
    fl.close()


BASE = dict(T_source=50.0, dTmin_source=5.0, T_steam=120.0, dTmin_steam=5.0, T_water_in=80.0, T_water_out=95.0,
            dTmin_water=5.0, dT_superheat=15.0, dTmin_ihx=5.0, T_ambient=15.0, preheat_subcool_share=1.0)


def py_solve(fl, **kw):
    cfg = ConventionalHPConfig(refrigerant="Propane", loss_UA_kWK=0.015, compressor=recip_r600a("Recip-M"), **(BASE | kw))
    return ConventionalHeatPump(cfg, props=fl).solve()


def test_gas_cooler_on_the_critical_isochore(propane):
    fl = propane
    assert 96 < fl.T_crit < 97.5 and 215 < fl.rho_crit < 225
    res = py_solve(fl)                                         # 125 C gas-cooler exit, above T_crit
    assert res.success, res.message
    r, st = res.report, res.states
    assert r["supercritical"] and r["rho_crit_kgm3"] == pytest.approx(fl.rho_crit)
    gc = fl.ref_TD(125.0, fl.rho_crit)
    assert res.P_H == pytest.approx(gc["P"], rel=1e-12) and res.P_H > fl.P_crit
    assert st["cond_bub"]["T"] == pytest.approx(125.0) and st["cond_bub"]["D"] == pytest.approx(fl.rho_crit)
    # no pressure drop: the map's high-side state is the gas-cooler exit itself
    assert st["cond_dew"]["T"] == pytest.approx(125.0, abs=1e-6) and st["cond_dew"]["D"] == pytest.approx(fl.rho_crit, rel=1e-6)
    assert r["steam_gen_min_approach_K"] == pytest.approx(5.0, abs=1e-6)
    assert st["subcool_out"]["T"] == pytest.approx(85.0)       # subcooler pinch at the cold end, as before
    # with a condenser pressure drop the compressor discharges higher, up the critical isochore
    dp = py_solve(fl, dP_K_cond=0.2)
    assert dp.success and dp.report["P_discharge_bar"] > r["P_discharge_bar"]
    assert dp.states["cond_dew"]["T"] > 125.0 and dp.states["cond_dew"]["D"] == pytest.approx(fl.rho_crit, rel=1e-6)
    assert dp.report["T_cond_bubble"] == pytest.approx(125.0)   # exit pinned


def test_browser_model_matches_python_supercritical(propane):
    fl = propane
    keys = ("COP", "P_cond_bar", "P_discharge_bar", "T_cond_dew", "T_cond_bubble", "T_discharge", "m_steam_kgs",
            "Q_cond_kW", "Q_subcool_kW", "W_electric_kW", "m_dot_kgs", "dT_subcool_cond", "T_subcool_out",
            "steam_gen_min_approach_K", "rho_crit_kgm3", "density_ratio", "density_ratio_beyond", "to_read", "tc_read",
            "w_swept_kJm3")
    for over in (dict(), dict(preheat_subcool_share=0.0), dict(dP_K_cond=0.2, dP_K_evap=0.05), dict(T_steam=100.0)):
        py = py_solve(fl, **over)
        js = fl.call(dict(op="solve", cfg=BASE | over, UA=0.015))
        assert py.success and js["ok"], (over, py.message, js.get("msg"))
        assert js["report"]["supercritical"] is True and py.report["supercritical"] is True
        for k in keys:
            assert js["report"][k] == pytest.approx(py.report[k], rel=1e-6, abs=1e-6), (over, k)
        assert [f["key"] for f in js["report"]["limit_flags"]] == [f["key"] for f in py.report["limit_flags"]]
        assert js["report"]["window_edge"] == py.report["window_edge"]


def test_continuous_through_the_critical_point(propane):
    """Bubble-point rule below T_crit, critical-isochore rule above: the two meet at the critical point."""
    fl = propane
    kw = dict(T_source=40.0, T_water_in=50.0, T_water_out=70.0)
    Ts_c = fl.T_crit - 5.0
    lo, hi = py_solve(fl, T_steam=Ts_c - 0.02, **kw), py_solve(fl, T_steam=Ts_c + 0.02, **kw)
    assert lo.success and hi.success
    assert not lo.report["supercritical"] and hi.report["supercritical"]
    assert hi.P_H == pytest.approx(lo.P_H, rel=2e-3)
    for k in ("COP", "Q_heat_kW", "W_electric_kW"):
        assert hi.report[k] == pytest.approx(lo.report[k], rel=5e-3), k
    # further up, the pressure keeps rising and the cycle still solves
    far = py_solve(fl, T_steam=110.0, **kw)
    assert far.success and far.P_H > hi.P_H


def test_supercritical_point_flags_the_pressure_rating(propane):
    """65 bar on a compressor rated (with R600a) to 28.4 bar: solved, with a red flag; the compressor map is
    read past its window (the critical density is beyond the R600a table)."""
    res = py_solve(propane)
    f = {x["key"]: x for x in res.report["limit_flags"]}
    assert "P_dis" in f and f["P_dis"]["value"] == pytest.approx(res.report["P_discharge_bar"]) and f["P_dis"]["limit"] == 28.4
    assert "T_dis" in f                                         # 177 C discharge
    assert res.report["window_edge"] is not None


@pytest.mark.parametrize("backend", ["refprop", "coolprop"])
def test_supercritical_solves_on_each_property_backend(backend):
    """rho_crit is kg/m3 on both backends (REFPROP's per-mole value was once scaled twice and every flash failed)."""
    from HeatPumpLab.cycles.ConventionalHeatPump.conventional_hp import FluidProperties
    try:
        props = FluidProperties("Propane", backend)
    except Exception as exc:                                   # backend not installed here
        pytest.skip(f"{backend} unavailable: {exc}")
    assert 215 < props.rho_crit < 225
    cfg = ConventionalHPConfig(refrigerant="Propane", compressor=recip_r600a("Recip-M"), **BASE)
    res = ConventionalHeatPump(cfg, props=props).solve()
    assert res.success, res.message
    assert res.report["supercritical"] and res.report["rho_crit_kgm3"] == pytest.approx(props.rho_crit)
