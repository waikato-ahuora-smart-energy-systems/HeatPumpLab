"""Assemble the Steam-generating heat pump app (apps/steam/index.html, CoolProp in the browser).

  uv run python apps/src/steam/build.py [--coolprop-js PATH] [--out DIR] [--mock]

Combines the dashboard's styles and chart renderers (dashboard_template.html), the
JavaScript cycle model (sghp_model.js), the CoolProp adapter and controls (sghp_live.js) and the
CoolProp Emscripten loader (coolprop.js, inlined). coolprop.wasm is copied next to index.html and
must be served with it (as an artifact supporting file, or by any static web server).
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

HERE = Path(__file__).resolve().parent
TEMPLATE = HERE / 'dashboard_template.html'
CP_DIR = HERE.parent / 'coolprop'

EXTRA_CSS = """
[hidden]{display:none!important}
.temp{gap:2px}
.slhead{display:flex;justify-content:space-between;align-items:baseline;gap:10px}
.slhead output{font-family:var(--f-data);font-size:14px;font-weight:500;color:var(--ink);font-variant-numeric:tabular-nums}
.controls{display:grid;grid-template-columns:minmax(0,1.2fr) minmax(0,1fr);gap:16px;align-items:start}   /* compressor and operating point left, refrigerant and settings right */
.controls > .ctl{align-content:start}
#h-set ~ .fields{grid-template-columns:repeat(4,minmax(0,1fr));align-items:end}
#h-set ~ .fields input{min-width:0;width:100%;box-sizing:border-box}
@media (max-width:1000px) and (min-width:861px),(max-width:480px){#h-set ~ .fields{grid-template-columns:repeat(2,minmax(0,1fr))}}
.ctl-col{display:grid;gap:16px;min-width:0}
.temps{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:6px 20px}
@media (max-width:560px){.temps{grid-template-columns:1fr}}
@media (max-width:860px){.controls{grid-template-columns:1fr}}   /* stacked on narrow screens: compressor and operating point first */
.ctl{background:var(--surface);border:1px solid var(--rule);border-radius:6px;padding:14px;display:grid;gap:10px;align-content:start;min-width:0}
.ctl h2{font-family:var(--f-display);font-weight:600;font-size:19px;margin:0}
.comps{display:grid;gap:8px}
.comp{display:grid;grid-template-columns:minmax(0,1.2fr) minmax(0,1fr) 74px 14px 28px;gap:8px;align-items:center}
@media (max-width:480px){.comp{grid-template-columns:minmax(0,1fr) 70px 14px 28px}.comp input[type=range]{display:none}}
select,input[type=number]{font:inherit;font-size:14px;color:var(--ink);background:var(--bg);border:1px solid var(--rule);border-radius:5px;padding:5px 7px;min-width:0;width:100%;box-sizing:border-box}
input[type=number]{font-family:var(--f-data);font-variant-numeric:tabular-nums}
input:disabled{opacity:.6}
select:focus-visible,input:focus-visible,.icon:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
.comp input[type=range]{width:100%;accent-color:var(--accent);min-width:0}
.unit{color:var(--ink-3);font-size:13px}
.pick{font:inherit;font-size:14px;display:flex;justify-content:space-between;align-items:center;gap:6px;min-width:0;width:100%;cursor:pointer;border:1px solid var(--rule);background:var(--bg);color:var(--ink);border-radius:5px;padding:5px 8px;text-align:left}
.pick span:first-child{overflow:hidden;text-overflow:ellipsis;white-space:nowrap}
.pick .caret{color:var(--ink-3)}
.pick[aria-expanded="true"]{border-color:var(--accent)}
.fluid-list{display:flex;flex-wrap:wrap;gap:6px;padding:8px;border:1px solid var(--rule);border-radius:6px;background:var(--bg)}
.fopt,.tbtn{font:inherit;font-size:13px;cursor:pointer;border:1px solid var(--rule);background:var(--surface);color:var(--ink-2);border-radius:999px;padding:4px 10px}
.fopt:hover:not(:disabled),.tbtn:hover{color:var(--ink);border-color:var(--accent)}
.fopt[aria-selected="true"],.tbtn[aria-pressed="true"]{background:var(--accent);border-color:var(--accent);color:var(--surface)}
.fopt:disabled{opacity:.4;cursor:default}
.toggle{display:inline-flex;gap:4px}
.wrap-toggle{flex-wrap:wrap}
/* cases and comparison */
.cases-head{display:flex;flex-wrap:wrap;justify-content:space-between;align-items:baseline;gap:6px 14px}
.cases-head h2,#compare h2{font-family:var(--f-display);font-weight:600;font-size:19px;margin:0}
#compare h3{font-family:var(--f-display);font-weight:600;font-size:15px;margin:10px 0 2px}
.hint.warn{color:var(--s2)}
#case-current{margin:0;color:var(--ink-2)}
.case-form{display:grid;grid-template-columns:minmax(0,1.3fr) minmax(0,1fr);gap:8px 12px;align-items:end}
.case-form label{display:grid;gap:3px;font-size:13px;color:var(--ink-2)}
.case-form .row{grid-column:1/-1}
@media (max-width:640px){.case-form{grid-template-columns:1fr}}
input[type=text]{font:inherit;font-size:14px;color:var(--ink);background:var(--bg);border:1px solid var(--rule);border-radius:5px;padding:5px 7px;min-width:0;width:100%;box-sizing:border-box}
#case-openbox,#case-linkbox{flex-wrap:nowrap}
#case-openbox input,#case-linkbox input{flex:1}
.btn:disabled{opacity:.5;cursor:default}
.case-list{display:grid;grid-template-columns:repeat(auto-fill,minmax(250px,1fr));gap:10px}
.case-card{border:1px solid var(--rule);border-left:4px solid var(--rule);border-radius:6px;padding:9px 10px;display:grid;gap:4px;background:var(--bg);align-content:start}
.case-card.on{border-left-color:var(--case)}
.case-card.loaded{box-shadow:inset 0 0 0 1px var(--accent)}
.case-top{display:grid;grid-template-columns:auto minmax(0,1fr) auto;gap:8px;align-items:start;cursor:pointer}
.case-top input{margin:3px 0 0;accent-color:var(--accent)}
.case-name{font-weight:600;font-size:14px;overflow-wrap:anywhere}
.case-cop{font-family:var(--f-data);font-size:13px;white-space:nowrap}
.case-sub,.case-note{font-size:12.5px;color:var(--ink-2);overflow-wrap:anywhere}
.case-meta{font-size:12px;color:var(--ink-3)}
.case-act{display:flex;flex-wrap:wrap;gap:6px;margin-top:2px}
.btn.sm{font-size:12px;padding:3px 8px}
.btn.danger{color:var(--s2);border-color:color-mix(in oklab,var(--s2) 45%,var(--rule))}
.cmp-wrap{overflow-x:auto}
table.cmp{min-width:0}
table.cmp th.ch,table.cmp td:nth-child(n+2){text-align:right}
table.cmp th.ch{vertical-align:bottom;min-width:120px;font-weight:600}
table.cmp td:nth-child(n+2){font-family:var(--f-data)}
table.cmp td small{color:var(--ink-3);font-size:11.5px;white-space:nowrap}
table.cmp td.best{background:color-mix(in oklab,var(--accent) 26%,transparent)}
table.cmp td.best span{font-weight:700}
table.cmp i.sw{display:inline-block;width:9px;height:9px;border-radius:2px;margin-right:6px}
.base-tag{display:inline-block;font-size:11px;font-weight:600;padding:1px 7px;margin-top:4px;border-radius:999px;background:var(--accent);color:var(--surface)}
.side-head .base-tag{grid-column:2/-1;justify-self:start;margin-top:0}
/* the columns share row tracks (subgrid), so every section lines up across cases whatever the name length or tags */
.side{display:grid;grid-auto-flow:column;grid-auto-columns:minmax(230px,1fr);grid-template-rows:repeat(8,auto);column-gap:12px;row-gap:0;padding-bottom:4px}
.side-col{grid-row:span 8;display:grid;grid-template-rows:subgrid;row-gap:8px;align-content:start;border:1px solid var(--rule);border-top:4px solid var(--case);border-radius:6px;padding:10px;background:var(--bg);min-width:0}
@supports not (grid-template-rows:subgrid){.side{grid-template-rows:none}.side-col{grid-row:auto;grid-template-rows:none}}
.side-head{display:grid;grid-template-columns:auto minmax(0,1fr) auto auto;gap:6px;align-items:start}
.grip{cursor:grab;touch-action:none;user-select:none;-webkit-user-select:none;color:var(--ink-3);font-size:16px;line-height:1.1;padding:0 2px}
.grip:hover{color:var(--ink)}
.icon.mv{width:22px;height:22px;font-size:14px}
.icon.mv:disabled{opacity:.35;cursor:default}
.side-col.dragging{opacity:.9;box-shadow:0 8px 22px rgba(0,0,0,.22);cursor:grabbing;position:relative;z-index:2}
.side-col h4{margin:0;font-size:14px;font-weight:600;overflow-wrap:anywhere}
.side-col .sub{font-size:12px;color:var(--ink-2);overflow-wrap:anywhere}
.side-cop{display:grid;justify-items:center;gap:4px}
.side-cop .v{font-family:var(--f-display);font-weight:600;font-size:40px;line-height:1;font-variant-numeric:tabular-nums}
.side-cop .stars svg{width:20px;height:20px}
.side-cop small{font-size:12px;color:var(--ink-3)}
.side-dials{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:2px 8px}
.side-dials .dial .dk{font-size:11.5px}
.side-kv{display:grid;grid-template-columns:auto auto;gap:3px 0;font-size:12.5px;margin:0}
.side-kv dt{padding-right:10px}
.side-kv dt{color:var(--ink-2)}
.side-kv dd{margin:0;text-align:right;font-family:var(--f-data);font-variant-numeric:tabular-nums}
.side-bar{display:flex;height:12px;border-radius:3px;overflow:hidden;background:var(--grid)}
.side-bar i{display:block;height:100%}
.side-barlab{font-size:11.5px;color:var(--ink-3)}
.side-col .mini svg{width:100%;height:auto}
.side-col .empty{font-size:13px;color:var(--ink-3)}
/* best value per metric: the same accent tint as the results table */
.side-kv dd.win,.win-chip{background:var(--accent);border-radius:3px;font-weight:700;color:var(--surface)}
.side-kv dd.win{padding:0 5px;margin-right:-5px;justify-self:end}
.side-kv .dlt{font-size:10.5px;font-weight:400;color:var(--ink-3);margin-left:2px}
.side-kv dd.win .dlt{color:inherit;opacity:.85}
.win-chip{padding:0 4px}
.side-cop.win .v{background:var(--accent);color:var(--surface);border-radius:6px;padding:2px 8px}
.side-dials .dial.win .dk{background:var(--accent);border-radius:3px;padding:0 5px;font-weight:700;color:var(--surface)}
.cmp-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:14px}
@media (max-width:860px){.cmp-grid{grid-template-columns:1fr}}
.pick:focus-visible,.fopt:focus-visible,.tbtn:focus-visible{outline:2px solid var(--accent);outline-offset:1px}
.icon{font:inherit;font-size:16px;line-height:1;cursor:pointer;border:1px solid var(--rule);background:var(--bg);color:var(--ink-2);border-radius:5px;width:28px;height:28px}
.row{display:flex;flex-wrap:wrap;gap:10px 14px;align-items:center}
.inline{display:inline-flex;gap:6px;align-items:center;font-size:13px;color:var(--ink-2)}
.inline select{width:auto}
.cmp-order{display:flex;flex-wrap:wrap;gap:6px 10px;align-items:center;margin:6px 0 4px}
.cmp-order .pick{width:auto;min-width:220px;max-width:100%}
.cmp-order .fluid-list{flex-basis:100%;box-sizing:border-box}
.hint{margin:0;color:var(--ink-3);font-size:13px}
.fields{display:grid;grid-template-columns:repeat(auto-fit,minmax(150px,1fr));gap:8px 12px}
.fields label{display:grid;gap:3px;font-size:12.5px;color:var(--ink-2)}
.fig.wide{grid-column:1/-1}
.fig .msg{font-size:14px;color:var(--ink)}
"""

MOCK = r"""
/* Mock CoolProp (test fluid) for layout checks only. */
var Module = (function () {
  const R = 8.314462618 / 58.12, cpv = 1.25, cpl = 2.40, L0 = 385.0, T0 = 273.15, P0 = 0.1035, rho_l = 580.0;
  const Psat = Tc => { const T = Tc + 273.15, L = L0 + (cpv - cpl) * (T - T0); return P0 * Math.exp((T * (L0 / T0 + (cpv - cpl) * Math.log(T / T0)) - L) / (R * T)); };
  const Tsat = P => SGHP.brent(T => Math.log(Psat(T) / P), -150, 250, 1e-12);
  const liq = (Tc, P) => { const T = Tc + 273.15; return { T: Tc, P, h: cpl * (T - T0), s: cpl * Math.log(T / T0), D: rho_l, q: -1 }; };
  const vap = (Tc, P) => { const T = Tc + 273.15; return { T: Tc, P, h: L0 + cpv * (T - T0), s: L0 / T0 + cpv * Math.log(T / T0) - R * Math.log(P / P0), D: P * 1000 / (R * T), q: -1 }; };
  const mix = (P, q) => { const Ts = Tsat(P), l = liq(Ts, P), v = vap(Ts, P); return { T: Ts, P, q, h: l.h + q * (v.h - l.h), s: l.s + q * (v.s - l.s), D: 1 / (1 / l.D + q * (1 / v.D - 1 / l.D)) }; };
  const from = (P, k, val) => { const l = mix(P, 0), v = mix(P, 1);
    if (val < l[k]) return liq(SGHP.brent(T => liq(T, P)[k] - val, -200, l.T + 1e-9, 1e-12), P);
    if (val > v[k]) return vap(SGHP.brent(T => vap(T, P)[k] - val, v.T - 1e-9, 600, 1e-12), P);
    return mix(P, (val - l[k]) / (v[k] - l[k])); };
  const ip = { PT_INPUTS: { value: 1 }, HmassP_INPUTS: { value: 2 }, PSmass_INPUTS: { value: 3 }, PQ_INPUTS: { value: 4 }, QT_INPUTS: { value: 5 } };
  class AS { update(p, a, b) { let s; const MPa = 1e6;
      if (p === ip.PT_INPUTS) { const P = a / MPa, T = b - 273.15; s = T > Tsat(P) ? vap(T, P) : liq(T, P); }
      else if (p === ip.HmassP_INPUTS) s = from(b / MPa, 'h', a / 1000);
      else if (p === ip.PSmass_INPUTS) s = from(a / MPa, 's', b / 1000);
      else if (p === ip.PQ_INPUTS) s = mix(a / MPa, b); else s = mix(Psat(b - 273.15), a);
      this.s = s; }
    T() { return this.s.T + 273.15; } p() { return this.s.P * 1e6; } hmass() { return this.s.h * 1000; } smass() { return this.s.s * 1000; }
    rhomass() { return this.s.D; } Q() { return this.s.q; } T_critical() { return 425.13; } p_critical() { return 3.796e6; }
    set_mass_fractions() {} set_mole_fractions() {} delete() {} }
  const M = { input_pairs: ip, factory: () => new AS(), get_global_param_string: () => "mock (test fluid)" };
  setTimeout(() => window.__cpReady(), 50);
  return M;
})();
"""


def build(out_dir: Path, coolprop_js: Path | None, mock: bool = False) -> Path:
    tpl = TEMPLATE.read_text()
    head = tpl[:tpl.index('</style>')] + EXTRA_CSS + '</style>\n'
    script = tpl[tpl.index('<script>') + len('<script>'):tpl.rindex('</script>')]
    consts = '\n'.join(line for line in script.splitlines()
                       if line.startswith(('const COMP ', 'const COMP_ORDER', 'const css', 'const NS', 'const el ',
                                           'const fmt', 'const key ')))
    shared = script[script.index('function niceTicks'):script.index('/* ---------- figures, table, page ---------- */')]
    shared = shared.replace('function copMap() {', 'function copMapGrid() {')
    shared += script[script.index('function figs(c)'):script.index('function chips(c)')]
    model = (HERE / 'sghp_model.js').read_text()
    live = (HERE / 'sghp_live.js').read_text() + '\n' + (HERE / 'sghp_cases.js').read_text()
    body = (HERE / 'body.html').read_text()
    if mock:
        loader = MOCK
    elif coolprop_js and coolprop_js.exists():
        loader = coolprop_js.read_text().replace('</script', '<\\/script')
        # expose the WebAssembly stack save/restore helpers so the page can undo the stack left behind
        # when a CoolProp C++ exception reaches JavaScript (see cpCall in sghp_live.js)
        anchor = 'Module["decrementExceptionRefcount"]=decrementExceptionRefcount;'
        if anchor in loader and 'var stackSave=' in loader:
            loader = loader.replace(anchor, anchor + 'Module["sghpStackSave"]=stackSave;Module["sghpStackRestore"]=stackRestore;', 1)
        else:
            print('warning: could not expose the CoolProp stack helpers; repeated CoolProp errors may abort the page')
    else:
        raise SystemExit(f"coolprop.js not found at {coolprop_js} (pass --coolprop-js, or --mock for a layout-only build)")
    # CoolProp 8 ships an ES-module Emscripten loader ("export default Module"): run it as an inline
    # module script that hands the factory to the page; older classic builds define a global Module.
    es_module = 'export default' in loader
    if es_module:
        loader_tag = ('<script type="module">\n' + re.sub(r'export\s+default\s+Module\s*;?', '', loader)
                      + '\nwindow.__cpFactoryReady(Module);\n</script>\n')
    else:
        loader_tag = '<script>\n' + loader + '\n</script>\n'
    preset = ("window.__cpReadyPromise = new Promise(r => { window.__cpReady = r; });\n"
              "window.__cpFactoryPromise = new Promise(r => { window.__cpFactoryReady = r; });\n"
              + ("" if es_module else "var Module = { locateFile: p => (p.endsWith('.wasm') ? 'coolprop.wasm' : p), "
                 "onRuntimeInitialized: () => window.__cpReady() };\n"))
    # a full document head: the charset keeps the degree signs and Greek letters intact on plain static
    # servers (GitHub Pages, python -m http.server), and the doctype keeps the page out of quirks mode
    skeleton = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
                '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n')
    html = (skeleton + head + body + '\n<script>\n' + model + '\n</script>\n<script>\n' + preset + '</script>\n'
            + loader_tag + '<script>\n' + consts + '\n' + live + '\n' + shared
            + '\nfunction copMap() { if (copMapLive()) copMapGrid(); }\nboot();\ncasesBoot();\n</script>\n')
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / 'index.html'
    path.write_text(html)
    wasm = (coolprop_js.parent / 'coolprop.wasm') if coolprop_js else None
    if wasm and wasm.exists() and not mock:
        shutil.copy(wasm, out_dir / 'coolprop.wasm')
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--coolprop-js', type=Path, default=CP_DIR / 'coolprop.js')
    ap.add_argument('--out', type=Path, default=HERE.parent.parent / 'steam')
    ap.add_argument('--mock', action='store_true', help='test fluid instead of CoolProp (layout checks)')
    a = ap.parse_args()
    print('->', build(a.out, a.coolprop_js, a.mock))


if __name__ == '__main__':
    main()
