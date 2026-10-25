/* Advanced dryer heat pumps: the REFPROP study (exact results) and the shared chart, cycle-diagram and
   process-flow renderers. The study data are fetched from the files published with the page. */
let D = null, CY = null, S5 = {};
let TEMPS = [150, 160, 170, 180, 190, 200], WS = [], WSA = [], GROUPS = [], OPTS = [];
const SHORT = {Propane: 'C3', Isobutane: 'iC4', Butane: 'nC4', Isopentane: 'iC5', Cyclopentane: 'cC5'};
const WLAB = new Proxy({}, {get: (_, w) => {
  const [mix, comp] = String(w).includes(':') ? String(w).split(':') : [null, String(w)];
  const frac = comp.split('-').map(v => +(parseFloat(v) * 100).toFixed(1)).join('/');
  return mix && mix !== 'Propane_Isopentane' ? `${frac} ${mix.split('_').map(n => SHORT[n] || n).join('/')}` : frac;
}});
/* cycle options: RHP-A and SLTHP with selectable composition; ternary RHP-A blends as their own options */
const OW = o => o.group ? (o.group.includes(state['w' + o.id]) ? state['w' + o.id] : o.group[0]) : (o.src === 'A' ? state.wA : state.wS);
function initData(d) {
  D = d; TEMPS = D.T; WS = D.W;
  /* RHP-A binary list: the shared Propane/Isopentane grid plus RHP-A-only higher-propane runs */
  WSA = WS.concat((D.WX || []).filter(k => k.startsWith('Propane_Isopentane:'))).sort((a, b) => parseFloat(a.split(':').pop()) - parseFloat(b.split(':').pop()));
  /* ternary blends are grouped by mixture: one card per mixture, with a composition switch when it has several */
  GROUPS = [];
  for (const k of (D.WX || [])) { const mix = k.split(':')[0]; if (mix === 'Propane_Isopentane') continue; let g = GROUPS.find(g => g.mix === mix); if (!g) GROUPS.push(g = {mix, keys: []}); g.keys.push(k); }
  OPTS = [{id: 'A', src: 'A', name: 'RHP-A', color: 'var(--sA)'}]
    .concat(GROUPS.map((g, i) => ({id: 'X' + (i + 1), src: 'A', group: g.keys, name: 'RHP-A ' + g.mix.split('_').map(n => SHORT[n] || n).join('/'), color: `var(--sX${Math.min(i + 1, 4)})`})))
    .concat([{id: 'S', src: 'S', name: 'SLTHP', color: 'var(--sS)'}]);
}
const PMAX = 70;  // slider top = no limit
/* the swept volume and volumetric efficiency scale the duties; with an air source the source air's outlet depends on
   the duty, so the live design re-solves */
function scaleChanged() { renderPFD(); renderSizeBars(); if (typeof LV !== 'undefined' && LV.src === 'air' && typeof changed === 'function') changed(); }
const vswTxt = () => (+state.vsw).toLocaleString('en-NZ', {maximumFractionDigits: state.vsw < 100 ? 1 : 0});
const state = {sizeBy: 'vsw', duty: 100, flow: 2, mode: 'live', cap: 27, a: '20', wA: '0.200-0.800', wS: '0.275-0.725', tol: '1', td: 200, cy: 'A', cyA: 'A', vsw: 145, nvMode: 'model', nvFix: 0.85};
const W = cyc => cyc === 'A' ? state.wA : state.wS;
try { Object.assign(state, JSON.parse(localStorage.getItem('dryerdash2') || '{}')); } catch (e) {}
state.mode = 'live';   // the page always opens on the live design
const save = () => { try { localStorage.setItem('dryerdash2', JSON.stringify(state)); } catch (e) {} };
const capv = c => c >= PMAX ? Infinity : c;
const key = (w, a, T) => `${w}|${a}|${T}`;

/* candidate designs per case: RHP-A grid records; SLTHP optimum per lift plus its pressure scan */
const CAND = {};
function cands(cyc, w, a, T) {
  const k = cyc + key(w, a, T);
  if (CAND[k]) return CAND[k];
  const out = [];
  if (cyc === 'A') for (const r of D.A[key(w, a, T)] || [])
    out.push({COP: r[2], VHC: r[3], P: r[4], Td: r[5], T1: r[6], lift: r[0], split: r[1], PL: r[7] ? r[4] / r[7] : null});
  else for (const r of D.S[key(w, a, T)] || []) {
    if (r[1] != null) out.push({COP: r[1], VHC: r[2], P: r[3], Td: r[4], T1: r[5], alpha: r[6], LMTD: r[7], lift: r[0], opt: true, PL: r[9] ?? null});
    for (const sc of r[8]) if (sc[1] != null) out.push({COP: sc[1], VHC: sc[2], P: sc[0], Td: sc[3], lift: r[0], opt: false, PL: r[9] ?? null, plApprox: true});
  }
  return CAND[k] = out;
}
/* best design under a limit; with a tolerance, the lowest-pressure design within tol of the best COP */
function best(cyc, w, a, T, cap) {
  const ok = cands(cyc, w, a, T).filter(c => c.P <= cap + 1e-9);
  if (!ok.length) return null;
  const top = Math.max(...ok.map(c => c.COP)), tol = +state.tol / 100;
  let b = null;
  for (const c of ok) if (c.COP >= top * (1 - tol) - 1e-12 && (!b || c.P < b.P - 1e-9 || (Math.abs(c.P - b.P) <= 1e-9 && c.COP > b.COP))) b = c;
  return b;
}
const bestA = (w, a, T, cap) => best('A', w, a, T, cap), bestS = (w, a, T, cap) => best('S', w, a, T, cap);
const BEST = {A: bestA, S: bestS};
/* best under the limit, else the unlimited optimum flagged as over the limit */
function bestOr(cyc, w, a, T, cap) {
  const b = BEST[cyc](w, a, T, cap);
  if (b) return b;
  const u = BEST[cyc](w, a, T, Infinity);
  return u ? {...u, over: true} : null;
}
function headline(cyc, w, a, cap) {
  const per = TEMPS.map(T => bestOr(cyc, w, a, T, cap));
  const n = per.filter(p => p && !p.over).length;
  const full = per.every(Boolean), over = per.some(p => p && p.over);
  return {per, n, mean: full ? per.reduce((s, p) => s + p.COP, 0) / TEMPS.length : null,
          vhc: full ? per.reduce((s, p) => s + p.VHC, 0) / TEMPS.length : null,
          pmax: full ? Math.max(...per.map(p => p.P)) : null, over,
          pl: full && per.every(p => p.PL) ? Math.min(...per.map(p => p.PL)) : null,
          opr: full && per.every(p => p.PL) ? Math.max(...per.map(p => p.P / p.PL)) : null,
          plApprox: per.some(p => p && p.plApprox)};
}
const pct = (x, y) => (x / y - 1) * 100;
const chip = (x, y) => { if (x == null || y == null) return '<span class="dchip muted">n/a</span>'; const p = pct(x, y); return `<span class="dchip ${p < 0 ? 'neg' : 'pos'}">${p >= 0 ? '+' : '−'}${Math.abs(p).toFixed(1)} %</span>`; };
const capTxt = c => c >= PMAX ? 'no limit' : c.toFixed(1) + ' bar';

/* ---------- controls ---------- */
const slider = document.getElementById('pcap'), pval = document.getElementById('pval');
slider.addEventListener('input', () => { state.cap = +slider.value; update(); });
const pnum = document.getElementById('pnum');
pnum.addEventListener('change', () => { state.cap = Math.min(Math.max(+pnum.value || 27, 18), PMAX); update(); });
const vsw = document.getElementById('vsw');
const VMIN = 20, VMAX = 15600, vswn = document.getElementById('vswn');
const toV = x => { const v = VMIN * Math.pow(VMAX / VMIN, x / 1000); return v < 200 ? Math.round(v) : v < 2000 ? Math.round(v / 5) * 5 : Math.round(v / 50) * 50; };
const toX = v => Math.round(1000 * Math.log(Math.min(Math.max(v, VMIN), VMAX) / VMIN) / Math.log(VMAX / VMIN));
/* size by the low-stage swept volume, by a target sink duty (heat to air) or by a target sink flow (dryer air), with
   the swept volume solved for it (the duties scale with the swept volume, so it is a ratio; a sink flow is the duty
   it takes to heat that air from the inlet to the dryer temperature; with an air source the design re-solves) */
function renderSize() {
  const live = state.mode === 'live', by = live ? (state.sizeBy || 'vsw') : 'vsw';
  const sc = document.getElementById('sizeCtl'); if (!sc) return;
  sc.hidden = !live;
  // the buttons are made once and only restyled, so a re-render (say from the duty box losing focus) mid-click
  // does not swallow the click
  const box = document.getElementById('segSize');
  if (!box.children.length) for (const [v, lab] of [['vsw', 'Swept volume'], ['duty', 'Sink duty'], ['flow', 'Sink flow']]) {
    const b = document.createElement('button'); b.type = 'button'; b.textContent = lab; b.dataset.v = v;
    b.onclick = () => { if (state.sizeBy === v) return; const cur = typeof sizeValues === 'function' ? sizeValues() : {}; if (v === 'duty' && cur.duty > 0) state.duty = cur.duty; if (v === 'flow' && cur.flow > 0) state.flow = cur.flow; state.sizeBy = v; save(); renderSize(); sizeChanged(); }; box.appendChild(b);
  }
  for (const b of box.children) b.setAttribute('aria-pressed', String((state.sizeBy || 'vsw') === b.dataset.v));
  renderSizeBars();
  const hint = document.getElementById('sizeHint');
  hint.hidden = by === 'vsw';
  hint.textContent = by === 'duty' ? 'Set the total sink heating duty under Operating point.' : by === 'flow' ? 'Set the dryer inlet air mass flow under Operating point.' : '';
  vsw.disabled = by !== 'vsw'; vswn.disabled = by !== 'vsw';
  vsw.value = toX(state.vsw); vswn.value = +(+state.vsw).toFixed(state.vsw < 100 ? 1 : 0);
}
function sizeChanged() {
  if (typeof applyDuty === 'function') applyDuty();
  if (typeof LV !== 'undefined' && LV.src === 'air' && typeof changed === 'function') { changed(); return; }
  if (typeof renderLive === 'function' && state.mode === 'live') renderLive();
  renderPFD();
  if (typeof casesChanged === 'function') casesChanged();
}
/* the sink duty and the dryer air flow as bars (log scales) under Operating point: the one the compressor is sized by
   is set there; the others are derived from the solved design (at the swept volume, or the other target) */
const logScale = (lo, hi, labels, minors) => ({
  to: v => Math.min(Math.max(Math.log(Math.max(+v, lo) / lo) / Math.log(hi / lo), 0), 1),
  from: u => { const v = lo * Math.pow(hi / lo, Math.min(Math.max(+u, 0), 1)); return +v.toPrecision(v < 1 ? 2 : 3); },
  labels, minors,
});
const SIZE_BARS = {
  duty: {id: 'lvDuty', unit: 'kW', lo: 1, hi: 10000, sc: logScale(1, 10000, [1, 10, 100, 1000, 10000], [2, 5, 20, 50, 200, 500, 2000, 5000])},
  flow: {id: 'lvFlow', unit: 'kg/s', lo: 0.05, hi: 200, sc: logScale(0.05, 200, [0.1, 1, 10, 100], [0.2, 0.5, 2, 5, 20, 50])},
};
const sizeFmt = v => v < 1 ? v.toFixed(3) : v < 10 ? v.toFixed(2) : v < 100 ? v.toFixed(1) : v.toFixed(0);
function sizeValues() {
  const out = {duty: null, flow: null};
  if (typeof LIVE !== 'undefined' && LIVE.res && LIVE.res.success && LIVE.cfg && typeof kwAt === 'function') {
    const kw = kwAt(LIVE.res, LIVE.cfg), dh = typeof airDh === 'function' ? airDh(LIVE.cfg) : null;
    if (kw && kw.air > 0) { out.duty = kw.air; if (dh) out.flow = kw.air / dh; }
  }
  if (state.sizeBy === 'duty') out.duty = state.duty;
  if (state.sizeBy === 'flow') out.flow = state.flow;
  return out;
}
function sizeTicks(id, sc, cur) {
  const box = document.getElementById(`${id}-ticks`); box.innerHTML = '';
  const put = (v, lab) => {
    const left = `calc(8px + (100% - 16px) * ${sc.to(v)})`, i = document.createElement('i'); i.style.left = left; box.append(i);
    if (lab) { const t = document.createElement('span'); t.style.left = left; t.textContent = v.toLocaleString('en-NZ'); if (cur != null && Math.abs(sc.to(v) - sc.to(cur)) < 1e-3) t.className = 'on'; box.append(t); }
  };
  for (const v of sc.minors) put(v, false);
  for (const v of sc.labels) put(v, true);
}
function renderSizeBars() {
  const by = state.sizeBy || 'vsw', vals = sizeValues(), from = {vsw: 'the swept volume', duty: 'the sink duty', flow: 'the air flow'}[by];
  for (const k of ['duty', 'flow']) {
    const b = SIZE_BARS[k], inp = document.getElementById(b.id), out = document.getElementById(`${b.id}-out`);
    if (!inp) continue;
    const active = by === k, v = vals[k];
    inp.disabled = !active;
    if (document.activeElement !== inp && v != null) inp.value = b.sc.to(v);
    out.textContent = v == null ? '–' : `${sizeFmt(v)} ${b.unit}${active ? '' : ' · derived'}`;
    out.classList.toggle('derived', !active); out.classList.toggle('editable', active);
    out.tabIndex = active ? 0 : -1; out.title = active ? 'Click to type a value' : `Derived from ${from}; choose this under Compressor, Size by, to set it`;
    sizeTicks(b.id, b.sc, v);
  }
}
for (const k of ['duty', 'flow']) {
  const b = SIZE_BARS[k], inp = document.getElementById(b.id), out = document.getElementById(`${b.id}-out`);
  const apply = v => { state[k] = Math.min(Math.max(v, b.lo), b.hi); save(); sizeChanged(); renderSize(); };
  inp.addEventListener('input', () => { out.textContent = `${sizeFmt(b.sc.from(inp.value))} ${b.unit}`; });
  inp.addEventListener('change', () => apply(b.sc.from(inp.value)));
  // the value can be typed while this bar sets the size: click it (or focus it and press Enter), type, Enter to apply
  const edit = () => {
    if (state.sizeBy !== k || out.hidden) return;
    const box = document.createElement('input');
    box.type = 'number'; box.step = 'any'; box.className = 'outedit'; box.value = +(+state[k]).toPrecision(4); box.title = `${b.lo} to ${b.hi}`;
    box.setAttribute('aria-label', `${out.previousElementSibling.textContent} (${b.lo} to ${b.hi} ${b.unit})`);
    out.hidden = true; out.after(box); box.focus(); box.select();
    let done = false;
    const finish = ok => {
      if (done) return; done = true;
      const v = +box.value; box.remove(); out.hidden = false;
      if (ok && box.value !== '' && isFinite(v) && v > 0) apply(v); else renderSizeBars();
      out.focus();
    };
    box.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); finish(true); } else if (e.key === 'Escape') { e.preventDefault(); finish(false); } });
    box.addEventListener('blur', () => finish(true));
  };
  out.addEventListener('click', edit);
  out.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); edit(); } });
}
vsw.addEventListener('input', () => { state.vsw = toV(+vsw.value); vswn.value = state.vsw; save(); scaleChanged(); });
vswn.addEventListener('change', () => { const v = Math.min(Math.max(+vswn.value || 145, VMIN), VMAX); state.vsw = v; vswn.value = v; vsw.value = toX(v); save(); scaleChanged(); });
const nvs = document.getElementById('nvs'), nvv = document.getElementById('nvv');
nvs.addEventListener('input', () => { state.nvFix = +nvs.value; state.nvMode = 'fixed'; renderNv(); save(); scaleChanged(); });
nvv.addEventListener('change', () => { state.nvFix = Math.min(Math.max(+nvv.value || 0.85, 0.3), 1); state.nvMode = 'fixed'; renderNv(); save(); scaleChanged(); });
function renderNv() {
  nvs.value = state.nvFix; nvv.value = state.nvFix.toFixed(2);
  nvs.style.opacity = state.nvMode === 'fixed' ? 1 : .45; nvv.style.opacity = state.nvMode === 'fixed' ? 1 : .45;
  const el = document.getElementById('segNv'); el.innerHTML = '';
  for (const [v, lab] of [['model', 'Model'], ['fixed', 'Fixed']]) {
    const b = document.createElement('button'); b.type = 'button'; b.textContent = lab; b.setAttribute('aria-pressed', String(state.nvMode === v));
    b.onclick = () => { state.nvMode = v; renderNv(); save(); scaleChanged(); }; el.appendChild(b);
  }
}
function seg(id, items, k) {
  const el = document.getElementById(id); el.innerHTML = '';
  for (const [v, lab] of items) {
    const b = document.createElement('button'); b.type = 'button'; b.textContent = lab;
    b.setAttribute('aria-pressed', String(state[k] === v)); b.onclick = () => { state[k] = v; update(); };
    el.appendChild(b);
  }
}

/* ---------- headline: one result per cycle at its own composition ---------- */
function renderHeadline(cap) {
  const H = {};
  for (const o of OPTS) H[o.id] = headline(o.src, OW(o), state.a, cap);
  document.getElementById('hlHead').textContent = `Headline · ${state.a} K approach · ${capTxt(state.cap)} · mean of 150–200 °C air`;
  const ref = H.A;
  document.getElementById('hlCards').innerHTML = OPTS.map(o => {
    const h = H[o.id];
    const sel = o.group ? (o.group.length > 1 ? `<div class="ctl"><span class="lab">${o.name.split(' ')[1]}, mass %</span><div class="seg" id="segW${o.id}"></div></div>` : `<div class="vs">${WLAB[o.group[0]]} (mass %)</div>`)
                        : `<div class="ctl"><span class="lab">Propane/Isopentane, mass %</span><div class="seg" id="segW${o.src}"></div></div>`;
    const vs = o.id === 'A' ? '<div class="vs">reference</div>' : `<div class="vs">vs RHP-A · COP ${chip(h.mean, ref.mean)} · VHC ${chip(h.vhc, ref.vhc)}</div>`;
    return `<div class="panel hcard"><div class="big"><span class="k"><i class="sw" style="background:${o.color}"></i>${o.name}</span>
      <span class="v ${h.over ? 'over' : ''}">${h.mean == null ? '—' : h.mean.toFixed(3)}</span>
      <span class="s ${h.over ? 'over' : ''}">${h.mean == null ? 'no feasible design at every temperature' : h.over ? `over the limit · only ${h.n} of 6 temperatures fit` : `COP, mean of 6 · ${Math.min(...h.per.map(p => p.COP)).toFixed(2)}–${Math.max(...h.per.map(p => p.COP)).toFixed(2)}`}</span>
      <div class="kv ${h.over ? 'over' : ''}"><span>VHC mean</span><b>${h.vhc == null ? '—' : h.vhc.toFixed(2)}</b><span>MJ/m³</span>
        <span>Max pressure</span><b>${h.pmax == null ? '—' : h.pmax.toFixed(1)}</b><span>bar</span>
        <span>LS suction pressure</span><b>${h.pl == null ? '—' : (h.plApprox ? '≈' : '') + h.pl.toFixed(2)}</b><span>bar, min</span>
        <span>Overall PR</span><b>${h.opr == null ? '—' : (h.plApprox ? '≈' : '') + h.opr.toFixed(1)}</b><span>max</span></div></div>${vs}${sel}</div>`;
  }).join('');
  seg('segWA', WSA.map(w => [w, WLAB[w]]), 'wA');
  for (const o of OPTS) if (o.group && o.group.length > 1) { if (!o.group.includes(state['w' + o.id])) state['w' + o.id] = o.group[0]; seg('segW' + o.id, o.group.map(k => [k, WLAB[k].split(' ')[0]]), 'w' + o.id); }
  seg('segWS', WS.map(w => [w, WLAB[w]]), 'wS');
}

/* ---------- chart engine ---------- */
const NS = 'http://www.w3.org/2000/svg', tip = document.getElementById('tip');
function el(tag, at, p) { const e = document.createElementNS(NS, tag); for (const k in at) e.setAttribute(k, at[k]); if (p) p.appendChild(e); return e; }
function ticks(lo, hi, n = 5) {
  const span = hi - lo || 1, raw = span / n, mag = Math.pow(10, Math.floor(Math.log10(raw)));
  const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= n) || 10 * mag, t = [];
  for (let v = Math.floor(lo / step) * step; v <= Math.ceil(hi / step) * step + step * 1e-6; v += step) t.push(+v.toFixed(10));
  return t;
}
/* series: [{name,color,dash,pts:[{x,y,hollow,note}]}]; opt: {xs:[lo,hi], xt:[ticks], xlabel, dec, vline, xfmt} */
function chart(id, series, opt) {
  const svg = document.getElementById(id); svg.innerHTML = '';
  const W = Math.max(280, Math.round((svg.parentElement.clientWidth || 560) - 36)), H = Math.round(Math.min(opt.h || 250, Math.max(190, W * 0.55)));
  const m = {l: 46, r: 14, t: 10, b: 34};
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const ys = series.flatMap(s => s.pts.map(p => p.y)).filter(v => v != null);
  if (!ys.length) { el('text', {x: W / 2, y: H / 2, 'text-anchor': 'middle', class: 'lbl'}, svg).textContent = 'no feasible designs under this limit'; return; }
  let lo = Math.min(...ys), hi = Math.max(...ys); const pad = (hi - lo) * .1 || Math.abs(hi) * .02 || .1; lo -= pad; hi += pad;
  const yt = ticks(lo, hi, 5); lo = yt[0]; hi = yt[yt.length - 1];
  const [x0, x1] = opt.xs, X = v => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r), Y = v => H - m.b - (v - lo) / (hi - lo) * (H - m.t - m.b);
  const g = el('g', {class: 'axis'}, svg);
  for (const t of yt) { el('line', {x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), class: 'gridline'}, g); el('text', {x: m.l - 6, y: Y(t) + 3.5, 'text-anchor': 'end'}, g).textContent = opt.yfmt ? opt.yfmt(t) : +t.toFixed(3); }
  for (const t of opt.xt) el('text', {x: X(t), y: H - m.b + 14, 'text-anchor': 'middle'}, g).textContent = opt.xfmt ? opt.xfmt(t) : t;
  el('line', {x1: m.l, x2: W - m.r, y1: H - m.b, y2: H - m.b, class: 'base'}, g);
  el('text', {x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', class: 'title'}, g).textContent = opt.xlabel;
  if (opt.vline != null) { const vx = X(Math.min(opt.vline, x1)); el('line', {x1: vx, x2: vx, y1: m.t, y2: H - m.b, stroke: 'var(--accent)', 'stroke-width': 1.2, 'stroke-dasharray': '4 3'}, svg); }
  if (opt.hline != null && opt.hline >= lo && opt.hline <= hi) { el('line', {x1: m.l, x2: W - m.r, y1: Y(opt.hline), y2: Y(opt.hline), stroke: 'var(--accent)', 'stroke-width': 1.2, 'stroke-dasharray': '4 3'}, svg); el('text', {x: W - m.r, y: Y(opt.hline) - 4, 'text-anchor': 'end', class: 'lbl'}, svg).textContent = opt.hlabel || ''; }
  for (const s of series) {
    let d = '', pen = false;
    for (const p of s.pts) { if (p.y == null) { pen = false; continue; } d += (pen ? 'L' : 'M') + X(p.x).toFixed(1) + ' ' + Y(p.y).toFixed(1); pen = true; }
    el('path', {d, fill: 'none', stroke: s.color, 'stroke-width': s.width || 2, 'stroke-dasharray': s.dash || '', 'stroke-linejoin': 'round', 'stroke-linecap': 'round', opacity: s.opacity || 1}, svg);
    if (s.markers !== false) for (const p of s.pts) if (p.y != null)
      el('circle', {cx: X(p.x), cy: Y(p.y), r: p.hollow ? 4.6 : 3.8, fill: p.hollow ? 'var(--panel)' : s.color, stroke: p.hollow ? 'var(--bad)' : 'var(--panel)', 'stroke-width': p.hollow ? 2 : 1.5}, svg);
  }
  const hit = el('rect', {x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent'}, svg);
  const cross = el('line', {y1: m.t, y2: H - m.b, stroke: 'var(--ink-3)', 'stroke-width': 1, 'stroke-dasharray': '2 3', visibility: 'hidden'}, svg);
  const xsAll = [...new Set(series.flatMap(s => s.pts.filter(p => p.y != null).map(p => p.x)))].sort((a, b) => a - b);
  const show = ev => {
    const r = svg.getBoundingClientRect(), vx = (ev.clientX - r.left) / r.width * W;
    let bx = null, bd = 1e9; for (const x of xsAll) { const dd = Math.abs(X(x) - vx); if (dd < bd) { bd = dd; bx = x; } }
    if (bx == null) return;
    cross.setAttribute('x1', X(bx)); cross.setAttribute('x2', X(bx)); cross.setAttribute('visibility', 'visible');
    let h = `<div class="tt">${opt.ttx(bx)}</div>`;
    for (const s of series) { if (s.tip === false) continue; const p = s.pts.find(q => q.x === bx); h += `<div class="tr"><i class="sw" style="background:${s.color}"></i><span>${s.name}</span><b>${p && p.y != null ? (p.hollow ? '(' : '') + p.y.toFixed(opt.dec) + (p.hollow ? ')' : '') : '—'}</b></div>`; if (p && p.note) h += `<div class="note">${p.note}</div>`; }
    tip.innerHTML = h; tip.hidden = false;
    tip.style.left = Math.min(ev.clientX + 14, innerWidth - tip.offsetWidth - 8) + 'px'; tip.style.top = Math.min(ev.clientY + 14, innerHeight - tip.offsetHeight - 8) + 'px';
  };
  hit.addEventListener('pointermove', show); hit.addEventListener('pointerdown', show);
  hit.addEventListener('pointerleave', () => { tip.hidden = true; cross.setAttribute('visibility', 'hidden'); });
}

/* ---------- per-temperature charts ---------- */
function perT(cyc, f, w = W(cyc)) {
  const cap = capv(state.cap);
  return TEMPS.map(T => {
    const b = bestOr(cyc, w, state.a, T, cap);
    if (!b) return {x: T, y: null, note: 'no design'};
    return {x: T, y: f(b), hollow: !!b.over, note: `lift ${b.lift} K · ${b.P.toFixed(1)} bar${b.over ? ' · over the limit' : cyc === 'S' && !b.opt ? ' · pressure scan point' : ''}`};
  });
}
function renderPerT() {
  const mk = f => OPTS.map(o => ({name: o.name, color: o.color, pts: perT(o.src, f, OW(o))}));
  const o = {xs: [145, 205], xt: TEMPS, xlabel: 'Dryer inlet air target (°C)', ttx: x => `${x} °C dryer inlet air`};
  const lg = document.getElementById('lgPerT');
  if (lg && !lg.dataset.done) { lg.dataset.done = 1; lg.insertAdjacentHTML('afterbegin', OPTS.filter(o => o.group).map(o => `<span><i class="lg-line" style="border-color:${o.color}"></i>${o.name}</span>`).join('')); }
  chart('cCOP', mk(b => b.COP), {...o, dec: 3});
  chart('cP', mk(b => b.P), {...o, dec: 1, hline: state.cap < PMAX ? state.cap : null, hlabel: 'limit'});
  chart('cTd', mk(b => b.Td), {...o, dec: 0});
  chart('cV', mk(b => b.VHC), {...o, dec: 2});
}

/* ---------- detail table ---------- */
function renderTable() {
  const cap = capv(state.cap);
  document.getElementById('tblTitle').textContent = `Detail · RHP-A ${WLAB[state.wA]} · SLTHP ${WLAB[state.wS]} · ${state.a} K approach · ${capTxt(state.cap)}`;
  let h = `<thead><tr><th></th><th colspan="7" class="gA">RHP-A</th><th colspan="8" class="gS">SLTHP</th></tr>
    <tr><th>Air</th><th>COP</th><th>marg.</th><th>lift/split</th><th>Pmax</th><th>Tdis</th><th>stage-1 out</th><th>VHC</th>
    <th>COP</th><th>Δ vs A</th><th>marg.</th><th>lift</th><th>Pmax</th><th>Tdis</th><th>RHX1 share</th><th>LMTD</th></tr></thead><tbody>`;
  const marg = (arr, i) => { if (i === 0 || !arr[i] || !arr[i - 1]) return '—'; const q0 = TEMPS[i - 1] - 15, q1 = TEMPS[i] - 15; return ((q1 - q0) / (q1 / arr[i].COP - q0 / arr[i - 1].COP)).toFixed(2); };
  const A = TEMPS.map(T => bestOr('A', state.wA, state.a, T, cap)), S = TEMPS.map(T => bestOr('S', state.wS, state.a, T, cap));
  TEMPS.forEach((T, i) => {
    const a = A[i], s = S[i];
    h += `<tr><td>${T} °C</td>` + (a ? `<td class="${a.over ? 'over' : ''}">${a.COP.toFixed(3)}</td><td>${marg(A, i)}</td><td>${a.lift}/${a.split.toFixed(2)}</td><td class="${a.over ? 'over' : ''}">${a.P.toFixed(1)}</td><td>${a.Td.toFixed(0)}</td><td>${a.T1.toFixed(0)}</td><td>${a.VHC.toFixed(2)}</td>` : '<td>none</td><td></td><td></td><td></td><td></td><td></td><td></td>')
      + (s ? `<td class="${s.over ? 'over' : ''}">${s.COP.toFixed(3)}</td><td>${chip(s.COP, a?.COP)}</td><td>${marg(S, i)}</td><td>${s.lift}</td><td class="${s.over ? 'over' : ''}">${s.P.toFixed(1)}</td><td>${s.Td.toFixed(0)}</td><td>${s.opt ? (s.alpha * 100).toFixed(0) + ' %' : '—'}</td><td>${s.opt && s.LMTD != null ? s.LMTD.toFixed(1) + ' K' : '—'}</td>` : '<td>none</td><td></td><td></td><td></td><td></td><td></td><td></td><td></td>') + '</tr>';
  });
  const hA = headline('A', state.wA, state.a, cap), hS = headline('S', state.wS, state.a, cap);
  h += `<tr><td><b>mean / max</b></td><td class="${hA.over ? 'over' : ''}"><b>${hA.mean?.toFixed(3) ?? '—'}</b></td><td></td><td></td><td class="${hA.over ? 'over' : ''}"><b>${hA.pmax?.toFixed(1) ?? '—'}</b></td><td></td><td></td><td class="${hA.over ? 'over' : ''}"><b>${hA.vhc?.toFixed(2) ?? '—'}</b></td><td class="${hS.over ? 'over' : ''}"><b>${hS.mean?.toFixed(3) ?? '—'}</b></td><td>${chip(hS.mean, hA.mean)}</td><td></td><td></td><td class="${hS.over ? 'over' : ''}"><b>${hS.pmax?.toFixed(1) ?? '—'}</b></td><td colspan="3"></td></tr>`;
  document.getElementById('tbl').innerHTML = h + '</tbody>';
}

/* ---------- cycle diagrams ---------- */
function xy(id, layers, opt) {
  const svg = document.getElementById(id); svg.innerHTML = '';
  const W = Math.max(260, Math.round((svg.parentElement.clientWidth || 400) - 36)), H = Math.round(Math.max(240, W * 0.85));
  const m = {l: 50, r: 10, t: 8, b: 36};
  svg.setAttribute('viewBox', `0 0 ${W} ${H}`);
  const all = layers.filter(L => L.fit !== false).flatMap(L => L.pts);
  if (!all.length) { el('text', {x: W / 2, y: H / 2, 'text-anchor': 'middle', class: 'lbl'}, svg).textContent = 'no diagram for this case'; return; }
  const fy = opt.ylog ? Math.log10 : v => v;
  let x0 = Math.min(...all.map(p => p[0])), x1 = Math.max(...all.map(p => p[0])), y0 = Math.min(...all.map(p => fy(p[1]))), y1 = Math.max(...all.map(p => fy(p[1])));
  const px = (x1 - x0) * .06, py = (y1 - y0) * .06; x0 -= px; x1 += px; y0 -= py; y1 += py;
  const xt = ticks(x0, x1, 5); x0 = xt[0]; x1 = xt[xt.length - 1];
  let yt;
  if (opt.ylog) { yt = [0.2, 0.5, 1, 2, 5, 10, 20, 50, 100].filter(v => Math.log10(v) >= y0 && Math.log10(v) <= y1); }
  else { yt = ticks(y0, y1, 6); y0 = yt[0]; y1 = yt[yt.length - 1]; }
  const X = v => m.l + (v - x0) / (x1 - x0) * (W - m.l - m.r), Y = v => H - m.b - (fy(v) - y0) / (y1 - y0) * (H - m.t - m.b);
  const g = el('g', {class: 'axis'}, svg);
  for (const t of yt) { el('line', {x1: m.l, x2: W - m.r, y1: Y(t), y2: Y(t), class: 'gridline'}, g); el('text', {x: m.l - 6, y: Y(t) + 3.5, 'text-anchor': 'end'}, g).textContent = t; }
  for (const t of xt) { el('line', {x1: X(t), x2: X(t), y1: m.t, y2: H - m.b, class: 'gridline'}, g); el('text', {x: X(t), y: H - m.b + 14, 'text-anchor': 'middle'}, g).textContent = +t.toFixed(3); }
  el('text', {x: (m.l + W - m.r) / 2, y: H - 4, 'text-anchor': 'middle', class: 'title'}, g).textContent = opt.xlabel;
  el('text', {x: 12, y: (m.t + H - m.b) / 2, 'text-anchor': 'middle', class: 'title', transform: `rotate(-90 12 ${(m.t + H - m.b) / 2})`}, g).textContent = opt.ylabel;
  const clip = 'c' + id; const cp = el('clipPath', {id: clip}, svg); el('rect', {x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b}, cp);
  const body = el('g', {'clip-path': `url(#${clip})`}, svg);
  for (const L of layers) {
    const d = L.pts.map((p, i) => (i ? 'L' : 'M') + X(p[0]).toFixed(1) + ' ' + Y(p[1]).toFixed(1)).join('');
    el('path', {d, fill: L.fill || 'none', stroke: L.color, 'stroke-width': L.width || 2, 'stroke-dasharray': L.dash || '', 'stroke-linejoin': 'round', 'stroke-linecap': 'round'}, body);
  }
  const pts = layers.filter(L => L.tip).flatMap(L => L.pts.map((p, i) => ({p, L, raw: L.raw[i]})));
  const hit = el('rect', {x: m.l, y: m.t, width: W - m.l - m.r, height: H - m.t - m.b, fill: 'transparent'}, svg);
  const dot = el('circle', {r: 4.5, fill: 'var(--panel)', 'stroke-width': 2, visibility: 'hidden'}, svg);
  const show = ev => {
    const r = svg.getBoundingClientRect(), vx = (ev.clientX - r.left) / r.width * W, vy = (ev.clientY - r.top) / r.height * H;
    let b = null, bd = 1e9; for (const q of pts) { const dd = (X(q.p[0]) - vx) ** 2 + (Y(q.p[1]) - vy) ** 2; if (dd < bd) { bd = dd; b = q; } }
    if (!b || bd > 900) { tip.hidden = true; dot.setAttribute('visibility', 'hidden'); return; }
    dot.setAttribute('cx', X(b.p[0])); dot.setAttribute('cy', Y(b.p[1])); dot.setAttribute('stroke', b.L.color); dot.setAttribute('visibility', 'visible');
    const [T, P, h, s_] = b.raw;
    tip.innerHTML = `<div class="tt">${b.L.name}</div><div class="tr"><i class="sw" style="background:${b.L.color}"></i><span>${b.L.comp}</span><b></b></div>
      <div class="note">T ${T.toFixed(1)} °C · P ${P.toFixed(2)} bar<br>h ${h.toFixed(1)} kJ/kg · s ${s_.toFixed(3)} kJ/kg·K</div>`;
    tip.hidden = false;
    tip.style.left = Math.min(ev.clientX + 14, innerWidth - tip.offsetWidth - 8) + 'px'; tip.style.top = Math.min(ev.clientY + 14, innerHeight - tip.offsetHeight - 8) + 'px';
  };
  hit.addEventListener('pointermove', show); hit.addEventListener('pointerdown', show);
  hit.addEventListener('pointerleave', () => { tip.hidden = true; dot.setAttribute('visibility', 'hidden'); });
}
/* the RHP-A fluid shown in the diagrams: one of the RHP-A headline options, at that card's composition */
const aOpt = () => OPTS.find(o => o.id === state.cyA && o.src === 'A') || OPTS[0];
const DW = k => k === 'A' ? OW(aOpt()) : state.wS;
/* component colours for the cycle diagrams */
const CCOL = {
  comp: ['Compressor', '#d62728'], gc: ['Gas cooler', '#ff7f0e'], rhx: ['RHX1 / recuperator', '#9467bd'],
  ihx: ['Internal heat exchanger', '#8c564b'], valve: ['Expansion valve', '#e377c2'], evap: ['Evaporator', '#1f77b4'],
  mix: ['Mixing', '#17becf'], turb: ['Turbine (expander)', '#e377c2'], casc: ['Cascade exchanger / steam generator', '#9467bd'], mvr: ['MVR / water injection', '#7f7f7f'], other: ['Other', '#bcbd22'],
};
const compGroup = n => /turbine|expander/i.test(n) ? 'turb' : /lthx/i.test(n) ? 'evap' : /preheater/i.test(n) ? 'rhx' : /cascade|steam generator/i.test(n) ? 'casc' : /mvr|injection/i.test(n) ? 'mvr' : /compressor/i.test(n) ? 'comp'
  : /gas cooler|condenser|subcooler/i.test(n) ? 'gc' : /rhx|recuperator/i.test(n) ? 'rhx'
  : /internal/i.test(n) ? 'ihx' : /valve|flash/i.test(n) ? 'valve' : /evaporator/i.test(n) ? 'evap' : /mix/i.test(n) ? 'mix' : 'other';
/* what the diagrams and the process flow show: the study case chosen below the diagrams, or the live design */
function cycleSrc() {
  if (state.mode === 'live') return typeof liveSrc === 'function' ? liveSrc() : null;
  const setName = state.cap >= PMAX ? 'none' : '27';
  const k = state.cy, w = DW(k), c = CY && CY.cases[`${w}|${state.a}|${state.td}`];
  const nm = k === 'A' ? aOpt().name : 'SLTHP', col = k === 'A' ? aOpt().color : 'var(--sS)';
  const d = c && c[`${k}|${setName}`];
  let info;
  if (!CY) info = `<div class="cyrow"><i class="sw" style="background:${col}"></i><b>${nm}</b><span class="f">loading the cycle diagrams…</span></div>`;
  else if (!d || d.error) info = `<div class="cyrow"><i class="sw" style="background:${col}"></i><b>${nm}</b><span class="f">no diagram${d && d.error ? ' (' + d.error.slice(0, 60) + ')' : ''}</span></div>`;
  else {
    const sel = bestOr(k, w, state.a, state.td, capv(state.cap));
    const same = sel && Math.abs(sel.COP - d.chosen_COP) < 5e-4 && Math.abs(sel.P - d.Pmax) < 0.05;
    const desc = k === 'A' ? `RHX1 lift ${d.lift} K · power split ${(+d.split).toFixed(2)} · Pmax ${d.Pmax.toFixed(1)} bar · COP ${d.COP.toFixed(3)}`
                           : `RHX1 lift ${d.lift} K · RHX1 share ${(d.alpha * 100).toFixed(0)} % · Pmax ${d.Pmax.toFixed(1)} bar · COP ${d.COP.toFixed(3)}`;
    info = `<div class="cyrow"><i class="sw" style="background:${col}"></i><b>${nm}</b><span class="f">${WLAB[w]} · ${desc}${same ? '' : ' <span class="flag">· differs from the design selected above (diagrams use the 1 % rule at ' + (setName === 'none' ? 'no limit' : '27 bar') + ')</span>'}</span></div>`;
  }
  return {k, name: nm, color: col, w, d, Tair: state.td, Tsrc: 40, dome: CY && CY.domes[w], info,
          title: k === 'A' ? `${nm} · RHX1 in series` : 'SLTHP · RHX1 parallel to the LT gas cooler'};
}
function renderCycles() {
  const src = cycleSrc(), info = document.getElementById('cyInfo');
  if (!src) { info.innerHTML = '<p class="lbl">Solving…</p>'; for (const id of ['dTs', 'dPh', 'dPr', 'dPr2', 'lgPr2']) document.getElementById(id).innerHTML = ''; document.getElementById('lgPr').innerHTML = ''; return; }
  info.innerHTML = src.info;
  if (src.generic) { genericDiagrams(src); return; }
  const d = src.d;
  const make = (ix, iy) => {
    const L = [], dm = src.dome;
    if (dm) L.push({pts: [...dm.bubble, ...dm.dew.slice().reverse()].map(p => [p[ix], p[iy]]), color: 'var(--dome)', width: 1.1, dash: '', fill: 'var(--dome-fill)', fit: false});
    if (d && d.paths) for (const path of d.paths) L.push({pts: path.pts.map(p => [p[ix], p[iy]]), raw: path.pts, color: CCOL[compGroup(path.name)][1], width: 2.2, dash: path.name === 'RHX1 (hot)' ? '6 4' : '', name: src.name, comp: path.name, tip: true});
    return L;
  };
  const used = new Set(d && d.paths ? d.paths.map(p => compGroup(p.name)) : []);
  document.getElementById('lgComp').innerHTML = Object.entries(CCOL).filter(([g]) => used.has(g))
    .map(([, [lab, col]]) => `<span><i class="lg-line" style="border-color:${col}"></i>${lab}</span>`).join('') +
    (src.k === 'S' ? '<span class="lbl">RHX1 hot branch dashed (runs parallel to the LT gas cooler)</span>' : '') +
    '<span><i class="lg-line" style="border-color:var(--dome)"></i>bubble / dew line</span>';
  xy('dTs', make(3, 0), {xlabel: 's (kJ/kg·K)', ylabel: 'T (°C)'});
  xy('dPh', make(2, 1), {xlabel: 'h (kJ/kg)', ylabel: 'P (bar, log)', ylog: true});
  renderProfiles(src);
}

/* heat exchange profiles, as the steam page: temperature against heat transferred (kW at the swept volume), each
   exchanger from its cold end; the other stream dashed, the heat transfer between the streams shaded */
const HXCOL = {gc: '#ff7f0e', rhx: '#9467bd', evap: '#1f77b4', ihx: '#8c564b'};
function profileData(src) {
  if (src.generic) return genericProfileData(src);
  const {k, d} = src;
  if (!d || !d.paths) return null;
  const u = energy(k, d, state.vsw, src.Tair, src.Tain ?? 15);
  if (!u) return null;
  const m = u.m, by = n => d.paths.filter(p => p.name === n), H = p => p[2], T = p => p[0];
  const Tain = src.Tain ?? 15, A = src.air;      // live: the air property table; study: constant cp
  const mAir = A ? u.kw.air / (A.h(src.Tair) - A.h(Tain)) : u.mAir;
  const airT = A ? (q => A.T(A.h(Tain) + q / mAir)) : (q => Tain + q / (mAir * CP_AIR));
  // a refrigerant path as [kW from x0, T], cold end first (hot streams run hot to cold, so reverse them)
  const prof = (pts, flow, x0 = 0, hot = true) => {
    const p = hot ? pts.slice().reverse() : pts.slice(), h0 = H(p[0]);
    return p.map(q => [x0 + flow * Math.abs(H(q) - h0), T(q)]);
  };
  const S = [];
  // gas coolers in series on the air, cold end first
  let gcs, flows;
  if (k === 'A') { gcs = by('Gas cooler').slice().reverse(); flows = gcs.map(() => m); }
  else { const a = d.alpha ?? 0, L = by('LT gas cooler')[0], Hh = by('HT gas cooler')[0]; gcs = [L, Hh].filter(Boolean); flows = [Hh ? (1 - a) * m : m, m]; }
  // each gas cooler against its own stretch of the air (the refrigerant steps between them where RHX1 sits)
  let x = 0;
  gcs.forEach((g, i) => {
    const p = prof(g.pts, flows[i], x), x1 = p[p.length - 1][0];
    const air = Array.from({length: 11}, (_, j) => { const q = x + (x1 - x) * j / 10; return [q, airT(q)]; });
    const nm = gcs.length === 1 ? 'Gas cooler' : i === 0 ? 'LT gas cooler' : i === gcs.length - 1 ? 'HT gas cooler' : 'MT gas cooler';
    S.push({key: 'gc', name: nm, col: HXCOL.gc, ref: p, other: air, otherName: 'Dryer air', dash: i === 0 && gcs.length > 1 ? '' : ''});
    x = x1;
  });
  // recuperators: hot side against the gas it heats (RHX1 the stage-1 suction, RHX2 the interstage gas)
  const rxs = k === 'A' ? [['Recuperator', 'RHX1', 'suction'], ['Recuperator 2', 'RHX2', 'interstage gas']] : [['RHX1', 'RHX1', 'suction']];
  for (const [pn, nm, what] of rxs) {
    const rc = k === 'A' ? by(pn).slice().sort((a, b) => a.pts[0][1] - b.pts[0][1]) : [by('RHX1 (cold)')[0], by('RHX1 (hot)')[0]];
    const cold = rc[0], hotP = rc.length > 1 ? rc[rc.length - 1] : null;
    if (!cold || !hotP) continue;
    const fh = k === 'A' ? m : (d.alpha ?? 0) * m;
    S.push({key: 'rhx', name: nm, col: nm === 'RHX2' ? '#c5b0d5' : HXCOL.rhx, ref: prof(hotP.pts, fh), other: prof(cold.pts, m, 0, false), otherName: `${nm} cold side (${what})`});
  }
  // evaporator and superheater against the source water
  const ev = by('Evaporator').flatMap(p => p.pts).sort((a, b) => H(a) - H(b));
  if (ev.length > 1) {
    const e = prof(ev, m, 0, false), qe = e[e.length - 1][0];
    const water = src.srcCurve ? Array.from({length: 17}, (_, i) => [qe * i / 16, src.srcCurve(i / 16)])
      : src.Twout != null ? [[0, src.Twout], [qe, src.Tsrc]] : [[0, src.Tsrc], [qe, src.Tsrc]];
    S.push({key: 'evap', name: d.scLow ? 'Gas heater (supercritical)' : 'Evaporator and superheater', short: d.scLow ? 'Gas heater' : null, col: HXCOL.evap, ref: e, other: water, otherName: src.srcCurve ? 'Source air' : src.Twout != null ? 'Source water' : 'Source water (isothermal, as the study)', cold: true});
    const ih = by('Internal heat exchanger');
    if (ih.length === 2) {
      const ic = ih.slice().sort((a, b) => a.pts[0][1] - b.pts[0][1]);
      S.push({key: 'ihx', name: 'Internal heat exchanger', col: HXCOL.ihx, ref: prof(ic[1].pts, m, qe), other: prof(ic[0].pts, m, qe, false), otherName: 'IHX cold side'});
    }
  }
  // approach: the smallest gap between the streams at the same heat
  const at = (line, q) => { for (let i = 1; i < line.length; i++) if (line[i][0] >= q - 1e-9) { const [a, b] = [line[i - 1], line[i]]; return b[0] === a[0] ? b[1] : a[1] + (b[1] - a[1]) * (q - a[0]) / (b[0] - a[0]); } return line[line.length - 1][1]; };
  for (const s of S) s.approach = Math.min(...s.ref.map(([q, t]) => (s.cold ? -1 : 1) * (t - at(s.other, q))));
  return {S, u};
}
/* two charts: the exchangers with the outside streams (dryer air, heat source) and the internal ones (recuperators,
   IHX, cascade exchangers, steam generator) */
const INTERNAL_HX = new Set(['rhx', 'ihx']);
function renderProfiles(src) {
  const P = profileData(src), panel2 = document.getElementById('dPr2Panel');
  const ext = P ? {...P, S: P.S.filter(s => !INTERNAL_HX.has(s.key))} : null;
  // internal exchangers each from 0 kW (they are not in series with one another)
  const zero = s => { const lo = Math.min(...s.ref.map(v => v[0]), ...s.other.map(v => v[0])); return lo ? {...s, ref: s.ref.map(v => [v[0] - lo, v[1]]), other: s.other.map(v => [v[0] - lo, v[1]])} : s; };
  const int = P ? {...P, S: P.S.filter(s => INTERNAL_HX.has(s.key)).map(zero)} : null;
  drawProfiles('dPr', 'lgPr', ext && ext.S.length ? ext : null);
  if (panel2) panel2.hidden = !(int && int.S.length);
  if (int && int.S.length) drawProfiles('dPr2', 'lgPr2', int, true); else { document.getElementById('dPr2').innerHTML = ''; document.getElementById('lgPr2').innerHTML = ''; }
}
function drawProfiles(svgId, lgId, P, internal = false) {
  // dashed: the other stream; internal exchangers draw it in the exchanger's own colour so each pair reads together
  const ocol = s => internal ? s.col : s.cold ? 'var(--accent)' : 'var(--ink-2)';
  if (P && internal) { const pal = ['#9467bd', '#8c564b', '#e377c2', '#17becf', '#bcbd22', '#d62728', '#1f77b4', '#7f7f7f']; P = {...P, S: P.S.map((x, i) => ({...x, col: pal[i % pal.length]}))}; }
  const svg = document.getElementById(svgId), lg = document.getElementById(lgId); svg.innerHTML = '';
  if (!P) { lg.innerHTML = ''; el('text', {x: 200, y: 120, 'text-anchor': 'middle', class: 'lbl'}, svg).textContent = 'no profiles for this case'; svg.setAttribute('viewBox', '0 0 400 240'); return; }
  const W = Math.max(260, Math.round((svg.parentElement.clientWidth || 400) - 36)), Hh = Math.round(Math.max(240, W * 0.72));
  const mg = {l: 50, r: 12, t: 10, b: 36};
  svg.setAttribute('viewBox', `0 0 ${W} ${Hh}`);
  const all = P.S.flatMap(s => [...s.ref, ...s.other]);
  let x1 = Math.max(...all.map(p => p[0])), y0 = Math.min(...all.map(p => p[1])), y1 = Math.max(...all.map(p => p[1]));
  const xt = ticks(0, x1 * 1.03, 5); x1 = xt[xt.length - 1];
  const yt = ticks(y0 - 5, y1 + 5, 6); y0 = yt[0]; y1 = yt[yt.length - 1];
  const X = v => mg.l + v / x1 * (W - mg.l - mg.r), Y = v => Hh - mg.b - (v - y0) / (y1 - y0) * (Hh - mg.t - mg.b);
  const g = el('g', {class: 'axis'}, svg);
  for (const t of yt) { el('line', {x1: mg.l, x2: W - mg.r, y1: Y(t), y2: Y(t), class: 'gridline'}, g); el('text', {x: mg.l - 6, y: Y(t) + 3.5, 'text-anchor': 'end'}, g).textContent = t; }
  for (const t of xt) { el('line', {x1: X(t), x2: X(t), y1: mg.t, y2: Hh - mg.b, class: 'gridline'}, g); el('text', {x: X(t), y: Hh - mg.b + 14, 'text-anchor': 'middle'}, g).textContent = +t.toFixed(3); }
  el('text', {x: (mg.l + W - mg.r) / 2, y: Hh - 4, 'text-anchor': 'middle', class: 'title'}, g).textContent = `Heat transferred (kW at ${vswTxt()} m³/h swept)`;
  el('text', {x: 12, y: (mg.t + Hh - mg.b) / 2, 'text-anchor': 'middle', class: 'title', transform: `rotate(-90 12 ${(mg.t + Hh - mg.b) / 2})`}, g).textContent = 'T (°C)';
  const path = pts => pts.map((p, i) => (i ? 'L' : 'M') + X(p[0]).toFixed(1) + ' ' + Y(p[1]).toFixed(1)).join('');
  for (const s of P.S) {   // shading first
    const top = s.cold ? s.other : s.ref, bot = s.cold ? s.ref : s.other;
    el('path', {d: path([...top, ...bot.slice().reverse()]) + 'Z', fill: s.col, 'fill-opacity': 0.14, stroke: 'none'}, svg);
  }
  const hov = [], placed = [];
  for (const s of P.S) {
    el('path', {d: path(s.other), fill: 'none', stroke: ocol(s), 'stroke-width': 1.5, 'stroke-dasharray': '5 3'}, svg);
    el('path', {d: path(s.ref), fill: 'none', stroke: s.col, 'stroke-width': 2.4, 'stroke-linecap': 'round', 'stroke-linejoin': 'round'}, svg);
    for (const p of s.ref) hov.push({p, name: s.name, col: s.col});
    for (const p of s.other) hov.push({p, name: internal ? `${s.name} cold side: ${s.otherName}` : s.otherName, col: ocol(s)});
    // one label per exchanger at its hot end (the gas coolers together), kept inside the plot
    const gcsAll = P.S.filter(q => q.key === 'gc'), lastGc = gcsAll[gcsAll.length - 1];
    if (s.key === 'gc' && s !== lastGc) continue;
    const e = s.ref[s.ref.length - 1], q0 = s.key === 'gc' ? 0 : s.ref[0][0];
    const ly = Y(e[1]) - 7 < mg.t + 10 ? Y(e[1]) + 16 : Y(e[1]) - 7;
    const txt = `${internal ? s.name : s.key === 'gc' ? (gcsAll.length > 1 ? 'Gas coolers' : 'Gas cooler') : s.key === 'evap' ? s.short || 'Evaporator' : s.key === 'ihx' ? 'IHX' : s.name} ${(e[0] - q0).toFixed(1)} kW`;
    const left = X(e[0]) - 6 - txt.length * 6.3 < mg.l, tw = txt.length * 6.3, lx = left ? mg.l + 4 : X(e[0]) - 6;
    // keep labels from sitting on one another (the stacked air exchangers can end near other labels)
    let yy = ly;
    const box = y => ({x0: left ? lx : lx - tw, x1: left ? lx + tw : lx, y0: y - 10, y1: y + 3});
    const hits = b => placed.some(q => b.x0 < q.x1 && q.x0 < b.x1 && b.y0 < q.y1 && q.y0 < b.y1);
    for (let k = 1; k < 12 && hits(box(yy)); k++) yy = ly + (k % 2 ? 1 : -1) * 13 * Math.ceil(k / 2);
    placed.push(box(yy));
    el('text', {x: lx, y: yy, 'text-anchor': left ? 'start' : 'end', class: 'lbl'}, svg).textContent = txt;
  }
  const hit = el('rect', {x: mg.l, y: mg.t, width: W - mg.l - mg.r, height: Hh - mg.t - mg.b, fill: 'transparent'}, svg);
  const dot = el('circle', {r: 4.5, fill: 'var(--panel)', 'stroke-width': 2, visibility: 'hidden'}, svg);
  const show = ev => {
    const r = svg.getBoundingClientRect(), vx = (ev.clientX - r.left) / r.width * W, vy = (ev.clientY - r.top) / r.height * Hh;
    let b = null, bd = 1e9; for (const q of hov) { const dd = (X(q.p[0]) - vx) ** 2 + (Y(q.p[1]) - vy) ** 2; if (dd < bd) { bd = dd; b = q; } }
    if (!b || bd > 900) { tip.hidden = true; dot.setAttribute('visibility', 'hidden'); return; }
    dot.setAttribute('cx', X(b.p[0])); dot.setAttribute('cy', Y(b.p[1])); dot.setAttribute('stroke', b.col); dot.setAttribute('visibility', 'visible');
    tip.innerHTML = `<div class="tt">${b.name}</div><div class="note">${b.p[0].toFixed(2)} kW · ${b.p[1].toFixed(1)} °C</div>`;
    tip.hidden = false;
    tip.style.left = Math.min(ev.clientX + 14, innerWidth - tip.offsetWidth - 8) + 'px'; tip.style.top = Math.min(ev.clientY + 14, innerHeight - tip.offsetHeight - 8) + 'px';
  };
  hit.addEventListener('pointermove', show); hit.addEventListener('pointerdown', show);
  hit.addEventListener('pointerleave', () => { tip.hidden = true; dot.setAttribute('visibility', 'hidden'); });
  // the dashed lines by name (the dryer air in grey, the heat source in blue), so it is clear which stream each is
  const others = new Map();
  for (const s of P.S) {
    const on = s.otherName || 'other stream', nm = internal && !/cold side/i.test(on) ? `${s.name} cold side: ${on}` : on.replace(/^./, c => c.toUpperCase());
    if (!others.has(nm)) others.set(nm, ocol(s));
  }
  lg.innerHTML = P.S.map(s => `<span><i class="lg-line" style="border-color:${s.col}"></i>${s.name} · approach ${s.approach.toFixed(1)} K</span>`).join('') +
    [...others].map(([nm, col]) => `<span><i class="lg-dash" style="border-color:${col}"></i>${nm}</span>`).join('');
}

/* ---------- process flow schematics ---------- */
const f0 = v => v == null ? '–' : v.toFixed(0), f1 = v => v == null ? '–' : v.toFixed(1);
function endpoints(d) {
  const by = n => d.paths.filter(p => p.name === n), first = p => p && p.pts[0], last = p => p && p.pts[p.pts.length - 1];
  return {by, first, last};
}
function arrow(x1, y1, x2, y2) { /* small arrowhead at (x2,y2) pointing along the segment */
  const a = Math.atan2(y2 - y1, x2 - x1), L = 7;
  const p = (t, r) => `${(x2 - L * Math.cos(a + t)).toFixed(1)},${(y2 - L * Math.sin(a + t)).toFixed(1)}`;
  return `<polygon points="${x2},${y2} ${p(0.45)} ${p(-0.45)}" fill="currentColor"/>`;
}
const T_ = (x, y, t, cls = '', anc = 'middle') => `<text x="${x}" y="${y}" text-anchor="${anc}" class="${cls}">${t}</text>`;
const comp = (x, y, lab) => `<polygon class="eq" points="${x - 18},${y - 26} ${x + 18},${y - 14} ${x + 18},${y + 14} ${x - 18},${y + 26}"/>${T_(x, y + 4, lab, 'cn')}`;
const box = (x, y, w, h, lab) => `<rect class="eq" x="${x}" y="${y}" width="${w}" height="${h}" rx="3"/>${T_(x + w / 2, y + h / 2 + 4, lab, 'cn')}`;
const valveSym = (x, y) => `<polygon class="eq" points="${x - 9},${y - 12} ${x + 9},${y - 12} ${x - 9},${y + 12} ${x + 9},${y + 12}"/>`;
/* duties and power from the stored states, scaled by the low-stage swept volume */
const MOTOR = 0.93, CP_AIR = 1.01;
function energy(k, d, Vsw, Tair, Tain = 15) {
  if (!d || !d.paths || !d.rho1) return null;
  const by = n => d.paths.filter(p => p.name === n), h0 = p => p.pts[0][2], h1 = p => p.pts[p.pts.length - 1][2];
  const comps = k === 'A' ? by('Compressor') : [by('Compressor 1')[0], by('Compressor 2')[0]];
  const c1 = comps[0], c2 = comps[1];
  const PR1 = c1.pts[c1.pts.length - 1][1] / c1.pts[0][1], nv = state.nvMode === 'fixed' ? state.nvFix : 1 - 0.08244 * Math.pow(Math.max(PR1 - 1, 0), 0.72773);
  const m = d.rho1 * nv * Vsw / 3600;   // kg/s
  let ht, lt, rx, ev;
  if (k === 'A') {
    const gcs = by('Gas cooler'), recs = d.paths.filter(p => /^Recuperator/.test(p.name));
    const PH = Math.max(...d.paths.map(p => p.pts[0][1]));
    ht = h0(gcs[0]) - h1(gcs[0]); lt = gcs.length > 1 ? gcs.slice(1).reduce((s, g) => s + h0(g) - h1(g), 0) : 0;
    rx = recs.filter(p => Math.abs(p.pts[0][1] - PH) < 1e-6 * PH + 1e-9).reduce((s, p) => s + h0(p) - h1(p), 0);   // hot sides
    ev = by('Evaporator').reduce((s, p) => s + h1(p) - h0(p), 0);
  } else {
    const a = d.alpha ?? 0, H = by('HT gas cooler')[0], L = by('LT gas cooler')[0], R = by('RHX1 (hot)')[0];
    ht = H ? h0(H) - h1(H) : 0; lt = (H ? 1 - a : 1) * (h0(L) - h1(L)); rx = R ? a * (h0(R) - h1(R)) : 0;
    ev = by('Evaporator').reduce((s, p) => s + h1(p) - h0(p), 0);
  }
  const w1 = h1(c1) - h0(c1), w2 = c2 ? h1(c2) - h0(c2) : 0, wx = by('Expander').reduce((s, p) => s + h0(p) - h1(p), 0);
  const kw = {ht: m * ht, lt: m * lt, rx: m * rx, ev: m * ev, w1: m * w1, w2: m * w2, ex: m * wx}, mot = d.motor ?? MOTOR;
  kw.el = (kw.w1 + kw.w2) / mot - kw.ex * mot; kw.air = kw.ht + kw.lt; kw.loss = kw.el - kw.w1 - kw.w2 + kw.ex;
  let V2 = null, PR2 = null, nv2 = null;
  if (c2) {
    PR2 = c2.pts[c2.pts.length - 1][1] / c2.pts[0][1];
    nv2 = state.nvMode === 'fixed' ? state.nvFix : 1 - 0.08244 * Math.pow(Math.max(PR2 - 1, 0), 0.72773);
    if (d.rho2) V2 = m / (d.rho2 * nv2) * 3600;
  }
  return {m, nv, PR1, PR2, nv2, V2, kw, rho: d.rho1, mAir: kw.air / (CP_AIR * (Tair - Tain)), cop: kw.air / kw.el};
}
const fk = v => v >= 100 ? Math.round(v).toLocaleString('en-NZ') : v.toFixed(1);
const Q = (x, y, t, cls, anc = 'middle') => `<text x="${x}" y="${y}" text-anchor="${anc}" class="q ${cls}">${t}</text>`;

/* saturation temperature at pressure P (bar) from the stored bubble or dew line, interpolated in ln P */
function satT(line, P) {
  if (!line || line.length < 2) return null;
  const pts = line.slice().sort((a, b) => a[1] - b[1]);
  if (P < pts[0][1] || P > pts[pts.length - 1][1]) return null;
  for (let i = 1; i < pts.length; i++) if (pts[i][1] >= P) {
    const [a, b] = [pts[i - 1], pts[i]], f = (Math.log(P) - Math.log(a[1])) / (Math.log(b[1]) - Math.log(a[1]));
    return a[0] + f * (b[0] - a[0]);
  }
  return null;
}
const evapSat = (dm, PL, x) => {
  if (!dm) return '';
  const tb = satT(dm.bubble, PL), td = satT(dm.dew, PL);
  return tb == null || td == null ? '' : T_(x, 225, `bubble ${f0(tb)} °C`, 'v', 'end') + T_(x, 238, `dew ${f0(td)} °C`, 'v', 'end');
};
function pfdRHP(d, Tair, u, dm, Tsrc = 40, Tain = 15, Twout = null) {
  const {by, first, last} = endpoints(d);
  const comps = by('Compressor'), gcs = by('Gas cooler'), recs = by('Recuperator'), ev = by('Evaporator');
  const rc = recs.slice().sort((a, b) => a.pts[0][1] - b.pts[0][1]), cold = rc[0], hot = rc[rc.length - 1];
  const ht = gcs[0], lt = gcs.length > 1 ? gcs[gcs.length - 1] : null;   // one gas cooler when there is no RHX1 lift
  const suc = cold ? last(cold) : first(comps[0]), s1 = last(comps[0]), dis = last(comps[1] || comps[0]), htOut = last(ht);
  const rxOut = hot && rc.length > 1 ? last(hot) : htOut, ltOut = lt ? last(lt) : htOut;
  const vIn = first(by('Expansion valve')[0] || by('Expander')[0]), eIn = first(ev[0]), sup = last(ev[ev.length - 1]);
  const qHT = ht.pts[0][2] - htOut[2], qLT = lt ? lt.pts[0][2] - ltOut[2] : 0, Tmid = Tain + (Tair - Tain) * qLT / (qHT + qLT);
  const PH = dis[1], PI = s1[1], PL = suc[1];
  const sat = evapSat(dm, PL, 505);
  const k = u && u.kw;
  const duty = !k ? '' : `
  ${Q(190, 47, `▲ ${fk(k.ht)} kW`, 'qh')}${Q(440, 47, `▲ ${fk(k.lt)} kW`, 'qh')}${Q(320, 188, `${fk(k.rx)} kW`, 'qw')}
  ${Q(440, 297, `▼ ${fk(k.ev)} kW`, 'ql')}${Q(34, 166, `${fk(k.el)} kW`, 'qw', 'end')}${T_(34, 178, 'electric', '', 'end')}${Q(150, 278, `${(u.m * 3600).toFixed(u.m * 3600 >= 100 ? 0 : 1)} kg/h · ${(u.m / u.rho * 3600).toFixed(u.m / u.rho * 3600 >= 100 ? 0 : 1)} m³/h`, 'ql')}`;
  return `<svg viewBox="-50 0 610 330" role="img" aria-label="RHP-A process flow">
  <g style="color:var(--hp)">
    <path class="ph" d="M60 112 V70 H150"/>${arrow(140, 70, 150, 70)}
    <path class="ph" d="M230 70 H290"/>${arrow(280, 70, 290, 70)}
    <path class="ph" d="M350 70 H400"/>${arrow(390, 70, 400, 70)}
    <path class="ph" d="M480 70 H515 V150"/>${arrow(515, 140, 515, 150)}
  </g>
  <g style="color:var(--lp)">
    <path class="pl" d="M515 174 V262 H480"/>${arrow(490, 262, 480, 262)}
    <path class="pl" d="M400 262 H350"/>${arrow(360, 262, 350, 262)}
    <path class="pl" d="M290 262 H60 V228"/>${arrow(60, 238, 60, 228)}
    <path class="pl" d="M60 176 V150" style="stroke:var(--hp)"/>
  </g>
  ${comp(60, 202, 'C1')}${comp(60, 128, 'C2')}
  ${box(150, 52, 80, 36, 'HT GC')}${box(400, 52, 80, 36, 'LT GC')}
  <rect class="eq" x="290" y="52" width="60" height="228" rx="3"/>${T_(320, 170, 'RHX1', 'cn')}
  ${box(400, 244, 80, 36, d.scLow ? 'Gas heater' : 'Evaporator')}${valveSym(515, 162)}${T_(501, 166, by('Expander').length ? 'expander' : 'valve', '', 'end')}
  <path class="air" d="M520 22 H130"/><g style="color:var(--ink-3)">${arrow(140, 22, 130, 22)}</g>
  ${T_(522, 18, `air ${f1(Tain)} °C`, '', 'end')}${T_(365, 18, `${f0(Tmid)} °C`, 'v')}${T_(132, 18, `${Tair} °C`, 'v', 'start')}
  <path class="wat" d="M496 304 H384"/><g style="color:var(--accent)">${arrow(394, 304, 384, 304)}</g>${T_(498, 318, Twout == null ? `water ${f1(Tsrc)} °C in` : `water ${f1(Tsrc)} → ${f1(Twout)} °C`, '', 'end')}
  ${T_(60, 62, `${f0(dis[0])} °C`, 'v')}${T_(260, 64, `${f0(htOut[0])}`, 'v')}${T_(375, 64, `${f0(rxOut[0])}`, 'v')}${T_(500, 64, `${f0(ltOut[0])}`, 'v')}
  ${T_(500, 256, `${f0(eIn[0])}`, 'v')}${T_(375, 256, `${f0(sup[0])}`, 'v')}${T_(150, 256, `${f0(suc[0])} °C`, 'v')}${T_(84, 168, `${f0(s1[0])} °C`, 'v', 'start')}
  ${T_(84, 182, `${f1(PI)} bar`, '', 'start')}${T_(160, 112, `P high ${f1(PH)} bar`, '', 'start')}${T_(240, 306, `P low ${f1(PL)} bar`, '', 'start')}${T_(34, 196, `overall PR`, '', 'end')}${Q(34, 208, f1(PH / PL), 'qw', 'end')}
  ${T_(320, 120, '▲ hot', '')}${T_(320, 230, 'cold ▼', '')}${duty}${sat}
</svg>`;
}
function pfdSLT(d, Tair, u, dm, Tsrc = 40, Tain = 15, Twout = null) {
  const {by, first, last} = endpoints(d);
  const c1 = by('Compressor 1')[0], c2 = by('Compressor 2')[0], ht = by('HT gas cooler')[0], lt = by('LT gas cooler')[0], rh = by('RHX1 (hot)')[0], rcld = by('RHX1 (cold)')[0], ev = by('Evaporator')[0], vv = by('Expansion valve')[0] || by('Expander')[0];
  const suc = last(rcld), s1 = last(c1), dis = last(c2), sp = ht ? last(ht) : dis, ltOut = last(lt), rxOut = rh && last(rh), vIn = first(vv), eIn = first(ev), sup = last(ev);
  const a = d.alpha ?? 0, qHT = ht ? ht.pts[0][2] - sp[2] : 0, qLT = (ht ? 1 - a : 1) * (lt.pts[0][2] - ltOut[2]), Tmid = Tain + (Tair - Tain) * qLT / (qHT + qLT);
  const PH = dis[1], PI = s1[1], PL = suc[1];
  const sat = evapSat(dm, PL, 505);
  const k = u && u.kw;
  const duty = !k ? '' : `
  ${Q(170, 47, `▲ ${fk(k.ht)} kW`, 'qh')}${Q(440, 17, `▲ ${fk(k.lt)} kW`, 'qh')}${Q(320, 214, `${fk(k.rx)} kW`, 'qw')}
  ${Q(440, 297, `▼ ${fk(k.ev)} kW`, 'ql')}${Q(34, 166, `${fk(k.el)} kW`, 'qw', 'end')}${T_(34, 178, 'electric', '', 'end')}${Q(150, 278, `${(u.m * 3600).toFixed(u.m * 3600 >= 100 ? 0 : 1)} kg/h · ${(u.m / u.rho * 3600).toFixed(u.m / u.rho * 3600 >= 100 ? 0 : 1)} m³/h`, 'ql')}`;
  return `<svg viewBox="-50 -18 610 348" role="img" aria-label="SLTHP process flow">
  <g style="color:var(--hp)">
    <path class="ph" d="M60 112 V70 H130"/>${arrow(120, 70, 130, 70)}
    <path class="ph" d="M210 70 H240"/>
    <path class="ph" d="M240 70 V40 H400"/>${arrow(390, 40, 400, 40)}
    <path class="ph" d="M240 70 V120 H290"/>${arrow(280, 120, 290, 120)}
    <path class="ph" d="M480 40 H515 V120"/>
    <path class="ph" d="M350 120 H515"/>${arrow(505, 120, 515, 120)}
    <path class="ph" d="M515 120 V150"/>${arrow(515, 140, 515, 150)}
  </g>
  <g style="color:var(--lp)">
    <path class="pl" d="M515 174 V262 H480"/>${arrow(490, 262, 480, 262)}
    <path class="pl" d="M400 262 H350"/>${arrow(360, 262, 350, 262)}
    <path class="pl" d="M290 262 H60 V228"/>${arrow(60, 238, 60, 228)}
    <path class="pl" d="M60 176 V150" style="stroke:var(--hp)"/>
  </g>
  <circle class="node" cx="240" cy="70" r="4"/><circle class="node" cx="515" cy="120" r="4"/>
  ${comp(60, 202, 'C1')}${comp(60, 128, 'C2')}
  ${box(130, 52, 80, 36, 'HT GC')}${box(400, 22, 80, 36, 'LT GC')}
  <rect class="eq" x="290" y="100" width="60" height="180" rx="3"/>${T_(320, 196, 'RHX1', 'cn')}
  ${box(400, 244, 80, 36, d.scLow ? 'Gas heater' : 'Evaporator')}${valveSym(515, 162)}${T_(501, 166, by('Expander').length ? 'expander' : 'valve', '', 'end')}
  <path class="air" d="M520 0 H110"/><g style="color:var(--ink-3)">${arrow(120, 0, 110, 0)}</g>
  ${T_(522, -5, `air ${f1(Tain)} °C`, '', 'end')}${T_(330, -5, `${f0(Tmid)} °C`, 'v')}${T_(112, -5, `${Tair} °C`, 'v', 'start')}
  <path class="wat" d="M496 304 H384"/><g style="color:var(--accent)">${arrow(394, 304, 384, 304)}</g>${T_(498, 318, Twout == null ? `water ${f1(Tsrc)} °C in` : `water ${f1(Tsrc)} → ${f1(Twout)} °C`, '', 'end')}
  ${T_(60, 62, `${f0(dis[0])} °C`, 'v')}${T_(232, 92, `${f0(sp[0])}`, 'v', 'end')}${T_(500, 34, `${f0(ltOut[0])}`, 'v')}
  ${T_(270, 34, `${((1 - a) * 100).toFixed(0)} %`, '')}${T_(266, 136, `${(a * 100).toFixed(0)} %`, '')}
  ${rxOut ? T_(372, 114, `${f0(rxOut[0])}`, 'v') : ''}${T_(509, 112, `mix ${f0(vIn[0])} °C`, 'v', 'end')}
  ${T_(500, 256, `${f0(eIn[0])}`, 'v')}${T_(375, 256, `${f0(sup[0])}`, 'v')}${T_(150, 256, `${f0(suc[0])} °C`, 'v')}${T_(84, 168, `${f0(s1[0])} °C`, 'v', 'start')}
  ${T_(84, 182, `${f1(PI)} bar`, '', 'start')}${T_(140, 112, `P high ${f1(PH)} bar`, '', 'start')}${T_(240, 306, `P low ${f1(PL)} bar`, '', 'start')}${T_(34, 196, `overall PR`, '', 'end')}${Q(34, 208, f1(PH / PL), 'qw', 'end')}
  ${T_(320, 140, '▲ hot', '')}${T_(320, 236, 'cold ▼', '')}${duty}${sat}
</svg>`;
}
/* RHP-B and RHP-C: the cycle schematic with a table of the states around the cycle */
function pfdStates(d, Tair, u, dm, Tsrc, Tain, Twout) {
  const rows = d.paths.filter(p => p.name !== 'Expansion valve' && p.pts.length && Math.abs(p.pts[p.pts.length - 1][2] - p.pts[0][2]) > 0.05).map(p => {
    const a = p.pts[0], b = p.pts[p.pts.length - 1], q = u ? u.m * Math.abs(b[2] - a[2]) : null;
    return `<tr><td>${p.name === 'Evaporator' && d.scLow ? 'Gas heater' : p.name}</td><td>${f1(a[0])} → ${f1(b[0])}</td><td>${a[1].toFixed(2)}${Math.abs(b[1] - a[1]) > 1e-3 ? ' → ' + b[1].toFixed(2) : ''}</td><td>${q == null ? '–' : fk(q)}</td></tr>`;
  }).join('');
  return `${typeof miniRHP === 'function' ? `<figure class="mini" style="max-width:460px">${miniRHP(d.variant, d.stages)}</figure>` : ''}
    <table class="sumt"><thead><tr><th>Step</th><th>T, °C</th><th>P, bar</th><th>kW</th></tr></thead><tbody>${rows}</tbody></table>`;
}
function renderPFD() {
  if (state.mode === 'live' && typeof lvHeadline === 'function') lvHeadline();
  { const ps = cycleSrc(); if (ps) renderProfiles(ps); }
  const host = document.getElementById('pfd'), src = cycleSrc();
  if (src && src.generic) { host.innerHTML = genericPFD(src); document.getElementById('vswSum').innerHTML = genericSummary(src); return; }
  const out = [], sums = [];
  if (src) {
    const {k, d, w, Tair} = src, fn = src.variant === 'B' || src.variant === 'C' || src.variant === 'H' ? pfdStates : k === 'A' ? pfdRHP : pfdSLT;
    const u = energy(k, d, state.vsw, Tair, src.Tain ?? 15);
    if (u && src.Twout != null && !src.srcAir) u.mWater = u.kw.ev / (4.18 * (src.Tsrc - src.Twout));
    if (u && src.srcAir) u.srcAir = src.srcAir;
    if (u && src.live && src.air) { u.mAir = u.kw.air / (src.air.h(Tair) - src.air.h(src.Tain ?? 15)); u.mSec = (src.phi || 0) * u.mAir; }
    if (u) sums.push([src.name, u, d, w]);
    let body;
    try { body = d && d.paths ? fn(d, Tair, u, src.dome, src.Tsrc, src.Tain ?? 15, src.Twout ?? null) : '<p class="lbl">no design</p>'; if (src.srcAir) body = body.replace(/>water (?=-?[0-9])/g, '>air '); } catch (e) { body = `<p class="lbl">schematic unavailable (${e.message})</p>`; }
    out.push(`<div class="panel" style="display:grid;gap:6px"><h3>${src.title}</h3><p class="lbl" style="margin:0">${src.label || WLAB[w]} · ${Tair} °C air · COP ${d ? d.COP.toFixed(3) : '–'}${src.live ? ' · CoolProp' : ' · REFPROP'}</p>${body}</div>`);
  }
  host.innerHTML = out.join('');
  const sum = document.getElementById('vswSum');
  if (!sums.length) { sum.innerHTML = src && src.d ? '<p class="lbl">Duties need the suction density: run <code>main_dryer_cycles.py --augment</code>.</p>' : ''; return; }
  const ref = null;
  const ABBR = { Propane: 'C3', Isobutane: 'iC4', Butane: 'nC4', Isopentane: 'iC5', PROPANE: 'C3', ISOBUTAN: 'iC4', BUTANE: 'nC4', IPENTANE: 'iC5' };
  const drawRows = (r, key) => {
    const s = S5[key];
    if (!s) return key ? r('Drawdown 5 °C, 1 bar', '—', 'run main_dryer_state_5C.py') : r('Drawdown 5 °C, 1 bar', '—', 'study blends only (REFPROP)');
    const bad = t => `<span style="color:var(--bad,#c0392b);font-weight:600">${t}</span>`;
    const gas = s.phase === 'vapour';
    const comp = v => v.map((f, i) => `${(f * 100).toFixed(0)} ${ABBR[s.names[i]] || s.names[i]}`).join(' / ');
    return [
      `<tr><td colspan="3" class="lbl" style="padding-top:8px">Drawdown after shutdown · 5 °C, 1 bar abs</td></tr>`,
      r('State', gas ? 'all vapour' : bad(s.phase === 'liquid' ? 'all liquid' : `${(s.q * 100).toFixed(0)} % vapour`), gas ? '' : bad(`${((1 - s.q) * 100).toFixed(0)} % liquid by mass`)),
      r('Dew point at 1 bar', s.Tdew.toFixed(1), `°C · ${s.Tdew < 5 ? `${(5 - s.Tdew).toFixed(1)} K margin` : bad(`${(s.Tdew - 5).toFixed(1)} K above 5 °C`)}`),
      r('Max all-vapour pressure at 5 °C', s.Pdew.toFixed(2), `bar abs (dew) · bubble ${s.Pbub.toFixed(2)}`),
      r(gas ? 'Gas density' : 'Overall / vapour density', gas ? s.D.toFixed(2) : `${s.D.toFixed(1)} / ${s.Dv.toFixed(2)}`, 'kg/m³'),
      s.x ? r('Liquid / vapour comp.', `${comp(s.x)}`, `${comp(s.y)} · mass %`) : '',
    ].join('');
  };
  const single = (r, d) => {   /* one compressor for the whole lift: stage-1 suction to final discharge pressure */
    const S = d && d.gam1;
    const cs = d.paths.filter(p => /^Compressor/.test(p.name));
    const PRo = cs[cs.length - 1].pts[cs[cs.length - 1].pts.length - 1][1] / cs[0].pts[0][1];
    if (!S) return r('Single stage, whole lift', `PR ${PRo.toFixed(1)}`, 'run main_dryer_cycles.py --gamma');
    const gm = (S[0] + S[1]) / 2;
    return [
      `<tr><td colspan="3" class="lbl" style="padding-top:8px">One compressor for the whole lift · overall PR ${PRo.toFixed(2)}</td></tr>`,
      r('γ, suction → isentropic discharge', `${S[0].toFixed(3)}→${S[1].toFixed(3)}`, `γm ${gm.toFixed(3)}`),
      r('Single-stage Vi = PR^(1/γm)', Math.pow(PRo, 1 / gm).toFixed(2), ''),
      r('Single-stage Vi, real-gas isentrope', S[2].toFixed(2), 'ρ(P_dis, s_suc) / ρ_suc'),
    ].join('');
  };
  const gamRows = (r, u, d) => {
    const G = d && d.gam;
    if (!G || !G.length) return d && d.live ? '' : r('γ = cp/cv · Vi', '—', 'run main_dryer_cycles.py --gamma');
    const PR = [u.PR1, u.PR2];
    const st = G.map((g, i) => g && PR[i] != null ? {g1: g[0], g2: g[1], gm: (g[0] + g[1]) / 2, vi: Math.pow(PR[i], 2 / (g[0] + g[1])), vis: g[3]} : null);
    const both = f => st.map(s => s ? f(s) : '–').join(' / ');
    const lab = st.length > 1 ? 'stage 1 / stage 2' : '';
    return [
      r('γ = cp/cv, suction → discharge', both(s => `${s.g1.toFixed(3)}→${s.g2.toFixed(3)}`), lab),
      r('Screw Vi = PR^(1/γm)', both(s => s.vi.toFixed(2)), `γm = mean of suction and discharge · ${both(s => s.gm.toFixed(3))}`),
      r('Vi, real-gas isentrope', both(s => s.vis.toFixed(2)), 'ρ(P_dis, s_suc) / ρ_suc'),
      single(r, d),
    ].join('');
  };
  sum.innerHTML = sums.map(([n, u, d, key]) => {
    const r = (lab, v, unit) => `<tr><td>${lab}</td><td>${v}</td><td class="lbl">${unit}</td></tr>`;
    return `<div class="panel sumpanel"><h3>${n} · ${vswTxt()} m³/h</h3><table class="sumt"><tbody>
      ${r('Heat to air', fk(u.kw.air), 'kW')}
      ${r('Electrical input', fk(u.kw.el), 'kW')}
      ${r('Motor losses', fk(u.kw.loss), 'kW')}
      ${u.kw.ex > 0 ? r('Expander work', fk(u.kw.ex), 'kW') : ''}
      ${r(u.srcAir ? 'Heat from source air' : 'Heat from water', fk(u.kw.ev), 'kW')}
      ${r('RHX1 duty', fk(u.kw.rx), 'kW')}
      ${r('COP', u.cop.toFixed(3), ref && ref[1] !== u ? chip(u.cop, ref[1].cop) : '')}
      ${r('Refrigerant flow', (u.m * 3.6).toFixed(2), 't/h')}
      ${r('Dryer air flow', (u.mAir * 3.6).toFixed(1), 't/h')}
      ${u.mSec ? r('Secondary air flow', (u.mSec * 3.6).toFixed(1), 't/h') : ''}
      ${u.mWater != null ? r('Source water flow', (u.mWater * 3.6).toFixed(1), 't/h') : ''}
      ${u.srcAir ? r('Source air flow / outlet', `${u.srcAir.m.toFixed(1)} kg/s / ${u.srcAir.Tout.toFixed(1)} °C`, `condensate ${(u.srcAir.cond * 3600).toFixed(1)} kg/h`) : ''}
      ${r('Stage-1 PR / ηvol', `${u.PR1.toFixed(2)} / ${u.nv.toFixed(2)}`, state.nvMode === 'fixed' ? 'fixed' : 'model')}
      ${r('Stage-1 suction density', u.rho.toFixed(2), 'kg/m³')}
      ${r('High-stage swept volume needed', u.V2 == null ? '—' : Math.round(u.V2).toLocaleString('en-NZ'), u.V2 == null ? 'run --augment' : `m³/h · LS/HS ${(state.vsw / u.V2).toFixed(2)}`)}
      ${r('Stage-2 PR / ηvol', u.PR2 == null ? '—' : `${u.PR2.toFixed(2)} / ${u.nv2.toFixed(2)}`, '')}
      ${gamRows(r, u, d)}
      ${drawRows(r, key)}
    </tbody></table></div>`;
  }).join('');
}

const STUDY_NOTE = 'One cycle at a time, coloured by component, with the bubble and dew lines of its fluid. Choose the RHP-A fluid here; its composition follows that fluid\'s headline card.';
const LIVE_NOTE = 'The live design, coloured by component, with the bubble and dew lines of its fluid from CoolProp.';
function setMode(m) {
  state.mode = m === 'live' ? 'live' : 'study';
  for (const b of document.querySelectorAll('#modeSwitch button')) b.setAttribute('aria-pressed', String(b.dataset.mode === state.mode));
  const study = state.mode === 'study';
  document.getElementById('study').hidden = !study;
  document.getElementById('live').hidden = study;
  document.getElementById('tblSect').hidden = !study;
  document.getElementById('cySel').hidden = !study;
  document.getElementById('status').hidden = study;
  dataChip();
  document.getElementById('solve-time').hidden = study || !document.getElementById('solve-time').textContent;
  document.getElementById('cyNote').textContent = study ? STUDY_NOTE : LIVE_NOTE;
  // live view in the steam page's arrangement: swept volume in the Compressor card, duties as the fourth diagram tile
  document.getElementById('cySect').classList.toggle('live-layout', !study);
  const ctl = document.getElementById('vswCtl'), aside = document.getElementById('vswAside'), sum = document.getElementById('vswSum');
  (study ? aside : document.getElementById('lvCompCard')).appendChild(ctl);
  aside.hidden = !study;
  (study ? document.querySelector('.pfdwrap') : document.getElementById('cyGrid')).appendChild(sum);
  update();
  if (!study && typeof liveEnter === 'function') liveEnter();
}
function update() {
  save();
  vsw.value = toX(state.vsw); vswn.value = state.vsw; renderNv(); renderSize();
  if (D && state.mode === 'study') {
    slider.value = state.cap; pnum.value = state.cap >= PMAX ? '' : +state.cap.toFixed(1); pnum.placeholder = 'none'; pval.textContent = state.cap >= PMAX ? 'no limit' : 'bar abs';
    seg('segA', (D.AP || ['20', '30']).map(a => [a, a + ' K']), 'a');
    seg('segTol', [['1', 'Lowest pressure within 1 %'], ['0', 'Highest COP']], 'tol');
    const cap = capv(state.cap);
    seg('segTd', TEMPS.map(T => [T, T + ' °C']), 'td');
    if (state.cy !== 'A' && state.cy !== 'S') state.cy = 'A';
    seg('segCy', [['A', 'RHP-A'], ['S', 'SLTHP']], 'cy');
    seg('segCyA', OPTS.filter(o => o.src === 'A').map(o => [o.id, o.id === 'A' ? 'C3/iC5' : o.name.replace('RHP-A ', '')]), 'cyA');
    renderHeadline(cap); renderPerT(); renderTable();
  }
  if (D) { renderCycles(); renderPFD(); }
}
/* the study-data chip shows in study mode; in live mode only if the study results (used to load designs and to
   compare) failed to load */
const DSTAT = (msg, cls = 'chip', live = false) => { const s = document.getElementById('dataStatus'); s.textContent = msg; s.className = cls; s.dataset.live = live ? '1' : ''; dataChip(); };
function dataChip() { const s = document.getElementById('dataStatus'); s.hidden = state.mode === 'live' && !s.dataset.live; }
async function getJSON(name) {
  const r = await fetch(name);
  if (!r.ok) throw new Error(`${name}: HTTP ${r.status}`);
  return r.json();
}
async function studyBoot() {
  for (const b of document.querySelectorAll('#modeSwitch button')) b.addEventListener('click', () => setMode(b.dataset.mode));
  setMode(state.mode);
  try {
    initData(await getJSON('dashboard_data.json'));
    DSTAT('REFPROP study results loaded · loading the cycle diagrams…');
    update();
    if (state.mode === 'live' && typeof liveEnter === 'function') liveEnter();
  } catch (e) { DSTAT(`The REFPROP study results could not load (${e.message}): study designs and comparisons are unavailable`, 'chip warn', true); return; }
  getJSON('state_5C_1bar.json').then(j => { S5 = j; renderPFD(); }).catch(() => {});
  try {
    CY = await getJSON('cycles_data.json');
    DSTAT('REFPROP study results loaded', 'chip info');
    if (state.mode === 'study') { renderCycles(); renderPFD(); }
  } catch (e) { DSTAT(`The cycle diagrams could not load (${e.message})`, 'chip warn'); }
}
let rz; addEventListener('resize', () => { clearTimeout(rz); rz = setTimeout(() => { if (D && state.mode === 'study') renderPerT(); if (D) renderCycles(); }, 150); });

