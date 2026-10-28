/* Advanced dryer heat pumps, live design: one RHP-A or SLTHP design solved here with CoolProp (dryer_model.js over
   the steam-generating heat pump page's CoolPropFluid), compared with the REFPROP study where it has the same case. */
const $ = id => document.getElementById(id);
const CPNAME = {Propane: 'Propane', Isobutane: 'IsoButane', Butane: 'n-Butane', Isopentane: 'Isopentane', Cyclopentane: 'Cyclopentane'};
const famOf = w => w.includes(':') ? w.split(':')[0] : 'Propane_Isopentane';
const fracOf = w => (w.includes(':') ? w.split(':')[1] : w).split('-').map(Number);
const famLabel = f => f.split('_').map(n => SHORT[n] || n).join('/');
const fracLabel = w => fracOf(w).map(v => +(v * 100).toFixed(1)).join('/');
const compsOf = w => famOf(w).split('_').map((n, i) => [CPNAME[n], fracOf(w)[i] * 100]);

const RHP = c => c === 'A' || c === 'B' || c === 'C' || c === 'H';      // RHP-A/B/C and the conventional cycle (H): one model, RHX1 and/or RHX2 or neither
// compressors in series: the conventional cycle and RHP-A can run one or two stages, the other recuperative cycles two
const stagesOf = c => c.cyc === 'H' ? (c.stagesH || 1) : c.cyc === 'A' ? (c.stagesA === 1 ? 1 : 2) : c.cyc === 'S' ? (c.stagesS === 1 ? 1 : 2) : 2;
const lifts = c => ({lift: c.cyc === 'B' || c.cyc === 'H' ? 0 : c.lift, lift2: c.cyc === 'A' || c.cyc === 'H' ? 0 : c.lift2, split: stagesOf(c) === 1 ? 0 : c.split, ihx_on: c.cyc !== 'H' || !!c.ihxH});
const MULTI = c => ['K', 'N', 'M', 'X', 'E', 'D', 'R', 'Y'].includes(c);    // multi-loop cycles (dryer_cascades.js, dryer_more.js)
const NOFLUID = c => c === 'R' || c === 'Y';    // Stirling (helium, not modelled) and the ammonia-water hybrid: no refrigerant to pick
const LV = {stagesS: 2, pmax: 300, expn: 'valve', etaX: 0.75, srcWm: 'out', srcMw: 10, stagesA: 2, dTs: 5, src: 'water', srcT: 20, srcM: 5, srcW: 7, srcHum: 'rh', srcRH: 50, srcDew: 9, airHum: 'rh', airWg: 6.4, airDew: 7.2, Patm: 101.325, fanModel: 5, fanEta: 0.8, fanPR: 1.5, fanSh: 5, fanMotor: 0.93, Tmt: 102, phi: 0, prMax: 20, stagesH: 1, ihxH: true, RH: 60, motor: 0.93, compModel: 5, etaC: 0.73, lift2: 40, cyc: 'A', comps: [['Butane', 100]], basis: 'mass', w: '', T: 170, a: 20, Ts: 40, Tain: 15, Two: 30, lift: 80, split: 0.45, off: 0, liftS: 80, share: 0.5, PH: 25, alpha: 0.3, PHg: null, PLg: null, gKey: '', init: false};
const LV_DEF = JSON.parse(JSON.stringify(LV));            // the defaults, for cases saved before a setting existed
try { const st = JSON.parse(localStorage.getItem('dryerlive1') || '{}'); Object.assign(LV, st); if (st.srcW != null && st.srcHum == null) LV.srcRH = null; } catch (e) {}   // saved before RH input: srcRH set from srcW once CoolProp is ready
for (const [m, k] of [['airHum', 'rh'], ['srcHum', 'rh']]) if (!['rh', 'w', 'dew'].includes(LV[m])) LV[m] = k;
if (!(+LV.Patm >= 70 && +LV.Patm <= 110)) LV.Patm = 101.325;
if (+LV.compModel === 0) { LV.compModel = 1; LV.etaC = 1; }   // the ideal model is now Constant at 1.00
for (const k of ['compModel', 'fanModel']) if (+LV[k] === 4) LV[k] = 5;   // the first Pareto front (4) was removed
if (!Array.isArray(LV.comps) || !LV.comps.length) LV.comps = LV.w ? compsOf(LV.w).map(([n, f]) => [n, +f.toFixed(6)]) : [['Butane', 100]];
if (LV.basis !== 'mole') LV.basis = 'mass';
// refrigerants of the multi-loop cycles, one list per loop (bottom first), and the parameters of each
const MULTI_DEF = {K: [[['CarbonDioxide', 100]], [['CarbonDioxide', 100]]], N: [[['IsoButane', 100]], [['Cyclopentane', 100]], [['Cyclopentane', 100]], [['Cyclopentane', 100]]], M: [[['Ammonia', 100]]], X: [[['CarbonDioxide', 100]]],
  E: [[['Cyclopentane', 100]]], D: [[['Cyclopentane', 100]]], R: [[['Helium', 100]]], Y: [[['Ammonia', 50], ['Water', 50]]]};
const PAR_DEF = {
  K: {v: 3, P_b: 150, P_t: 150, r: 0.43, gc_off: 0, T_ev_b: 10, T_ev_t: 20, T_mt: 102, air_ratio: 0.5, dT_cgc: 5, dT_ihx_b: 5, dT_ihx_t: 5, Pmax: 150},
  N: {n: 2, T_cond: [100, 150, 180], cond_air: [true, true, true], sub_air: [true, true, true, true], dT_casc: 5, superheat: 5, ihx: [false, false, false, false], dT_ihx: 5},
  M: {T_cond: 90, dT_sg: 5, n_mvr: 9, inject: true, sub_air: true, ihx: false, dT_ihx: 5, split_fill: false, split: Array(19).fill(0)},
  X: {arch: 1, P_L: 80, P_H: 165, P_H_min: false, x: 0.37, eps_ihx: 0.85, eps_pre: 0.85, eps_lthx: 0.85, eta_c: 0.85, eta_t: 0.85, eta_m: 0.99, dP_hp: 50, dP_lp: 10, Pmax: 300},
  E: {arr: 0, T_int_auto: true, T_int: 100, dT_econ: 5, superheat_inj: 5, superheat: 5, Pmax: 60},
  D: {n: 3, T_cond_auto: true, T_cond: [85, 120, 145, 160], superheat: 5, Pmax: 60},
  R: {n: 1, bounds_auto: true, bounds: [60, 100, 140], eta_carnot: 0.37, dT_hot: 20, dT_cold: 5},
  Y: {w_rich: 0.5, q_des: 0.4, Pmax: 50, T_dis_max: 170, wet: true, limit_dis: false, eta_pump: 0.6},
};
const ECON_ARR = ['flash', 'economiser', 'parallel', 'single'];
if (!LV.mainComps) LV.mainComps = LV.comps;
LV.multi = Object.assign(JSON.parse(JSON.stringify(MULTI_DEF)), LV.multi || {});
LV.par = Object.assign(JSON.parse(JSON.stringify(PAR_DEF)), LV.par || {});
const savedKv = LV.par.K && LV.par.K.v;
for (const k in PAR_DEF) LV.par[k] = Object.assign(JSON.parse(JSON.stringify(PAR_DEF[k])), LV.par[k]);
if (savedKv !== PAR_DEF.K.v) LV.par.K = JSON.parse(JSON.stringify(PAR_DEF.K));   // the TT layout changed: start from the new defaults
LV.loopSel = Object.assign({K: 1, N: 0, M: 0, X: 0, E: 0, D: 0, R: 0, Y: 0}, LV.loopSel || {});
function linkComps() {
  if (MULTI(LV.cyc)) { const L = LV.multi[LV.cyc]; LV.loopSel[LV.cyc] = Math.min(LV.loopSel[LV.cyc], loopCount(LV.cyc) - 1); LV.comps = L[LV.loopSel[LV.cyc]]; }
  else LV.comps = LV.mainComps;
}
const loopCount = c => c === 'K' ? 2 : c === 'N' ? LV.par.N.n : 1;
const loopNames = c => c === 'K' ? ['Bottom', 'Top'] : c === 'N' ? Array.from({length: LV.par.N.n}, (_, i) => `Level ${i + 1}`) : c === 'X' || c === 'R' ? ['Working fluid'] : c === 'E' || c === 'D' ? ['Refrigerant'] : c === 'Y' ? ['Solution'] : ['Bottom'];
// the working fluid of loop i as text: the hybrid's solution is set by its ammonia fraction under Design
const loopFluidText = (c, v, i) => v === 'Y' ? `ammonia–water, ${(100 * c.par.Y.w_rich).toFixed(1)} % NH₃ (mass)` : v === 'R' ? 'helium (not modelled)' : fluidLabel({comps: c.multi[v][i], basis: c.basis});
const lvSave = () => { try { localStorage.setItem('dryerlive1', JSON.stringify(LV)); } catch (e) {} };
let LIVE = {res: null, cfg: null, fl: null, flKey: '', dome: null, air: null, busy: false};
let cpPromise = null;

/* ---------------------------------------------------------------- fluids offered: the study's */
function liveKeys(cyc) { return cyc === 'S' ? WS.slice() : WS.concat((D.WX || []).filter(w => !WS.includes(w))); }
function families(cyc) {
  const f = [...new Set(liveKeys(cyc).map(famOf))];
  return f.sort((a, b) => (a === 'Propane_Isopentane' ? -1 : b === 'Propane_Isopentane' ? 1 : 0));
}
const keysOfFam = (cyc, fam) => liveKeys(cyc).filter(w => famOf(w) === fam).sort((a, b) => fracOf(a)[0] - fracOf(b)[0]);
function nearestKey(cyc, fam, prop) {
  const ks = keysOfFam(cyc, fam);
  return ks.reduce((b, w) => Math.abs(fracOf(w)[0] - prop) < Math.abs(fracOf(b)[0] - prop) ? w : b, ks[0]);
}

/* ---------------------------------------------------------------- the refrigerant: any CoolProp mixture, as on the steam page;
   LV.w is the study blend it matches (mass shares within 0.01 %), or '' */
const active = c => c.comps.filter(x => x[1] > 0);
const fluidKey = c => active(c).map(x => `${x[0]}:${(+x[1]).toFixed(6)}`).join('&') + '|' + c.basis;
const fluidLabel = c => { const a = active(c); return a.length === 1 ? LABEL(a[0][0]) : `${a.map(x => LABEL(x[0])).join('/')} ${a.map(x => +(+x[1]).toFixed(2)).join('/')}${c.basis === 'mole' ? ' (mole)' : ''}`; };
const canon = list => list.filter(x => x[1] > 1e-9).map(x => `${x[0]}:${(x[1] / 100).toFixed(4)}`).sort().join('&');
let BLENDS = null;   // canonical composition -> study key
function studyKeyOf(c) {
  if (!D || c.basis !== 'mass') return '';
  if (!BLENDS) { BLENDS = new Map(); for (const w of liveKeys('A')) BLENDS.set(canon(compsOf(w)), w); }
  return BLENDS.get(canon(c.comps)) || '';
}
const syncW = () => { LV.w = studyKeyOf(LV); };
const blendLabel = w => `${famLabel(famOf(w))} ${fracLabel(w)}`;
function setComps(list) { LV.mainComps = list.map(x => [x[0], x[1]]); LV.basis = 'mass'; linkComps(); syncW(); }

/* ---------------------------------------------------------------- study rows for the live case */
const onGrid = (c = LV) => TEMPS.includes(c.T) && (D.AP || []).includes(String(c.a)) && c.Ts === 40 && c.Tain === 15;
// the MVR fans of HP + MVR: their own efficiency model, pressure ratio, exit superheat and motor
const fanCfg = o => ({PR_mvr: +o.fanPR, sh_inject: +o.fanSh, fan_model: +o.fanModel, fan_eta: +o.fanEta, fan_motor: +o.fanMotor});
const cfgKey = o => JSON.stringify([waterFlowSrc(o) ? ['wf', o.srcMw, state.vsw, state.nvMode, state.nvFix] : null, o.src === 'air' ? [o.srcT, o.srcM, o.srcHum, o.srcRH, o.srcW, o.srcDew, state.vsw, state.nvMode, state.nvFix] : null, o.cyc === 'M' ? fanCfg(o) : null, o.motor, o.compModel, o.compModel === 1 ? o.etaC : null, o.cyc, MULTI(o.cyc) ? [o.multi[o.cyc].slice(0, loopCount(o.cyc)), o.par[o.cyc], o.basis] : fluidKey(o), o.T, o.a, srcDT(o), o.Ts, o.Tain, o.RH, o.airHum, o.airWg, o.airDew, o.Patm, o.Two, o.stagesH, o.stagesA, o.ihxH, o.prMax, o.Tmt, o.phi, o.lift, o.split, o.liftS, o.share, o.PH, o.alpha]);
const snapT = T => TEMPS.reduce((b, t) => Math.abs(t - T) < Math.abs(b - T) ? t : b, TEMPS[0]);
const snapA = a => (D.AP || ['20']).map(Number).reduce((b, x) => Math.abs(x - a) < Math.abs(b - a) ? x : b, 20);
const rowsA = (w, a, T) => D.A[`${w}|${a}|${T}`] || [];
const rowsS = (w, a, T) => (D.S[`${w}|${a}|${T}`] || []).filter(r => r[1] != null);
/* the study design under a rule: 'cop' highest COP, 'p1' lowest pressure within 1 % of it */
function studyPick(cyc, w, a, T, rule) {
  if (cyc !== 'A' && cyc !== 'S') return null;
  const rows = cyc === 'A' ? rowsA(w, a, T) : rowsS(w, a, T);
  if (!rows.length) return null;
  const cop = r => cyc === 'A' ? r[2] : r[1], P = r => cyc === 'A' ? r[4] : r[3];
  const top = Math.max(...rows.map(cop));
  if (rule === 'cop') return rows.find(r => cop(r) === top);
  return rows.filter(r => cop(r) >= top * 0.99 - 1e-12).reduce((b, r) => (P(r) < P(b) - 1e-9 || (Math.abs(P(r) - P(b)) <= 1e-9 && cop(r) > cop(b)) ? r : b));
}
/* gas-cooler and suction pressure guesses (MPa) for RHP-A: the study's nearest design in this case, else the last solve */
function guessA() {
  const T = snapT(LV.T), a = snapA(LV.a), rows = rowsA(LV.w, a, T);
  const near = rows.length ? rows.reduce((b, r) => {
    const d = x => Math.abs(x[0] - LV.lift) + 100 * Math.abs(x[1] - LV.split);
    return d(r) < d(b) ? r : b;
  }) : null;
  const fromRow = r => ({PH: r[4] / 10, PL: r[7] ? r[4] / r[7] / 10 : null});
  if (LV.gKey === fluidKey(LV) && LV.PHg) return {PH: LV.PHg, PL: LV.PLg};
  return near ? fromRow(near) : {PH: null, PL: null};
}
function loadStudy(rule) {
  syncW();
  if (!LV.w) { $('lvLoadHint').textContent = 'This refrigerant is not one of the study blends.'; return; }
  const T = snapT(LV.T), a = snapA(LV.a), r = studyPick(LV.cyc, LV.w, a, T, rule);
  if (!r && LV.cyc !== 'A' && LV.cyc !== 'S') { $('lvLoadHint').textContent = 'The REFPROP study covered RHP-A and the SLTHP only.'; return; }
  if (!r) { $('lvLoadHint').textContent = `The study has no ${LV.cyc === 'A' ? 'RHP-A' : 'SLTHP'} design for this blend${LV.cyc === 'S' ? ' (the SLTHP was studied with propane/isopentane 15/85 to 30/70)' : ''}.`; return; }
  LV.T = T; LV.a = a; LV.Ts = 40; LV.Tain = 15;
  if (LV.cyc === 'A') { LV.lift = r[0]; LV.split = r[1]; LV.stagesA = 2; }
  else { LV.liftS = r[0]; LV.share = 0.5; LV.PH = r[3]; LV.alpha = r[6]; }
  LV.init = true; renderLiveControls(); changed();
}

/* ---------------------------------------------------------------- controls */
function lvSeg(id, items, cur, on) {
  const el = $(id); el.innerHTML = '';
  for (const [v, lab] of items) {
    const b = document.createElement('button'); b.type = 'button'; b.textContent = lab; b.setAttribute('aria-pressed', String(cur === v));
    b.onclick = () => on(v); el.appendChild(b);
  }
}
/* a non-linear slider: the input runs 0..1 and the value is max·u³, so the low flows get most of the travel */
const cubeScale = (max, labels, minors) => ({
  to: v => Math.cbrt(Math.max(+v, 0) / max),
  from: u => { const v = max * Math.pow(Math.min(Math.max(+u, 0), 1), 3), q = v < 2 ? 0.1 : v < 20 ? 0.5 : 1; return Math.max(0.1, Math.round(v / q) * q); },
  labels, minors,
});
const CUBE200 = cubeScale(200, [0, 1, 5, 20, 50, 200], [2, 10, 30, 40, 60, 70, 80, 90, 100, 125, 150, 175]);
const CUBE2000 = cubeScale(2000, [0, 1, 10, 100, 500, 2000], [2, 5, 20, 50, 200, 300, 400, 700, 1000, 1500]);
function ticksScaled(id, sc) {
  const box = $(`${id}-ticks`), cur = sc.to(LV[SLIDERS[id][0]]); box.innerHTML = '';
  const put = (v, lab) => {
    const left = `calc(8px + (100% - 16px) * ${sc.to(v)})`, i = document.createElement('i'); i.style.left = left; box.append(i);
    if (lab) { const s = document.createElement('span'); s.style.left = left; s.textContent = v; if (Math.abs(sc.to(v) - cur) < 1e-3) s.className = 'on'; box.append(s); }
  };
  for (const v of sc.minors) put(v, false);
  for (const v of sc.labels) put(v, true);
}
function ticksFor(id, minor, major, dec = 0) {
  const inp = $(id), box = $(`${id}-ticks`); box.innerHTML = '';
  const lo = +inp.min, hi = +inp.max, cur = +inp.value, n = Math.round((hi - lo) / minor);
  const first = Math.ceil(lo / minor - 1e-9) * minor;
  for (let k = 0; ; k++) {
    const v = +(first + k * minor).toFixed(6); if (v > hi + 1e-9 || k > n + 1) break;
    const left = `calc(8px + (100% - 16px) * ${(v - lo) / (hi - lo)})`;
    const i = document.createElement('i'); i.style.left = left; box.append(i);
    if (Math.abs(Math.round(v / major) * major - v) < 1e-9) { const s = document.createElement('span'); s.style.left = left; s.textContent = v.toFixed(dec); if (Math.abs(v - cur) < 1e-9) s.className = 'on'; box.append(s); }
  }
}
/* humidity given three ways: relative humidity (%), moisture content (g/kg dry air) or dew point (C), at the stream's
   temperature and the atmospheric pressure. humW gives the moisture content (kg/kg) for the chosen way; the read-out
   beside the active slider shows the other two, or that the air is saturated (more moisture than it can hold). */
const patm = (o = LV) => 1000 * (+o.Patm > 0 ? +o.Patm : 101.325);
const HUM = {air: {mode: 'airHum', T: 'Tain', rh: 'RH', w: 'airWg', dew: 'airDew', ids: {rh: 'lvRH', w: 'lvAirWg', dew: 'lvAirDew'}},
             src: {mode: 'srcHum', T: 'srcT', rh: 'srcRH', w: 'srcW', dew: 'srcDew', ids: {rh: 'lvSrcRH', w: 'lvSrcW', dew: 'lvSrcDew'}}};
function wSat(T, P) { try { const w = cpCall(() => CP.HAPropsSI('W', 'T', T + 273.15, 'P', P, 'R', 1)); return w > 0 && isFinite(w) ? w : Infinity; } catch (e) { return Infinity; } }
function humW(o, side) {
  const h = HUM[side], T = +o[h.T], P = patm(o), m = o[h.mode];
  if (m === 'w') return Math.max(+o[h.w] || 0, 0) / 1000;
  if (m === 'dew') { const w = wSat(+o[h.dew], P); if (!isFinite(w)) throw new Error(`dew point ${o[h.dew]} °C is above boiling at ${(P / 1000).toFixed(1)} kPa`); return w; }
  const rh = Math.min(Math.max(+o[h.rh] || 0, 0), 100) / 100;
  return rh > 0 ? cpCall(() => CP.HAPropsSI('W', 'T', T + 273.15, 'P', P, 'R', rh)) : 0;
}
/* the other two measures for a moisture content W (kg/kg) at T and P */
function humEquiv(T, P, W) {
  const Ws = wSat(T, P), out = {sat: W > Ws + 1e-7, Ws};
  const Wc = Math.min(W, Ws);
  try { out.rh = W <= 0 ? 0 : 100 * cpCall(() => CP.HAPropsSI('R', 'T', T + 273.15, 'P', P, 'W', Wc)); } catch (e) {}
  try { out.dew = W <= 0 ? null : cpCall(() => CP.HAPropsSI('D', 'T', T + 273.15, 'P', P, 'W', Wc)) - 273.15; } catch (e) {}
  return out;
}
function humNote(side, o = LV) {
  if (typeof CP === 'undefined' || !CP) return '';
  if (side === 'src' && o.src !== 'air') return '';
  try {
    const h = HUM[side], T = +o[h.T], P = patm(o), W = humW(o, side), m = o[h.mode], e = humEquiv(T, P, W);
    if (e.sat) return ` · saturated: ${(1000 * e.Ws).toFixed(1)} g/kg max at ${T.toFixed(1)} °C`;
    const parts = [];
    if (m !== 'rh' && e.rh != null) parts.push(`${e.rh.toFixed(0)} %`);
    if (m !== 'w') parts.push(`${(1000 * W).toFixed(1)} g/kg`);
    if (m !== 'dew') parts.push(W > 0 && e.dew != null ? `dew ${e.dew.toFixed(1)} °C` : 'dry');
    return parts.length ? ` · ${parts.join(' · ')}` : '';
  } catch (err) { return ` · ${cpError(err)}`; }
}
function refreshHumOuts() {
  for (const id of ['lvRH', 'lvAirWg', 'lvAirDew', 'lvSrcRH', 'lvSrcW', 'lvSrcDew']) if ($(`${id}-out`)) sliderOut(id);
  // the g/kg scales end at saturation, which moves with the air temperature and pressure (and is known once CoolProp loads)
  for (const id of ['lvAirWg', 'lvSrcW']) if ($(id) && document.activeElement !== $(id)) { const sc = SLIDERS[id][5]; $(id).value = sc.to(LV[SLIDERS[id][0]]); ticksScaled(id, sc); }
}
/* the g/kg sliders run from dry air (0) to saturation at the stream's temperature and the atmospheric pressure, on a
   square-root scale so the low moisture contents get most of the travel; the end of the track is saturated air */
function satScale(side) {
  const h = HUM[side], fallback = side === 'src' ? 100 : 50, r1 = x => Math.round(10 * x) / 10;
  const max = () => {
    if (typeof CP === 'undefined' || !CP) return fallback;            // until CoolProp loads
    const w = 1000 * wSat(+LV[h.T], patm(LV));
    return w > 0 ? Math.min(w, 1000) : fallback;                       // near boiling the air could hold almost any amount
  };
  return {
    to: v => Math.sqrt(Math.min(Math.max(+v, 0) / max(), 1)),
    from: u => {
      const m = max(), v = m * Math.pow(Math.min(Math.max(+u, 0), 1), 2), q = m < 50 ? 0.1 : m < 200 ? 0.5 : 1;
      return Math.min(+(Math.round(v / q) * q).toFixed(1), r1(m));
    },
    get labels() {             // 0, round values spread along the track, and the saturation value at its end
      const m = max(), out = [0];
      for (const x of [0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500]) {
        const pos = Math.sqrt(x / m), last = Math.sqrt(out[out.length - 1] / m);
        if (pos - last >= 0.16 && 1 - pos >= 0.16) out.push(x);
      }
      out.push(r1(m));
      return out;
    },
    get minors() { const m = max(); return [0.25, 0.5, 0.75].map(f => r1(m * f * f)); },
  };
}
/* switching the way humidity is given keeps the same moisture: the newly shown value is set from the current one */
function setHumMode(side, m) {
  const h = HUM[side];
  if (LV[h.mode] === m) return;
  if (typeof CP !== 'undefined' && CP) try {
    const T = +LV[h.T], P = patm(LV), W = humW(LV, side), e = humEquiv(T, P, W), Wc = Math.min(W, e.Ws);
    if (m === 'rh') LV[h.rh] = e.sat ? 100 : Math.round(e.rh);
    if (m === 'w') LV[h.w] = +(1000 * Wc).toFixed(1);
    if (m === 'dew' && e.dew != null) LV[h.dew] = +e.dew.toFixed(1);
  } catch (e) {}
  LV[h.mode] = m; renderLiveControls(); changed();
}
const SLIDERS = {
  lvT: ['T', 5, 20, 0, v => `${+(+v).toFixed(1)} °C`], lvTain: ['Tain', 2.5, 10, 0, v => `${+(+v).toFixed(1)} °C`], lvRH: ['RH', 5, 25, 0, v => `${+(+v).toFixed(0)} %${humNote('air')}`], lvAirWg: ['airWg', 1, 10, 0, v => `${+(+v).toFixed(1)} g/kg${humNote('air')}`, satScale('air')], lvAirDew: ['airDew', 5, 10, 0, v => `${+(+v).toFixed(1)} °C${humNote('air')}`], lvPatm: ['Patm', 1, 10, 0, v => `${(+v).toFixed(1)} kPa`], lvDTs: ['dTs', 2.5, 10, 0, v => `${+(+v).toFixed(1)} K`], lvTs: ['Ts', 5, 10, 0, v => `${+(+v).toFixed(1)} °C`], lvTwo: ['Two', 5, 10, 0, v => `${+(+v).toFixed(1)} °C`],
  lva: ['a', 2.5, 10, 0, v => `${+(+v).toFixed(1)} K`],
  lvLift: ['lift', 5, 25, 0, v => `${v} K`], lvLift2: ['lift2', 5, 25, 0, v => `${v} K`], lvSplit: ['split', 0.05, 0.1, 2, v => v.toFixed(2)],
  lvLiftS: ['liftS', 10, 20, 0, v => `${v} K`], lvShare: ['share', 0.05, 0.1, 2, v => v.toFixed(2)],
  lvEtaC: ['etaC', 0.01, 0.05, 2, v => (+v).toFixed(2)],
  lvTmt: ['Tmt', 10, 40, 0, v => `${+(+v).toFixed(1)} °C`], lvPhi: ['phi', 0.1, 0.5, 1, v => (+v).toFixed(2)],
  lvPrMax: ['prMax', 1, 5, 0, v => `${+(+v).toFixed(0)}`],
  lvMotor: ['motor', 0.01, 0.05, 2, v => (+v).toFixed(2)], lvPmax: ['pmax', 25, 50, 0, v => `${+(+v).toFixed(0)} bar`], lvEtaX: ['etaX', 0.05, 0.1, 1, v => (+v).toFixed(2)],
  lvSrcT: ['srcT', 5, 20, 0, v => `${+(+v).toFixed(1)} °C`], lvSrcM: ['srcM', 1, 5, 0, v => `${+(+v).toFixed(+v < 10 ? 1 : 0)} kg/s`, CUBE200], lvSrcMw: ['srcMw', 1, 5, 0, v => `${+(+v).toFixed(+v < 10 ? 1 : 0)} kg/s`, CUBE2000],
  lvSrcRH: ['srcRH', 5, 25, 0, v => v == null ? '–' : `${+(+v).toFixed(0)} %${humNote('src')}`], lvSrcW: ['srcW', 5, 20, 0, v => `${+(+v).toFixed(1)} g/kg${humNote('src')}`, satScale('src')], lvSrcDew: ['srcDew', 5, 20, 0, v => `${+(+v).toFixed(1)} °C${humNote('src')}`],
  lvFanEta: ['fanEta', 0.01, 0.05, 2, v => (+v).toFixed(2)], lvFanPR: ['fanPR', 0.05, 0.25, 2, v => (+v).toFixed(2)],
  lvFanSh: ['fanSh', 1, 5, 0, v => `${+(+v).toFixed(0)} K`], lvFanMotor: ['fanMotor', 0.01, 0.05, 2, v => (+v).toFixed(2)],
};
function sliderOut(id) { const [k, , , , f] = SLIDERS[id]; $(`${id}-out`).textContent = f(LV[k]); }
function renderLiveControls() {
  if (!D) return;
  linkComps(); syncW();
  LV.etaC = Math.min(Math.max(+LV.etaC || 0.73, 0.5), 1);
  $('lvEtaCRow').hidden = LV.compModel !== 1;
  lvSeg('lvComp', DRYER.COMP_MODELS.map((n, i) => [i, n]).slice(1).filter(([, n]) => n), LV.compModel, i => { LV.compModel = i; renderLiveControls(); changed(); });
  LV.motor = Math.min(Math.max(+LV.motor || 0.93, 0.8), 1);
  LV.expn = LV.expn === 'exp' ? 'exp' : 'valve'; LV.etaX = Math.min(Math.max(+LV.etaX || 0.75, 0.3), 1);
  const brayton = LV.cyc === 'X', noExp = ['X', 'E', 'D', 'R', 'Y'].includes(LV.cyc);
  if (LV.cyc === 'S') { lvSeg('lvStagesS', [[1, 'Single stage'], [2, 'Two stages']], stagesOf(LV), v => { LV.stagesS = v; renderLiveControls(); changed(); }); $('lvShareRow').hidden = stagesOf(LV) === 1; }
  LV.pmax = Math.min(Math.max(+LV.pmax || 300, 50), 400); $('lvPmaxRow').hidden = MULTI(LV.cyc);
  lvSeg('lvExp', [['valve', 'Expansion valve'], ['exp', 'Expander']], LV.expn, v => { LV.expn = v; renderLiveControls(); changed(); });
  $('lvExpRow').hidden = noExp; $('lvEtaXRow').hidden = noExp || LV.expn !== 'exp';
  $('lvExpHint').textContent = brayton ? 'The sCO₂ Brayton cycle always expands through its turbine; its efficiency is set under Design.'
    : LV.cyc === 'R' ? 'Stirling units have no compressor or valve: their performance is a fraction of the Carnot COP, set under Design.'
    : noExp ? 'This cycle expands through valves (no expander option).'
    : LV.expn === 'exp' ? `Every expansion valve becomes an expander of this isentropic efficiency; its work is returned as electricity at the motor efficiency. A two-phase outlet (the usual case below the critical point) is flagged: two-phase expanders are less proven than gas turbines.`
    : 'Isenthalpic throttling through an expansion valve.';
  const air = LV.src === 'air';
  lvSeg('lvSrc', [['water', 'Water'], ['air', 'Air']], air ? 'air' : 'water', v => { LV.src = v; renderLiveControls(); changed(); });
  const wFlow = !air && LV.srcWm === 'flow';
  $('lvTsRow').hidden = air; $('lvTwoRow').hidden = air || wFlow; $('lvSrcMwRow').hidden = !wFlow; $('lvSrcWmRow').hidden = air;
  lvSeg('lvSrcWm', [['out', 'Set water outlet'], ['flow', 'Set water flow']], wFlow ? 'flow' : 'out', v => { LV.srcWm = v; renderLiveControls(); changed(); });
  for (const e of document.querySelectorAll('.src-air')) e.hidden = !air;
  const HM = [['rh', 'RH %'], ['w', 'g/kg'], ['dew', 'Dew °C']];
  lvSeg('lvAirHum', HM, LV.airHum, m => setHumMode('air', m));
  lvSeg('lvSrcHum', HM, LV.srcHum, m => setHumMode('src', m));
  for (const side of ['air', 'src']) for (const [m, id] of Object.entries(HUM[side].ids)) $(`${id}Row`).hidden = (side === 'src' && !air) || LV[HUM[side].mode] !== m;
  $('lvSrcHint').textContent = air
    ? `The source air cools through the evaporator to the temperature at which its enthalpy drop carries the evaporator duty at the swept volume; below its dew point its moisture condenses. The evaporator and superheater stay ${srcDT(LV)} K below the air (counterflow); the superheater outlet is ${srcDT(LV)} K below the air inlet.`
    : `The evaporator pressure is the highest at which the evaporator and superheater stay ${srcDT(LV)} K below the source water (counterflow); the superheater outlet is ${srcDT(LV)} K below the water inlet.${LV.srcWm === 'flow' ? ' The water flow is held constant; its outlet temperature follows from the evaporator duty.' : ''}`;
  $('lvFanSec').hidden = LV.cyc !== 'M';  $('lvFanEtaRow').hidden = +LV.fanModel !== 1;
  lvSeg('lvFanModel', DRYER.COMP_MODELS.map((n, i) => [i, n]).slice(1).filter(([, n]) => n), +LV.fanModel, i => { LV.fanModel = i; renderLiveControls(); changed(); });
  const mu = MULTI(LV.cyc);
  $('lvLoopRow').hidden = !mu || loopCount(LV.cyc) < 2; $('lvLoopHint').hidden = LV.cyc !== 'M' && !NOFLUID(LV.cyc);
  for (const e of [$('lvComps'), $('lvAdd').parentElement]) e.hidden = NOFLUID(LV.cyc);
  if (mu) lvSeg('lvLoop', loopNames(LV.cyc).map((n, i) => [i, n]), LV.loopSel[LV.cyc], i => { LV.loopSel[LV.cyc] = i; openPicker = -1; renderLiveControls(); });
  $('lvLoopHint').textContent = LV.cyc === 'R' ? 'Stirling units run on helium; the model does not need its properties.' : LV.cyc === 'Y' ? 'The hybrid cycle runs on an ammonia–water solution (IAPWS G4-01 properties, solved in the page); set its ammonia fraction under Design.' : 'The HP + MVR cycle\'s steam loop is water.';
  lvSeg('lvCyc', [['H', 'Conventional'], ['A', 'RHP-A'], ['B', 'RHP-B'], ['C', 'RHP-C'], ['S', 'SLTHP'], ['K', 'TT cascade'], ['N', 'n-cascade'], ['M', 'HP + MVR'], ['X', 'sCO₂ Brayton'], ['E', 'Flash / economiser'], ['D', 'n-pressure'], ['Y', 'NH₃–H₂O hybrid'], ['R', 'Stirling']], LV.cyc, v => { LV.cyc = v; renderLiveControls(); changed(); });
  $('lvMini').innerHTML = RHP(LV.cyc) ? miniRHP(LV.cyc, stagesOf(LV)) : mu ? miniMulti(LV.cyc) : miniSchematic(LV.cyc);
  compRows();
  for (const b of document.querySelectorAll('.tbtn[data-basis]')) b.setAttribute('aria-pressed', String(b.dataset.basis === LV.basis));
  const studied = LV.w && (LV.cyc === 'A' || (LV.cyc === 'S' && WS.includes(LV.w)));
  if (LV.Two > LV.Ts - 1) LV.Two = LV.Ts - 1;
  // the 2/3 rule on the total recuperator lift
  const liftMax = Math.floor((LV.T - LV.Tain) * 2 / 3);
  const l1 = LV.cyc === 'B' ? 0 : LV.lift, l2 = LV.cyc === 'A' ? 0 : LV.lift2;
  if (l1 + l2 > liftMax) { if (LV.cyc === 'B') LV.lift2 = liftMax; else if (LV.cyc === 'A') LV.lift = liftMax; else LV.lift2 = Math.max(liftMax - LV.lift, 0); }
  $('lvLift').max = LV.cyc === 'C' ? liftMax - Math.min(LV.lift2, liftMax) : liftMax; $('lvLift2').max = LV.cyc === 'C' ? liftMax - Math.min(LV.lift, liftMax) : liftMax;
  LV.lift = Math.min(LV.lift, +$('lvLift').max); LV.lift2 = Math.min(LV.lift2, +$('lvLift2').max);
  for (const id in SLIDERS) { const [k, minor, major, dec, , sc] = SLIDERS[id]; $(id).value = sc ? sc.to(LV[k]) : LV[k]; if (sc) ticksScaled(id, sc); else ticksFor(id, minor, major, dec); sliderOut(id); }
  $('lvAlpha').value = +LV.alpha.toFixed(4);
  $('lvDesA').hidden = !RHP(LV.cyc); $('lvDesS').hidden = LV.cyc !== 'S'; $('lvDesX').hidden = !mu;
  if (mu) renderMultiFields();
  document.querySelector('.row.load').hidden = mu; $('lvLoadHint').hidden = mu;
  $('lvLiftRow').hidden = LV.cyc === 'B' || LV.cyc === 'H'; $('lvLift2Row').hidden = LV.cyc === 'A' || LV.cyc === 'H';
  $('lvStagesRow').hidden = LV.cyc !== 'H' && LV.cyc !== 'A'; $('lvSplitRow').hidden = stagesOf(LV) === 1;
  $('lvIhxRow').hidden = LV.cyc !== 'H';
  if (LV.cyc === 'H' || LV.cyc === 'A') lvSeg('lvStages', [[1, 'Single stage'], [2, 'Two stages']], stagesOf(LV), v => { if (LV.cyc === 'H') LV.stagesH = v; else LV.stagesA = v; renderLiveControls(); changed(); });
  if (LV.cyc === 'H') {
    lvSeg('lvIhx', [[0, 'None'], [1, 'IHX']], LV.ihxH ? 1 : 0, v => { LV.ihxH = !!v; renderLiveControls(); changed(); }); }
  $('lvDesAHint').textContent = (LV.cyc === 'H'
    ? `No recuperators: the ${LV.stagesH === 1 ? 'compressor discharges' : 'two compressors in series discharge'} straight into the gas coolers, and the evaporator and superheater take the source water.${LV.ihxH ? ' The IHX heats the suction with the gas-cooler outlet when that is hotter than the source water.' : ''} The gas-cooler pressure is solved (the lowest at which the gas-cooler pinch reaches the approach) and the gas-cooler outlet offset is optimised for COP; an offset also keeps compression dry.`
    : 'The gas-cooler pressure is solved (the lowest at which the gas-cooler pinch reaches the approach) and the gas-cooler outlet offset is optimised for COP: the gas-cooler outlet sits at air inlet + approach + offset, and an offset keeps stage 1 out of wet compression.')
    + ` When the heater outlet (source inlet less the source approach) is above the refrigerant\'s critical temperature, as for CO₂, the low side has no dew point to fix its pressure: the evaporator becomes a gas heater, and its pressure is chosen for COP, no higher than the source pinch allows, with the gas cooler at ${pmaxOf(LV)} bar at most.`;
  const T = snapT(LV.T), a = snapA(LV.a), snapped = T !== LV.T || a !== LV.a || LV.Ts !== 40 || LV.Tain !== 15;
  $('lvLoadHint').textContent = !(LV.cyc === 'A' || LV.cyc === 'S') ? 'The REFPROP study covered RHP-A and the SLTHP only.' : !studied ? 'Loading a study design needs a study blend.'
    : `Loads the study's ${LV.cyc === 'A' ? 'RHX1 lift and stage-2 power fraction' : 'RHX1 lift and RHX1 flow share'} for ${blendLabel(LV.w)} at ${T} °C and ${a} K${snapped ? ' (the nearest study case: air inlet 15 °C, source water in 40 °C)' : ''}.`;
  for (const b of document.querySelectorAll('[data-load]')) b.disabled = !studied;
  lvTitle();
  lvSave();
}
/* Fluid choice uses page buttons rather than a native drop-down, which some embedded viewers do not open. */
let openPicker = -1, blendsOpen = false;
function balance() {
  const c = LV.comps;
  if (c.length === 1) { c[0][1] = 100; return; }
  let others = c.slice(0, -1).reduce((a, x) => a + x[1], 0);
  if (others > 100) { const f = 100 / others; c.slice(0, -1).forEach(x => x[1] *= f); others = 100; }
  c[c.length - 1][1] = 100 - others;
}
function fluidChanged() { balance(); if (!MULTI(LV.cyc)) LV.mainComps = LV.comps; syncW(); renderLiveControls(); changed(); }
function compRows() {
  const box = $('lvComps'); box.innerHTML = '';
  const comps = LV.comps;
  comps.forEach((c, i) => {
    const last = i === comps.length - 1;
    const row = document.createElement('div'); row.className = 'comp';
    const sel = document.createElement('button'); sel.type = 'button'; sel.className = 'pick';
    sel.setAttribute('aria-label', `Component ${i + 1}: ${LABEL(c[0])}. Change`); sel.setAttribute('aria-expanded', String(openPicker === i));
    sel.innerHTML = `<span>${LABEL(c[0])}</span><span class="caret" aria-hidden="true">▾</span>`;
    sel.addEventListener('click', () => { openPicker = openPicker === i ? -1 : i; compRows(); });
    const rng = document.createElement('input'); rng.type = 'range'; rng.min = 0; rng.max = 100; rng.step = 0.5; rng.value = c[1];
    rng.setAttribute('aria-label', `Component ${i + 1} percent`);
    const num = document.createElement('input'); num.type = 'number'; num.min = 0; num.max = 100; num.step = 0.5; num.value = +(+c[1]).toFixed(2);
    num.setAttribute('aria-label', `Component ${i + 1} percent`);
    const pct = document.createElement('span'); pct.className = 'unit'; pct.textContent = '%';
    if (last) { rng.disabled = num.disabled = comps.length > 1; num.title = rng.title = 'Balance: 100 % minus the other components'; }
    const set = v => { c[1] = Math.max(0, Math.min(100, +v || 0)); fluidChanged(); };
    rng.addEventListener('input', () => { num.value = rng.value; });
    rng.addEventListener('change', () => set(rng.value));
    num.addEventListener('change', () => set(num.value));
    const rm = document.createElement('button'); rm.type = 'button'; rm.className = 'icon'; rm.textContent = '×'; rm.title = 'Remove component';
    rm.setAttribute('aria-label', `Remove ${LABEL(c[0])}`); rm.hidden = comps.length === 1;
    rm.addEventListener('click', () => { comps.splice(i, 1); openPicker = -1; fluidChanged(); });
    row.append(sel, rng, num, pct, rm); box.append(row);
    if (openPicker === i) {
      const list = document.createElement('div'); list.className = 'fluid-list'; list.setAttribute('role', 'listbox');
      list.setAttribute('aria-label', `Choose component ${i + 1}`);
      const used = new Set(comps.map((x, k) => (k === i ? null : x[0])));
      for (const n of FLUIDS) {
        const b = document.createElement('button'); b.type = 'button'; b.className = 'fopt'; b.textContent = LABEL(n);
        b.setAttribute('role', 'option'); b.setAttribute('aria-selected', String(n === c[0]));
        if (used.has(n)) { b.disabled = true; b.title = 'Already in the mixture'; }
        b.addEventListener('click', () => { c[0] = n; openPicker = -1; fluidChanged(); });
        list.append(b);
      }
      box.append(list);
    }
  });
  $('lvAdd').hidden = comps.length >= 3;
}
function bindLive() {
  $('lvAdd').addEventListener('click', () => {
    const used = new Set(LV.comps.map(c => c[0])); const next = FLUIDS.find(n => !used.has(n));
    LV.comps.push([next, 0]); const n = LV.comps.length; LV.comps.forEach(c => c[1] = 100 / n);
    openPicker = n - 1; fluidChanged();
  });
  for (const b of document.querySelectorAll('.tbtn[data-basis]')) b.addEventListener('click', () => { LV.basis = b.dataset.basis; fluidChanged(); });
  document.addEventListener('keydown', e => { if (e.key === 'Escape' && openPicker >= 0) { openPicker = -1; compRows(); } });
  for (const id in SLIDERS) {
    const [k, , , , , sc] = SLIDERS[id], val = () => sc ? sc.from($(id).value) : +$(id).value;
    $(id).addEventListener('input', () => { LV[k] = val(); sliderOut(id); });
    $(id).addEventListener('change', () => { LV[k] = val(); renderLiveControls(); changed(); });
    // the value beside each slider can be typed: click it (or focus it and press Enter), type, Enter to apply
    const out = $(`${id}-out`);
    out.tabIndex = 0; out.classList.add('editable'); out.title = 'Click to type a value';
    const edit = () => {
      if (out.hidden) return;
      const inp = document.createElement('input'), lo = sc ? sc.from(0) : +$(id).min, hi = sc ? sc.from(1) : +$(id).max;
      inp.type = 'number'; inp.step = 'any'; inp.className = 'outedit'; inp.value = +(+LV[k]).toFixed(4); inp.setAttribute('aria-label', `${out.previousElementSibling ? out.previousElementSibling.textContent : k} (${lo} to ${hi})`);
      inp.title = `${lo} to ${hi}`;
      out.hidden = true; out.after(inp); inp.focus(); inp.select();
      let done = false;
      const finish = ok => {
        if (done) return; done = true;
        const v = +inp.value; inp.remove(); out.hidden = false;
        if (ok && inp.value !== '' && isFinite(v)) { LV[k] = Math.min(Math.max(v, lo), hi); renderLiveControls(); changed(); } else sliderOut(id);
        out.focus();
      };
      inp.addEventListener('keydown', e => { if (e.key === 'Enter') { e.preventDefault(); finish(true); } else if (e.key === 'Escape') { e.preventDefault(); finish(false); } });
      inp.addEventListener('blur', () => finish(true));
    };
    out.addEventListener('click', edit);
    out.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); edit(); } });
  }
  const num = (id, k, lo, hi) => $(id).addEventListener('change', () => { const v = +$(id).value; if (isFinite(v)) LV[k] = Math.min(Math.max(v, lo), hi); renderLiveControls(); changed(); });
  num('lvAlpha', 'alpha', 0, 0.95);
  for (const b of document.querySelectorAll('[data-load]')) b.addEventListener('click', () => loadStudy(b.dataset.load));
}

/* ---------------------------------------------------------------- CoolProp and solving */
function ensureCoolProp() {
  if (cpPromise) return cpPromise;
  STATUS('Loading CoolProp…'); $('status').className = 'chip';
  cpPromise = (async () => {
    CP = await loadCoolProp();
    try { CP_VERSION = CP.get_global_param_string('version'); } catch (_) { CP_VERSION = ''; }
    airFor(LV);
    if (LV.srcRH == null) { try { srcAirModel(LV); } catch (e) { LV.srcRH = 50; } lvSave(); }
    refreshHumOuts();
    STATUS(`CoolProp ${CP_VERSION} ready`); $('status').className = 'chip info';
  })().catch(e => { $('status').className = 'chip warn'; STATUS(`CoolProp could not load: ${cpError(e)}`); cpPromise = null; throw e; });
  return cpPromise;
}
/* humid air at the inlet: absolute humidity from the inlet temperature and relative humidity (CoolProp HAPropsSI,
   101.325 kPa); enthalpy per kg of dry air, so every air flow is per kg of dry air */
const airKey = (o = LV) => `${(+o.Tain).toFixed(2)}|${patm(o)}|${o.airHum}|${(+o.RH).toFixed(2)}|${(+o.airWg).toFixed(3)}|${(+o.airDew).toFixed(3)}`;
/* secondary air: primary and secondary air are heated together from the inlet to the secondary target, where the
   secondary leaves; the primary goes on to the dryer target. Per kg of primary air the pair is one stream whose
   enthalpy rises (1 + phi) times as fast below the secondary target, so every model takes it as its air (counterflow
   along the air path) and its air flows are kg of primary (dryer) air. The TT cascade models the split itself. */
const secOn = o => +o.phi > 0 && +o.Tmt > +o.Tain && +o.Tmt < +o.T;
function airFor(o) {
  const base = airBase(o);
  if (!secOn(o)) { LIVE.air = base; LIVE.airSec = null; return base; }
  const k = `${airKey(o)}|${(+o.Tmt).toFixed(2)}|${(+o.phi).toFixed(4)}`;
  if (LIVE.airSec && LIVE.airSec.key === k) { LIVE.air = LIVE.airSec.air; return LIVE.air; }
  const comp = DRYER.airSecondary(base, +o.Tain, +o.Tmt, +o.phi, +o.T);
  LIVE.airSec = {key: k, air: comp}; LIVE.air = comp;
  return comp;
}
// share of the heat that goes to the secondary air (per kg of primary air), from the humid-air enthalpy
function secShare(o) {
  return LIVE.airBase ? DRYER.secondaryShare(LIVE.airBase, +o.Tain, +o.Tmt, +o.T, +o.phi) : 0;
}
function airBase(o) {
  const k = airKey(o);
  if (LIVE.airBase && LIVE.airKey === k) return LIVE.airBase;
  const P = patm(o), Win = humW(o, 'air'), W = Math.min(Win, wSat(+o.Tain, P));   // above saturation: taken as saturated
  LIVE.airBase = DRYER.airModel(T => cpCall(() => CP.HAPropsSI('H', 'T', T + 273.15, 'P', P, 'W', W)) / 1000);
  LIVE.airBase.W = W; LIVE.airW = W; LIVE.airKey = k; LIVE.airSec = null;
  refreshHumOuts();
  return LIVE.airBase;
}
/* ---------------------------------------------------------------- the heat source: water, or ambient/exhaust air
   Water: in at Ts, out at Two (linear). Air: in at srcT with moisture srcW (g/kg dry air) and mass flow srcM (kg/s,
   moist); its outlet depends on the evaporator duty at the swept volume, so the design is solved on the outlet
   temperature until the duty cools the air to it. Below the dew point the air is saturated and the water vapour
   it loses condenses (leaving as liquid at the air temperature). */
function srcAirModel(o) {
  if (o.srcHum === 'rh' && o.srcRH == null) {    // saved before the RH input: carry the old moisture over as RH
    const e = humEquiv(+o.srcT, patm(o), Math.max(+o.srcW || 0, 0) / 1000); o.srcRH = e.sat ? 100 : Math.round(e.rh || 0);
    if (o !== LV && LV.srcRH == null) LV.srcRH = o.srcRH;
  }
  const key = `${(+o.srcT).toFixed(3)}|${patm(o)}|${o.srcHum}|${(+o.srcRH).toFixed(2)}|${(+o.srcW).toFixed(4)}|${(+o.srcDew).toFixed(3)}`;
  if (LIVE.srcAir && LIVE.srcAir.key === key) return LIVE.srcAir;
  const P = patm(o), Tin = +o.srcT;
  const Wsat = T => { try { const w = cpCall(() => CP.HAPropsSI('W', 'T', T + 273.15, 'P', P, 'R', 1)); return w > 0 && isFinite(w) ? w : Infinity; } catch (e) { return Infinity; } };
  // more moisture than the air can hold at its temperature would be fog: the air is taken as saturated
  const Win = humW(o, 'src'), WsIn = Wsat(Tin), W = Math.min(Win, WsIn), capped = Win > WsIn + 1e-7;
  const hAt = T => {
    let Ws = Infinity;          // saturation moisture: unbounded as the air nears boiling at P
    try { Ws = cpCall(() => CP.HAPropsSI('W', 'T', T + 273.15, 'P', P, 'R', 1)); } catch (e) {}
    if (!(Ws > 0) || !isFinite(Ws)) Ws = Infinity;
    if (W <= Ws) return cpCall(() => CP.HAPropsSI('H', 'T', T + 273.15, 'P', P, 'W', W)) / 1000;
    return cpCall(() => CP.HAPropsSI('H', 'T', T + 273.15, 'P', P, 'R', 1)) / 1000 + (W - Ws) * 4.186 * T;
  };
  const Ts = [], Hs = [];
  for (let T = Tin - 60; T <= Tin + 1e-9; T += 0.5) { Ts.push(T); Hs.push(hAt(T)); }
  if (Ts[Ts.length - 1] < Tin) { Ts.push(Tin); Hs.push(hAt(Tin)); }
  const interp = (xs, ys, x) => {
    if (x <= xs[0]) return ys[0] + (x - xs[0]) * (ys[1] - ys[0]) / (xs[1] - xs[0]);
    let i = 1; while (i < xs.length - 1 && xs[i] < x) i++;
    return ys[i - 1] + (x - xs[i - 1]) * (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1]);
  };
  let Td = null; try { Td = capped ? Tin : cpCall(() => CP.HAPropsSI('D', 'T', Tin + 273.15, 'P', P, 'W', W)) - 273.15; } catch (e) {}
  LIVE.srcAir = {key, W, Win, capped, Tin, Td, h: T => interp(Ts, Hs, T), T: h => interp(Hs, Ts, h), Wout: T => Math.min(W, Wsat(T))};
  return LIVE.srcAir;
}
/* solver inputs from a configuration: the dryer air, approach, motor and the source (water line or air curve) */
/* the source-side minimum approach (K): evaporator and superheater to the water or air, 0-40 K, 5 K before it was an input */
const srcDT = o => Math.min(Math.max(o && o.dTs != null && isFinite(+o.dTs) ? +o.dTs : 5, 0), 40);
const pmaxOf = o => Math.min(Math.max(+o.pmax || 300, 50), 400);   // bar: the high-side limit of the RHP cycles and the SLTHP
function baseOf(cfg) {
  const dTs = srcDT(cfg);
  const b = {T_air_tar: cfg.T, dT_GC_min: cfg.a, T_source_in: cfg.Ts, T_air_inlet: cfg.Tain, T_source_out: cfg.Two, motor_eff: cfg.motor, dT_src: dTs, dT_ev: dTs,
    expander: cfg.expn === 'exp' && cfg.cyc !== 'X' ? {eta: Math.min(Math.max(+cfg.etaX || 0.75, 0.3), 1)} : null,
    P_max: MULTI(cfg.cyc) ? undefined : pmaxOf(cfg) / 10};
  if (cfg.src === 'air' && cfg.srcTo != null) {
    const A = srcAirModel(cfg), hi = A.h(+cfg.srcT), ho = A.h(+cfg.srcTo);
    Object.assign(b, {T_source_in: +cfg.srcT, T_source_out: +cfg.srcTo, source_name: 'source air', source_curve: f => A.T(ho + f * (hi - ho))});
  }
  return b;
}
/* the source, for the headline: water flow and temperatures, or the source air's flow, outlet and condensate */
function srcOf(c, r) {           // for the diagrams: the source air curve and its result
  const a = r && r.success && r.srcAir;
  return a ? {srcAir: a, srcCurve: baseOf({...c, srcTo: a.Tout}).source_curve} : {};
}
function srcLine(c, r, kw, water) {
  const a = r && r.srcAir;
  if (a) return `${hfmt(kw, 1)} kW from ${a.m.toFixed(1)} kg/s of source air, ${hfmt(a.Tin, 1)} → ${hfmt(a.Tout, 1)} °C, ${(1000 * a.W).toFixed(1)} → ${(1000 * a.Wout).toFixed(1)} g/kg${a.cond > 1e-6 ? ` (${(a.cond * 3600).toFixed(1)} kg/h condensate)` : ''}`;
  return `${hfmt(kw, 1)} kW from ${water != null ? `${(water * 3.6).toFixed(1)} t/h of ` : ''}source water, ${hfmt(c.Ts, 1)} → ${hfmt(c.Two, 1)} °C`;
}
/* duties (kW) of any cycle's result at the page's swept volume */
function kwAt(r, cfg) {
  if (!r || !r.success) return null;
  if (MULTI(cfg.cyc)) { const u = energyG(r, cfg); return u ? u.kw : null; }
  const rep = r.report, d = {paths: r.paths, rho1: rep.rho_suction_1, rho2: rep.rho_suction_2, motor: cfg.motor, alpha: cfg.cyc === 'S' ? rep.alpha_RHX1 : undefined};
  const u = energy(RHP(cfg.cyc) ? 'A' : 'S', d, state.vsw, cfg.T, cfg.Tain);
  return u ? u.kw : null;
}
/* heat (kJ) per kg of dryer air from the inlet to the dryer temperature: humid air, with any secondary air carried
   per kg of primary air. The page's air is swapped in for the call and restored after it. */
function airDh(cfg) {
  if (!cfg || typeof CP === 'undefined' || !CP) return null;
  const keep = {air: LIVE.air, airBase: LIVE.airBase, airKey: LIVE.airKey, airSec: LIVE.airSec, airW: LIVE.airW};
  try { const A = airFor(cfg), dh = A.h(+cfg.T) - A.h(+cfg.Tain); return dh > 0 ? dh : null; }
  catch (e) { return null; }
  finally { Object.assign(LIVE, keep); }
}
/* the heat to air (kW) the sizing asks for: the sink duty, or the sink flow times its heat per kg; null at a fixed
   swept volume */
function targetDuty(cfg) {
  if (state.mode !== 'live') return null;
  if (state.sizeBy === 'duty') return state.duty;
  if (state.sizeBy === 'flow') { const dh = airDh(cfg); return dh ? state.flow * dh : null; }
  return null;
}
/* kW drawn from the source: at the swept volume, or scaled to the target sink duty or flow */
function srcKW(r, cfg) {
  const kw = kwAt(r, cfg); if (!kw) return null;
  const t = targetDuty(cfg);
  return t != null && kw.air > 0 ? kw.ev * t / kw.air : kw.ev;
}
/* sizing by sink duty or flow: set the swept volume that gives the target heat to air for the current design */
function applyDuty() {
  if (state.sizeBy === 'vsw' || state.mode !== 'live' || !LIVE.res || !LIVE.res.success || !LIVE.cfg) return;
  const kw = kwAt(LIVE.res, LIVE.cfg), t = targetDuty(LIVE.cfg);
  if (kw && kw.air > 0 && t != null) state.vsw = state.vsw * t / kw.air;
  if (typeof renderSize === 'function') renderSize();
}
/* solve with the chosen source; for air, a secant on the air outlet temperature (the duty falls as the outlet,
   and with it the evaporating temperature, falls, so there is one crossing) */
const WATER_SRC = {h: T => 4.18 * T, T: h => h / 4.18, W: 0, Wout: () => 0, capped: false, Td: null};
const waterFlowSrc = cfg => cfg.src !== 'air' && cfg.srcWm === 'flow';
function withSource(cfg, solveOne) {
  const wf = waterFlowSrc(cfg);
  if (cfg.src !== 'air' && !wf) return solveOne(cfg);
  const A = wf ? WATER_SRC : srcAirModel(cfg), Tin = wf ? +cfg.Ts : +cfg.srcT, mDry = wf ? Math.max(+cfg.srcMw, 1e-6) : Math.max(+cfg.srcM, 1e-6) / (1 + A.W), hin = A.h(Tin);
  const lo = Tin - 55, hi = Tin - 0.3;
  const run = To => {
    const c = wf ? {...cfg, Ts: Tin, Two: To} : {...cfg, Ts: Tin, Two: To, srcTo: To}, r = solveOne(c), Q = srcKW(r, c);
    const Tc = Q != null ? A.T(hin - Q / mDry) : null;
    return {To, r, c, Q, g: Tc == null ? null : Tc - To};
  };
  let a = run(Math.min(Math.max(LIVE.srcToGuess ?? Tin - 6, lo), hi));
  if (a.g == null) for (const To of [Tin - 3, Tin - 10, Tin - 18]) { a = run(Math.min(Math.max(To, lo), hi)); if (a.g != null) break; }
  if (a.g == null) { const r = a.r || {success: false}; if (!r.success) r.message = `${r.message || 'no solution'} (source ${wf ? 'water' : 'air'} ${Tin} °C, no outlet found)`; return r; }
  // the duty depends only weakly on the outlet, so the first step goes to the outlet the duty gives (a fixed-point
  // step), then secant; a warm start from the last solve usually lands within the tolerance at once
  const tol = 0.05;
  if (Math.abs(a.g) <= tol) return finishSrc(a, A, cfg);
  let b = run(Math.min(Math.max(a.To + a.g, lo), hi));
  if (Math.abs(a.To + a.g - b.To) < 1e-6 && b.g != null && Math.abs(b.g) > tol) b = run(Math.min(Math.max(a.To + (a.g > 0 ? 2 : -2), lo), hi));
  let best = Math.abs(a.g) < Math.abs(b.g ?? Infinity) ? a : b;
  for (let k = 0; k < 14 && Math.abs(best.g) > tol; k++) {
    let x;
    if (b.g == null) x = (a.To + b.To) / 2;
    else if (Math.abs(b.g - a.g) < 1e-9) x = b.To + (b.g > 0 ? 1 : -1);
    else x = b.To - b.g * (b.To - a.To) / (b.g - a.g);
    x = Math.min(Math.max(x, lo), hi);
    if (Math.abs(x - b.To) < 1e-4) break;
    const n = run(x);
    if (n.g == null) { b = {...b, To: (b.To + x) / 2}; continue; }
    a = b.g == null ? a : b; b = n;
    if (Math.abs(n.g) < Math.abs(best.g)) best = n;
  }
  if (cfg.cyc === 'X' && Math.abs(best.g) > tol) {
    const lim = srcLimited(cfg, solveOne, A, Tin, mDry, hin, best, wf);
    if (lim) return finishSrc(lim, A, cfg);
  }
  return finishSrc(best, A, cfg);
}
/* sCO2 Brayton with a fixed source flow: the LTHX inlet is the turbine exhaust, set by the cycle, so when the source
   is too small (or the approach too large) for the duty the approach allows, no source outlet balances. The source is
   then cooled to the turbine exhaust plus the source approach (the cold-end pinch), and the LTHX duty per kg of CO2 is
   capped at what that source gives; the cap is solved with the cycle, whose CO2 flow follows the sizing. */
function srcLimited(cfg, solveOne, A, Tin, mDry, hin, start, wf) {
  if (!start.r || !start.r.success || !start.r.report) return null;
  const dTs = srcDT(cfg);
  // the largest LTHX duty per kg CO2 (q) the source can sustain: with the cap q the cycle solves, the source is
  // cooled to the turbine exhaust plus the approach, and what it gives there covers the CO2 flow the sizing sets.
  // Bisection on q: a cap that fails to solve (too little heat) or that the source cannot cover moves the bracket.
  const at = q => {
    let Tt = start.r.report.T_turb_out, run = null;
    for (let k = 0; k < 3; k++) {                       // the turbine exhaust moves a little with the cap
      const To = Math.min(Tin - 0.3, Tt + dTs);
      const c = wf ? {...cfg, Ts: Tin, Two: To, _qcap: q} : {...cfg, Ts: Tin, Two: To, srcTo: To, _qcap: q};
      const r = solveOne(c);
      // a colder outlet can start the cycle on its collapsed branch (the LTHX cold end briefly too warm): keep the
      // last good outlet, whose approach is then at or above the one set
      if (!r.success) return run || {ok: false, r, c, To};
      const Q = srcKW(r, c), m = Q / r.q_source;
      run = {ok: true, r, c, To, Q, spare: mDry * (hin - A.h(To)) - m * r.q_source};
      if (Math.abs(r.report.T_turb_out - Tt) < 0.02) break;
      Tt = r.report.T_turb_out;
    }
    return run;
  };
  let lo = 0, hi = start.r.q_source, best = null, last = null;
  for (let k = 0; k < 22 && hi - lo > 1e-3 * Math.max(hi, 1); k++) {
    const q = (lo + hi) / 2, t = at(q); last = t;
    if (!t.ok) { lo = q; continue; }                      // too little heat for the target: the cap must rise
    if (t.spare >= 0) { lo = q; best = t; } else hi = q;  // the source covers it: try more
  }
  if (!best) {
    const r = (last && last.r) || start.r, To = last ? last.To : Tin - 0.3, Qs = mDry * (hin - A.h(To));
    const f = {...r, success: false, message: `${r.success ? 'no design' : r.message}. The source limits the LTHX: ${wf ? 'the water' : 'the source air'} (${(wf ? +cfg.srcMw : +cfg.srcM).toFixed(1)} kg/s) can only be cooled to the turbine exhaust plus the ${dTs} K approach (${To.toFixed(1)} °C), giving ${hfmt(Qs, 1)} kW, too little for this design; more source flow, a smaller approach or a higher discharge pressure gives more`};
    return {To, r: f, c: last ? last.c : cfg, Q: null, g: 0};
  }
  const out = {To: best.To, r: best.r, c: best.c, Q: best.Q, g: A.T(hin - best.Q / mDry) - best.To};
  out.r.srcLimited = true;
  out.r.warnings = [...(out.r.warnings || []), `the source limits the LTHX: ${wf ? 'the water' : 'the source air'} is cooled to the turbine exhaust plus the ${dTs} K approach (${out.To.toFixed(1)} °C) and gives ${hfmt(out.Q, 1)} kW, less than the approach alone would allow; more source flow or a smaller approach gives more`];
  return out;
}
function finishSrc(best, A, cfg) {
  const r = best.r;
  if (!r.success) return r;
  LIVE.srcToGuess = best.To;
  if (waterFlowSrc(cfg)) {           // source water at a fixed flow: its outlet is the solved temperature
    r.warnings = r.warnings || [];
    if (Math.abs(best.g) > 0.2) r.warnings.push(`source water outlet not converged (${best.g.toFixed(2)} K)`);
    r.cfgSolved = {Ts: +cfg.Ts, Two: best.To};
    return r;
  }
  const Tin = +cfg.srcT, mDry = Math.max(+cfg.srcM, 1e-6) / (1 + A.W);
  const Wout = A.Wout(best.To);
  r.srcAir = {Tin, Tout: best.To, Td: A.Td, W: A.W, Wout, m: +cfg.srcM, mDry, kw: best.Q, cond: mDry * Math.max(A.W - Wout, 0), resid: best.g};
  r.warnings = r.warnings || [];
  if (A.capped) r.warnings.push(`source air moisture ${(1000 * A.Win).toFixed(1)} g/kg is more than air at ${Tin.toFixed(1)} °C can hold; taken as saturated at ${(1000 * A.W).toFixed(1)} g/kg`);
  if (Math.abs(best.g) > 0.2) r.warnings.push(`source air outlet not converged (${best.g.toFixed(2)} K)`);
  if (best.To < 0 && r.srcAir.cond > 0) r.warnings.push(`source air leaves at ${best.To.toFixed(1)} °C: its moisture would freeze on the evaporator (defrost not modelled)`);
  r.cfgSolved = {Ts: Tin, Two: best.To, srcTo: best.To};
  return r;
}
let lvTimer = null;
/* automatic solving can be paused (a per-viewer choice kept in this browser): changes are then saved but not
   solved until Solve now or resuming */
let AUTO = (() => { try { return localStorage.getItem('dryerAutoSolve') !== 'off'; } catch (_) { return true; } })(), PENDING = false;
function renderAuto() {
  const b = $('autoBtn'); b.textContent = AUTO ? 'Pause auto-solve' : 'Resume auto-solve'; b.setAttribute('aria-pressed', String(!AUTO));
  $('solveBtn').hidden = AUTO; $('staleChip').hidden = AUTO || !PENDING;
  if (!AUTO && $('status').textContent === 'Loading CoolProp…') $('status').textContent = 'Auto-solve paused: CoolProp loads on the first solve';
}
function solveOnDemand() { PENDING = false; renderAuto(); clearTimeout(lvTimer); ensureCoolProp().then(solveNow, () => {}); }
$('autoBtn').addEventListener('click', () => {
  AUTO = !AUTO; try { localStorage.setItem('dryerAutoSolve', AUTO ? 'on' : 'off'); } catch (_) {}
  if (AUTO && PENDING) solveOnDemand(); else renderAuto();
});
$('solveBtn').addEventListener('click', solveOnDemand);
renderAuto();
function changed() {
  lvSave(); clearTimeout(lvTimer);
  if (typeof casesChanged === 'function') casesChanged();
  if (!AUTO) { PENDING = true; renderAuto(); return; }
  lvTimer = setTimeout(() => { ensureCoolProp().then(solveNow, () => {}); }, 120);
}
function solveNow() {
  if (state.mode !== 'live') return;
  syncW();
  const first = LIVE.flKey !== fluidKey(LV);
  STATUS(MULTI(LV.cyc) ? 'Solving…' : first ? `Preparing ${fluidLabel(LV)} and solving…` : 'Solving…'); $('status').className = 'chip';
  LIVE.busy = true; renderLive();
  setTimeout(() => {          // let the status paint before the solve blocks the page
    const t0 = performance.now(), cfg = JSON.parse(JSON.stringify(LV)), cfgIn = JSON.parse(JSON.stringify(LV)), fk = fluidKey(cfg);
    airFor(cfg);
    DRYER.setCompModel(cfg.compModel, cfg.etaC);
    let res;
    try {
      if (MULTI(cfg.cyc)) { res = withSource(cfg, solveMulti); throw {done: true}; }
      if (LIVE.flKey !== fk) {
        LIVE.flKey = ''; LIVE.dome = null;
        if (LIVE.fl && LIVE.fl.destroy) try { LIVE.fl.destroy(); } catch (_) {}
        LIVE.fl = new CoolPropFluid(CP, active(cfg), cfg.basis); LIVE.flKey = fk;
        try { const dm = SGHP.dome(LIVE.fl, -50, 40); LIVE.dome = {bubble: dm.liq, dew: dm.vap}; } catch (e) { LIVE.dome = null; }
      }
      res = withSource(cfg, c => {
        const base = baseOf(c);
        let r;
        if (RHP(c.cyc)) {
          r = solveA(c, base);
          if (r.success) { LV.PHg = r.P_H; LV.PLg = r.P_L; LV.gKey = fk; LV.off = r.offset; }
        } else {
          r = solveS(c, base);
          if (r.success) LV.PH = r.P_H * 10;
        }
        if (r.success) r.paths = DRYER.paths(LIVE.fl, r);
        return r;
      });
    } catch (e) { if (!(e && e.done)) res = {success: false, message: cpError(e)}; }
    res.ms = performance.now() - t0;
    if (res.cfgSolved) Object.assign(cfg, res.cfgSolved);
    LIVE.res = res; LIVE.cfg = cfg; LIVE.cfgIn = cfgIn; LIVE.busy = false;
    STATUS(`CoolProp ${CP_VERSION} ready`); $('status').className = 'chip info';
    const st = $('solve-time'); st.textContent = `${res.success ? 'Solved' : 'Stopped'} in ${(res.ms / 1000).toFixed(1)} s`; st.hidden = false;
    lvSave(); applyDuty(); renderLive(); renderCycles(); renderPFD(); showAutoSplits(); refreshHumOuts(); renderSizeBars();
    if (typeof casesChanged === 'function') casesChanged();
  }, 30);
}
/* solve a set of inputs away from the page (saved cases, comparison, sweeps): the same solvers as the live design,
   with the page's fluid and air swapped in for the call and restored after it */
function solveCfg(cfg) {
  const keep = {fl: LIVE.fl, air: LIVE.air, airBase: LIVE.airBase, airKey: LIVE.airKey, airSec: LIVE.airSec};
  const t0 = performance.now();
  let res;
  try {
    airFor(cfg);
    DRYER.setCompModel(cfg.compModel, cfg.etaC);
    if (MULTI(cfg.cyc)) res = withSource(cfg, solveMulti);
    else {
      const fl = fluidFor((cfg.mainComps || cfg.comps).filter(x => x[1] > 0), cfg.basis);
      LIVE.fl = fl;
      res = withSource(cfg, c => {
        const base = baseOf(c), r = RHP(c.cyc) ? solveA(c, base) : solveS(c, base);
        if (r.success) r.paths = DRYER.paths(fl, r);
        return r;
      });
      res.dome = fl.domeL || null;
    }
    if (res.success) res.secShare = secOn(cfg) ? DRYER.secondaryShare(LIVE.airBase, +cfg.Tain, +cfg.Tmt, +cfg.T, +cfg.phi) : 0;
  } catch (e) { res = {success: false, message: cpError(e)}; }
  finally {
    Object.assign(LIVE, keep);
    DRYER.setCompModel(LV.compModel, LV.etaC);
    refreshHumOuts();
  }
  res.ms = performance.now() - t0;
  return res;
}
/* RHP-A with the gas-cooler outlet offset optimised for COP: the smallest offset that keeps stage 1 dry (0 when it
   already is), then a golden-section search above it if COP still rises with the offset */
/* when the stage pressure-ratio limit is what stops a design, solve it once without the limit (at the offsets that
   lower the ratio most) and say what it would need */
function prNeed(cfg, base) {
  let best = null;
  for (const off of [50, 25, 0]) {
    const r = DRYER.rhpa(LIVE.fl, LIVE.air, {...base, ...lifts(cfg), max_stage_pr: 1e9, gc_out_offset: off});
    if (r.success) { const pr = Math.max(r.report.PR_1, r.report.PR_2); if (!best || pr < best.pr) best = {pr, r, off}; }
  }
  if (!best || best.pr <= (cfg.prMax || 20)) return '';
  const p = best.r.report;
  return `. Without the limit the design needs a stage pressure ratio of ${best.pr.toFixed(1)} (${p.P_evap_bar.toFixed(2)} → ${p.Pmax_bar.toFixed(1)} bar at a ${best.off} K offset, ηis ${p.eta_is_1.toFixed(2)}, COP ${best.r.COP.toFixed(2)}); raise the max stage pressure ratio in the Compressor card to accept it`;
}
function solveA(cfg, base) {
  const g = guessA();
  let PH = g.PH, PL = g.PL, n = 0;
  const run = off => {
    n++;
    const r = DRYER.rhpa(LIVE.fl, LIVE.air, {...base, ...lifts(cfg), max_stage_pr: cfg.prMax || 20, gc_out_offset: off, P_H_guess: PH, PL_guess: PL});
    if (r.success) { PH = r.P_H; PL = r.P_L; }
    r.offset = off; return r;
  };
  const wet = r => !r.success && /wet compression|suction is not gas|stage pressure ratio above|no gas-cooler pressure reaches/.test(r.message);   // a raised outlet can cure all four: the last when the pressure limit binds
  let best = run(0);
  if (!best.success) {
    if (!wet(best)) return best;
    let lo = 0, rh = null, hi = null;
    for (const o of [5, 10, 20, 30, 50]) { const r = run(o); if (r.success) { hi = o; rh = r; break; } if (!wet(r)) return r; lo = o; }
    if (!rh && /no gas-cooler pressure reaches/.test(best.message)) return Object.assign(best, {message: `${best.message}; no gas-cooler outlet offset up to 50 K meets the approach within the ${pmaxOf(cfg)} bar max pressure: raise the max pressure or the approach`});
    if (!rh) return Object.assign(best, {message: `${best.message}; no gas-cooler outlet offset up to 50 K keeps the compression dry with every stage pressure ratio at ${cfg.prMax || 20} or below${prNeed(cfg, base)}${cfg.cyc === 'H' && !cfg.ihxH ? '. Without an IHX the suction leaves the superheater at the water inlet less the source approach: a dry-expanding fluid such as isopentane then compresses into the two-phase region, so turn on the IHX or choose a fluid that compresses dry' : ''}`});
    while (hi - lo > 0.5) { const m = (lo + hi) / 2, r = run(m); if (r.success) { hi = m; rh = r; } else if (wet(r)) lo = m; else break; }
    best = rh;
  }
  const o0 = best.offset, up = run(o0 + 1);
  if (up.success && up.COP > best.COP) {
    let a = o0, b = Math.min(o0 + 20, 50), gr = (Math.sqrt(5) - 1) / 2;
    let x1 = b - gr * (b - a), x2 = a + gr * (b - a), r1 = run(x1), r2 = run(x2);
    const f = r => r.success ? r.COP : -Infinity;
    for (let k = 0; k < 6; k++) {
      if (isFinite(f(r1)) && isFinite(f(r2)) && Math.abs(f(r1) - f(r2)) < 1e-3 * Math.max(f(r1), f(r2))) break;   // COP flat to 0.1 %: the offset hardly matters
      if (f(r1) > f(r2)) { b = x2; x2 = x1; r2 = r1; x1 = b - gr * (b - a); r1 = run(x1); }
      else { a = x1; x1 = x2; r1 = r2; x2 = a + gr * (b - a); r2 = run(x2); }
    }
    for (const r of [up, r1, r2]) if (f(r) > best.COP) best = r;
  }
  best.offset_evals = n;
  return best;
}
/* SLTHP with the gas-cooler pressure solved: the lowest at which every approach is met (the study's optima sit on
   their approach constraints), searched from the last or loaded pressure */
/* SLTHP at a low-side cap (MPa, or null for the pinch-set low side): the lowest gas-cooler pressure at which every
   approach is met (the study's optima sit on their approach constraints), searched from P0 */
function solveSAt(cfg, base, cap, P0, split) {
  let n = 0, last = null;
  const run = P => { n++; const r = DRYER.slthp(LIVE.fl, LIVE.air, {...base, lift: cfg.liftS, split, P_H: P, alpha: cfg.alpha, ...(cap ? {P_L_max: cap} : {})}); if (r.alpha_needed) last = r; return r; };
  const margin = r => r.success ? r.report.min_approach_margin_K : (r.margin ?? NaN);
  const Pmax = pmaxOf(cfg) / 10;
  let P = Math.min(Math.max(P0, 0.5), Pmax), r = run(P), m = margin(r);
  let lo = null, hi = null, rHi = null;
  if (m >= 0) { hi = P; rHi = r; for (let k = 0; k < 80; k++) { const Pn = P / 1.05; const rn = run(Pn), mn = margin(rn); if (!(mn >= 0)) { lo = Pn; break; } P = Pn; hi = P; rHi = rn; if (P < 0.3) break; } }
  else { lo = P; for (let k = 0; k < 80; k++) { const Pn = Math.min(P * 1.05, Pmax); if (Pn <= P) break; const rn = run(Pn), mn = margin(rn); P = Pn; if (mn >= 0) { hi = P; rHi = rn; break; } lo = P; } }
  if (!rHi) {
    const a = last && last.alpha_needed && last.alpha_needed < 1 ? `; the RHX1 flow share must be at least ${last.alpha_needed.toFixed(2)} for the RHX1 branch to carry its duty` : '';
    return Object.assign(r, {message: `no gas-cooler pressure up to ${pmaxOf(cfg)} bar (the max pressure) meets the approaches (${r.message})${a}`, evals: n});
  }
  if (lo != null) for (let k = 0; k < 30 && hi - lo > 1e-4; k++) { const mid = Math.sqrt(lo * hi), rm = run(mid); if (margin(rm) >= 0) { hi = mid; rHi = rm; } else lo = mid; }
  rHi.evals = n; return rHi;
}
/* SLTHP. A subcritical low side evaporates at the highest pressure the source pinch allows. Above the critical point
   (CO2 with a warm source) the gas heater's pressure has no dew point to bound it, and the pinch-set pressure leaves
   too small a pressure ratio, so the low-side pressure is chosen for COP, as for the RHP cycles: a cap on it (the
   source pinch may lower it further) scanned in ln P over 0.4 P_crit to the max pressure / 1.2, then golden section,
   with the gas-cooler pressure solved at each */
function solveS(cfg, base) {
  const split = stagesOf(cfg) === 1 ? 0 : cfg.share, fl = LIVE.fl;
  let sc = false;
  try { sc = isFinite(fl.P_crit) && DRYER.scLow(fl, base.T_source_in - base.dT_src); } catch (_) { sc = false; }
  const guess = +cfg.PH > 0 ? cfg.PH / 10 : 2.5;
  if (!sc) return solveSAt(cfg, base, null, Math.min(guess, 9), split);
  let PHg = Math.max(guess, 1.5 * fl.P_crit), n = 0, best = null, lastFail = null;
  const memo = new Map();
  const f = x => {
    const k = x.toFixed(6); if (memo.has(k)) return memo.get(k);
    const r = solveSAt(cfg, base, Math.exp(x), PHg, split); n += r.evals || 0;
    if (r.success) { PHg = r.P_H; if (!best || r.COP > best.COP) best = r; } else lastFail = r;
    const v = r.success ? r.COP : -Infinity; memo.set(k, v); return v;
  };
  const lo = Math.log(0.4 * fl.P_crit), hi = Math.log(pmaxOf(cfg) / 10 / 1.2);
  const xs = Array.from({length: 7}, (_, i) => hi - (hi - lo) * i / 6), fs = xs.map(f);   // from the top: warm starts walk down
  if (!best) return Object.assign(lastFail, {message: `${lastFail.message} (supercritical low side: no gas-heater pressure from ${(Math.exp(lo) * 10).toFixed(0)} to ${(Math.exp(hi) * 10).toFixed(0)} bar solves)`, evals: n});
  const i = fs.indexOf(Math.max(...fs));
  let a = xs[Math.min(i + 1, xs.length - 1)], b = xs[Math.max(i - 1, 0)];
  const g = (Math.sqrt(5) - 1) / 2;
  let x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = f(x1), f2 = f(x2);
  for (let k = 0; k < 7; k++) {
    if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = f(x1); }
    else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = f(x2); }
  }
  best.evals = n; best.report.P_L_by = 'COP'; return best;
}
function liveEnter() {
  if (!D) return;
  if (!LV.init) {
    // first visit: RHP-A with n-butane heating 15 C air to 170 C with a 20 K approach from 40 C source water (butane
    // is not a study blend, so the design starts from the page's defaults rather than a loaded study design)
    setComps([['Butane', 100]]); LV.T = 170; LV.a = 20; LV.cyc = 'A'; LV.Ts = 40; LV.Tain = 15;
    LV.init = true; renderLiveControls(); changed(); return;
  }
  renderLiveControls();
  if (!LIVE.res || cfgKey(LIVE.cfg) !== cfgKey(LV)) changed(); else { renderLive(); renderCycles(); renderPFD(); }
}

/* ---------------------------------------------------------------- results */
const fx = (v, d) => v == null || !isFinite(v) ? '–' : v.toFixed(d);
function liveD() {
  const r = LIVE.res, c = LIVE.cfg;
  if (!r || !r.success) return null;
  const rep = r.report;
  return {live: true, COP: r.COP, chosen_COP: r.COP, Pmax: r.P_H * 10, lift: RHP(c.cyc) ? c.lift : c.liftS, split: RHP(c.cyc) ? c.split : c.share, variant: c.cyc, stages: stagesOf(c),
          motor: c.motor, alpha: c.cyc === 'S' ? rep.alpha_RHX1 : undefined, paths: r.paths, rho1: rep.rho_suction_1, rho2: rep.rho_suction_2, scLow: !!rep.low_side_above_Pcrit};
}
function liveSrc() {
  if (LIVE.cfg && MULTI(LIVE.cfg.cyc) || (!LIVE.cfg && MULTI(LV.cyc))) return multiSrc();
  const c = LIVE.cfg || LV, v = c.cyc, k = RHP(v) ? 'A' : 'S', nm = v === 'H' ? 'Conventional' : RHP(v) ? `RHP-${v}` : 'SLTHP', col = k === 'A' ? 'var(--sA)' : 'var(--sS)';
  const d = liveD(), r = LIVE.res;
  const desc = !r ? (LIVE.busy ? 'solving…' : 'waiting for CoolProp…') : !r.success ? `<span class="flag">${r.message}</span>`
    : k === 'A' ? `${v !== 'B' && v !== 'H' ? `RHX1 lift ${c.lift} K · ` : ''}${v !== 'A' && v !== 'H' ? `RHX2 lift ${c.lift2} K · ` : ''}${stagesOf(c) === 1 ? 'single stage' : `stage-2 power fraction ${c.split.toFixed(2)}`} · offset ${fx(r.offset, 1)} K · Pmax ${d.Pmax.toFixed(1)} bar · COP ${r.COP.toFixed(3)}`
    : `RHX1 lift ${c.liftS} K · RHX1 share ${(c.alpha * 100).toFixed(1)} % · Pmax ${d.Pmax.toFixed(1)} bar · COP ${r.COP.toFixed(3)}`;
  return {k, name: nm, color: col, w: c.w, label: fluidLabel(c), d, Tair: c.T, Tsrc: c.Ts, Tain: c.Tain, Twout: c.Two, ...srcOf(c, r), air: LIVE.air, phi: secOn(c) ? +c.phi : 0, dome: LIVE.flKey === fluidKey(c) ? LIVE.dome : null, live: true,
          variant: v, stages: stagesOf(c), title: k === 'A' ? `${nm} · ${v === 'H' ? `${(c.stagesH || 1) === 1 ? 'single-stage' : 'two-stage'}${c.ihxH ? ' with IHX' : ', no IHX'}` : v === 'A' ? (stagesOf(c) === 1 ? 'single stage, RHX1 on the suction' : 'RHX1 on the stage-1 suction') : v === 'B' ? 'RHX2 on the interstage gas' : 'RHX1 and RHX2'}` : 'SLTHP · RHX1 parallel to the LT gas cooler',
          info: `<div class="cyrow"><i class="sw" style="background:${col}"></i><b>${nm}</b><span class="f">${fluidLabel(c)} · ${c.T} °C air · ${c.a} K approach · ${desc}</span></div>`};
}
function renderLive() {
  lvTitle();
  if (LIVE.cfg && MULTI(LIVE.cfg.cyc)) return renderLiveMulti();
  const r = LIVE.res, c = LIVE.cfg, chips = [];
  if (LIVE.busy) chips.push('<span class="chip">Solving…</span>');
  if (r && !r.success) chips.push(`<span class="chip bad">No solution: ${r.message}</span>`);
  if (r && r.success) for (const w of r.warnings || []) {
    const m = /short of the minimum by ([\d.]+) K/.exec(w), small = m && +m[1] < 0.2;   // rounding of a stored design
    chips.push(`<span class="chip ${small || /flow share raised/.test(w) ? 'info' : 'warn'}">${small ? `Approach ${m[1]} K short of the minimum (within the rounding of the study's stored design)` : w}</span>`);
  }
  if (r && r.success && RHP(c.cyc) && r.report.ihx_on) chips.push('<span class="chip info">Internal heat exchanger on (source water below the gas-cooler outlet)</span>');
  if (c && +c.phi > 0 && !secOn(c)) chips.push(`<span class="chip warn">Secondary air ignored: its target (${c.Tmt} °C) must lie between the air inlet and the dryer inlet target</span>`);
  if (r && r.success && r.report.low_side_supercritical) {
    const p = r.report, q = p.evap_in_quality;
    chips.push(`<span class="chip info">${p.low_side_above_Pcrit ? `Supercritical low side: ${p.P_evap_bar.toFixed(1)} bar, above the critical pressure (${fx(p.P_crit_bar, 1)} bar); the evaporator is a gas heater` : `Heater outlet above the critical temperature: the low side (${p.P_evap_bar.toFixed(1)} bar) evaporates and then heats the gas`}, ${p.T_evap_inlet.toFixed(1)} °C${q > 0 && q < 1 ? ` (quality ${q.toFixed(2)})` : ''} → ${p.T_superheat_out.toFixed(1)} °C. Its pressure ${p.P_L_by === 'source pinch' ? 'is the highest the source pinch allows (also the COP optimum)' : 'is chosen for COP, below the pressure the source pinch allows'}</span>`);
  }
  if (r && r.success && r.report.T_dew_evap > c.Ts - srcDT(c) - 0.05) chips.push(`<span class="chip info">Evaporator dew point held at the superheater outlet (${(c.Ts - srcDT(c)).toFixed(0)} °C); inlet ${r.report.T_evap_inlet.toFixed(1)} °C</span>`);
  $('lvChips').innerHTML = chips.join('');
  lvHeadline();
  const figs = $('lvFigs');
  if (!r || !r.success) { figs.innerHTML = ''; return; }
  const p = r.report, fig = (k, v, u, sub = '') => `<div class="fig"><span class="k">${k}</span><span class="v">${v}<small>${u}</small></span>${sub ? `<span class="sub">${sub}</span>` : ''}</div>`;
  const one = stagesOf(c) === 1;
  const cells = one ? [
    fig('COP', r.COP.toFixed(3), '', 'heat to air / electricity'),
    fig('VHC', p.VHC.toFixed(2), 'MJ/m³', 'per m³ of suction'),
    fig('Gas-cooler pressure', p.Pmax_bar.toFixed(2), 'bar', 'solved for the approach'),
    fig('Suction pressure', p.P_evap_bar.toFixed(2), 'bar', p.T_dew_evap == null ? `supercritical (P_crit ${fx(p.P_crit_bar, 1)} bar)` : `dew ${p.T_dew_evap.toFixed(1)} °C`),
    fig('Pressure ratio', p.Overall_PR.toFixed(1), '', `ηis ${p.eta_is_1.toFixed(3)}`),
    fig('Discharge', p.T_comp_discharge.toFixed(1), '°C', `suction ${p.T_suction.toFixed(1)} °C`),
    fig('Heat to air', p.q_air_kJkg.toFixed(1), 'kJ/kg', 'gas coolers'),
    fig(r.srcAir ? 'Heat from source air' : 'Heat from water', p.q_evap_kJkg.toFixed(1), 'kJ/kg', r.srcAir ? `air ${r.srcAir.Tin.toFixed(1)} → ${r.srcAir.Tout.toFixed(1)} °C · ${(r.srcAir.cond * 3600).toFixed(1)} kg/h condensate` : `evaporator and superheater`),
  ] : [
    fig('COP', r.COP.toFixed(3), '', 'heat to air / electricity'),
    fig('VHC', p.VHC.toFixed(2), 'MJ/m³', 'per m³ of stage-1 + stage-2 suction'),
    fig('Gas-cooler pressure', p.Pmax_bar.toFixed(2), 'bar', 'solved for the approach'),
    fig('Suction pressure', p.P_evap_bar.toFixed(2), 'bar', p.T_dew_evap == null ? `supercritical (P_crit ${fx(p.P_crit_bar, 1)} bar)` : `dew ${p.T_dew_evap.toFixed(1)} °C`),
    fig('Intermediate pressure', p.P_int_bar.toFixed(2), 'bar', `PR ${p.PR_1.toFixed(2)} / ${p.PR_2.toFixed(2)}`),
    fig('Overall PR', p.Overall_PR.toFixed(1), '', `ηis ${p.eta_is_1.toFixed(3)} / ${p.eta_is_2.toFixed(3)}`),
    fig('Stage-2 discharge', p.T_comp_discharge.toFixed(1), '°C', `stage 1 out ${p.Comp_1_T_out.toFixed(1)} °C`),
    fig('Heat to air', p.q_air_kJkg.toFixed(1), 'kJ/kg', `RHX1 ${p.q_rhx1_kJkg.toFixed(1)}${p.q_rhx2_kJkg ? ` · RHX2 ${p.q_rhx2_kJkg.toFixed(1)}` : ''} kJ/kg`),
    fig(r.srcAir ? 'Heat from source air' : 'Heat from water', p.q_evap_kJkg.toFixed(1), 'kJ/kg', r.srcAir ? `air ${r.srcAir.Tin.toFixed(1)} → ${r.srcAir.Tout.toFixed(1)} °C · ${(r.srcAir.cond * 3600).toFixed(1)} kg/h condensate` : `evaporator and superheater`),
  ];
  if (RHP(c.cyc)) {
    cells.push(fig('Gas-cooler pinch', p.gc_minimum_approach.toFixed(2), 'K', `approach ${c.a} K`));
    cells.push(fig('GC outlet offset', r.offset.toFixed(1), 'K', `optimised · ${r.offset_evals} solves`));
  }
  else cells.push(fig('Tightest approach margin', p.min_approach_margin_K.toFixed(2), 'K', p.rhx1_LMTD_K != null ? `RHX1 LMTD ${p.rhx1_LMTD_K.toFixed(1)} K` : ''));
  cells.push(fig('Electricity', p.w_el_kJkg.toFixed(1), 'kJ/kg', `per kg of refrigerant, motor ${c.motor.toFixed(2)}${p.expander ? ' · net of the expander' : ''}`));
  if (p.expander) cells.push(expFig(p, fig));
  if (!one) cells.push(fig('Stage-2 power share', (100 * p.power_frac_2).toFixed(1), '%', 'of the shaft power'));
  if (p.evap_pinch_K != null) cells.push(fig(p.low_side_above_Pcrit ? 'Gas-heater pinch' : 'Evaporator pinch', p.evap_pinch_K.toFixed(2), 'K', `to the ${r.srcAir ? 'source air' : 'source water'} · ${p.low_side_above_Pcrit ? 'heater' : 'evaporator'} in ${p.T_evap_inlet.toFixed(1)} °C`));
  if (secOn(c)) cells.push(fig('Secondary air', (100 * secShare(c)).toFixed(1), '%', `of the heat · ${(+c.phi).toFixed(2)} kg per kg primary to ${c.Tmt} °C`));
  cells.push(fig('Air per refrigerant', p.air_per_ref_kgkg.toFixed(2), 'kg/kg', `${secOn(c) ? 'primary ' : ''}dry air${secOn(c) ? ` · total ${(p.air_per_ref_kgkg * (1 + +c.phi)).toFixed(2)}` : ''}${p.water_per_ref_kgkg != null ? ` · water ${p.water_per_ref_kgkg.toFixed(2)} kg/kg` : ''}`));
  cells.push(fig('Solve time', (r.ms / 1000).toFixed(1), 's', `CoolProp ${CP_VERSION}`));
  figs.innerHTML = cells.join('');
}

/* ---------------------------------------------------------------- headline: COP with stars and four dials, as the steam page */
const hfmt = (v, d) => v == null || !isFinite(v) ? '–' : v.toFixed(d);
const niceTop = (v, step) => Math.max(step, Math.ceil(v * 1.15 / step) * step);
function starRowL(n, id) {
  let out = '';
  for (let i = 0; i < 5; i++) {
    const f = Math.max(0, Math.min(1, n - i)), cid = `${id}-${i}`;
    const P = 'M10 1.5 L12.6 7.1 L18.6 7.7 L14 11.8 L15.4 17.8 L10 14.7 L4.6 17.8 L6 11.8 L1.4 7.7 L7.4 7.1 Z';
    out += `<svg viewBox="0 0 20 20" aria-hidden="true"><defs><clipPath id="${cid}"><rect x="0" y="0" width="${20 * f}" height="20"/></clipPath></defs>
      <path d="${P}" fill="var(--grid)" stroke="var(--dome)" stroke-width="1"/>
      <path d="${P}" fill="var(--accent)" stroke="var(--accent)" stroke-width="1" clip-path="url(#${cid})"/></svg>`;
  }
  return out;
}
function dialL(value, max, colorVar, unit, dec = 1) {
  const W = 120, H = 92, cx = 60, cy = 58, R = 44, a0 = 210, a1 = -30;
  const pt = (deg, rr = R) => [cx + rr * Math.cos(deg * Math.PI / 180), cy - rr * Math.sin(deg * Math.PI / 180)];
  const arc = (from, to) => { const [x0, y0] = pt(from), [x1, y1] = pt(to); return `M${x0.toFixed(2)} ${y0.toFixed(2)} A${R} ${R} 0 ${from - to > 180 ? 1 : 0} 1 ${x1.toFixed(2)} ${y1.toFixed(2)}`; };
  const f = Math.max(0, Math.min(1, value / max)), av = a0 - (a0 - a1) * f;
  const [nx, ny] = pt(av, R - 10);
  let ticks = '';
  for (let k = 0; k <= 4; k++) { const d = a0 - (a0 - a1) * k / 4, [x0, y0] = pt(d, R + 4), [x1, y1] = pt(d, R + 9);
    ticks += `<line x1="${x0.toFixed(1)}" y1="${y0.toFixed(1)}" x2="${x1.toFixed(1)}" y2="${y1.toFixed(1)}" stroke="var(--dome)" stroke-width="1"/>`; }
  const [lx0, ly0] = pt(a0, R + 6), [lx1, ly1] = pt(a1, R + 6);
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${hfmt(value, dec)} ${unit} of ${max}">
    <path d="${arc(a0, a1)}" fill="none" stroke="var(--grid)" stroke-width="9" stroke-linecap="round"/>
    ${f > 0.002 ? `<path d="${arc(a0, av)}" fill="none" stroke="var(${colorVar})" stroke-width="9" stroke-linecap="round"/>` : ''}
    ${ticks}
    <line x1="${cx}" y1="${cy}" x2="${nx.toFixed(2)}" y2="${ny.toFixed(2)}" stroke="var(--ink)" stroke-width="2.2" stroke-linecap="round"/>
    <circle cx="${cx}" cy="${cy}" r="4" fill="var(--ink)"/>
    <text x="${cx}" y="${cy + 24}" text-anchor="middle" style="font-family:var(--f-data);font-size:15px;fill:var(--ink);font-weight:500">${hfmt(value, dec)}</text>
    <text x="${lx0.toFixed(1)}" y="${(ly0 + 12).toFixed(1)}" text-anchor="middle" style="font-family:var(--f-data);font-size:9.5px;fill:var(--ink-3)">0</text>
    <text x="${lx1.toFixed(1)}" y="${(ly1 + 12).toFixed(1)}" text-anchor="middle" style="font-family:var(--f-data);font-size:9.5px;fill:var(--ink-3)">${max}</text>
  </svg>`;
}
function lvHeadline() {
  if (LIVE.cfg && MULTI(LIVE.cfg.cyc)) return lvHeadlineMulti();
  const host = $('lvHeadline'), r = LIVE.res, c = LIVE.cfg, d = liveD();
  if (!r || !r.success || !d) { host.hidden = true; return; }
  const u = energy(RHP(c.cyc) ? 'A' : 'S', d, state.vsw, c.T, c.Tain);
  if (!u) { host.hidden = true; return; }
  host.hidden = false;
  const heat = u.kw.air, el = u.kw.el, kwTop = niceTop(Math.max(heat, el), 10), p = r.report;
  const n = Math.max(0, Math.min(5, r.COP));
  const water = c.Two != null ? u.kw.ev / (4.18 * (c.Ts - c.Two)) : null;
  host.innerHTML = `<div class="hl-txt hl-cop"><span class="hl-k">Heating COP</span>
      <span class="hl-v">${hfmt(r.COP, 2)}</span>
      <span class="stars" role="img" aria-label="${n.toFixed(2)} of 5 stars, one star per unit of COP">${starRowL(n, 'lvst')}</span></div>
    <div class="dial">${dialL(heat, kwTop, '--s4', 'kW heat to air')}<span class="dk">Heat to air, kW</span></div>
    <div class="dial">${dialL(el, kwTop, '--s2', 'kW electrical')}<span class="dk">Electrical, kW</span></div>
    <div class="dial">${dialL(p.Pmax_bar, niceTop(p.Pmax_bar, 10), '--s6', 'bar maximum pressure')}<span class="dk">Max pressure, bar</span></div>
    <div class="dial">${dialL(p.T_comp_discharge, niceTop(p.T_comp_discharge, 50), '--s5', '°C discharge', 0)}<span class="dk">Discharge, °C</span></div>
    <span class="hl-s">${vswTxt()} m³/h low-stage swept volume · dryer air ${(u.mAir * 3.6).toFixed(1)} t/h, ${hfmt(c.Tain, 1)} → ${hfmt(c.T, 1)} °C · ${srcLine(c, r, u.kw.ev, water)}</span>`;
}

/* ---------------------------------------------------------------- a small schematic of the chosen cycle (Design card) */
function miniSchematic(cyc) {
  const ar = (x, y, dir) => { const d = {r: [-7, -4, -7, 4], l: [7, -4, 7, 4], u: [-4, 7, 4, 7], d: [-4, -7, 4, -7]}[dir];
    return `<polygon points="${x},${y} ${x + d[0]},${y + d[1]} ${x + d[2]},${y + d[3]}" fill="currentColor"/>`; };
  const T = (x, y, s, cls = 'cn', a = 'middle') => `<text x="${x}" y="${y}" text-anchor="${a}" class="${cls}">${s}</text>`;
  const comp = (x, y, l) => `<polygon class="eq" points="${x - 13},${y - 18} ${x + 13},${y - 9} ${x + 13},${y + 9} ${x - 13},${y + 18}"/>${T(x, y + 4, l)}`;
  const box = (x, y, w, h, l) => `<rect class="eq" x="${x}" y="${y}" width="${w}" height="${h}" rx="3"/>${T(x + w / 2, y + h / 2 + 4, l)}`;
  const valve = (x, y) => `<polygon class="eq" points="${x - 7},${y - 9} ${x + 7},${y - 9} ${x - 7},${y + 9} ${x + 7},${y + 9}"/>`;
  const air = `<path class="air" d="M372 8 H86"/><g style="color:var(--ink-3)">${ar(86, 8, 'l')}</g>${T(376, 12, 'air in', '', 'start')}${T(80, 12, 'dryer', '', 'end')}`;
  const water = `<path class="wat" d="M336 172 H262"/><g style="color:var(--accent)">${ar(262, 172, 'l')}</g>${T(340, 176, 'water', '', 'start')}`;
  const low = `<g style="color:var(--lp)"><path class="pl" d="M356 104 V150 H336"/>${ar(336, 150, 'l')}<path class="pl" d="M262 150 H232"/>${ar(232, 150, 'l')}<path class="pl" d="M188 150 H40 V134"/>${ar(40, 134, 'u')}</g>
    <path class="ph" d="M40 98 V78" style="stroke:var(--hp)"/>`;
  let hi, rhx, cap;
  if (cyc === 'A') {
    hi = `<g style="color:var(--hp)"><path class="ph" d="M40 42 V30 H86"/>${ar(86, 30, 'r')}<path class="ph" d="M146 30 H188"/>${ar(188, 30, 'r')}
      <path class="ph" d="M232 30 H262"/>${ar(262, 30, 'r')}<path class="ph" d="M322 30 H356 V86"/>${ar(356, 86, 'd')}</g>`;
    rhx = `<rect class="eq" x="188" y="18" width="44" height="144" rx="3"/>${T(210, 94, 'RHX1')}${T(210, 50, 'hot', 'lbl2')}${T(210, 140, 'cold', 'lbl2')}`;
    cap = 'RHP-A: RHX1 in series between the HT and LT gas coolers; it heats the suction gas (lift) before stage 1.';
  } else {
    hi = `<g style="color:var(--hp)"><path class="ph" d="M40 42 V30 H86"/>${ar(86, 30, 'r')}<path class="ph" d="M146 30 H166"/>
      <path class="ph" d="M166 30 V22 H262"/>${ar(262, 22, 'r')}<path class="ph" d="M322 22 H356 V58"/>
      <path class="ph" d="M166 30 V58 H188"/>${ar(188, 58, 'r')}<path class="ph" d="M232 58 H356"/>${ar(352, 58, 'r')}<path class="ph" d="M356 58 V86"/>${ar(356, 86, 'd')}</g>
      <circle class="node" cx="166" cy="30" r="3"/><circle class="node" cx="356" cy="58" r="3"/>${T(178, 18, '1 − α', 'lbl2', 'start')}${T(172, 72, 'α', 'lbl2', 'start')}`;
    rhx = `<rect class="eq" x="188" y="46" width="44" height="116" rx="3"/>${T(210, 108, 'RHX1')}${T(210, 74, 'hot', 'lbl2')}${T(210, 142, 'cold', 'lbl2')}`;
    cap = 'SLTHP: after the HT gas cooler the flow splits; RHX1 runs parallel to the LT gas cooler and heats the suction gas.';
  }
  const ltY = cyc === 'A' ? 18 : 10;
  return `<svg viewBox="0 0 420 186" role="img" aria-label="${cap}">${air}${water}${hi}${low}
    ${comp(40, 116, 'C1')}${comp(40, 60, 'C2')}${box(86, 18, 60, 24, 'HT GC')}${box(262, ltY, 60, 24, 'LT GC')}${rhx}
    ${valve(356, 95)}${box(262, 138, 74, 24, 'Evaporator')}</svg><figcaption>${cap}</figcaption>`;
}

/* RHP-A/B/C: RHX1 heats the stage-1 suction, RHX2 the interstage gas; their hot sides sit between the gas coolers */
function miniRHP(v, stages = LV.stagesH || 1) {
  if (v === 'H') return miniConv(stages, (LIVE.cfg && LIVE.cfg.cyc === 'H' ? LIVE.cfg : LV).ihxH);
  const r1 = v !== 'B', r2 = v !== 'A';
  const ar = (x, y, dir) => { const d = {r: [-7, -4, -7, 4], l: [7, -4, 7, 4], u: [-4, 7, 4, 7], d: [-4, -7, 4, -7]}[dir];
    return `<polygon points="${x},${y} ${x + d[0]},${y + d[1]} ${x + d[2]},${y + d[3]}" fill="currentColor"/>`; };
  const T = (x, y, s, cls = 'cn', a = 'middle') => `<text x="${x}" y="${y}" text-anchor="${a}" class="${cls}">${s}</text>`;
  const comp = (x, y, l) => `<polygon class="eq" points="${x - 13},${y - 18} ${x + 13},${y - 9} ${x + 13},${y + 9} ${x - 13},${y + 18}"/>${T(x, y + 4, l)}`;
  const box = (x, y, w, h, l) => `<rect class="eq" x="${x}" y="${y}" width="${w}" height="${h}" rx="3"/>${T(x + w / 2, y + h / 2 + 4, l)}`;
  let rx1x = 0, rx2x = 0;
  const hp = [], lp = [], eq = [];
  // top line, left to right: C2 -> [GC3, RHX2 hot] -> [GC2, RHX1 hot] -> GC1 -> valve
  let x = 70; const Y = 24;
  const single = v === 'A' && stages === 1;          // one compressor: C1 discharges straight to the gas coolers
  hp.push(single ? `M40 132 V${Y} H${x}` : `M40 ${32} V${Y} H${x}`);
  const gcs = r1 && r2 ? ['HT GC', 'MT GC', 'LT GC'] : ['HT GC', 'LT GC'];
  const blocks = [];
  if (r2) blocks.push('rhx2');
  if (r1) blocks.push('rhx1');
  let gi = 0;
  const putGC = () => { eq.push(box(x, Y - 12, 52, 24, gcs[gi++])); hp.push(`M${x + 52} ${Y} H${x + 70}`); x += 70; };
  for (const bk of blocks) {
    putGC();
    const bw = 38;
    if (bk === 'rhx2') { eq.push(`<rect class="eq" x="${x}" y="${Y - 14}" width="${bw}" height="100" rx="3"/>${T(x + bw / 2, Y + 44, 'RHX2')}`); rx2x = x; }
    else { eq.push(`<rect class="eq" x="${x}" y="${Y - 14}" width="${bw}" height="178" rx="3"/>${T(x + bw / 2, Y + 84, 'RHX1')}`); rx1x = x; }
    hp.push(`M${x + bw} ${Y} H${x + bw + 18}`); x += bw + 18;
  }
  putGC();
  const xv = x + 6;
  hp.push(`M${x - 18} ${Y} H${xv} V96`);
  eq.push(`<polygon class="eq" points="${xv - 7},96 ${xv + 7},96 ${xv - 7},114 ${xv + 7},114"/>`);
  // low side: valve -> evaporator -> [RHX1 cold] -> C1
  const ey = 176, ex = xv - 80;
  lp.push(`M${xv} 114 V${ey} H${ex + 66}`);
  eq.push(box(ex, ey - 12, 66, 24, 'Evaporator'));
  if (r1) { lp.push(`M${ex} ${ey} H${rx1x + 38}`, `M${rx1x} ${ey} H40 V168`); }
  else lp.push(`M${ex} ${ey} H40 V168`);
  // interstage: C1 -> [RHX2 cold] -> C2
  const mid = [];
  if (r2) mid.push(`M40 132 V112 H${rx2x}`, `M${rx2x} 88 H52 V68`);
  else if (!single) mid.push('M40 132 V68');
  const air = `<path class="air" d="M${xv + 14} 6 H86"/><g style="color:var(--ink-3)">${ar(86, 6, 'l')}</g>${T(80, 10, 'dryer', '', 'end')}`;
  const water = `<path class="wat" d="M${ex + 66} 198 H${ex}"/><g style="color:var(--accent)">${ar(ex, 198, 'l')}</g>${T(ex + 70, 202, 'water', '', 'start')}`;
  const cap = v === 'A' ? `RHP-A${single ? ', single stage' : ''}: RHX1 heats the ${single ? 'compressor' : 'stage-1'} suction; its hot side sits between the HT and LT gas coolers.`
    : v === 'B' ? 'RHP-B: RHX2 heats the interstage gas before stage 2; an IHX on the low side when the source water is cooler than the gas-cooler outlet.'
    : 'RHP-C: RHX1 heats the stage-1 suction and RHX2 the interstage gas; three gas coolers on the air.';
  return `<svg viewBox="0 0 ${xv + 40} 210" role="img" aria-label="${cap}">${air}${water}
    <g style="color:var(--hp)">${hp.map(d => `<path class="ph" d="${d}"/>`).join('')}${mid.map(d => `<path class="ph" d="${d}" style="stroke:var(--hp)"/>`).join('')}</g>
    <g style="color:var(--lp)">${lp.map(d => `<path class="pl" d="${d}"/>`).join('')}</g>
    ${single ? comp(40, 150, 'C') : comp(40, 150, 'C1') + comp(40, 50, 'C2')}${eq.join('')}</svg><figcaption>${cap}</figcaption>`;
}

function miniConv(stages, ihx) {
  const ar = (x, y, dir) => { const d = {r: [-7, -4, -7, 4], l: [7, -4, 7, 4], u: [-4, 7, 4, 7], d: [-4, -7, 4, -7]}[dir];
    return `<polygon points="${x},${y} ${x + d[0]},${y + d[1]} ${x + d[2]},${y + d[3]}" fill="currentColor"/>`; };
  const T = (x, y, s, cls = 'cn', a = 'middle') => `<text x="${x}" y="${y}" text-anchor="${a}" class="${cls}">${s}</text>`;
  const comp = (x, y, l) => `<polygon class="eq" points="${x - 13},${y - 18} ${x + 13},${y - 9} ${x + 13},${y + 9} ${x - 13},${y + 18}"/>${T(x, y + 4, l)}`;
  const box = (x, y, w, h, l) => `<rect class="eq" x="${x}" y="${y}" width="${w}" height="${h}" rx="3"/>${T(x + w / 2, y + h / 2 + 4, l)}`;
  // top: compressor -> gas cooler [-> IHX hot] -> valve; bottom: valve -> evaporator [-> IHX cold] -> compressor
  const Y = 24, ey = 176, gx = 90, gw = 120, ix = 240, iw = 38, xv = ihx ? 380 : 320, ex = ihx ? 296 : 220, ew = ihx ? 66 : 80;
  const hp = (stages === 1 ? ['M40 82 V24 H' + gx] : ['M40 132 V68', 'M40 32 V24 H' + gx]).concat(ihx ? [`M${gx + gw} ${Y} H${ix}`, `M${ix + iw} ${Y} H${xv} V96`] : [`M${gx + gw} ${Y} H${xv} V96`]);
  const lp = [`M${xv} 114 V${ey} H${ex + ew}`].concat(ihx ? [`M${ex} ${ey} H${ix + iw}`, `M${ix} ${ey} H40 V${stages === 1 ? 118 : 168}`] : [`M${ex} ${ey} H40 V${stages === 1 ? 118 : 168}`]);
  const comps = stages === 1 ? comp(40, 100, 'C') : comp(40, 150, 'C1') + comp(40, 50, 'C2');
  const cap = (stages === 1 ? 'Conventional single-stage heat pump: compressor, gas cooler(s) on the air, valve, and evaporator and superheater on the source water'
    : 'Conventional two-stage heat pump: two compressors in series without intercooling or recuperation, gas cooler(s) on the air, valve, evaporator and superheater on the water') + (ihx ? '; the IHX heats the suction with the gas-cooler outlet.' : '.');
  const air = `<path class="air" d="M${xv + 14} 6 H${gx + 6}"/><g style="color:var(--ink-3)">${ar(gx + 6, 6, 'l')}</g>${T(gx, 10, 'dryer', '', 'end')}`;
  const water = `<path class="wat" d="M${ex + ew} 198 H${ex}"/><g style="color:var(--accent)">${ar(ex, 198, 'l')}</g>${T(ex + ew + 4, 202, 'water', '', 'start')}`;
  return `<svg viewBox="0 0 ${xv + 40} 210" role="img" aria-label="${cap}">${air}${water}
    <g style="color:var(--hp)">${hp.map(d => `<path class="ph" d="${d}"/>`).join('')}</g>
    <g style="color:var(--lp)">${lp.map(d => `<path class="pl" d="${d}"/>`).join('')}</g>
    ${comps}${box(gx, Y - 12, gw, 24, 'Gas cooler')}${ihx ? `<rect class="eq" x="${ix}" y="${Y - 14}" width="${iw}" height="178" rx="3"/>${T(ix + iw / 2, Y + 84, 'IHX')}` : ''}
    <polygon class="eq" points="${xv - 7},96 ${xv + 7},96 ${xv - 7},114 ${xv + 7},114"/>${box(ex, ey - 12, ew, 24, ihx ? 'Evap.' : 'Evaporator')}</svg><figcaption>${cap}</figcaption>`;
}

/* ================================================================== multi-loop cycles: TT cascade, n-cascade, HP + MVR */
const MULTI_NAME = {K: 'TT cascade', N: 'n-cascade', M: 'HP + MVR', X: 'sCO₂ Brayton', E: 'Flash / economiser', D: 'n-pressure', R: 'Stirling', Y: 'NH₃–H₂O hybrid'};
const MULTI_TITLE = {K: 'Transcritical–transcritical cascade (Kong et al.)', N: 'n-level cascade (Hoang et al.)', M: 'HP + MVR: heat pump and MVR steam', X: 'sCO₂ reverse Brayton (Dalal et al.)',
  E: 'Two-stage heat pump: flash tank, economiser or parallel compression', D: 'Multi-pressure condensation (after Dai et al.)', R: 'Reverse-Stirling units (semi-empirical)', Y: 'Ammonia–water hybrid absorption–compression (Osenbrück)'};
function renderMultiFields() {
  const c = LV.cyc, p = LV.par[c], box = $('lvDesX');
  const num = (key, lab, step, title = '', idx = null) => {
    const v = idx == null ? p[key] : p[key][idx];
    return `<label title="${title}">${lab}<input type="number" step="${step}" data-par="${key}"${idx == null ? '' : ` data-idx="${idx}"`} value="${v}"></label>`;
  };
  const chk = (key, lab, idx = null) => { const v = idx == null ? p[key] : p[key][idx]; return `<label class="inline"><input type="checkbox" data-par="${key}"${idx == null ? '' : ` data-idx="${idx}"`}${v ? ' checked' : ''}> ${lab}</label>`; };
  let h = '';
  if (c === 'K') {
    h = `<div class="fields">${num('P_b', 'Bottom pressure, bar', 'any', 'Bottom-cycle discharge; below the critical pressure the bottom cycle runs subcritical')}${num('P_t', 'Top pressure, bar', 'any')}${num('r', 'Bottom/top mass flow', 'any', 'Bottom-cycle refrigerant flow per kg of top-cycle flow')}${num('gc_off', 'Top GC outlet raise, K', 'any', 'Raises the top gas-cooler outlet above the coldest the approach allows, which gives the top IHX more heat (Kong: the outlet temperature is an optimisation variable)')}
      ${num('T_ev_b', 'Bottom evaporating, °C', 'any')}${num('T_ev_t', 'Top evaporating, °C', 'any')}
      ${num('dT_cgc', 'Cascade GC approach, K', 'any')}${num('dT_ihx_b', 'Bottom IHX approach, K', 'any')}${num('dT_ihx_t', 'Top IHX approach, K', 'any')}</div>
      <div class="row optrow">${num('Pmax', 'Max pressure, bar', 'any', 'Upper bound on both discharge pressures for the optimiser')}${optBtn()}</div>
      <p class="hint">Optimise adjusts the bottom and top pressures, the bottom/top mass flow and the top gas-cooler outlet raise for the highest COP, up to the max pressure.</p>
      <p class="hint">The bottom cycle rejects all its heat in the cascade gas cooler, superheating the top suction as far as the CGC pinch allows. The top discharge heats the primary air in the HT gas cooler, then the primary and secondary air in the MT gas cooler; both keep the gas-cooler approach. After Kong (2026, PhD thesis, Waikato), cycle TT-4 with the additional gas cooler (Figure 58): about a third of the air goes to the fluidised bed at 102 °C (secondary/primary air ratio 0.5 in the Operating point). Kong's reference: ηis 0.70 with 0.96 mechanical efficiency (Constant 0.70 and motor 0.96 here), 150 bar, 200 °C air, COP 2.55.</p>`;
  } else if (c === 'N') {
    h = `<div class="pickrow"><span class="lab">Levels</span><div class="seg" id="lvNlev"></div></div><div class="fields">`;
    for (let k = 0; k < p.n - 1; k++) h += num('T_cond', `Level ${k + 1} condensing, °C`, 'any', '', k);
    h += `${num('dT_casc', 'Cascade approach, K', 'any')}${num('superheat', 'Min suction superheat, K', 'any', 'Raised by the IHX where a dry fluid needs more for a dry discharge')}${num('dT_ihx', 'IHX approach, K', 'any', 'A level with its IHX on heats its suction to the liquid entering the IHX less this approach')}</div><div class="levels">`;
    for (let k = 0; k < p.n; k++) h += `<div class="row"><span class="lab lvlab">Level ${k + 1}</span>${k < p.n - 1 ? chk('cond_air', 'condense to air', k) : '<span class="hint">condenses to air</span>'}${chk('sub_air', 'subcool to air', k)}${chk('ihx', 'IHX', k)}</div>`;
    h += `</div><div class="row optrow">${optBtn()}</div><p class="hint">Optimise sets the lower levels' condensing temperatures for the highest COP, keeping every air approach.</p><p class="hint">Air passes the levels in series. A lower level that condenses into the air does so until the air leaves it at its condensing temperature less the approach; the rest drives the next level. The top level's condensing temperature is solved for the air approach.</p>`;
  } else if (c === 'X') {
    h = `<div class="pickrow"><span class="lab">Architecture</span><div class="seg" id="lvXarch"></div></div>
      <div class="fields">${num('P_L', 'Suction pressure, bar', 'any', 'Compressor 1 suction, the low side of the cycle: 80 bar in the paper, above the 73.8 bar critical pressure; Optimise also sets it')}${num('P_H', 'Discharge pressure, bar', 'any', 'Compressor discharge; compressor 2 must reach the dryer target plus the gas-cooler approach')}
      ${num('x', 'Secondary split', '0.01', 'Share of the CO₂ sent to compressor 1, whose hot discharge preheats the main stream before compressor 2')}
      ${num('eps_ihx', 'IHX effectiveness', '0.01')}${num('eps_pre', 'Preheater effectiveness', '0.01')}
      ${num('eta_c', 'Compressor ηis', '0.01', 'Both compressors; the Efficiency model under Compressor does not apply to this cycle')}${num('eta_t', 'Turbine ηis', '0.01')}${num('eta_m', 'Mechanical efficiency', '0.001', 'Compressors take work / this, the turbine gives back work × this')}
      ${num('dP_hp', 'High-side ΔP per exchanger, kPa', 'any')}${num('dP_lp', 'Low-side ΔP per exchanger, kPa', 'any')}</div>
      <div class="row">${chk('P_H_min', 'use the lowest discharge pressure that reaches the target')}</div>
      <div class="row optrow">${num('Pmax', 'Max pressure, bar', 'any', 'Upper bound for the optimiser and the lowest-pressure search')}${optBtn()}</div>
      <p class="hint">Optimise sets the suction pressure, the discharge pressure and the split for the highest COP. It starts from the values above (80 bar suction by default) and may take the suction below the 73.8 bar critical pressure (a transcritical low side) or up to where the discharge pressure meets the max pressure.${p.arch === 2 ? ' Architecture 2 is slow: a few minutes.' : ''}</p>
      <p class="hint">${p.arch === 2 ? 'Architecture 2 (with Held’s split stream): the secondary stream, after preheating the main stream, heats the cooler part of the air in HTHX 2 and rejoins before the turbine; the main stream heats the hotter part in HTHX 1, then gives heat back in the IHX.' : 'Architecture 1: compressor 2 heats the air in the HTHX; the secondary stream, after preheating the main stream, rejoins it after the HTHX and both give heat back in the IHX before the turbine.'} The LTHX (gas heater) keeps the source approach and the HTHX the gas-cooler approach, both under Operating point. The whole cycle is supercritical; the turbine (expander) replaces the valve. After Dalal, Klöppel, Rose and Walmsley (Euro sCO₂ 2027).</p>`;
  } else if (c === 'E') {
    const arr = ECON_ARR[p.arr] || 'flash';
    h = `<div class="pickrow"><span class="lab">Arrangement</span><div class="seg" id="lvEarr"></div></div>
      ${arr === 'single' ? '' : `<div class="row">${chk('T_int_auto', 'intermediate pressure at the geometric mean of the evaporating and high pressures')}</div>`}
      <div class="fields">${arr === 'single' || p.T_int_auto ? '' : num('T_int', 'Intermediate saturation, °C', 'any', 'Saturation (dew) temperature at the intermediate pressure; Optimise sets it')}
      ${num('superheat', 'Min suction superheat, K', 'any', 'Raised where an overhanging fluid needs more for a discharge 5 K above the dew point')}
      ${arr === 'economiser' ? num('dT_econ', 'Economiser approach, K', 'any', 'The main liquid leaves the economiser this far above the saturation temperature at the intermediate pressure') : ''}
      ${arr === 'economiser' || arr === 'parallel' ? num('superheat_inj', 'Injected vapour superheat, K', 'any', arr === 'parallel' ? 'The parallel compressor suction is heated at least this far above saturation by the high-pressure liquid (more if a dry discharge needs it)' : 'The injected vapour leaves the economiser this far above saturation') : ''}</div>
      <div class="row optrow">${num('Pmax', 'Max pressure, bar', 'any', 'Highest high-side pressure searched')}${arr === 'single' ? '' : optBtn()}</div>
      ${arr === 'single' ? '' : '<p class="hint">Optimise sets the intermediate saturation temperature for the highest COP.</p>'}
      <p class="hint">${{flash: 'Flash tank: the condenser liquid flashes at the intermediate pressure; the flash vapour joins the low-stage discharge at the high-stage suction and the saturated liquid expands to the evaporator.',
        economiser: 'Closed economiser (vapour injection): part of the condenser liquid expands to the intermediate pressure and evaporates against the rest, which is subcooled; the vapour joins the low-stage discharge.',
        parallel: 'Parallel compression: as the flash tank, but the flash vapour has its own compressor to the high pressure and the low stage compresses straight to it.',
        single: 'Single stage, for comparison: one compressor and one valve.'}[arr]} The high pressure is the lowest that heats the air to the target with the approach (subcritical or transcritical, from the fluid). After Arpagaus et al. (2018) and Adamson et al. (2022).</p>`;
  } else if (c === 'D') {
    h = `<div class="pickrow"><span class="lab">Condensing pressures</span><div class="seg" id="lvDn"></div></div>
      <div class="row">${chk('T_cond_auto', 'lower condensing temperatures in equal steps from the evaporator to the air target')}</div>
      <div class="fields">${p.T_cond_auto ? '' : Array.from({length: p.n - 1}, (_, k) => num('T_cond', `Condensing ${k + 1}, °C`, 'any', k ? '' : 'Saturation (dew) temperatures of the lower pressures, rising; Optimise sets them', k)).join('')}
      ${num('superheat', 'Min suction superheat, K', 'any', 'Raised where an overhanging fluid needs more for a discharge 5 K above the dew point')}</div>
      <div class="row optrow">${num('Pmax', 'Max pressure, bar', 'any', 'Highest top pressure searched')}${optBtn()}</div>
      <p class="hint">Optimise sets the lower condensing temperatures for the highest COP.</p>
      <p class="hint">Compressor 1 feeds condenser 1, which heats the air from the inlet; at each pressure part of the flow goes on through the next compressor to the next condenser, up to the top one, which heats the air to the target. Each condenser's liquid is throttled into the one below. The flows follow from the air balances, the air between condensers is the hottest the condensers allow, and the top pressure is the lowest that reaches the target. Two pressures is the dual-pressure cycle of Dai et al. (2020).</p>`;
  } else if (c === 'R') {
    h = `<div class="pickrow"><span class="lab">Units in series</span><div class="seg" id="lvRn"></div></div>
      ${p.n > 1 ? `<div class="row">${chk('bounds_auto', 'equal air temperature steps between the units')}</div>` : ''}
      <div class="fields">${p.n > 1 && !p.bounds_auto ? Array.from({length: p.n - 1}, (_, k) => num('bounds', `Air after unit ${k + 1}, °C`, 'any', '', k)).join('') : ''}
      ${num('eta_carnot', 'Carnot fraction', '0.01', 'COP = this × T_hot / (T_hot − T_cold), in K; 0.37 from the RFC Maasdam case (COP 2 for 155 °C steam from about 75 °C waste heat)')}
      ${num('dT_hot', 'Hot-side approach, K', 'any', 'Each unit\'s hot side sits this far above the air it delivers (through an intermediate water or steam loop)')}${num('dT_cold', 'Cold-side approach, K', 'any', 'The cold side sits this far below the source outlet')}</div>
      ${p.n > 1 ? `<div class="row optrow">${optBtn()}</div><p class="hint">Optimise sets the air temperatures between the units for the highest COP.</p>` : ''}
      <p class="hint">A semi-empirical block, not a physical Stirling model: industrial reverse-Stirling heat pumps (helium, about 400–720 kW) publish little performance data, so each unit's COP is a fraction of the Carnot COP between its hot and cold sides (Høeg et al. 2026; Häggqvist et al. 2023). Treat results as indicative; above about 200 °C hot side there is no field data. Size the units by sink duty or air flow (Size), as there is no compressor swept volume.</p>`;
  } else if (c === 'Y') {
    h = `<div class="fields">${num('w_rich', 'NH₃ mass fraction', '0.01', 'Ammonia mass fraction of the circulating solution; Optimise sets it')}${num('q_des', 'Degassing (vapour fraction)', '0.01', 'Vapour mass fraction boiled off in the desorber; Optimise sets it')}
      ${num('T_dis_max', 'Discharge limit, °C', 'any', 'Ammonia compressors and oil: 170 °C (Jensen et al. 2015)')}${num('eta_pump', 'Solution pump η', '0.01')}</div>
      <div class="row">${chk('wet', 'wet compression: inject lean liquid at the suction to hold the discharge limit')}${p.wet ? '' : chk('limit_dis', 'fail designs above the discharge limit')}</div>
      <div class="row optrow">${num('Pmax', 'Max pressure, bar', 'any', 'Component pressure class: Jensen et al. compare 28, 50 and 140 bar')}${optBtn()}</div>
      <p class="hint">Optimise sets the ammonia fraction and the degassing for the highest COP. A solve takes a few seconds; Optimise up to a few minutes.</p>
      <p class="hint">The solution boils partly in the desorber on the source and separates; the ammonia-rich vapour is compressed (with lean liquid injected to hold the discharge limit) and the lean liquid pumped, and both are absorbed again in the absorber, which heats the air with a large temperature glide. The absorber outlet is the coldest that keeps the air approach and the high pressure the lowest that heats the air to the target. Properties: the IAPWS G4-01 ammonia–water formulation (Tillner-Roth and Friend 1998), solved in the page. After Jensen et al. (2015).</p>`;
  } else {
    h = `<div class="fields">${num('T_cond', 'Bottom condensing, °C', 'any')}${num('dT_sg', 'Steam generator approach, K', 'any')}
      ${num('n_mvr', 'MVR stages', '1', 'Each at the fan pressure ratio (MVR fans card)')}</div>${splitFields(p)}
      <div class="row">${chk('inject', 'water injection between stages')}${chk('sub_air', 'bottom liquid subcooled to the air')}${chk('ihx', 'bottom IHX')}</div>
      ${p.ihx ? `<div class="fields">${num('dT_ihx', 'Bottom IHX approach, K', 'any', 'The IHX heats the bottom suction to the liquid entering it (after any subcooling to the air) less this approach')}</div>` : ''}
      <div class="row optrow">${optBtn()}</div><p class="hint">Optimise sets the bottom condensing temperature and the number of MVR stages for the highest COP, keeping both air approaches, with the process splits ${p.split_fill ? 'sized to the approach at every point' : 'held at the values above'}.</p>
      <p class="hint">The bottom condenser heats the air to its condensing temperature less the approach and raises steam with the rest; the steam (under vacuum below 100 °C) is compressed stage by stage and desuperheats and condenses into the air to the target.</p>`;
  }
  box.innerHTML = h;
  if (c === 'N') lvSeg('lvNlev', [[2, '2'], [3, '3'], [4, '4']], p.n, v => { p.n = v; renderLiveControls(); changed(); });
  for (const inp of box.querySelectorAll('[data-par]')) inp.addEventListener('change', () => {
    const k = inp.dataset.par, i = inp.dataset.idx, v = inp.type === 'checkbox' ? inp.checked : +inp.value;
    if (inp.type !== 'checkbox' && !isFinite(v)) return;
    if (i == null) p[k] = v; else p[k][+i] = k === 'split' ? Math.min(Math.max(v, 0), 1) : v;
    if (LIVE.optMsg) LIVE.optMsg[c] = '';
    if (inp.type === 'checkbox') { if ((k === 'ihx' && c === 'M') || ['T_int_auto', 'bounds_auto', 'wet', 'T_cond_auto'].includes(k)) renderMultiFields(); $('lvMini').innerHTML = miniMulti(LV.cyc); }   // the schematic shows the options
    changed();
  });
  const ob = $('lvOpt'); if (ob) ob.addEventListener('click', () => optimise(c));
  if (c === 'M') lvSeg('lvSplitMode', [[false, 'Manual'], [true, 'Auto (to the approach)']], !!p.split_fill, v => { p.split_fill = v; renderLiveControls(); changed(); });
  if (c === 'X') lvSeg('lvXarch', [[1, '1: split + preheater'], [2, '2: with split-stream HTHXs']], +p.arch === 2 ? 2 : 1, v => { p.arch = v; renderLiveControls(); changed(); });
  if (c === 'E') lvSeg('lvEarr', [[0, 'Flash tank'], [1, 'Economiser'], [2, 'Parallel'], [3, 'Single stage']], +p.arr || 0, v => { p.arr = v; if (LIVE.optMsg) LIVE.optMsg.E = ''; renderLiveControls(); changed(); });
  if (c === 'D') lvSeg('lvDn', [[2, '2'], [3, '3'], [4, '4'], [5, '5']], p.n, v => { p.n = v; if (LIVE.optMsg) LIVE.optMsg.D = ''; renderLiveControls(); changed(); });
  if (c === 'R') lvSeg('lvRn', [[1, '1'], [2, '2'], [3, '3'], [4, '4']], p.n, v => { p.n = v; if (LIVE.optMsg) LIVE.optMsg.R = ''; renderLiveControls(); changed(); });
}
/* HP + MVR: after each fan but the last, the fraction of its exhaust that goes to the process (condensing into the
   air at that pressure); the rest goes on to the next fan. Auto sizes each to heat the air to its saturation less
   the approach, and shows the values of the last solve. */
function splitFields(p) {
  const n = Math.max(1, Math.round(p.n_mvr)), auto = !!p.split_fill, rep = LIVE.res && LIVE.res.cycle === 'M' && LIVE.res.report.splits;
  let h = `<div class="pickrow"><span class="lab" title="Exhaust after each MVR fan: the fraction to the process (condensing into the air at that pressure); the rest goes on to the next fan">Process split</span><div class="seg" id="lvSplitMode"></div></div>`;
  if (n < 2) return h + '<p class="hint">One fan: all its exhaust goes to the process.</p>';
  h += '<div class="fields splits">';
  for (let k = 0; k < n - 1; k++) {
    const v = auto ? (rep && rep.length === n - 1 ? (+rep[k]).toFixed(3) : '') : p.split[k];
    h += `<label title="Fraction of fan ${k + 1}'s exhaust to the process">Fan ${k + 1} → process<input type="number" min="0" max="1" step="0.01" data-par="split" data-idx="${k}" value="${v}"${auto ? ' disabled' : ''}></label>`;
  }
  return h + `<label>Fan ${n} → process<input type="number" value="1" disabled></label></div>`;
}
function showAutoSplits() {          // Auto split: show the fractions of the last solve in the (disabled) fields
  if (LV.cyc !== 'M' || !LV.par.M.split_fill) return;
  const sp = LIVE.res && LIVE.res.success && LIVE.res.cycle === 'M' ? LIVE.res.report.splits : null;
  for (const inp of document.querySelectorAll('#lvDesX [data-par="split"]')) inp.value = sp && sp[+inp.dataset.idx] != null ? (+sp[+inp.dataset.idx]).toFixed(3) : '';
}
function optBtn() {
  const c = LV.cyc, busy = LIVE.opt === c;
  return `<button type="button" class="tbtn" id="lvOpt"${LIVE.opt ? ' disabled' : ''}>Optimise for COP</button><span class="hint" id="lvOptMsg">${busy ? 'Optimising…' : (LIVE.optMsg && LIVE.optMsg[c]) || ''}</span>`;
}
/* Multi-loop cycles: maximise the COP over each cycle's design variables, a few solves per frame so the page keeps
   painting. TT: the two pressures, the bottom/top flow and the top GC outlet raise; n-level cascade: the lower
   levels' condensing temperatures; booster: the bottom condensing temperature and the MVR stage count. */
function optimise(c) {
  if (LIVE.opt || !CP || !LIVE.air || !MULTI(c)) return;
  airFor(LV);
  const cfg = JSON.parse(JSON.stringify(LV)), p = cfg.par[c];
  if (cfg.src === 'air') {          // air source: optimise at the air outlet of the last solve; the page then re-solves it
    const To = LIVE.res && LIVE.res.srcAir ? LIVE.res.srcAir.Tout : cfg.srcT - 6;
    Object.assign(cfg, {Ts: cfg.srcT, Two: To, srcTo: To});
  } else if (cfg.srcWm === 'flow') {          // fixed water flow: at the water outlet of the last solve
    cfg.Two = LIVE.cfg && LIVE.cfg.src !== 'air' && LIVE.cfg.Two != null ? LIVE.cfg.Two : cfg.Ts - 10;
  }
  const base = baseOf(cfg);
  DRYER.setCompModel(cfg.compModel, cfg.etaC);
  let gen, apply;
  if (c === 'K') {
    const fls = cfg.multi.K.slice(0, 2).map(l => fluidFor(l, cfg.basis));
    gen = DRYCASC.kongOptimiseSteps(fls[0], fls[1], LIVE.airBase, {...base, ...p, ...ttAir(cfg), P_b: p.P_b / 10, P_t: p.P_t / 10, eta: null}, {Pmax: (p.Pmax || 150) / 10});
    apply = (o, P) => {
      P.P_b = +(o.P_b * 10).toFixed(1); P.P_t = +(o.P_t * 10).toFixed(1); P.r = +o.r.toFixed(3); P.gc_off = +(+o.gc_off).toFixed(1);
      return [o.atBound.P_b && 'bottom pressure', o.atBound.P_t && 'top pressure', o.atBound.r && 'mass-flow ratio'];
    };
  } else if (c === 'N') {
    const fls = cfg.multi.N.slice(0, p.n).map(l => fluidFor(l, cfg.basis));
    gen = DRYCASC.cascadeOptimiseSteps(fls, LIVE.air, {...base, ...p});
    apply = (o, P) => { o.T_cond.forEach((v, k) => { P.T_cond[k] = +(+v).toFixed(1); }); return o.atBound.map((b, k) => b && `level ${k + 1} condensing`); };
  } else if (c === 'X') {
    const fls = cfg.multi.X.slice(0, 1).map(l => fluidFor(l, cfg.basis));
    gen = DRYCASC.braytonOptimiseSteps(fls[0], LIVE.air, {...base, ...braytonCfg(p)});
    apply = (o, P) => { if (o.P_L != null) P.P_L = +(o.P_L * 10).toFixed(1); P.P_H = +(o.P_H * 10).toFixed(1); P.x = +o.x.toFixed(3); P.P_H_min = false; return [o.atBound.P_L && 'suction pressure', o.atBound.P_H && 'discharge pressure', o.atBound.x && 'split', o.budgetHit && 'solve budget (run it again to continue)']; };
  } else if (c === 'E' || c === 'D') {
    const fl = fluidFor(cfg.multi[c][0], cfg.basis);
    if (c === 'E') {
      gen = DRYMORE.economiserOptimiseSteps(fl, LIVE.air, {...base, ...econCfg(p)});
      apply = (o, P) => { if (o.T_int != null) { P.T_int = +(+o.T_int).toFixed(1); P.T_int_auto = false; } return [o.atBound && 'intermediate temperature']; };
    } else {
      gen = DRYMORE.multiPressureOptimiseSteps(fl, LIVE.air, {...base, ...multiCfg(p)});
      apply = (o, P) => { o.T_cond.forEach((v, k) => { P.T_cond[k] = +(+v).toFixed(1); }); P.T_cond_auto = false; return o.atBound.map((b, k) => b && `condensing ${k + 1}`); };
    }
  } else if (c === 'R') {
    gen = DRYMORE.stirlingOptimiseSteps(LIVE.air, {...base, ...stirCfg(p)});
    apply = (o, P) => { o.boundaries.forEach((v, k) => { P.bounds[k] = +(+v).toFixed(1); }); if (o.boundaries.length) P.bounds_auto = false; return []; };
  } else if (c === 'Y') {
    gen = DRYMORE.hybridOptimiseSteps(LIVE.air, {...base, ...hybCfg(p)});
    apply = (o, P) => { P.w_rich = +(+o.w_rich).toFixed(3); P.q_des = +(+o.q_des).toFixed(3); return [o.atBound.w_rich && 'ammonia fraction', o.atBound.q_des && 'degassing']; };
  } else {
    const fls = cfg.multi.M.slice(0, 1).map(l => fluidFor(l, cfg.basis));
    gen = DRYCASC.boosterOptimiseSteps(fls[0], fluidFor([['Water', 100]], 'mass'), LIVE.air, {...base, ...p, ...fanCfg(cfg)});
    apply = (o, P) => { P.T_cond = +(+o.T_cond).toFixed(1); P.n_mvr = o.n_mvr; if (!P.split_fill && o.splits) o.splits.forEach((v, k) => { P.split[k] = v; }); return [o.atBound.T_cond && 'bottom condensing', o.atBound.n_mvr && 'MVR stages']; };
  }
  LIVE.opt = c; LIVE.optMsg = LIVE.optMsg || {}; LIVE.optMsg[c] = '';
  const say = t => { const m = $('lvOptMsg'); if (m && LV.cyc === c) m.textContent = t; };
  const btn = $('lvOpt'); if (btn) btn.disabled = true; say('Optimising…');
  const t0 = performance.now();
  const tick = () => {
    let s;
    try { const tEnd = performance.now() + 60; do { s = gen.next(); } while (!s.done && performance.now() < tEnd); }
    catch (e) { s = {done: true, value: {success: false, message: String(e && e.message || e)}}; }
    if (!s.done) { say(`Optimising… ${s.value.nEval} solves · COP ${isFinite(s.value.best) ? s.value.best.toFixed(3) : '–'}`); setTimeout(tick, 0); return; }
    LIVE.opt = null;
    const o = s.value;
    if (!o.success) { LIVE.optMsg[c] = `Optimiser: ${o.message}`; renderLiveControls(); return; }
    const bound = apply(o, LV.par[c]).filter(Boolean);
    const st = isFinite(o.start) ? o.start.toFixed(3) + (o.startOk === false ? ' (broke an approach)' : '') : '–';
    LIVE.optMsg[c] = `COP ${st} → ${o.COP.toFixed(3)} in ${o.nEval} solves, ${((performance.now() - t0) / 1000).toFixed(1)} s` + (bound.length ? ` · at its bound: ${bound.join(', ')}` : '');
    renderLiveControls(); changed();
  };
  setTimeout(tick, 0);
}
function fluidFor(list, basis) {
  const key = list.map(x => `${x[0]}:${(+x[1]).toFixed(6)}`).join('&') + '|' + basis;
  LIVE.flc = LIVE.flc || {};
  if (!LIVE.flc[key]) {
    const f = new CoolPropFluid(CP, list.filter(x => x[1] > 0), basis);
    try { const dm = SGHP.dome(f, -50, 40); f.domeL = {bubble: dm.liq, dew: dm.vap}; } catch (e) { f.domeL = null; }
    LIVE.flc[key] = f;
  }
  return LIVE.flc[key];
}
// sCO2 Brayton inputs in the page's units (bar, kPa) to the model's (MPa)
const braytonCfg = p => ({...p, P_L: p.P_L / 10, P_H: p.P_H / 10, Pmax: p.Pmax / 10, dP_hp: p.dP_hp / 1000, dP_lp: p.dP_lp / 1000});
// the new cycles' design inputs (dryer_more.js) from the page's
const econCfg = p => ({arrangement: ECON_ARR[p.arr] || 'flash', T_int: p.T_int_auto ? null : +p.T_int, dT_econ: +p.dT_econ, superheat_inj: +p.superheat_inj, superheat: +p.superheat, P_max_bar: +p.Pmax});
const multiCfg = p => ({n_pressures: p.n, T_cond: p.T_cond_auto ? null : p.T_cond.slice(0, p.n - 1).map(Number), superheat: +p.superheat, P_max_bar: +p.Pmax});
const stirCfg = p => ({n_units: p.n, boundaries: p.n > 1 && !p.bounds_auto ? p.bounds.slice(0, p.n - 1).map(Number) : null, eta_carnot: +p.eta_carnot, dT_hot: +p.dT_hot, dT_cold: +p.dT_cold});
const hybCfg = p => ({w_rich: +p.w_rich, q_des: +p.q_des, P_max_bar: +p.Pmax, T_dis_max: +p.T_dis_max, wet_compression: !!p.wet, limit_discharge: !p.wet && !!p.limit_dis, eta_pump: +p.eta_pump});
const ttAir = cfg => ({air_ratio: secOn(cfg) ? +cfg.phi : 0, T_mt: secOn(cfg) ? +cfg.Tmt : (cfg.Tain + cfg.T) / 2});
function solveMulti(cfg) {
  const c = cfg.cyc, p = cfg.par[c];
  const base = baseOf(cfg);
  const fls = NOFLUID(c) ? [] : cfg.multi[c].slice(0, loopCount(c)).map(l => fluidFor(l, cfg.basis));
  let r;
  if (c === 'K') r = DRYCASC.kong(fls[0], fls[1], LIVE.airBase, {...base, ...p, ...ttAir(cfg), P_b: p.P_b / 10, P_t: p.P_t / 10, eta: null});
  else if (c === 'N') r = DRYCASC.cascade(fls, LIVE.air, {...base, ...p});
  else if (c === 'X') r = DRYCASC.brayton(fls[0], LIVE.air, {...base, ...braytonCfg(p), q_lthx_max: cfg._qcap, _quick: cfg._qcap != null});   // the source-limit search skips the lowest-pressure search on failures
  else if (c === 'E') r = DRYMORE.economiser(fls[0], LIVE.air, {...base, ...econCfg(p)});
  else if (c === 'D') r = DRYMORE.multiPressure(fls[0], LIVE.air, {...base, ...multiCfg(p)});
  else if (c === 'R') r = DRYMORE.stirling(LIVE.air, {...base, ...stirCfg(p)});
  else if (c === 'Y') r = DRYMORE.hybrid(LIVE.air, {...base, ...hybCfg(p)});
  else r = DRYCASC.booster(fls[0], fluidFor([['Water', 100]], 'mass'), LIVE.air, {...base, ...p, ...fanCfg(cfg)});
  r.fluids = fls; if (c === 'M') r.fluids = [fls[0], fluidFor([['Water', 100]], 'mass')];
  return r;
}
const loopFluid = (r, L) => L.fluid === 'water' ? r.fluids[r.fluids.length - 1] : r.fluids[L.fluid];
/* duties at the swept volume: the reference compressor (the bottom loop's) sets the mass flows */
function energyG(r, c) {
  if (!r || !r.success) return null;
  // Stirling: no compressor, so the volume is the dryer air's at the inlet (dry-air density) and the basis kg of air
  const nv = r.ref.air ? 1 : state.nvMode === 'fixed' ? state.nvFix : DRYER.etaVol(r.ref.PR);
  const rho = r.ref.air ? (+c.Patm || 101.325) * 1000 / (287.05 * (+c.Tain + 273.15)) : r.ref.rho;
  const m = rho * nv * state.vsw / 3600 / r.ref.flow;        // kg/s of the result's basis flow
  const wEl = r.report && isFinite(r.report.w_el_kJkg) ? r.report.w_el_kJkg : r.w_shaft / (c.motor ?? 0.93);   // the fans' own motor in HP + MVR
  const kw = {air: m * r.q_air, el: m * wEl, ev: m * r.q_source};
  kw.ex = m * (r.report && r.report.w_exp_kJkg || 0);
  kw.loss = kw.el - m * r.w_shaft + kw.ex;
  return {m, nv, kw};
}
function multiSrc() {
  const c = LIVE.cfg || LV, v = c.cyc, r = LIVE.res && LIVE.cfg && MULTI(LIVE.cfg.cyc) ? LIVE.res : null;
  const desc = !r ? (LIVE.busy ? 'solving…' : 'waiting for CoolProp…') : !r.success ? `<span class="flag">${r.message}</span>` : `COP ${r.COP.toFixed(3)}${isFinite(r.Pmax_bar) ? ` · Pmax ${r.Pmax_bar.toFixed(1)} bar` : ''}`;
  const fl = loopNames(v).map((n, i) => `${n} ${loopFluidText(c, v, i)}`).join(' · ');
  return {generic: true, k: v, name: MULTI_NAME[v], color: 'var(--sA)', res: r, cfg: c, Tair: c.T, Tsrc: c.Ts, Tain: c.Tain, Twout: c.Two, ...srcOf(c, r), live: true,
          title: MULTI_TITLE[v], label: fl,
          info: `<div class="cyrow"><i class="sw" style="background:var(--sA)"></i><b>${MULTI_NAME[v]}</b><span class="f">${fl} · ${c.T} °C air · ${desc}</span></div>`};
}
function renderLiveMulti() {
  const r = LIVE.res, c = LIVE.cfg, chips = [];
  if (LIVE.busy) chips.push('<span class="chip">Solving…</span>');
  if (r && !r.success) chips.push(`<span class="chip bad">No solution: ${r.message}</span>`);
  if (r && r.success) for (const w of r.warnings || []) chips.push(`<span class="chip ${/vacuum/.test(w) ? 'info' : 'warn'}">${w}</span>`);
  $('lvChips').innerHTML = chips.join('');
  lvHeadlineMulti();
  const figs = $('lvFigs');
  if (!r || !r.success) { figs.innerHTML = ''; return; }
  const p = r.report, fig = (k, v, u, sub = '') => `<div class="fig"><span class="k">${k}</span><span class="v">${v}<small>${u}</small></span>${sub ? `<span class="sub">${sub}</span>` : ''}</div>`;
  const cells = c.cyc === 'R' ? [fig('COP', r.COP.toFixed(3), '', 'heat to air / electricity'), fig('Hottest unit', r.T_dis_max.toFixed(1), '°C', 'hot side'),
      fig('Cold side', p.T_cold.toFixed(1), '°C', 'every unit'), fig('Carnot fraction', p.units[0].eta_carnot.toFixed(2), '', p.eta_model)]
    : [fig('COP', r.COP.toFixed(3), '', 'heat to air / electricity'), fig('VHC', r.VHC.toFixed(2), 'MJ/m³', 'per m³ of suction gas, all compressors'),
    fig('Max pressure', r.Pmax_bar.toFixed(1), 'bar', ''), fig('Hottest discharge', r.T_dis_max.toFixed(1), '°C', '')];
  if (r.report && r.report.expander) cells.push(expFig(r.report, fig));
  if (c.cyc === 'K') cells.push(
    fig('Top cycle', `${p.P_evap_top_bar.toFixed(1)} → ${p.P_top_bar.toFixed(0)}`, 'bar', `PR ${p.PR_top.toFixed(2)} · ηis ${p.eta_top.toFixed(2)} · discharge ${p.T_dis_top.toFixed(0)} °C`),
    fig('Bottom cycle', `${p.P_evap_bot_bar.toFixed(1)} → ${p.P_bot_bar.toFixed(0)}`, 'bar', `PR ${p.PR_bot.toFixed(2)} · ${p.bottom_supercritical ? 'transcritical' : 'subcritical'} · discharge ${p.T_dis_bot.toFixed(0)} °C`),
    fig('Top suction', p.T_suc_top.toFixed(1), '°C', `cascade GC ${p.q_cgc.toFixed(1)} kJ/kg`),
    ...(p.air_ratio > 0 ? [fig('Secondary air', (100 * p.share_sec).toFixed(1), '%', `of the heat · ${p.air_ratio.toFixed(2)} kg per kg primary to ${c.Tmt} °C`)] : []),
    fig('HT gas cooler', p.gap_HT.toFixed(1), 'K', `approach · ${p.q_HT.toFixed(0)} kJ/kg · gas out ${p.T_gc_mid.toFixed(1)} °C`),
    fig('MT gas cooler', p.gap_MT.toFixed(1), 'K', `approach · ${p.q_MT.toFixed(0)} kJ/kg · gas out ${p.T_gc_out.toFixed(1)} °C`),
    fig('Bottom after CGC', p.T_bot_cgc_out.toFixed(1), '°C', `valve inlet ${p.T_bot_valve.toFixed(1)} °C`),
    fig(r.srcAir ? 'Air between evaporators' : 'Water between evaporators', p.T_water_mid.toFixed(1), '°C', `series, top then bottom · approaches ${p.gap_evap_top.toFixed(1)} / ${p.gap_evap_bot.toFixed(1)} K`),
    fig('IHX duties', `${p.q_ihx_top.toFixed(0)} / ${p.q_ihx_bot.toFixed(0)}`, 'kJ/kg', 'top / bottom, per kg top'),
    fig('Bottom/top flow', p.r.toFixed(3), 'kg/kg', `refrigerant mass flow · GC outlet raise ${p.gc_off.toFixed(1)} K`));
  else if (c.cyc === 'N') p.levels.forEach((l, k) => cells.push(fig(`Level ${k + 1}`, `${l.Te.toFixed(0)} → ${l.Tc.toFixed(1)}`, '°C', `${l.Pe_bar.toFixed(2)} → ${l.Pc_bar.toFixed(1)} bar · PR ${l.PR.toFixed(2)} · ${(100 * l.q_air / r.q_air).toFixed(0)} % of the heat · superheat ${l.superheat.toFixed(0)} K${k < p.levels.length - 1 ? ` · ${(100 * l.f_air).toFixed(0)} % condensing to air` : ''}`)));
  else if (c.cyc === 'E') cells.push(
    fig('Pressures', `${p.P_evap_bar.toFixed(2)}${p.P_int_bar != null ? ` → ${p.P_int_bar.toFixed(1)}` : ''} → ${p.P_high_bar.toFixed(1)}`, 'bar', `evaporating ${p.T_evap.toFixed(1)} °C${p.T_int != null ? ` · intermediate ${p.T_int.toFixed(1)} °C` : ''}${p.supercritical ? ' · transcritical' : ''}`),
    fig('Compression', p.PR_high != null ? `${p.PR_low.toFixed(2)} / ${p.PR_high.toFixed(2)}` : p.PR_low.toFixed(2), 'PR', `${p.arrangement === 'parallel' ? 'low stage / parallel' : p.PR_high != null ? 'low / high stage' : 'single stage'} · discharge ${p.T_dis_low.toFixed(0)}${p.T_dis_high != null ? ` / ${p.T_dis_high.toFixed(0)}` : p.T_dis_parallel != null ? ` / ${p.T_dis_parallel.toFixed(0)}` : ''} °C`),
    ...(p.arrangement === 'single' ? [] : [fig(p.arrangement === 'economiser' ? 'Injected' : 'Flash vapour', (100 * p.flows.injected).toFixed(1), '%', `of the high-pressure flow · ${(100 * p.flows.low).toFixed(1)} % through the evaporator`)]),
    fig('Condenser', p.gap_air.toFixed(1), 'K', `air approach · suction superheat ${p.superheat.toFixed(1)} K`));
  else if (c.cyc === 'D') cells.push(
    fig('Pressures', [p.P_evap_bar.toFixed(2), ...p.P_cond_bar.map(v => v.toFixed(1))].join(' → '), 'bar', `evaporating ${p.T_evap.toFixed(1)} °C · condensing ${p.T_cond.map(v => v.toFixed(1)).join(' / ')} °C${p.supercritical ? ' · top transcritical' : ''}`),
    fig('Air between condensers', p.T_air_between.map(v => v.toFixed(0)).join(' / '), '°C', `heat split ${p.q_air_levels.map(q => (100 * q / p.q_air).toFixed(0)).join(' / ')} %`),
    fig('Compressor flows', p.flows.map(v => v.toFixed(2)).join(' / '), 'kg/kg', `of the evaporator flow · discharge ${p.T_dis.map(v => v.toFixed(0)).join(' / ')} °C`),
    fig('Air approaches', p.gaps.map(v => v.toFixed(1)).join(' / '), 'K', `PR ${p.PR.map(v => v.toFixed(2)).join(' / ')}`));
  else if (c.cyc === 'R') p.units.forEach((u, k) => cells.push(fig(`Unit ${k + 1}`, u.COP.toFixed(2), 'COP', `air ${u.T_air_in.toFixed(0)} → ${u.T_air_out.toFixed(0)} °C · hot side ${u.T_hot.toFixed(0)} °C · Carnot ${u.carnot_COP.toFixed(2)} · ${(100 * u.q_air / r.q_air).toFixed(0)} % of the heat`)));
  else if (c.cyc === 'Y') cells.push(
    fig('Pressures', `${p.P_low_bar.toFixed(1)} → ${p.P_high_bar.toFixed(1)}`, 'bar', `PR ${p.PR.toFixed(2)} · discharge ${p.T_dis.toFixed(0)} °C`),
    fig('Solution', (100 * p.w_rich).toFixed(1), '% NH₃', `vapour ${(100 * p.w_vapour).toFixed(1)} % · lean liquid ${(100 * p.w_lean).toFixed(1)} % · degassing ${p.q_des.toFixed(2)}`),
    fig('Absorber', `${p.T_abs_in.toFixed(0)} → ${p.T_abs_out.toFixed(1)}`, '°C', `air approach ${p.gap_air.toFixed(1)} K · glide ${p.T_dew_high.toFixed(0)} → ${p.T_bubble_high.toFixed(0)} °C`),
    fig('Desorber', `${p.T_des_in.toFixed(1)} → ${p.T_des_out.toFixed(1)}`, '°C', `source approach ${p.gap_source.toFixed(1)} K`),
    fig('Wet compression', p.wet_compression ? (100 * p.m_injected).toFixed(1) : 'off', p.wet_compression ? '%' : '', p.wet_compression ? 'lean liquid injected, per kg of solution' : 'dry compression of the vapour'),
    fig('Work split', (100 * p.w_pump / (p.w_comp + p.w_pump)).toFixed(1), '% pump', `compressor ηis ${p.eta_comp.toFixed(3)}`));
  else if (c.cyc === 'X') cells.push(
    fig('Pressures', `${p.P_L_bar.toFixed(0)} → ${p.P_H_bar.toFixed(1)}`, 'bar', `PR ${p.PR.toFixed(2)} · architecture ${p.arch}`),
    fig('Shaft COP', p.COP_shaft.toFixed(3), '', 'heat to air / net shaft work, as in the paper'),
    fig('Compressor discharge', `${p.T_c1_out.toFixed(0)} / ${p.T_c2_out.toFixed(0)}`, '°C', `1 / 2 · split ${p.x.toFixed(3)} · ηis ${p.eta_c}`),
    fig('Suction', p.T_suction.toFixed(1), '°C', `compressor 2 inlet ${p.T_c2_in.toFixed(1)} °C after the preheater`),
    fig('Turbine', `${p.T_turb_in.toFixed(0)} → ${p.T_turb_out.toFixed(0)}`, '°C', `gives back ${p.w_turb.toFixed(1)} kJ/kg, ${(100 * p.w_turb / (p.w_c1 + p.w_c2)).toFixed(0)} % of the compression`),
    fig(p.arch === 2 ? 'HTHX 1' : 'HTHX', p.gap_hthx.toFixed(1), 'K', `approach · CO₂ out ${p.T_hthx_out.toFixed(1)} °C`),
    ...(p.arch === 2 ? [fig('HTHX 2', p.gap_hthx2 != null ? p.gap_hthx2.toFixed(1) : '–', 'K', `approach · air to ${p.T_air_mid.toFixed(1)} °C · CO₂ out ${p.T_hthx2_out.toFixed(1)} °C`)] : []),
    fig('IHX / preheater', `${p.q_ihx.toFixed(0)} / ${p.q_pre.toFixed(0)}`, 'kJ/kg', `effectiveness ${p.eps.ihx} / ${p.eps.pre}`),
    fig('LTHX', p.gap_lthx.toFixed(1), 'K', `approach to the source (set ${p.dT_src} K) · CO₂ ${p.T_turb_out.toFixed(1)} → ${p.T_lthx_out.toFixed(1)} °C${r.srcLimited ? ' · duty limited by the source flow' : ''}`));
  else cells.push(
    fig('Bottom cycle', `${p.P_evap_bar.toFixed(1)} → ${p.P_cond_bar.toFixed(1)}`, 'bar', `PR ${p.PR_bottom.toFixed(2)} · discharge ${p.T_dis_bottom.toFixed(0)} °C`),
    fig('Steam raised', p.T_steam.toFixed(1), '°C', `${p.P_steam_bar.toFixed(3)} bar${p.P_steam_bar < 1.013 ? ' (vacuum)' : ''}`),
    fig('After MVR', p.P_final_bar.toFixed(2), 'bar', `saturation ${p.T_sat_final.toFixed(1)} °C · steam ${p.T_steam_final.toFixed(0)} °C`),
    fig('Steam per kg bottom', p.steam_per_kg.toFixed(3), 'kg', `${p.steam_final_per_kg.toFixed(3)} kg after injection`),
    fig('Work split', (100 * p.w_mvr / (p.w_mvr + p.w_bottom)).toFixed(0), '% MVR', `fan ηis ${p.eta_mvr.toFixed(3)} at PR ${(+p.PR_mvr).toFixed(2)}`),
    fig('Process steam', p.stages.filter(s => s.to_process > 1e-9).length, `of ${p.stages.length} fans`, p.stages.filter(s => s.to_process > 1e-9).map((s, k) => `${s.T_sat.toFixed(0)} °C`).join(' · ')),
    fig('Air after the bottom cycle', p.T_air_mid.toFixed(1), '°C', `approach ${p.gap_bottom.toFixed(1)} K · steam ${p.gap_steam.toFixed(1)} K${p.sub_raise > 0.05 ? ` · subcooler outlet raised ${p.sub_raise.toFixed(1)} K` : ''}${p.t1_drop > 0.05 ? ` · ${p.t1_drop.toFixed(1)} K below condensing − approach, to keep the subcooler approach` : ''}`));
  if (c.cyc !== 'K' && secOn(c)) cells.push(fig('Secondary air', (100 * secShare(c)).toFixed(1), '%', `of the heat · ${(+c.phi).toFixed(2)} kg per kg primary to ${c.Tmt} °C`));
  cells.push(fig('Solve time', (r.ms / 1000).toFixed(1), 's', `CoolProp ${CP_VERSION}`));
  while (cells.length % 4) cells.push(`<div class="fig"></div>`);
  figs.innerHTML = cells.join('');
}
function lvHeadlineMulti() {
  const host = $('lvHeadline'), r = LIVE.res, c = LIVE.cfg;
  if (!r || !r.success || !c) { host.hidden = true; return; }
  const u = energyG(r, c);
  host.hidden = false;
  const heat = u.kw.air, el = u.kw.el, kwTop = niceTop(Math.max(heat, el), 10), n = Math.max(0, Math.min(5, r.COP));
  const water = c.Two != null ? u.kw.ev / (4.18 * (c.Ts - c.Two)) : null;
  const air = (r.air || []).map(a => `${a.name} ${hfmt(u.m * a.q, 1)} kW${a.to ? ` to ${hfmt(a.to, 0)} °C` : ''}`).join(' · ');
  host.innerHTML = `<div class="hl-txt hl-cop"><span class="hl-k">Heating COP</span><span class="hl-v">${hfmt(r.COP, 2)}</span>
      <span class="stars" role="img" aria-label="${n.toFixed(2)} of 5 stars, one star per unit of COP">${starRowL(n, 'lvst')}</span></div>
    <div class="dial">${dialL(heat, kwTop, '--s4', 'kW heat to air')}<span class="dk">Heat to air, kW</span></div>
    <div class="dial">${dialL(el, kwTop, '--s2', 'kW electrical')}<span class="dk">Electrical, kW</span></div>
    ${isFinite(r.Pmax_bar) ? `<div class="dial">${dialL(r.Pmax_bar, niceTop(r.Pmax_bar, 10), '--s6', 'bar maximum pressure')}<span class="dk">Max pressure, bar</span></div>` : ''}
    <div class="dial">${dialL(r.T_dis_max, niceTop(r.T_dis_max, 50), '--s5', '°C discharge', 0)}<span class="dk">${c.cyc === 'R' ? 'Hot side' : 'Discharge'}, °C</span></div>
    <span class="hl-s">${vswTxt()} m³/h ${c.cyc === 'R' ? 'of dryer air at the inlet' : `swept by ${c.cyc === 'X' ? 'both compressors' : ['E', 'D', 'Y'].includes(c.cyc) ? 'all compressors' : 'the bottom compressor'}`}${air ? ` · ${air}` : ''} · ${srcLine(c, r, u.kw.ev, water)}</span>`;
}
/* diagrams of every loop: each loop's own dome, the loops told apart by dash */
function genericDiagrams(src) {
  const r = src.res;
  if (!r || !r.success) { for (const id of ['dTs', 'dPh']) xy(id, [], {}); document.getElementById('lgComp').innerHTML = ''; return; }
  const dashes = ['', '6 4', '2 3', '9 3 2 3'];
  const make = (ix, iy) => {
    const L = [], seen = new Set();
    for (const lp of r.loops) {
      const f = loopFluid(r, lp);
      if (f && f.domeL && !seen.has(f)) { seen.add(f); L.push({pts: [...f.domeL.bubble, ...f.domeL.dew.slice().reverse()].map(p => [p[ix], p[iy]]), color: 'var(--dome)', width: 1.1, fill: 'var(--dome-fill)', fit: false}); }
    }
    r.loops.forEach((lp, li) => { for (const path of lp.paths) if (path.pts.length > 1) L.push({pts: path.pts.map(p => [p[ix], p[iy]]), raw: path.pts, color: CCOL[compGroup(path.name)][1], width: 2.2, dash: dashes[li % 4], name: lp.name, comp: path.name, tip: true}); });
    return L;
  };
  xy('dTs', make(3, 0), {xlabel: 's (kJ/kg·K)', ylabel: 'T (°C)'});
  xy('dPh', make(2, 1), {xlabel: 'h (kJ/kg)', ylabel: 'P (bar, log)', ylog: true});
  const used = new Set(r.loops.flatMap(l => l.paths.map(p => compGroup(p.name))));
  document.getElementById('lgComp').innerHTML = Object.entries(CCOL).filter(([g]) => used.has(g)).map(([, [lab, col]]) => `<span><i class="lg-line" style="border-color:${col}"></i>${lab}</span>`).join('') +
    r.loops.map((l, i) => `<span><i class="lg-line" style="border-color:var(--ink-2);border-top-style:${i ? (i === 1 ? 'dashed' : 'dotted') : 'solid'}"></i>${l.name}</span>`).join('') + '<span><i class="lg-line" style="border-color:var(--dome)"></i>bubble / dew lines</span>';
}
const GHXCOL = {ihx: '#8c564b', gc: '#ff7f0e', gc2: '#e3a33b', gc3: '#d62728', gc4: '#bcbd22', gc5: '#7f7f7f', rhx: '#9467bd', evap: '#1f77b4', evap2: '#17becf'};
function genericProfileData(src) {
  const r = src.res; if (!r || !r.success) return null;
  const u = energyG(r, src.cfg), m = u.m;
  const sc = pts => pts.map(p => [p[0] * m, p[1]]);
  const S = r.hx.map(h => h.coldRef
    ? {key: h.key, name: h.name, col: GHXCOL[h.key] || '#1f77b4', ref: sc(h.cold), other: sc(h.hot), otherName: h.otherName, cold: true}
    : {key: h.key === 'gc' ? 'gcx' : h.key, name: h.name, col: GHXCOL[h.key] || '#ff7f0e', ref: sc(h.hot), other: sc(h.cold), otherName: h.otherName});
  // the dryer air passes its exchangers in series: stack them along the heat axis, coldest air first, so the air
  // reads as one continuous line (each starts where the one before it left the air)
  const airHx = S.filter(x => !x.cold && /air/i.test(x.otherName || '') && x.other.length);
  airHx.sort((p, q) => Math.min(...p.other.map(v => v[1])) - Math.min(...q.other.map(v => v[1])));
  let q0 = 0;
  for (const x of airHx) {
    const lo = Math.min(...x.ref.map(v => v[0]), ...x.other.map(v => v[0])), hi = Math.max(...x.ref.map(v => v[0]), ...x.other.map(v => v[0]));
    const d = q0 - lo;
    if (Math.abs(d) > 1e-9) { x.ref = x.ref.map(v => [v[0] + d, v[1]]); x.other = x.other.map(v => [v[0] + d, v[1]]); }
    q0 += hi - lo;
  }
  const at = (line, q) => { for (let i = 1; i < line.length; i++) if (line[i][0] >= q - 1e-9) { const [a, b] = [line[i - 1], line[i]]; return b[0] === a[0] ? b[1] : a[1] + (b[1] - a[1]) * (q - a[0]) / (b[0] - a[0]); } return line[line.length - 1][1]; };
  for (const s of S) s.approach = Math.min(...s.ref.map(([q, T]) => (s.cold ? -1 : 1) * (T - at(s.other, q))));
  return {S, u};
}
function genericPFD(src) {
  const r = src.res;
  if (!r || !r.success) return `<div class="panel"><h3>${src.title}</h3><p class="lbl">no design</p></div>`;
  const u = energyG(r, src.cfg), m = u.m;
  const rows = r.loops.map(lp => `<tr><td colspan="4" class="lbl" style="padding-top:8px">${lp.name} · ${hfmt(lp.flow * m * 3600, 0)} kg/h</td></tr>` + lp.paths.filter(p => p.name !== 'Expansion valve' && p.pts.length > 1).map(p => {
    const a = p.pts[0], b = p.pts[p.pts.length - 1], q = lp.flow * m * Math.abs(b[2] - a[2]);
    return `<tr><td>${p.name}</td><td>${f1(a[0])} → ${f1(b[0])}</td><td>${a[1].toFixed(2)}${Math.abs(b[1] - a[1]) > 1e-3 ? ' → ' + b[1].toFixed(2) : ''}</td><td>${fk(q)}</td></tr>`;
  }).join('')).join('');
  return `<div class="panel" style="display:grid;gap:6px"><h3>${src.title}</h3><p class="lbl" style="margin:0">${src.label} · ${src.Tair} °C air · COP ${r.COP.toFixed(3)} · CoolProp</p>
    <figure class="mini" style="max-width:520px">${miniMulti(src.k, src.cfg)}</figure>
    <table class="sumt"><thead><tr><th>Step</th><th>T, °C</th><th>P, bar</th><th>kW</th></tr></thead><tbody>${rows}</tbody></table></div>`;
}
/* expander result: work per kg of the reference flow and the outlet state (two-phase quality where it flashes) */
function expFig(p, fig) {
  const q = [].concat(p.exp_out_quality).filter(v => v != null);
  return fig('Expander work', p.w_exp_kJkg.toFixed(1), 'kJ/kg', q.length ? `two-phase outlet · ${q.map(v => (100 * v).toFixed(0)).join(' / ')} % vapour` : 'single-phase outlet');
}
function genericSummary(src) {
  const r = src.res; if (!r || !r.success) return '';
  const u = energyG(r, src.cfg), row = (lab, v, unit) => `<tr><td>${lab}</td><td>${v}</td><td class="lbl">${unit}</td></tr>`;
  const sa = r.srcAir, water = src.Twout != null && !sa ? u.kw.ev / (4.18 * (src.Tsrc - src.Twout)) : null;
  return `<div class="panel sumpanel"><h3>${src.name} · ${vswTxt()} m³/h</h3><table class="sumt"><tbody>
    ${row('Heat to air', fk(u.kw.air), 'kW')}${(r.air || []).map(a => row(`&nbsp;&nbsp;${a.name}`, fk(u.m * a.q), 'kW')).join('')}
    ${row('Electrical input', fk(u.kw.el), 'kW')}${row('Motor losses', fk(u.kw.loss), 'kW')}${u.kw.ex > 0 ? row('Expander work', fk(u.kw.ex), 'kW') : ''}${row(sa ? 'Heat from source air' : 'Heat from water', fk(u.kw.ev), 'kW')}
    ${row('COP', r.COP.toFixed(3), '')}${r.loops.map(l => row(`${l.name} flow`, (l.flow * u.m * 3.6).toFixed(2), 't/h')).join('')}
    ${water != null ? row('Source water flow', (water * 3.6).toFixed(1), 't/h') : ''}
    ${sa ? row('Source air flow / outlet', `${sa.m.toFixed(1)} kg/s / ${sa.Tout.toFixed(1)} °C`, `condensate ${(sa.cond * 3600).toFixed(1)} kg/h`) : ''}
    ${row('Bottom compressor PR / ηvol', `${r.ref.PR.toFixed(2)} / ${u.nv.toFixed(2)}`, state.nvMode === 'fixed' ? 'fixed' : 'model')}
  </tbody></table></div>`;
}
/* small schematics of the multi-loop cycles */
function miniMulti(v, cfg = LV) {
  const ar = (x, y, dir) => { const d = {r: [-7, -4, -7, 4], l: [7, -4, 7, 4], u: [-4, 7, 4, 7], d: [-4, -7, 4, -7]}[dir];
    return `<polygon points="${x},${y} ${x + d[0]},${y + d[1]} ${x + d[2]},${y + d[3]}" fill="currentColor"/>`; };
  const T = (x, y, s, cls = 'cn', a = 'middle') => `<text x="${x}" y="${y}" text-anchor="${a}" class="${cls}">${s}</text>`;
  const comp = (x, y, l) => `<polygon class="eq" points="${x - 13},${y - 16} ${x + 13},${y - 8} ${x + 13},${y + 8} ${x - 13},${y + 16}"/>${T(x, y + 4, l)}`;
  const box = (x, y, w, h, l) => `<rect class="eq" x="${x}" y="${y}" width="${w}" height="${h}" rx="3"/>${T(x + w / 2, y + h / 2 + 4, l)}`;
  const valve = (x, y) => `<polygon class="eq" points="${x - 6},${y - 8} ${x + 6},${y - 8} ${x - 6},${y + 8} ${x + 6},${y + 8}"/>`;
  const hp = d => `<path class="ph" d="${d}"/>`, lp = d => `<path class="pl" d="${d}"/>`;
  let body = '', W = 440, H = 220, cap = '';
  if (v === 'K') {
    // top loop (upper half): HT GC then MT GC on the air; bottom loop (lower half) feeds only the CGC
    body = `<g transform="translate(0 30)"><g style="color:var(--hp)">${hp('M60 34 V22 H200')}${ar(200, 22, 'r')}${hp('M260 22 H300')}${ar(300, 22, 'r')}${hp('M360 22 H396 V58')}${hp('M60 150 V140 H170')}${ar(170, 140, 'r')}${hp('M214 140 H396 V170')}</g>
      <g style="color:var(--lp)">${lp('M396 74 V96 H330')}${ar(330, 96, 'l')}${lp('M270 96 H214')}${lp('M170 96 H60 V66')}${ar(60, 66, 'u')}${lp('M396 186 V204 H330')}${ar(330, 204, 'l')}${lp('M270 204 H60 V182')}${ar(60, 182, 'u')}</g>
      ${comp(60, 50, 'C top')}${comp(60, 166, 'C bot')}${box(200, 10, 60, 24, 'HT GC')}${box(300, 10, 60, 24, 'MT GC')}
      <rect class="eq" x="170" y="84" width="44" height="68" rx="3"/>${T(192, 122, 'CGC')}${valve(396, 66)}${valve(396, 178)}
      ${box(270, 84, 60, 24, 'Evap top')}${box(270, 192, 60, 24, 'Evap bot')}
      <path class="wat" d="M300 48 V84"/><path class="wat" d="M300 108 V192"/><path class="wat" d="M300 216 V232 H236"/><g style="color:var(--accent)">${ar(300, 84, 'd')}${ar(300, 192, 'd')}${ar(236, 232, 'l')}</g>${T(306, 56, 'water', 'lbl2', 'start')}
      <path class="air" d="M436 -6 H180"/><path class="air" d="M290 -6 V-26"/>${T(436, -12, 'primary + secondary', 'lbl2', 'end')}${T(182, -12, 'primary', 'lbl2', 'start')}${T(284, -18, 'secondary', 'lbl2', 'end')}</g>`;
    H = 270;
    cap = 'TT cascade: the bottom cycle only superheats the top suction in the cascade gas cooler (CGC); the top cycle heats primary and secondary air together in the MT gas cooler, then the primary air alone in the HT gas cooler; an IHX in each cycle (not drawn); the source water passes the top evaporator, then the bottom one.';
  } else if (v === 'N') {
    const n = (cfg.par || LV.par).N.n, rowH = 72; H = 30 + n * rowH;
    for (let k = 0; k < n; k++) {
      const y = H - 30 - k * rowH, top = k === n - 1;
      body += `${comp(60, y - 18, `L${k + 1}`)}<g style="color:var(--hp)">${hp(`M60 ${y - 34} V${y - 44} H${top ? 262 : 166}`)}${ar(top ? 262 : 166, y - 44, 'r')}</g>`;
      // cascade exchanger k→k+1: level k's discharge condenses (at y − 44) against level k+1's evaporating flow,
      // whose return line runs at y − rowH + 6; the box spans just those two lines, between the levels
      body += top ? box(262, y - 56, 70, 24, 'Air cond.') : `<rect class="eq" x="166" y="${y - rowH - 3}" width="52" height="${rowH - 30}" rx="3"/>${T(192, y - 50, `C${k + 1}→${k + 2}`)}`;
      body += `<g style="color:var(--lp)">${lp(`M${top ? 332 : 218} ${y - 44} H396 V${y - 18}`)}${lp(`M396 ${y - 2} V${y + 6} H${k === 0 ? 330 : 214}`)}${lp(`M${k === 0 ? 270 : 170} ${y + 6} H60 V${y}`)}${ar(60, y, 'u')}</g>${valve(396, y - 10)}`;
      if (k === 0) body += box(262, y - 6, 76, 24, 'Evaporator');
      if (!top && (cfg.par || LV.par).N.sub_air[k]) body += T(232, y - 48, 'subcool + air', 'lbl2', 'start');
      // IHX: hot side on the liquid line before the valve, cold side on the suction line, joined by a dashed tie
      if (((cfg.par || LV.par).N.ihx || [])[k]) {
        const xs = k === 0 ? 210 : 100, xh = 342;
        body += `<path d="M${xh + 18} ${y - 34} V${y - 22} M${xs + 15} ${y - 4} V${y - 22} H${xh + 18}" fill="none" stroke="var(--ink-3)" stroke-width="1" stroke-dasharray="2 2"/>`
          + box(xh, y - 54, 36, 20, 'IHX') + box(xs, y - 4, 30, 20, 'IHX');
      }
    }
    body += `<path class="air" d="M430 ${H - 4} V8"/>${T(426, 16, 'air', 'lbl2', 'end')}`;
    cap = `${n}-level cascade: each level evaporates in the condenser of the one below (level 1 on the source water); the air passes the levels in series and the top level condenses into it${((cfg.par || LV.par).N.ihx || []).slice(0, n).some(Boolean) ? '; an IHX heats a level\'s suction with its liquid before the valve' : ''}.`;
  } else if (v === 'X') {
    // sCO2 reverse Brayton: C2 -> HTHX (air) -> IHX hot -> turbine -> LTHX (source) -> IHX cold -> split: C1 (secondary)
    // -> preheater hot -> [arch 2: HTHX 2] -> rejoins; main -> preheater cold -> C2
    const a2 = +((cfg.par || LV.par).X.arch) === 2;
    H = 236;
    body = `<g transform="translate(0 14)"><g style="color:var(--hp)">${hp('M60 34 V22 H140')}${ar(140, 22, 'r')}${hp('M200 22 H330')}${ar(330, 22, 'r')}${hp('M370 22 H404 V38')}
        ${hp('M250 180 V196 H158 V122')}${ar(158, 122, 'u')}
        ${a2 ? `${hp('M158 98 V80 H216')}${ar(216, 80, 'r')}${hp('M268 80 H300 V2 H392 V22')}${ar(392, 22, 'd')}` : `${hp('M158 98 V80 H300 V22')}${ar(300, 26, 'u')}`}</g>
      <g style="color:var(--lp)">${lp('M404 74 V190 H382')}${ar(382, 190, 'l')}${lp('M354 178 V126')}${ar(354, 126, 'u')}${lp('M330 110 H250')}${lp('M250 110 H190')}${ar(190, 110, 'l')}${lp('M130 110 H60 V70')}${ar(60, 70, 'u')}${lp('M250 110 V148')}${ar(250, 148, 'd')}</g>
      ${comp(60, 52, 'C2')}${comp(250, 164, 'C1')}${box(140, 10, 60, 24, a2 ? 'HTHX 1' : 'HTHX')}${box(130, 98, 60, 24, 'Preheater')}
      <rect class="eq" x="330" y="10" width="40" height="112" rx="3"/>${T(350, 70, 'IHX')}${box(330, 178, 52, 24, 'LTHX')}
      <polygon class="eq" points="${404 - 13},${40} ${404 + 13},${48} ${404 + 13},${64} ${404 - 13},${72}"/>${T(404, 60, 'T')}
      ${a2 ? `${box(216, 68, 52, 24, 'HTHX 2')}<path class="air" d="M242 60 V48"/>${T(246, 52, 'cool air', 'lbl2', 'start')}` : ''}
      <path class="air" d="M300 -8 H150"/><g style="color:var(--ink-3)">${ar(150, -8, 'l')}</g>${T(304, -4, 'air', 'lbl2', 'start')}
      <path class="wat" d="M382 214 H330"/><g style="color:var(--accent)">${ar(330, 214, 'l')}</g>${T(386, 218, 'source', 'lbl2', 'start')}
      ${T(262, 106, 'split', 'lbl2', 'start')}</g>`;
    cap = a2 ? 'sCO₂ Brayton, architecture 2: the secondary stream (C1) preheats the main stream before C2, then heats the cooler air in HTHX 2 and rejoins before the turbine (T); the main stream heats the hotter air in HTHX 1 and gives heat back to the suction in the IHX; the LTHX takes heat from the source.'
      : 'sCO₂ Brayton, architecture 1: the secondary stream (C1) preheats the main stream before C2; C2 heats the air in the HTHX; the two streams rejoin, give heat back to the suction in the IHX and expand in the turbine (T); the LTHX takes heat from the source.';
  } else if (v === 'E') {
    const arr = ECON_ARR[(cfg.par || LV.par).E.arr] || 'flash', single = arr === 'single', par = arr === 'parallel', eco = arr === 'economiser';
    H = 236;
    body = `<g transform="translate(0 8)"><g style="color:var(--hp)">${hp(`M60 ${single ? 150 : 44} V22 H210`)}${ar(210, 22, 'r')}${hp('M290 22 H396 V60')}
        ${single ? '' : hp('M60 150 V112')}${par ? `${hp('M190 92 V22')}${ar(190, 30, 'u')}` : ''}</g>
      <g style="color:var(--lp)">${single ? lp('M396 76 V204 H330') : `${lp('M396 76 V96')}${ar(396, 96, 'd')}${lp(`M${eco ? 396 : 396} 132 V170`)}${lp('M396 186 V204 H330')}`}${ar(330, 204, 'l')}${lp('M270 204 H60 V182')}${ar(60, 182, 'u')}
        ${single ? '' : `${lp(`M370 110 H${par ? 204 : 74}`)}${ar(par ? 204 : 74, 110, 'l')}${par ? '' : lp('M60 112 V76')}`}</g>
      ${single ? comp(60, 166, 'C') : `${comp(60, 166, 'C1')}${comp(60, 60, 'C2')}`}${par ? comp(190, 108, 'C3') : ''}
      ${box(210, 10, 80, 24, single ? 'Condenser' : 'Condenser')}${valve(396, 68)}${single ? '' : box(370, 96, 52, 36, eco ? 'Econ.' : 'Flash')}${valve(396, 178)}${box(270, 192, 60, 24, 'Evaporator')}
      ${single || par ? '' : T(70, 106, 'mix', 'lbl2', 'start')}${single ? '' : T(300, 104, eco ? 'injection' : 'flash vapour', 'lbl2', 'middle')}
      <path class="air" d="M300 -2 H200"/><g style="color:var(--ink-3)">${ar(200, -2, 'l')}</g>${T(304, 2, 'air', 'lbl2', 'start')}
      <path class="wat" d="M330 222 H270"/><g style="color:var(--accent)">${ar(270, 222, 'l')}</g>${T(334, 226, 'source', 'lbl2', 'start')}</g>`;
    cap = {flash: 'Flash tank: the condenser liquid flashes at the intermediate pressure; the flash vapour joins the low-stage (C1) discharge at the high-stage (C2) suction and the saturated liquid expands to the evaporator.',
      economiser: 'Closed economiser: part of the condenser liquid expands to the intermediate pressure and evaporates against the rest (subcooling it); the vapour is injected at the high-stage (C2) suction.',
      parallel: 'Parallel compression: the flash vapour has its own compressor (C3) to the high pressure; the low stage (C1) compresses the evaporator vapour straight to it.',
      single: 'Single stage: one compressor, condenser, valve and evaporator, for comparison.'}[arr];
  } else if (v === 'D') {
    // n levels from the bottom: compressor k discharges at y_k into condenser k (and on up to compressor k + 1);
    // each condenser's liquid drops through a valve into the one below; the bottom one feeds the evaporator
    const n = (cfg.par || LV.par).D.n, rowH = 46, base = 56 + n * rowH;
    H = base + 46;
    const y = k => base - 24 - k * rowH;
    for (let k = 0; k < n; k++) {
      body += `<g style="color:var(--hp)">${hp(`M60 ${y(k) - 2} V${y(k) - 16} H210`)}${ar(210, y(k) - 16, 'r')}${k < n - 1 ? hp(`M60 ${y(k) - 16} V${y(k + 1) + 16}`) : ''}</g>`;
      body += `${comp(60, y(k) + 14, `C${k + 1}`)}${box(210, y(k) - 28, 70, 24, `Cond. ${k + 1}`)}`;
      if (k > 0) body += `<g style="color:var(--lp)">${lp(`M280 ${y(k) - 16} H330 V${y(k - 1) - 16} H282`)}${ar(282, y(k - 1) - 16, 'l')}</g>${valve(330, (y(k) + y(k - 1)) / 2 - 16)}`;
    }
    body += `<g style="color:var(--lp)">${lp(`M280 ${y(0) - 10} H396 V${base + 2}`)}${lp(`M396 ${base + 18} V${base + 26} H330`)}${ar(330, base + 26, 'l')}${lp(`M270 ${base + 26} H60 V${y(0) + 32}`)}${ar(60, y(0) + 32, 'u')}</g>
      ${valve(396, base + 10)}${box(270, base + 14, 60, 24, 'Evaporator')}
      <path class="air" d="M190 ${y(0)} V${y(n - 1) - 34}"/><g style="color:var(--ink-3)">${ar(190, y(n - 1) - 34, 'u')}</g>${T(186, y(n - 1) - 30, 'air', 'lbl2', 'end')}`;
    cap = `${n}-pressure condensation: compressor 1 feeds condenser 1, which heats the air first; at each pressure part of the flow goes on to the next compressor and condenser, and each condenser's liquid is throttled into the one below.${n === 2 ? ' This is the dual-pressure cycle of Dai et al. (2020).' : ''}`;
  } else if (v === 'R') {
    const n = (cfg.par || LV.par).R.n, w = 300 / n;
    H = 170;
    body = `${Array.from({length: n}, (_, k) => box(70 + k * w, 50, w - 20, 50, `Unit ${k + 1}`)).join('')}
      <path class="air" d="M40 30 H410"/><g style="color:var(--ink-3)">${ar(410, 30, 'r')}</g>${T(40, 22, 'air in', 'lbl2', 'start')}${T(410, 22, 'to the dryer', 'lbl2', 'end')}
      <path class="wat" d="M410 128 H40"/><g style="color:var(--accent)">${ar(40, 128, 'l')}</g>${T(410, 146, 'source (every unit\'s cold side)', 'lbl2', 'end')}
      ${Array.from({length: n}, (_, k) => `<path class="ph" d="M${60 + (k + 0.5) * w} 50 V30"/><path class="pl" d="M${60 + (k + 0.5) * w} 100 V128"/>`).join('')}`;
    cap = `Reverse-Stirling: ${n > 1 ? `${n} units in series on the air, each heating one band of it` : 'one unit heats the air'} at a hot side above the air it delivers, taking heat from the source at its cold side; each unit's COP is a fraction of its Carnot COP.`;
  } else if (v === 'Y') {
    H = 236;
    body = `<g transform="translate(0 8)"><g style="color:var(--hp)">${hp('M60 150 V22 H220')}${ar(220, 22, 'r')}${hp('M150 160 V70 H220')}${ar(220, 70, 'r')}</g>
      <g style="color:var(--lp)">${lp('M300 46 H396 V170')}${lp('M396 186 V204 H330')}${ar(330, 204, 'l')}${lp('M270 204 H230')}${ar(230, 204, 'l')}${lp('M190 204 H60 V182')}${ar(60, 182, 'u')}${lp('M190 204 H150 V176')}${ar(150, 176, 'u')}
        <path class="pl" style="stroke-dasharray:3 3" d="M150 120 H74"/>${ar(74, 120, 'l')}</g>
      ${comp(60, 166, 'C')}<circle class="eq" cx="150" cy="168" r="10"/>${T(150, 172, 'P')}${box(220, 10, 80, 72, 'Absorber')}${valve(396, 178)}${box(270, 192, 60, 24, 'Desorber')}${box(190, 192, 40, 24, 'Sep.')}
      ${T(152, 116, 'wet', 'lbl2', 'start')}${T(66, 196, 'vapour', 'lbl2', 'start')}${T(160, 186, 'lean liquid', 'lbl2', 'start')}
      <path class="air" d="M340 -4 H230"/><g style="color:var(--ink-3)">${ar(230, -4, 'l')}</g>${T(344, 0, 'air', 'lbl2', 'start')}
      <path class="wat" d="M330 222 H270"/><g style="color:var(--accent)">${ar(270, 222, 'l')}</g>${T(334, 226, 'source', 'lbl2', 'start')}</g>`;
    cap = 'Ammonia–water hybrid: the solution boils partly in the desorber on the source; the separator sends the ammonia-rich vapour to the compressor (C, with lean liquid injected to hold the discharge limit) and the lean liquid to the pump (P); both are absorbed again in the absorber, which heats the air, and the solution returns through the valve.';
  } else {
    body = `<g style="color:var(--hp)">${hp('M60 150 V140 H170')}${ar(170, 140, 'r')}${hp('M230 140 H300')}${ar(300, 140, 'r')}${hp('M360 140 H396 V170')}</g>
      <g style="color:var(--lp)">${lp('M396 186 V204 H330')}${ar(330, 204, 'l')}${lp('M270 204 H60 V182')}${ar(60, 182, 'u')}</g>
      ${comp(60, 166, 'C')}${box(170, 128, 60, 24, 'Steam gen')}${box(300, 128, 60, 24, 'Air cond.')}${valve(396, 178)}${box(270, 192, 60, 24, 'Evaporator')}
      <g style="color:#7f7f7f">${hp('M200 128 V60 H230').replace('class="ph"', 'class="ph" style="stroke:#7f7f7f"')}</g>
      ${[0, 1, 2].map(i => comp(250 + i * 40, 60, i === 2 ? 'n' : String(i + 1))).join('')}${box(360, 30, 64, 24, 'Steam cond.')}
      <path class="ph" style="stroke:#7f7f7f" d="M330 52 V42 H360"/><path class="air" d="M430 8 H244"/>${T(240, 12, 'air', 'lbl2', 'end')}`;
    cap = `HP + MVR: the heat pump heats the air first and raises (vacuum) steam with the rest of its condenser heat; MVR fans at a pressure ratio of ${(+(cfg.fanPR ?? 1.5)).toFixed(2)} each compress the steam. After each fan a share of the exhaust can go to the process, condensing into the air at that pressure; the last fan's steam heats the air to the target.`;
    if ((cfg.par || LV.par).M.ihx) body += `<path d="M376 150 V160 M165 194 V160 H376" fill="none" stroke="var(--ink-3)" stroke-width="1" stroke-dasharray="2 2"/>${box(362, 130, 28, 20, 'IHX')}${box(150, 194, 30, 20, 'IHX')}`;
    body += `<g style="color:#7f7f7f">${[0, 1].map(i => `<path class="ph" style="stroke:#7f7f7f;stroke-dasharray:3 3" d="M${270 + i * 40} 43 V12"/>${ar(270 + i * 40, 12, 'u')}`).join('')}</g>${T(290, 32, 'split', 'lbl2')}`;
  }
  return `<svg viewBox="0 0 ${W} ${H}" role="img" aria-label="${cap}">${body}</svg><figcaption>${cap}</figcaption>`;
}

/* full name of the cycle and its refrigerant(s), above the headline */
const FULL_NAME = {
  H: 'Conventional heat pump', H1: 'Conventional single-stage heat pump', H2: 'Conventional two-stage heat pump',
  A: 'Recuperative heat pump A (RHX1 on the stage-1 suction)', B: 'Recuperative heat pump B (RHX2 on the interstage gas, with IHX)',
  C: 'Recuperative heat pump C (RHX1 and RHX2)', S: 'Split-stream heat pump (SLTHP: RHX1 parallel to the LT gas cooler)',
  K: 'Transcritical–transcritical cascade heat pump (Kong et al.)', N: 'cascade heat pump (Hoang et al.)', M: 'Heat pump with MVR steam', X: 'Supercritical CO₂ reverse Brayton heat pump (Dalal et al.)',
  E: 'Two-stage heat pump', D: 'condensation heat pump', R: 'Reverse-Stirling heat pump units (semi-empirical)', Y: 'Ammonia–water hybrid absorption–compression heat pump',
};
const ECON_NAME = {flash: 'with a flash tank', economiser: 'with a closed economiser (vapour injection)', parallel: 'with parallel compression', single: 'single stage, for comparison'};
function lvTitle() {
  const c = LIVE.cfg && LIVE.cfg.cyc === LV.cyc ? LIVE.cfg : LV, v = c.cyc;
  const name = v === 'N' ? `${c.par.N.n}-level ${FULL_NAME.N}` : v === 'D' ? `${c.par.D.n}-pressure ${FULL_NAME.D}` : v === 'E' ? `${FULL_NAME.E} ${ECON_NAME[ECON_ARR[c.par.E.arr] || 'flash']}` : v === 'R' ? `Reverse-Stirling heat pump, ${c.par.R.n > 1 ? `${c.par.R.n} units in series` : 'one unit'} (semi-empirical)` : v === 'H' ? FULL_NAME['H' + (c.stagesH || 1)] + (c.ihxH ? ' with IHX' : '') : FULL_NAME[v];
  const nameOut = v === 'A' && stagesOf(c) === 1 ? name.replace('stage-1 suction', 'suction') + ', single stage' : name;
  const fl = MULTI(v) ? loopNames(v).map((n, i) => `${n.toLowerCase()} ${loopFluidText(c, v, i)}`).join(', ') + (v === 'M' ? ', steam (water)' : '')
    : fluidLabel({comps: c.mainComps || c.comps, basis: c.basis}) + ((c.mainComps || c.comps).filter(x => x[1] > 0).length > 1 && c.basis !== 'mole' ? ' (mass)' : '');
  $('lvTitle').innerHTML = `${nameOut}: <span>${fl}</span>`;
}

