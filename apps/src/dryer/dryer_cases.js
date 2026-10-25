/* ================================================================ cases (after the steam-generating heat pump page)
   Save the live design's inputs as named cases, keep them in a library shared by everyone who opens the page (the
   artifact's db capability; this browser only when that is unavailable), send one as a link, and compare up to four:
   results side by side against a baseline, the inputs that differ, overlaid T-s cycles, and a sweep of one operating
   variable. Every case is re-solved here with CoolProp; a stored case keeps a summary of its results for the list
   and for drift checks. */

const CASE_SCHEMA = 1;
const MODEL_VERSION = "dryer-js 2026-10-13";                 // bump when the cycle models change results
const CASE_COLORS = ["--s1", "--s2", "--s3", "--s4"];
const MAX_COMPARE = 4;
const PAGE_URL = typeof location !== "undefined" && /(\.github\.io|^localhost|^127\.0\.0\.1)$/.test(location.hostname)
  ? location.origin + location.pathname : "https://waikato-ahuora-smart-energy-systems.github.io/HeatPumpLab/dryer/";   // where case links point: this page on GitHub Pages or localhost, else the public GitHub Pages app
const LOCAL_KEY = "dryer-cases-local";
const CASE_KEYS = ["cyc", "mainComps", "basis", "multi", "par", "compModel", "etaC", "motor", "prMax", "T", "Tain", "RH", "Tmt", "phi",
  "a", "Ts", "Two", "lift", "lift2", "split", "liftS", "share", "alpha", "stagesH", "stagesA", "stagesS", "ihxH", "fanModel", "fanEta", "fanPR", "fanSh", "fanMotor", "src", "srcT", "srcM", "srcW", "srcHum", "srcRH", "srcDew", "airHum", "airWg", "airDew", "Patm", "dTs", "srcWm", "srcMw", "expn", "etaX", "pmax"];
const CYC_SHORT = { H: "Conventional", A: "RHP-A", B: "RHP-B", C: "RHP-C", S: "SLTHP", K: "TT cascade", N: "n-cascade", M: "HP + MVR", X: "sCO2 Brayton" };
// numeric inputs: [key, lo, hi]
const CASE_NUM = [["compModel", 1, 5], ["etaC", 0.5, 1], ["motor", 0.8, 1], ["prMax", 5, 50], ["T", 100, 220], ["Tain", 0, 40], ["RH", 0, 100],
  ["Tmt", 30, 200], ["phi", 0, 2], ["a", 0, 40], ["Ts", 0, 80], ["Two", 0, 80], ["lift", 0, 140], ["lift2", 0, 140], ["split", 0.3, 0.7],
  ["liftS", 0, 140], ["share", 0.3, 0.7], ["alpha", 0, 0.95], ["stagesH", 1, 2], ["stagesA", 1, 2], ["stagesS", 1, 2],
  ["fanModel", 1, 5], ["fanEta", 0.5, 1], ["fanPR", 1.1, 2.5], ["fanSh", 0, 50], ["fanMotor", 0.8, 1], ["srcT", 0, 100], ["srcM", 0.1, 200], ["srcW", 0, 1000], ["srcRH", 0, 100], ["srcDew", -30, 100], ["airWg", 0, 1000], ["airDew", -30, 40], ["Patm", 70, 110], ["dTs", 0, 40], ["srcMw", 0.1, 2000], ["etaX", 0.3, 1], ["pmax", 50, 400]];
const SWEEPS = {
  T: { label: "Dryer inlet air target, °C", lo: 120, hi: 220, step: 10 },
  Two: { label: "Source water out, °C", lo: 10, hi: 35, step: 5 },
  phi: { label: "Secondary/primary air", lo: 0, hi: 1, step: 0.1 },
};
const SWEEP_Y = {
  COP: { label: "Heating COP", d: 3, f: s => s.COP },
  Pmax: { label: "Max pressure, bar", d: 1, f: s => s.Pmax },
  VHC: { label: "VHC, MJ/m³", d: 2, f: s => s.VHC },
  Tdis: { label: "Hottest discharge, °C", d: 1, f: s => s.Tdis },
};
// results rows: [key, label, decimals, better: 1 higher, -1 lower, 0 neither]
const RESULT_ROWS = [
  ["COP", "Heating COP", 3, 1], ["VHC", "VHC, MJ/m³", 2, 1], ["Pmax", "Max pressure, bar", 1, -1], ["Psuc", "Lowest suction pressure, bar", 2, 0],
  ["Tdis", "Hottest discharge, °C", 1, -1], ["qair", "Heat to air, kJ/kg ref.", 1, 0], ["wel", "Electricity, kJ/kg ref.", 1, -1],
  ["sec", "Secondary air, % of heat", 1, 0], ["Qkw", "Heat to air, kW", 1, 1], ["Wkw", "Electrical, kW", 1, -1], ["vsw", "Low-stage swept volume, m³/h", 0, -1],
];
const fixed = (v, d) => v == null || !isFinite(v) ? "–" : (+v).toFixed(d);
const humText = (m, rh, w, dew) => m === "w" ? `${fixed(w, 1)} g/kg` : m === "dew" ? `dew ${fixed(dew, 1)} °C` : `${fixed(rh, 0)} % RH`;
const INPUT_ROWS = [
  ["Cycle", i => cycleName(i)], ["Refrigerant", i => caseFluids(i)], ["Dryer inlet air, °C", i => i.T], ["Secondary air target, °C", i => i.Tmt],
  ["Secondary/primary air", i => fixed(i.phi, 2)], ["Air inlet, °C", i => i.Tain], ["Air inlet humidity", i => humText(i.airHum, i.RH, i.airWg, i.airDew)], ["Atmospheric pressure, kPa", i => fixed(i.Patm == null ? 101.325 : i.Patm, 1)],
  ["Source approach, K", i => fixed(i.dTs ?? 5, 1)], ["Heat source", i => i.src === "air" ? `air ${i.srcT} °C · ${fixed(i.srcM, 1)} kg/s · ${humText(i.srcHum || "w", i.srcRH, i.srcW, i.srcDew)}` : (i.srcWm === "flow" ? `water ${i.Ts} °C · ${fixed(i.srcMw, 1)} kg/s` : `water ${i.Ts} → ${i.Two} °C`)], ["Gas-cooler approach, K", i => i.a],
  ["Compressor model", i => i.compModel === 1 ? `Constant ${fixed(i.etaC, 2)}` : DRYER.COMP_MODELS[i.compModel]], ["Motor efficiency", i => fixed(i.motor, 2)],
  ["Max stage pressure ratio", i => i.prMax],
  ["Sized by", i => !i.size ? "–" : i.size.by === "duty" ? `sink duty ${fixed(i.size.value, 1)} kW` : i.size.by === "flow" ? `air flow ${fixed(i.size.value, 2)} kg/s` : `swept volume ${fixed(i.size.value, 0)} m³/h`],
  ["MVR fans", i => i.cyc !== "M" ? "–" : `${+i.fanModel === 1 ? `Constant ${fixed(i.fanEta, 2)}` : DRYER.COMP_MODELS[i.fanModel]} · PR ${fixed(i.fanPR, 2)} · exit +${i.fanSh} K · motor ${fixed(i.fanMotor, 2)}`],
  ["Design", i => designText(i)],
];

const CS = {
  store: null, db: null, user: null, me: null, canWrite: null, downloads: null,
  list: [], compare: [], colors: {}, loaded: null, fromLink: null,
  solved: new Map(),                 // id -> {json, res, sum}
  sweep: { x: "T", y: "COP", data: null, key: "", running: false },
  names: {}, nameAuto: true,
};

/* ---------------------------------------------------------------- inputs */
const clone = o => JSON.parse(JSON.stringify(o));
function caseInputs() {
  const o = {}; for (const k of CASE_KEYS) o[k] = clone(LV[k]);
  o.size = pageSize();
  return o;
}
/* the page's sizing saved with a case: Size by and the value it sets (the swept volume, sink duty or air flow);
   the other two are derived from the design, so they are not stored */
const SIZE_KEY = { vsw: ["vsw", 20, 15600], duty: ["duty", 1, 10000], flow: ["flow", 0.05, 200] };
const pageSize = () => { const by = SIZE_KEY[state.sizeBy] ? state.sizeBy : "vsw"; return { by, value: +state[SIZE_KEY[by][0]] }; };
function cleanSize(raw) {           // cases saved before sizing was stored keep the page's current sizing
  if (!raw || typeof raw !== "object" || !SIZE_KEY[raw.by] || raw.value == null || !isFinite(+raw.value)) return pageSize();
  const [, lo, hi] = SIZE_KEY[raw.by];
  return { by: raw.by, value: Math.min(hi, Math.max(lo, +raw.value)) };
}
const inputsJSON = inp => JSON.stringify(inp);
const okName = n => typeof n === "string" && /^[A-Za-z0-9 \-]{1,40}$/.test(n);
function cleanComps(list) {
  const cs = Array.isArray(list) ? list.filter(c => Array.isArray(c) && okName(c[0]) && isFinite(+c[1])).slice(0, 4) : [];
  if (!cs.length) throw new Error("no refrigerant in this case");
  return cs.map(c => [c[0], Math.min(100, Math.max(0, +c[1]))]);
}
function cleanInputs(raw) {           // anything from the store or a link is untrusted: known fields, in range, defaults otherwise
  if (!raw || typeof raw !== "object") throw new Error("no inputs in this case");
  const o = clone(LV_DEF);
  o.multi = clone(MULTI_DEF); o.par = clone(PAR_DEF);
  if (!CYC_SHORT[raw.cyc]) throw new Error("unknown cycle");
  o.cyc = raw.cyc;
  o.src = raw.src === "air" ? "air" : "water";
  o.srcWm = raw.srcWm === "flow" ? "flow" : "out";
  o.expn = raw.expn === "exp" ? "exp" : "valve";
  o.mainComps = cleanComps(raw.mainComps || raw.comps);
  o.basis = raw.basis === "mole" ? "mole" : "mass";
  const HMS = ["rh", "w", "dew"];
  o.airHum = HMS.includes(raw.airHum) ? raw.airHum : "rh";
  o.srcHum = HMS.includes(raw.srcHum) ? raw.srcHum : (raw.srcW != null ? "w" : "rh");   // cases saved before the RH input gave g/kg
  for (const [k, lo, hi] of CASE_NUM) { const x = +raw[k]; if (raw[k] != null && isFinite(x)) o[k] = Math.min(hi, Math.max(lo, x)); }
  o.compModel = Math.round(o.compModel); o.stagesH = Math.round(o.stagesH);
  for (const k of ['compModel', 'fanModel']) if (Math.round(o[k]) === 4) o[k] = 5;   // the first Pareto front (4) was removed o.stagesA = Math.round(o.stagesA); o.stagesS = Math.round(o.stagesS);
  o.ihxH = raw.ihxH == null ? o.ihxH : !!raw.ihxH;
  if (raw.multi && typeof raw.multi === "object") for (const c of Object.keys(MULTI_DEF)) if (Array.isArray(raw.multi[c]))
    o.multi[c] = MULTI_DEF[c].map((d, i) => { try { return cleanComps(raw.multi[c][i]); } catch (_) { return d; } });
  if (raw.par && typeof raw.par === "object") for (const c of Object.keys(PAR_DEF)) {
    const src = raw.par[c]; if (!src || typeof src !== "object") continue;
    for (const k of Object.keys(PAR_DEF[c])) {
      const d = PAR_DEF[c][k], v = src[k];
      if (Array.isArray(d)) { if (Array.isArray(v)) o.par[c][k] = d.map((x, i) => typeof x === "boolean" ? (v[i] == null ? x : !!v[i]) : (isFinite(+v[i]) ? +v[i] : x)); }
      else if (typeof d === "boolean") { if (v != null) o.par[c][k] = !!v; }
      else if (v != null && isFinite(+v)) o.par[c][k] = +v;
    }
  }
  o.par.N.n = Math.min(4, Math.max(2, Math.round(o.par.N.n))); o.par.M.n_mvr = Math.min(20, Math.max(1, Math.round(o.par.M.n_mvr))); o.par.M.split = o.par.M.split.map(v => Math.min(1, Math.max(0, v)));
  const out = {}; for (const k of CASE_KEYS) out[k] = o[k];
  out.size = cleanSize(raw.size);
  return out;
}
function applyInputs(inp) {
  for (const k of CASE_KEYS) LV[k] = clone(inp[k]);
  if (inp.size) { state.sizeBy = inp.size.by; state[SIZE_KEY[inp.size.by][0]] = inp.size.value; save(); renderSize(); }
  LV.PHg = null; LV.PLg = null; LV.gKey = '';
  openPicker = -1;
  renderLiveControls(); changed();
}
const cycleName = i => i.cyc === "A" && i.stagesA === 1 ? "RHP-A, single-stage" : i.cyc === "H" ? `Conventional, ${i.stagesH === 2 ? "two" : "single"}-stage${i.ihxH ? " + IHX" : ""}` : i.cyc === "N" ? `${i.par.N.n}-level cascade` : CYC_SHORT[i.cyc];
function caseFluids(i) {
  const lab = cs => fluidLabel({ comps: cs, basis: i.basis });
  if (i.cyc === "K") return `${lab(i.multi.K[0])} / ${lab(i.multi.K[1])}`;
  if (i.cyc === "N") return i.multi.N.slice(0, i.par.N.n).map(lab).join(" / ");
  if (i.cyc === "M") return `${lab(i.multi.M[0])} + steam`;
  if (i.cyc === "X") return lab(i.multi.X[0]);
  return lab(i.mainComps);
}
const designText = i => designText0(i) + (i.expn === "exp" && i.cyc !== "X" ? ` · expander ηis ${fixed(i.etaX, 2)}` : "");
function designText0(i) {
  const p = i.par;
  switch (i.cyc) {
    case "H": return i.stagesH === 2 ? `stage-2 fraction ${fixed(i.split, 2)}` : "single stage";
    case "A": return i.stagesA === 1 ? `RHX1 ${i.lift} K · single stage` : `RHX1 ${i.lift} K · split ${fixed(i.split, 2)}`;
    case "B": return `RHX2 ${i.lift2} K · split ${fixed(i.split, 2)}`;
    case "C": return `RHX1 ${i.lift} K · RHX2 ${i.lift2} K · split ${fixed(i.split, 2)}`;
    case "S": return `RHX1 ${i.liftS} K · ${i.stagesS === 1 ? "single stage" : `share ${fixed(i.share, 2)}`} · α ${fixed(i.alpha, 3)}`;
    case "K": return `${p.K.P_b}/${p.K.P_t} bar · r ${p.K.r} · evap ${p.K.T_ev_b}/${p.K.T_ev_t} °C`;
    case "N": return `condensing ${p.N.T_cond.slice(0, p.N.n - 1).join("/")} °C`;
    case "X": return `arch ${p.X.arch} · ${p.X.P_L}/${p.X.P_H_min ? "lowest" : p.X.P_H} bar · split ${fixed(p.X.x, 2)}`;
    case "M": return `bottom ${p.M.T_cond} °C · ${p.M.n_mvr} MVR stages · split ${p.M.split_fill ? "auto" : p.M.split.slice(0, Math.max(0, p.M.n_mvr - 1)).map(v => fixed(v, 2)).join("/") || "–"}`;
  }
  return "";
}
const caseSummary = i => `${cycleName(i)} · ${caseFluids(i)} · ${i.T} °C air${+i.phi > 0 ? ` (+${fixed(i.phi, 2)} at ${i.Tmt} °C)` : ""} · ${i.src === "air" ? `air source ${i.srcT} °C, ${fixed(i.srcM, 1)} kg/s` : `water ${i.Ts} → ${i.Two} °C`}`;
const suggestName = i => `${CYC_SHORT[i.cyc]} · ${caseFluids(i)} · ${i.T} °C`;
function summarize(r, i) {
  if (!r || !r.success) return { ok: false, msg: r && r.message };
  const multi = MULTI(i.cyc), p = r.report || {};
  const o = multi ? {
    ok: true, COP: r.COP, VHC: r.VHC, Pmax: r.Pmax_bar, Tdis: r.T_dis_max, qair: r.q_air, wel: r.w_shaft / (i.motor || 0.93),
    Psuc: i.cyc === "X" ? p.P_L_bar : i.cyc === "K" ? Math.min(p.P_evap_bot_bar, p.P_evap_top_bar) : i.cyc === "N" ? (p.levels && p.levels[0] ? p.levels[0].Pe_bar : null) : p.P_evap_bar,
  } : {
    ok: true, COP: r.COP, VHC: p.VHC, Pmax: p.Pmax_bar, Psuc: p.P_evap_bar, Tdis: p.T_comp_discharge, qair: p.q_air_kJkg, wel: p.w_el_kJkg,
  };
  o.sec = 100 * (i.cyc === "K" ? (p.share_sec || 0) : (r.secShare || 0));
  // duties at the page's low-stage swept volume (the reference compressor's suction density and volumetric efficiency)
  if (multi) { const u = energyG(r, i); if (u) { o.Qkw = u.kw.air; o.Wkw = u.kw.el; } }
  else if (p.rho_suction_1) {
    const nv = state.nvMode === "fixed" ? state.nvFix : DRYER.etaVol(p.PR_1 || 1), m = p.rho_suction_1 * nv * state.vsw / 3600;
    o.Qkw = m * o.qair; o.Wkw = m * o.wel;
  }
  for (const k in o) if (typeof o[k] === "number" && !isFinite(o[k])) o[k] = null;
  o.warn = Array.isArray(r.warnings) ? r.warnings.length : 0; o.warnings = Array.isArray(r.warnings) ? r.warnings.slice(0, 5) : [];
  // duties follow the page's sizing: at its swept volume, or scaled to the target sink duty or sink flow (the swept
  // volume each case would need is then reported; a sink flow takes each case's own air heating per kg)
  if (o.Qkw > 0 && state.vsw > 0) {
    const q1 = o.Qkw / state.vsw, w1 = o.Wkw / state.vsw, dh = airDh(i);
    delete o.Qkw; delete o.Wkw;
    const tgt = () => state.mode !== "live" ? null : state.sizeBy === "duty" ? state.duty : state.sizeBy === "flow" && dh ? state.flow * dh : null;
    Object.defineProperty(o, "Qkw", { enumerable: true, get: () => tgt() ?? q1 * state.vsw });
    Object.defineProperty(o, "Wkw", { enumerable: true, get: () => (tgt() != null ? tgt() / q1 : state.vsw) * w1 });
    Object.defineProperty(o, "vsw", { enumerable: true, get: () => tgt() != null ? tgt() / q1 : state.vsw });
  }
  return o;
}

/* ---------------------------------------------------------------- storage */
async function casesBoot() {
  casesUI();
  const hash = (location.hash || "").match(/case=([A-Za-z0-9_-]+)/);
  if (hash) openCode(hash[1], "link");
  let db = null, user = null, downloads = null;
  try {
    if (window.claude && typeof window.claude.use === "function") [db, user, downloads] = await Promise.all([claude.use("db"), claude.use("user"), claude.use("downloads")]);
  } catch (_) { db = null; }
  CS.downloads = downloads;
  $("case-export").hidden = !downloads;
  if (db) {
    CS.db = db; CS.user = user; CS.store = "shared";
    try { CS.me = user ? await user.id() : null; } catch (_) { CS.me = null; }
    try { CS.canWrite = user ? await user.can("data.write") : null; } catch (_) { CS.canWrite = null; }
    db.collection("cases").orderBy("createdAt", "desc").limit(300).onSnapshot(snap => {
      CS.list = snap.docs.map(d => Object.assign({ id: d.id }, d.data()));
      casesChanged();
    }, err => caseStatus(`The shared case library stopped updating (${err.code}). Reload the page to reconnect.`, "warn"));
  } else {
    CS.store = "local";
    try { CS.list = JSON.parse(localStorage.getItem(LOCAL_KEY) || "[]"); } catch (_) { CS.list = []; }
    casesChanged();
  }
}
function saveLocal() { try { localStorage.setItem(LOCAL_KEY, JSON.stringify(CS.list)); } catch (_) {} }
async function storeAdd(doc) {
  if (CS.store === "shared") { const ref = await CS.db.collection("cases").add(doc); return ref.id; }
  const id = "l" + Date.now().toString(36) + Math.random().toString(36).slice(2, 6);
  CS.list.unshift(Object.assign({ id }, doc)); saveLocal(); casesChanged(); return id;
}
async function storeUpdate(id, patch) {
  if (CS.store === "shared") return CS.db.doc(`cases/${id}`).update(patch);
  const c = CS.list.find(x => x.id === id); if (c) Object.assign(c, patch); saveLocal(); casesChanged();
}
async function storeDelete(id) {
  if (CS.store === "shared") return CS.db.doc(`cases/${id}`).delete();
  CS.list = CS.list.filter(x => x.id !== id); saveLocal(); casesChanged();
}
function writeError(e) {
  const code = e && e.code;
  if (code === "invalid_argument" || code === "permission_denied") { CS.canWrite = false; casesChanged(); return "You can view the shared cases but not change them (ask the page owner for Contributor access). Links still work."; }
  if (code === "quota_exceeded") return "The case library is full: delete some cases first.";
  return `Could not save (${code || (e && e.message) || e}).`;
}

/* ---------------------------------------------------------------- actions */
function currentResults() {
  const c = LIVE.cfg;
  // compare the inputs as given (an air-source solve writes its source temperatures into LIVE.cfg)
  if (!c || !LIVE.res || inputsJSON(cleanInputs(LIVE.cfgIn || c)) !== inputsJSON(cleanInputs(LV))) return { ok: false, stale: true };
  if (LIVE.res.success) LIVE.res.secShare = secOn(c) ? secShare(c) : 0;
  return summarize(LIVE.res, c);
}
async function saveCase(asNew) {
  const inp = cleanInputs(caseInputs()), name = ($("case-name").value.trim() || suggestName(inp)).slice(0, 120);
  const note = $("case-note").value.trim().slice(0, 500), now = new Date().toISOString();
  let res = currentResults();
  if (res.stale) {
    caseStatus("Solving the case before saving…");
    await new Promise(r => setTimeout(r, 30));          // let the status paint
    try { await ensureCoolProp(); const r = solveCfg(inp); res = summarize(r, inp); } catch (e) { res = { ok: false, msg: cpError(e) }; }
  }
  try {
    if (!asNew && CS.loaded && CS.list.some(c => c.id === CS.loaded.id)) {
      await storeUpdate(CS.loaded.id, { name, note, inputs: inp, results: res, model: MODEL_VERSION, coolprop: CP_VERSION, updatedAt: now });
      CS.loaded = { id: CS.loaded.id, name, json: inputsJSON(inp) };
      caseStatus(`Updated “${name}”.`);
    } else {
      const doc = { v: CASE_SCHEMA, name, note, inputs: inp, results: res, model: MODEL_VERSION, coolprop: CP_VERSION, createdBy: CS.me, createdAt: now, updatedAt: now };
      const id = await storeAdd(doc);
      CS.loaded = { id, name, json: inputsJSON(inp) }; CS.fromLink = null;
      caseStatus(`Saved “${name}”${CS.store === "local" ? " in this browser" : ""}.`);
    }
    setCaseName(name, inp);
    CS.solved.delete(CS.loaded.id);
    casesChanged();
  } catch (e) { caseStatus(writeError(e), "warn"); }
}
function loadCase(id) {
  const c = CS.list.find(x => x.id === id); if (!c) return;
  let inp; try { inp = cleanInputs(c.inputs); } catch (e) { caseStatus(`Could not open “${c.name}”: ${e.message}`, "warn"); return; }
  applyInputs(inp);
  CS.loaded = { id, name: c.name, json: inputsJSON(inp) }; CS.fromLink = null;
  setCaseName(c.name || "", inp); $("case-note").value = c.note || "";
  caseStatus(`Loaded “${c.name}”.`);
  casesChanged();
}
async function deleteCase(id) {
  const c = CS.list.find(x => x.id === id); if (!c) return;
  const btn = document.querySelector(`[data-del="${id}"]`);
  if (btn && btn.dataset.armed !== "1") { btn.dataset.armed = "1"; btn.textContent = "Confirm delete"; setTimeout(() => { if (btn.isConnected) { btn.dataset.armed = ""; btn.textContent = "Delete"; } }, 4000); return; }
  try {
    await storeDelete(id);
    CS.compare = CS.compare.filter(x => x !== id); delete CS.colors[id];
    if (CS.loaded && CS.loaded.id === id) CS.loaded = null;
    caseStatus(`Deleted “${c.name}”.`); casesChanged();
  } catch (e) { caseStatus(writeError(e), "warn"); }
}
function toggleCompare(id, on) {
  if (on && !CS.compare.includes(id)) {
    if (CS.compare.length >= MAX_COMPARE) { caseStatus(`Compare up to ${MAX_COMPARE} cases: untick one first.`, "warn"); casesChanged(); return; }
    CS.compare.push(id);
    const used = new Set(Object.values(CS.colors)); CS.colors[id] = [0, 1, 2, 3].find(k => !used.has(k));
  } else if (!on) { CS.compare = CS.compare.filter(x => x !== id); delete CS.colors[id]; }
  casesChanged();
}
const caseColor = id => CASE_COLORS[CS.colors[id] != null ? CS.colors[id] : 0];
function moveCompare(id, to) {
  const from = CS.compare.indexOf(id); if (from < 0) return;
  to = Math.max(0, Math.min(CS.compare.length - 1, to)); if (to === from) return;
  CS.compare.splice(from, 1); CS.compare.splice(to, 0, id); setOrderBy(""); casesChanged();
}

/* "Order by": keep the comparison sorted by name or by a result, best first (higher first where neither is better);
   cases not solved yet, or that failed, go last in their current order. A hand re-order (drag or ‹ ›) clears it. */
const ORDER_KEY = "dryer-compare-order";
CS.orderBy = "";
function orderOptions() {
  return [["", "Your order (drag)"], ["name", "Name"]].concat(RESULT_ROWS.map(([k, lab, , better]) =>
    [k, `${lab}${better < 0 ? " (lowest first)" : better > 0 ? " (highest first)" : ""}`]));
}
function setOrderBy(k) {
  CS.orderBy = orderOptions().some(o => o[0] === k) ? k : "";
  const o = orderOptions().find(x => x[0] === CS.orderBy), val = $("cmp-order-val");
  if (val) val.textContent = o[1];
  for (const b of document.querySelectorAll("#cmp-order-list .fopt")) b.setAttribute("aria-selected", String(b.dataset.value === CS.orderBy));
  try { localStorage.setItem(ORDER_KEY, CS.orderBy); } catch (_) {}
}
function sortedCompare(cs, key, valueOf) {
  if (!key) return cs.map(c => c.id);
  if (key === "name") return cs.slice().sort((a, b) => String(a.name).localeCompare(String(b.name), undefined, { numeric: true })).map(c => c.id);
  const row = RESULT_ROWS.find(r => r[0] === key), sgn = row && row[3] < 0 ? 1 : -1;
  const v = cs.map((c, i) => { const x = valueOf(c); return { id: c.id, i, x: x != null && isFinite(x) ? x : null }; });
  return v.sort((a, b) => (a.x == null) - (b.x == null) || (a.x == null ? a.i - b.i : sgn * (a.x - b.x) || a.i - b.i)).map(o => o.id);
}
function applyOrder(cs) {
  if (!CS.orderBy) return false;
  const order = sortedCompare(cs, CS.orderBy, c => { const s = solvedResults(c); return s && s.ok ? s[CS.orderBy] : null; });
  const next = order.concat(CS.compare.filter(id => !order.includes(id)));
  if (next.join() === CS.compare.join()) return false;
  CS.compare = next; return true;
}
function startDrag(ev, col) {
  if (ev.button > 0) return;
  ev.preventDefault();
  const host = $("cmp-side");
  col.classList.add("dragging");
  const move = e => {
    const others = [...host.children].filter(x => x !== col);
    const next = others.find(x => { const r = x.getBoundingClientRect(); return e.clientX < r.left + r.width / 2; });
    if (next !== col.nextSibling) host.insertBefore(col, next || null);
    const wrap = host.parentElement, r = wrap.getBoundingClientRect();
    if (e.clientX > r.right - 40) wrap.scrollLeft += 12; else if (e.clientX < r.left + 40) wrap.scrollLeft -= 12;
  };
  const end = () => {
    window.removeEventListener("pointermove", move); window.removeEventListener("pointerup", end); window.removeEventListener("pointercancel", end);
    col.classList.remove("dragging");
    const order = [...host.children].map(x => x.dataset.id).filter(Boolean);
    if (order.join() !== CS.compare.join()) { CS.compare = order.concat(CS.compare.filter(x => !order.includes(x))); setOrderBy(""); casesChanged(); }
  };
  window.addEventListener("pointermove", move); window.addEventListener("pointerup", end); window.addEventListener("pointercancel", end);
}

/* links: the inputs, base64url JSON after #case= */
const b64url = s => btoa(unescape(encodeURIComponent(s))).replace(/\+/g, "-").replace(/\//g, "_").replace(/=+$/, "");
const unb64url = s => decodeURIComponent(escape(atob(s.replace(/-/g, "+").replace(/_/g, "/") + "===".slice((s.length + 3) % 4))));
const caseCode = (name, inp) => b64url(JSON.stringify({ v: CASE_SCHEMA, n: name, i: inp }));
const caseLink = (name, inp) => `${PAGE_URL}#case=${caseCode(name, inp)}`;
function openCode(text, how) {
  const m = String(text).trim().match(/(?:case=)?([A-Za-z0-9_-]{20,})\s*$/);
  if (!m) { caseStatus("That doesn't look like a case link.", "warn"); return false; }
  let o, inp;
  try { o = JSON.parse(unb64url(m[1])); inp = cleanInputs(o.i); }
  catch (e) { caseStatus(`Could not read that case link${e && e.message ? `: ${e.message}` : ""}.`, "warn"); return false; }
  const name = typeof o.n === "string" ? o.n.slice(0, 120) : suggestName(inp);
  const go = () => {
    if (state.mode !== 'live' && typeof setMode === 'function') setMode('live');
    applyInputs(inp);
    CS.loaded = null; CS.fromLink = { name, json: inputsJSON(inp) };
    setCaseName(name, inp); $("case-note").value = "";
    caseStatus(`Opened “${name}” from a ${how}. Save it to add it to the cases.`);
    casesChanged();
  };
  if (D) go(); else { const t = setInterval(() => { if (D) { clearInterval(t); go(); } }, 200); }
  return true;
}
async function copyText(text, label) {
  const box = $("case-linkbox"), inp = $("case-linktext");
  inp.value = text; box.hidden = false; inp.select();
  try { await navigator.clipboard.writeText(text); caseStatus(`${label} copied.`); }
  catch (_) { caseStatus(`${label} is selected below: press Ctrl+C (⌘C) to copy it.`); }
}

/* ---------------------------------------------------------------- export */
const csvCell = v => { const s = v == null ? "" : String(v); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
async function exportCompare() {
  const cs = compareCases(); if (!cs.length || !CS.downloads) return;
  const rows = [["", ...cs.map(c => c.name)]];
  for (const [lab, f] of INPUT_ROWS) rows.push([lab, ...cs.map(c => f(c.inputs))]);
  for (const [k, lab, d] of RESULT_ROWS) rows.push([lab, ...cs.map(c => { const r = solvedResults(c); return r && r.ok && r[k] != null ? (+r[k]).toFixed(d + 2) : ""; })]);
  rows.push(["Case link", ...cs.map(c => caseLink(c.name, c.inputs))]);
  const csv = rows.map(r => r.map(csvCell).join(",")).join("\n");
  const json = JSON.stringify({ schema: CASE_SCHEMA, model: MODEL_VERSION, coolprop: CP_VERSION, exported: new Date().toISOString(),
    cases: cs.map(c => ({ name: c.name, note: c.note || "", inputs: c.inputs, results: solvedResults(c) })) }, null, 2);
  try {
    await CS.downloads.save({ filename: "dryer_case_comparison.csv", data: csv });
    await CS.downloads.save({ filename: "dryer_cases.json", data: json });
  } catch (e) { if (e && e.code !== "declined") caseStatus(`Export failed (${e.code || e}).`, "warn"); }
}

/* ---------------------------------------------------------------- page */
function setCaseName(name, inp) { CS.nameAuto = !name || name === suggestName(inp); $("case-name").value = CS.nameAuto ? "" : name; }
function caseStatus(text, cls) { const s = $("case-status"); s.textContent = text || ""; s.className = "hint" + (cls ? " " + cls : ""); }
function casesUI() {
  const ord = $("cmp-order"), list = $("cmp-order-list");
  const openList = on => { list.hidden = !on; ord.setAttribute("aria-expanded", String(on)); };
  for (const [k, lab] of orderOptions()) {
    const b = document.createElement("button"); b.type = "button"; b.className = "fopt"; b.dataset.value = k; b.textContent = lab;
    b.setAttribute("role", "option");
    b.addEventListener("click", () => { setOrderBy(k); openList(false); ord.focus(); casesChanged(); });
    list.append(b);
  }
  let saved = ""; try { saved = localStorage.getItem(ORDER_KEY) || ""; } catch (_) {}
  setOrderBy(saved);
  ord.addEventListener("click", () => openList(list.hidden));
  document.addEventListener("keydown", e => { if (e.key === "Escape" && !list.hidden) { openList(false); ord.focus(); } });
  document.addEventListener("pointerdown", e => { if (!list.hidden && !e.target.closest(".cmp-order")) openList(false); });
  $("case-name").addEventListener("input", () => { CS.nameAuto = !$("case-name").value.trim(); });
  $("case-save").addEventListener("click", () => saveCase(false));
  $("case-saveas").addEventListener("click", () => saveCase(true));
  $("case-copy").addEventListener("click", () => { const inp = cleanInputs(caseInputs()); copyText(caseLink($("case-name").value.trim() || suggestName(inp), inp), "Link to the current inputs"); });
  $("case-open").addEventListener("click", () => { const b = $("case-openbox"); b.hidden = !b.hidden; if (!b.hidden) $("case-opentext").focus(); });
  $("case-opengo").addEventListener("click", () => { if (openCode($("case-opentext").value, "link")) { $("case-opentext").value = ""; $("case-openbox").hidden = true; } });
  $("case-export").addEventListener("click", exportCompare);
  $("case-none").addEventListener("click", () => { CS.compare = []; CS.colors = {}; casesChanged(); });
  const group = (id, opts, cur, pick) => {
    const box = $(id); box.innerHTML = "";
    for (const [k, s] of Object.entries(opts)) {
      const b = document.createElement("button"); b.type = "button"; b.className = "tbtn"; b.dataset.v = k;
      b.textContent = s.label.replace(/, [^,]*$/, ""); b.setAttribute("aria-pressed", String(k === cur));
      b.addEventListener("click", () => { for (const o of box.children) o.setAttribute("aria-pressed", String(o === b)); pick(k); });
      box.append(b);
    }
  };
  group("sweep-x", SWEEPS, CS.sweep.x, k => { CS.sweep.x = k; CS.sweep.data = null; CS.sweep.key = ""; drawSweep(); });
  group("sweep-y", SWEEP_Y, CS.sweep.y, k => { CS.sweep.y = k; drawSweep(); });
  $("sweep-run").addEventListener("click", () => { if (CS.sweep.running) CS.sweep.key = ""; else runSweep(); });
}
let casesDirty = false;
function casesChanged() { if (casesDirty) return; casesDirty = true; requestAnimationFrame(() => { casesDirty = false; casesRender(); }); }
function casesHeader() {
  if (!$("case-current")) return;
  let j = ""; try { j = inputsJSON(cleanInputs(caseInputs())); } catch (_) {}
  const cur = $("case-current");
  if (CS.loaded) cur.textContent = j === CS.loaded.json ? `Showing “${CS.loaded.name}”` : `Showing “${CS.loaded.name}”, modified`;
  else if (CS.fromLink) cur.textContent = j === CS.fromLink.json ? `Showing “${CS.fromLink.name}” from a link (not saved)` : `Opened “${CS.fromLink.name}” from a link, modified`;
  else cur.textContent = "The current design is not saved as a case";
  const own = CS.loaded && CS.list.find(c => c.id === CS.loaded.id);
  $("case-save").textContent = own ? "Update case" : "Save case";
  $("case-saveas").hidden = !own;
  const ro = CS.store === "shared" && CS.canWrite === false;
  $("case-save").disabled = $("case-saveas").disabled = ro;
  try { $("case-name").placeholder = suggestName(cleanInputs(caseInputs())); } catch (_) {}
  if (CS.nameAuto && document.activeElement !== $("case-name")) $("case-name").value = "";
  $("case-store").textContent = CS.store === "shared" ? (ro ? "Shared cases, view only" : "Shared with everyone who can open this page")
    : CS.store === "local" ? "Saved in this browser only" : "Loading cases…";
}

let renderSeq = 0;
async function casesRender() {
  const seq = ++renderSeq;
  casesHeader();
  const ids = [...new Set(CS.list.map(c => c.createdBy).filter(Boolean))];
  if (CS.user && ids.length) { try { const ps = await CS.user.profiles(ids); for (const id of ids) CS.names[id] = ps[id] ? (ps[id].isMe ? "you" : ps[id].name) : ""; } catch (_) {} }
  let editor = false;
  if (CS.user) { try { editor = await CS.user.canEdit(); } catch (_) {} }
  if (seq !== renderSeq) return;
  CS.compare = CS.compare.filter(id => CS.list.some(c => c.id === id)).slice(0, MAX_COMPARE);
  const box = $("case-list"); box.innerHTML = "";
  if (!CS.list.length) { const p = document.createElement("p"); p.className = "hint"; p.textContent = CS.store ? "No cases yet. Set the design, then Save case." : "Loading cases…"; box.append(p); }
  for (const c of CS.list) {
    let inp; try { inp = cleanInputs(c.inputs); } catch (_) { continue; }
    const ci = CS.compare.indexOf(c.id), r = c.results || {};
    const card = document.createElement("div"); card.className = "case-card" + (ci >= 0 ? " on" : "") + (CS.loaded && CS.loaded.id === c.id ? " loaded" : "");
    if (ci >= 0) card.style.setProperty("--case", `var(${caseColor(c.id)})`);
    const top = document.createElement("label"); top.className = "case-top";
    const cb = document.createElement("input"); cb.type = "checkbox"; cb.checked = ci >= 0; cb.setAttribute("aria-label", `Compare ${c.name}`);
    cb.addEventListener("change", () => toggleCompare(c.id, cb.checked));
    const nm = document.createElement("span"); nm.className = "case-name"; nm.textContent = c.name || "Untitled case";
    const st = caseLight(c), cop = document.createElement("span"); cop.className = "case-cop";
    const dot = document.createElement("i"); dot.className = `tl tl-${st.light}`; dot.setAttribute("aria-hidden", "true");
    cop.append(dot, st.text); cop.title = st.tip; cop.setAttribute("aria-label", st.text ? `${st.label}: ${st.text}` : st.label);
    top.append(cb, nm, cop);
    const sub = document.createElement("div"); sub.className = "case-sub"; sub.textContent = caseSummary(inp);
    const meta = document.createElement("div"); meta.className = "case-meta";
    const who = c.createdBy ? (CS.names[c.createdBy] || "someone") : (CS.store === "local" ? "this browser" : "someone");
    meta.textContent = `${who} · ${c.updatedAt ? new Date(c.updatedAt).toLocaleDateString(undefined, { day: "numeric", month: "short", year: "numeric" }) : ""}${c.model && c.model !== MODEL_VERSION ? " · older model" : ""}`;
    card.append(top, sub);
    if (c.note) { const n = document.createElement("div"); n.className = "case-note"; n.textContent = c.note; card.append(n); }
    card.append(meta);
    const act = document.createElement("div"); act.className = "case-act";
    const b = (t, fn, attrs = {}) => { const x = document.createElement("button"); x.type = "button"; x.className = "btn sm"; x.textContent = t; for (const k in attrs) x.setAttribute(k, attrs[k]); x.addEventListener("click", fn); act.append(x); return x; };
    b("Load", () => loadCase(c.id));
    if (!solvedResults(c)) b("Check", () => checkCase(c.id), { title: "Solve this case now with the current model and update its status" });
    b("Copy link", () => copyText(caseLink(c.name, inp), `Link to “${c.name}”`));
    const mayDelete = CS.store === "local" || (c.createdBy && c.createdBy === CS.me) || editor;
    if (mayDelete && !(CS.store === "shared" && CS.canWrite === false)) b("Delete", () => deleteCase(c.id), { "data-del": c.id });
    card.append(act);
    box.append(card);
  }
  $("case-none").hidden = !CS.compare.length;
  $("case-count").textContent = CS.compare.length ? `${CS.compare.length} of ${MAX_COMPARE} ticked to compare` : `Tick up to ${MAX_COMPARE} cases to compare them`;
  renderCompare();
}

/* ---------------------------------------------------------------- comparison */
const compareCases = () => CS.compare.map(id => CS.list.find(c => c.id === id)).filter(Boolean)
  .map(c => { try { return Object.assign({}, c, { inputs: cleanInputs(c.inputs) }); } catch (_) { return null; } }).filter(Boolean);
/* traffic light for a case: green solved, amber solved with warnings (or not checked against the current model),
   red no design. A case solved this session wins over the results stored when it was saved. */
function caseLight(c) {
  const fresh = solvedResults(c), r = fresh || c.results || null;
  if (!r || r.stale) return { light: "a", label: "Not checked", text: "", tip: "Not checked: no results stored for this case. Check solves it now" };
  if (!r.ok) return { light: "r", label: "No design", text: "no design", tip: r.msg ? String(r.msg).slice(0, 300) : "The design did not solve" };
  const older = !fresh && c.model && c.model !== MODEL_VERSION;
  if (r.warn > 0) return { light: "a", label: "Solved with warnings", text: `COP ${fixed(r.COP, 2)} · ${r.warn} warning${r.warn > 1 ? "s" : ""}`, tip: (r.warnings || []).join("\n") || "Solved with warnings" };
  if (older) return { light: "a", label: "Solved with an older model", text: `COP ${fixed(r.COP, 2)} · older model`, tip: "Saved with an earlier version of the cycle models: Check re-solves it" };
  return { light: "g", label: "Solved", text: `COP ${fixed(r.COP, 2)}`, tip: fresh ? "Solved with the current model" : "Solved when saved" };
}
/* solve one case now; keep its status for this session and, where the case may be changed here, store the results */
function checkCase(id) {
  const c = CS.list.find(x => x.id === id); if (!c) return;
  caseStatus(`Solving “${c.name}”…`);
  setTimeout(async () => {
    try {
      await ensureCoolProp();
      const inp = cleanInputs(c.inputs), res = solveCfg(inp), sum = summarize(res, inp);
      CS.solved.set(c.id, { json: inputsJSON(c.inputs), res, sum });
      const mine = CS.store === "local" || (c.createdBy && c.createdBy === CS.me);
      if (mine && CS.canWrite !== false) { try { await storeUpdate(c.id, { results: sum, model: MODEL_VERSION, coolprop: CP_VERSION }); } catch (_) {} }
      caseStatus(sum.ok ? `“${c.name}”: COP ${fixed(sum.COP, 3)}${sum.warn ? `, ${sum.warn} warning${sum.warn > 1 ? "s" : ""}` : ""}.` : `“${c.name}”: no design (${sum.msg || "did not solve"}).`, sum.ok ? undefined : "warn");
    } catch (e) { caseStatus(`Could not solve “${c.name}”: ${cpError(e)}`, "warn"); }
    casesChanged();
  }, 30);
}
function solvedResults(c) { const s = CS.solved.get(c.id); return s && s.json === inputsJSON(c.inputs) ? s.sum : null; }
let compareJob = 0;
function renderCompare() {
  const host = $("compare"), cs = compareCases();
  host.hidden = cs.length < 2;
  $("case-export").disabled = cs.length < 1;
  if (cs.length < 2) return;
  if (applyOrder(cs)) { casesChanged(); return; }
  const job = ++compareJob, todo = cs.filter(c => !solvedResults(c));
  compareTables(cs); sideBySide(cs);
  $("cmp-status").textContent = todo.length ? `Solving ${todo.length} case${todo.length > 1 ? "s" : ""}…` : "";
  const step = () => {
    if (job !== compareJob) return;
    if (!CP || !LIVE.airBase && !LIVE.air) { ensureCoolProp().then(() => setTimeout(step, 50), () => {}); return; }
    const c = todo.shift();
    if (c) {
      const res = solveCfg(c.inputs);
      CS.solved.set(c.id, { json: inputsJSON(c.inputs), res, sum: summarize(res, c.inputs) });
      $("cmp-status").textContent = todo.length ? `Solving ${todo.length} more…` : "";
      compareTables(cs); sideBySide(cs);
      setTimeout(step, 30); return;
    }
    if (applyOrder(cs)) { casesChanged(); return; }
    compareTables(cs); sideBySide(cs); overlayTs(cs); drawSweep();
  };
  setTimeout(step, 30);
}
function compareTables(cs) {
  const swatch = c => `<i class="sw" style="background:var(${caseColor(c.id)})"></i>`;
  const msgs = cs.map(c => { const s = CS.solved.get(c.id); return s && s.json === inputsJSON(c.inputs) && !s.sum.ok ? `No design for ${c.name}: ${s.sum.msg}` : null; }).filter(Boolean);
  const drift = cs.filter(c => { const r = solvedResults(c), s = c.results; return r && r.ok && s && s.ok && Math.abs(r.COP - s.COP) > 1e-3 * Math.abs(s.COP); }).map(c => c.name);
  const notes = [...msgs, ...(drift.length ? [`Re-solved results differ from when these were saved (model or CoolProp changed): ${drift.join(", ")}`] : [])];
  $("cmp-notes").innerHTML = notes.map(n => `<span class="chip warn">${escapeHTML(n)}</span>`).join("");
  // results, one column per case; the left-most is the baseline
  const sums = cs.map(solvedResults), base = sums[0];
  const head = `<tr><th>Result</th>${cs.map((c, i) => `<th class="ch"><span class="cmphd">${swatch(c)}<span>${escapeHTML(c.name)}</span></span>
    ${i === 0 ? '<span class="cmpmv"><span class="base-tag">baseline</span></span>' : ""}</th>`).join("")}</tr>`;
  let body = "";
  for (const [k, lab, d, better] of RESULT_ROWS) {
    const v = sums.map(s => s && s.ok ? s[k] : null), fin = v.filter(x => x != null && isFinite(x));
    const best = better && fin.length > 1 && Math.max(...fin) - Math.min(...fin) > 1e-9 ? (better > 0 ? Math.max(...fin) : Math.min(...fin)) : null;
    body += `<tr><td>${lab}</td>${v.map((x, i) => {
      if (!sums[i]) return `<td class="lbl">solving…</td>`;
      if (x == null) return "<td>–</td>";
      const b = base && base.ok ? base[k] : null, pc = i > 0 && b != null && Math.abs(b) > 1e-12 ? 100 * (x - b) / b : null;
      return `<td class="${best != null && Math.abs(x - best) < 1e-9 ? "best" : ""}"><span>${fixed(x, d)}</span>${pc != null ? ` <small>(${Math.abs(pc) < 0.05 ? "=" : `${pc > 0 ? "+" : "−"}${Math.abs(pc).toFixed(1)} %`})</small>` : ""}</td>`;
    }).join("")}</tr>`;
  }
  $("cmp-results").innerHTML = `<table class="cmp">${head}${body}</table>`;
  // inputs that differ
  let rows = "", same = [];
  for (const [lab, f] of INPUT_ROWS) {
    const vals = cs.map(c => String(f(c.inputs)));
    if (vals.every(v => v === vals[0])) { same.push(lab.replace(/, .*$/, "")); continue; }
    rows += `<tr><td>${lab}</td>${vals.map(v => `<td>${escapeHTML(v)}</td>`).join("")}</tr>`;
  }
  $("cmp-inputs").innerHTML = rows ? `<table class="cmp"><tr><th>Input</th>${cs.map(c => `<th class="ch">${swatch(c)}${escapeHTML(c.name)}</th>`).join("")}</tr>${rows}</table>` : `<p class="hint">These cases have identical inputs.</p>`;
  $("cmp-same").textContent = same.length && rows ? `Same in every case: ${same.join(", ")}.` : "";
}
const escapeHTML = s => String(s).replace(/[&<>"']/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]));
const cssVar = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim() || "#888";
/* side-by-side dashboard: one column per case, the headline dials on scales shared by all the cases, key figures
   against the baseline (the left-most case), the heat and power split, and the case's own T-s cycle */
function sideBySide(cs) {
  const host = $("cmp-side"); host.innerHTML = "";
  const st = cs.map(c => { const s = CS.solved.get(c.id); return s && s.json === inputsJSON(c.inputs) ? s : null; });
  const sums = st.map(s => s && s.sum.ok ? s.sum : null), ok = sums.filter(Boolean);
  const top = (vals, step) => niceTop(Math.max(...vals.filter(v => v != null && isFinite(v)), step), step);
  const kwTop = ok.length ? top(ok.flatMap(s => [s.Qkw, s.Wkw]), 10) : 10, pTop = ok.length ? top(ok.map(s => s.Pmax), 10) : 10;
  const tTop = ok.length ? top(ok.map(s => s.Tdis), 50) : 50, powTop = ok.length ? Math.max(...ok.map(s => (s.Qkw || 0) + (s.Wkw || 0))) || 1 : 1;
  const mark = (k, better) => {
    const v = sums.map(s => s ? s[k] : null), fin = v.filter(x => x != null && isFinite(x));
    if (!better || fin.length < 2) return () => "";
    const hi = Math.max(...fin), lo = Math.min(...fin), tol = 1e-9 * Math.max(1, Math.abs(hi));
    if (hi - lo <= tol) return () => "";
    const best = better > 0 ? hi : lo;
    return i => (v[i] != null && Math.abs(v[i] - best) <= tol ? "win" : "");
  };
  const base = sums[0];
  const delta = (v, b, i) => {
    if (i === 0 || b == null || v == null || !isFinite(v) || !isFinite(b) || Math.abs(b) < 1e-12) return "";
    const pc = 100 * (v - b) / b;
    return Math.abs(pc) < 0.05 ? " (=)" : ` (${pc > 0 ? "+" : "−"}${Math.abs(pc).toFixed(1)} %)`;
  };
  const KV = [["VHC", "VHC, MJ/m³", 2, 1], ["Psuc", "Lowest suction, bar", 2, 0], ["Pmax", "Max pressure, bar", 1, -1], ["Tdis", "Hottest discharge, °C", 1, -1],
    ["qair", "Heat to air, kJ/kg ref.", 1, 0], ["wel", "Electricity, kJ/kg ref.", 1, -1], ["sec", "Secondary air, % of heat", 1, 0]].map(([k, lab, d, b]) => ({ k, lab, d, cls: mark(k, b) }));
  const M_ = { COP: mark("COP", 1), Q: mark("Qkw", 1), W: mark("Wkw", -1), P: mark("Pmax", -1), T: mark("Tdis", -1) };
  cs.forEach((c, i) => {
    const col = document.createElement("div");
    col.className = "side-col"; col.style.setProperty("--case", `var(${caseColor(c.id)})`); col.dataset.id = c.id;
    const head = document.createElement("div"); head.className = "side-head";
    const grip = document.createElement("span"); grip.className = "grip"; grip.textContent = "⠿"; grip.title = "Drag to re-order"; grip.setAttribute("aria-hidden", "true");
    const h = document.createElement("h4"); h.textContent = c.name;
    const mv = (d, lab, sym) => { const b = document.createElement("button"); b.type = "button"; b.className = "icon mv"; b.textContent = sym; b.setAttribute("aria-label", `Move ${c.name} ${lab}`); b.disabled = (d < 0 && i === 0) || (d > 0 && i === cs.length - 1); b.addEventListener("click", () => moveCompare(c.id, i + d)); return b; };
    head.append(grip, h, mv(-1, "left", "‹"), mv(1, "right", "›"));
    if (i === 0) { const bt = document.createElement("span"); bt.className = "base-tag"; bt.textContent = "baseline"; head.append(bt); }
    grip.addEventListener("pointerdown", ev => startDrag(ev, col));
    const sub = document.createElement("div"); sub.className = "sub"; sub.textContent = caseSummary(c.inputs);
    col.append(head, sub);
    const s = st[i], r = sums[i];
    if (!s || !r) {
      const e = document.createElement("div"); e.className = "empty"; e.textContent = !s ? "Solving…" : `No design: ${s.sum.msg || ""}`;
      col.append(e); host.append(col); return;
    }
    const n = Math.max(0, Math.min(5, r.COP));
    const cop = document.createElement("div"); cop.className = "side-cop " + M_.COP(i);
    cop.innerHTML = `<span class="v">${fixed(r.COP, 2)}</span><span class="stars" role="img" aria-label="${n.toFixed(2)} of 5 stars">${starRowL(n, `sst-${c.id}`)}</span><small>heating COP${delta(r.COP, base ? base.COP : null, i).replace(/^ \((.*)\)$/, " · $1 vs baseline")}</small>`;
    const dials = document.createElement("div"); dials.className = "side-dials";
    dials.innerHTML = `<div class="dial ${M_.Q(i)}">${dialL(r.Qkw, kwTop, "--s4", "kW heat to air")}<span class="dk">Heat to air, kW</span></div>
      <div class="dial ${M_.W(i)}">${dialL(r.Wkw, kwTop, "--s2", "kW electrical")}<span class="dk">Electrical, kW</span></div>
      <div class="dial ${M_.P(i)}">${dialL(r.Pmax, pTop, "--s6", "bar maximum pressure")}<span class="dk">Max pressure, bar</span></div>
      <div class="dial ${M_.T(i)}">${dialL(r.Tdis, tTop, "--s5", "°C discharge", 0)}<span class="dk">Discharge, °C</span></div>`;
    const bar = document.createElement("div"); bar.className = "side-bar"; bar.setAttribute("role", "img");
    const qs = (r.Qkw || 0) * (r.sec || 0) / 100, qp = (r.Qkw || 0) - qs;
    bar.setAttribute("aria-label", `heat to primary air ${fixed(qp, 1)} kW, secondary air ${fixed(qs, 1)} kW, electrical ${fixed(r.Wkw, 1)} kW`);
    for (const [v, cv] of [[qp, "--s4"], [qs, "--s3"], [r.Wkw || 0, "--s2"]]) { if (!(v > 0)) continue; const x = document.createElement("i"); x.style.width = `${(100 * v / powTop).toFixed(2)}%`; x.style.background = `var(${cv})`; bar.append(x); }
    const bl = document.createElement("div"); bl.className = "side-barlab";
    bl.textContent = `${fixed(r.Qkw, 1)} kW to air${r.sec > 0 ? ` (${fixed(r.sec, 0)} % secondary)` : ""} · ${fixed(r.Wkw, 1)} kW electrical at ${fixed(r.vsw ?? state.vsw, 0)} m³/h`;
    const kv = document.createElement("dl"); kv.className = "side-kv";
    for (const m of KV) {
      const dt = document.createElement("dt"), dd = document.createElement("dd"), cl = m.cls(i);
      dt.textContent = m.lab; dd.textContent = fixed(r[m.k], m.d);
      const dl = delta(r[m.k], base ? base[m.k] : null, i);
      if (dl) { const sm = document.createElement("small"); sm.className = "dlt"; sm.textContent = dl; dd.append(sm); }
      if (cl) dd.className = cl;
      kv.append(dt, dd);
    }
    const mini = document.createElement("div"); mini.className = "mini";
    mini.innerHTML = `<svg id="side-ts-${c.id}"></svg>`;
    col.append(cop, dials, bar, bl, kv, mini);
    host.append(col);
  });
  // draw the small cycles once every column is in place, so each is sized to its own column
  cs.forEach((c, i) => { const s = st[i]; if (!s || !sums[i]) return; const layers = tsLayers(c, s.res); if (layers.length) xy(`side-ts-${c.id}`, layers, { xlabel: "s, kJ/(kg·K)", ylabel: "T, °C" }); });
}
function tsLayers(c, res) {
  if (!res || !res.success) return [];
  const col = cssVar(caseColor(c.id)), paths = res.paths || (res.loops || []).flatMap(L => L.paths), out = [];
  for (const p of paths) {
    const pts = p.pts.filter(q => isFinite(q[3]) && isFinite(q[0]));
    if (pts.length < 2) continue;
    out.push({ name: c.name, comp: p.name, color: col, pts: pts.map(q => [q[3], q[0]]), raw: pts, tip: true, dash: /Compressor|Stage/.test(p.name) ? "6 3" : "" });
  }
  return out;
}
// every compared case's cycle on one T-s diagram, in its colour (both loops of the multi-loop cycles)
function overlayTs(cs) {
  const layers = [];
  for (const c of cs) { const s = CS.solved.get(c.id); if (s) layers.push(...tsLayers(c, s.res)); }
  if (!layers.length) { $("cmp-ts").innerHTML = ""; return; }
  xy("cmp-ts", layers, { xlabel: "Specific entropy, kJ/(kg·K)", ylabel: "Temperature, °C" });
}

/* sweep: one operating variable across its range for every compared case, a point per frame (started on request:
   each point is a full design solve) */
function runSweep() {
  const cs = compareCases(); if (cs.length < 2 || !CP) return;
  const sx = SWEEPS[CS.sweep.x], key = CS.sweep.x + "|" + cs.map(c => c.id + inputsJSON(c.inputs)).sort().join("|");
  CS.sweep.key = key; CS.sweep.running = true; $("sweep-run").textContent = "Stop";
  CS.sweep.data = cs.map(c => ({ id: c.id, name: c.name, pts: [] }));
  const xs = []; for (let x = sx.lo; x <= sx.hi + 1e-9; x += sx.step) xs.push(+x.toFixed(6));
  const jobs = []; cs.forEach((c, i) => xs.forEach(x => jobs.push([i, c, x])));
  let n = 0; const total = jobs.length;
  const done = () => { CS.sweep.running = false; $("sweep-run").textContent = "Run the sweep"; };
  const step = () => {
    if (CS.sweep.key !== key) { done(); $("sweep-status").textContent = "Sweep stopped."; return; }
    const [i, c, x] = jobs.shift(); n++;
    const inp = Object.assign(clone(c.inputs), { [CS.sweep.x]: x });
    const res = solveCfg(inp);
    CS.sweep.data[i].pts.push([x, summarize(res, inp)]);
    $("sweep-status").textContent = jobs.length ? `Solving ${n} of ${total} points…` : `${CS.sweep.data.reduce((a, d) => a + d.pts.filter(p => p[1].ok).length, 0)} of ${total} points solved`;
    drawSweep();
    if (jobs.length) setTimeout(step, 20); else done();
  };
  setTimeout(step, 20);
}
function drawSweep() {
  const d = CS.sweep.data, svg = $("cmp-sweep");
  if (!d) { svg.innerHTML = ""; $("sweep-status").textContent = "Each point is a full design solve, so the sweep runs on request."; return; }
  const sx = SWEEPS[CS.sweep.x], sy = SWEEP_Y[CS.sweep.y];
  const series = d.filter(s => CS.compare.includes(s.id)).map(s => ({ name: s.name, color: cssVar(caseColor(s.id)), pts: s.pts.map(([x, r]) => ({ x, y: r.ok ? sy.f(r) : null })) }));
  const xt = []; for (let x = sx.lo; x <= sx.hi + 1e-9; x += sx.step * (sx.hi - sx.lo > 8 * sx.step ? 2 : 1)) xt.push(+x.toFixed(6));
  chart("cmp-sweep", series, { xs: [sx.lo, sx.hi], xt, xlabel: sx.label, dec: sy.d, ttx: x => `${sx.label.replace(/, .*$/, "")} ${x}`, h: 300 });
}
