/* ================================================================ cases
   Save the page's inputs as named cases, keep them in a library shared by everyone who opens the
   page (the artifact's db capability; this browser only when that is unavailable), send one as a
   link, and compare up to four: a results table against a baseline, the inputs that differ,
   overlaid T-s and P-h cycles, and a sweep of one operating variable. Everything is re-solved here
   with CoolProp; a stored case keeps a summary of its results for the list and for drift checks. */

const CASE_SCHEMA = 1;
const MODEL_VERSION = "sghp-js 2026-10-11";                 // bump when the cycle model changes results
const CASE_COLORS = ["--s1", "--s2", "--s3", "--s4"];
const MAX_COMPARE = 4;
const PAGE_URL = typeof location !== "undefined" && /(\.github\.io|^localhost|^127\.0\.0\.1)$/.test(location.hostname)
  ? location.origin + location.pathname : "https://waikato-ahuora-smart-energy-systems.github.io/HeatPumpLab/steam/";   // where case links point: this page on GitHub Pages or localhost, else the public GitHub Pages app
const LOCAL_KEY = "sghp-cases-local";
const SWEEPS = {
  T_steam: { label: "Steam temperature, °C", lo: 90, hi: 140, step: 5 },
  T_source: { label: "Heat-source temperature, °C", lo: 20, hi: 80, step: 5 },
  compressor_Hz: { label: "Compressor speed, Hz", lo: 30, hi: 70, step: 5 },
};
const SWEEP_Y = {
  COP: { label: "Heating COP", d: 2, f: r => r.COP },
  Q_heat_kW: { label: "Heat delivered, kW", d: 1, f: r => r.Q_heat_kW },
  W_electric_kW: { label: "Electrical power, kW", d: 1, f: r => r.W_electric_kW },
  P_cond_bar: { label: "Max pressure, bar", d: 1, f: r => r.P_discharge_bar ?? r.P_cond_bar },
  steam_share: { label: "Steam share of heat, %", d: 0, f: r => 100 * r.Q_cond_kW / r.Q_heat_kW },
};
// results table rows: [key, label, decimals, better: 1 higher, -1 lower, 0 neither]
const RESULT_ROWS = [
  ["COP", "Heating COP", 3, 1], ["Q_heat_kW", "Heat delivered, kW", 1, 1], ["Q_cond_kW", "Steam, kW", 1, 1],
  ["Q_subcool_kW", "Hot water, kW", 1, 1], ["W_electric_kW", "Electrical, kW", 1, -1], ["P_discharge_bar", "Max pressure, bar", 2, -1],
  ["P_evap_bar", "Low pressure, bar", 2, 0], ["PR", "Pressure ratio", 2, -1], ["T_discharge", "Discharge, °C", 1, -1],
  ["eta_is_overall", "Isentropic efficiency", 3, 1], ["eta_vol", "Volumetric efficiency", 3, 1],
  ["m_steam_kgh", "Steam, kg/h", 1, 1], ["m_dot_kgh", "Mass flow, kg/h", 0, 0], ["m_water_kgh", "Hot-water flow, kg/h", 0, 0], ["rho_suction_kgm3", "Suction gas density, kg/m³", 2, 0],
  ["V_suction_m3h", "Suction flow, m³/h", 1, -1], ["V_discharge_m3h", "Discharge flow, m³/h", 1, -1],
  ["VHC", "VHC, MJ/m³", 2, 1],
];
const INPUT_ROWS = [
  ["fluid", "Refrigerant"], ["compressor", "Compressor"], ["T_source", "Heat source, °C"], ["T_steam", "Steam, °C"],
  ["T_water_in", "Hot water in, °C"], ["T_water_out", "Hot water out, °C"], ["compressor_Hz", "Speed, Hz"],
  ["dTmin_source", "ΔT source, K"], ["dTmin_steam", "ΔT steam, K"], ["dTmin_water", "ΔT hot water, K"],
  ["dT_superheat", "Suction superheat, K"], ["preheat_subcool_pct", "Feedwater preheat from subcooling, %"], ["dTmin_ihx", "ΔT IHX, K"], ["T_ambient", "Ambient, °C"], ["UA_W_K", "Compressor heat loss UA, W/K"],
  ["UA_cond_W_K", "Condenser heat loss UA, W/K"], ["UA_subcool_W_K", "Subcooler heat loss UA, W/K"], ["UA_evap_W_K", "Evaporator heat loss UA, W/K"],
  ["dPK_cond", "Condenser pressure-loss K, kPa/(m³/h)¹·⁸"], ["dPK_subcool", "Subcooler pressure-loss K, kPa/(m³/h)¹·⁸"], ["dPK_evap", "Evaporator pressure-loss K, kPa/(m³/h)¹·⁸"],
];

const CS = {
  orderBy: "",          // "" (your own order) | "name" | a RESULT_ROWS key: the comparison is kept sorted, best first
  store: null,          // "shared" | "local" | null (loading)
  db: null, user: null, me: null, canWrite: null, downloads: null,
  list: [],             // [{id, ...doc}]
  compare: [],          // ids, in display order (drag to re-order)
  colors: {},           // id -> colour index, fixed while the case is ticked
  loaded: null,         // {id, name, json} of the case last loaded or saved
  fromLink: null,       // a case opened from a link and not yet saved
  solved: new Map(),    // id -> {json, case, dome, label}
  fluids: new Map(),    // fluid key -> CoolPropFluid (compare and sweeps)
  sweep: { x: "T_steam", y: "COP", data: null, key: "", running: false },
  names: {},
  nameAuto: true,       // the name field follows the inputs until someone types a name
};

/* ---------------------------------------------------------------- inputs */
function caseInputs() {
  const v = id => +$(id).value;
  return {
    comps: comps.filter(c => c[1] > 0).map(c => [c[0], +c[1].toFixed(4)]), basis: $("frac-basis").value,
    compressor: COMPRESSOR.name,
    T_source: v("T_source"), T_steam: v("T_steam"), T_water_in: v("T_water_in"), T_water_out: v("T_water_out"),
    compressor_Hz: v("hz"), dTmin_source: v("dTmin_source"), dTmin_steam: v("dTmin_steam"), dTmin_water: v("dTmin_water"),
    dT_superheat: v("dT_superheat"), dTmin_ihx: v("dTmin_ihx"), T_ambient: v("T_ambient"), UA_W_K: v("UA_manual"),
    UA_cond_W_K: v("UA_cond"), UA_subcool_W_K: v("UA_subcool"), UA_evap_W_K: v("UA_evap"), preheat_subcool_pct: v("preheat_sub"),
    dPK_cond: v("dPK_cond"), dPK_subcool: v("dPK_subcool"), dPK_evap: v("dPK_evap"),
  };
}
const inputsJSON = inp => JSON.stringify(inp);
function cleanInputs(raw) {           // anything from the store or a link is untrusted: keep known fields, in range
  if (!raw || typeof raw !== "object") throw new Error("no inputs in this case");
  const opt = (k, lo, hi) => { const x = +raw[k]; return isFinite(x) ? Math.min(hi, Math.max(lo, x)) : 0; };
  const num = (k, lo, hi) => { const x = +raw[k]; if (!isFinite(x)) throw new Error(`${k} is missing`); return Math.min(hi, Math.max(lo, x)); };
  const cs = Array.isArray(raw.comps) ? raw.comps.filter(c => Array.isArray(c) && FLUIDS.includes(c[0])).slice(0, 3) : [];
  if (!cs.length) throw new Error("no known refrigerant in this case");
  return {
    comps: cs.map(c => [c[0], Math.min(100, Math.max(0, +c[1] || 0))]), basis: raw.basis === "mole" ? "mole" : "mass",
    compressor: SGHP.COEFFS[raw.compressor] ? raw.compressor : "Recip-M",
    T_source: num("T_source", 20, 80), T_steam: num("T_steam", 90, 140), T_water_in: num("T_water_in", 20, 120),
    T_water_out: num("T_water_out", 30, 135), compressor_Hz: num("compressor_Hz", 30, 70),
    dTmin_source: num("dTmin_source", 0, 50), dTmin_steam: num("dTmin_steam", 0, 50), dTmin_water: num("dTmin_water", 0, 50),
    dT_superheat: num("dT_superheat", 0, 80), dTmin_ihx: num("dTmin_ihx", 0, 50), T_ambient: num("T_ambient", -40, 60),
    UA_W_K: num("UA_W_K", 0, 1000),
    // added later: cases saved before have no exchanger losses
    UA_cond_W_K: opt("UA_cond_W_K", 0, 1000), UA_subcool_W_K: opt("UA_subcool_W_K", 0, 1000), UA_evap_W_K: opt("UA_evap_W_K", 0, 1000),
    preheat_subcool_pct: opt("preheat_subcool_pct", 0, 100),
    dPK_cond: opt("dPK_cond", 0, 1e4), dPK_subcool: opt("dPK_subcool", 0, 1e4), dPK_evap: opt("dPK_evap", 0, 1e4),
  };
}
function applyInputs(inp) {
  comps.length = 0; for (const c of inp.comps) comps.push([c[0], c[1]]);
  balance(); openPicker = -1;
  $("frac-basis").value = inp.basis;
  for (const o of document.querySelectorAll(".tbtn[data-basis]")) o.setAttribute("aria-pressed", String(o.dataset.basis === inp.basis));
  setCompressor(inp.compressor);
  const set = (id, x) => { $(id).value = x; };
  set("T_source", inp.T_source); set("T_steam", inp.T_steam); set("T_water_in", inp.T_water_in); set("T_water_out", inp.T_water_out);
  set("hz", inp.compressor_Hz); set("dTmin_source", inp.dTmin_source); set("dTmin_steam", inp.dTmin_steam);
  set("dTmin_water", inp.dTmin_water); set("dT_superheat", inp.dT_superheat); set("dTmin_ihx", inp.dTmin_ihx);
  set("T_ambient", inp.T_ambient); set("UA_manual", inp.UA_W_K);
  set("UA_cond", inp.UA_cond_W_K); set("UA_subcool", inp.UA_subcool_W_K); set("UA_evap", inp.UA_evap_W_K); set("preheat_sub", inp.preheat_subcool_pct);
  set("dPK_cond", inp.dPK_cond); set("dPK_subcool", inp.dPK_subcool); set("dPK_evap", inp.dPK_evap);
  compRows(); tempTicks(); hzTicks(); preheatTicks(); superheatTicks(); changed("fluid");
}
const fluidLabel = inp => inp.comps.length === 1 ? LABEL(inp.comps[0][0])
  : `${inp.comps.map(c => LABEL(c[0])).join("/")} ${inp.comps.map(c => +(+c[1]).toFixed(1)).join("/")} (${inp.basis})`;
const fluidKey = inp => inp.comps.map(c => c[0]).join("&") + "|" + inp.comps.map(c => c[1] / 100).join(";") + "|" + inp.basis;
const caseSummary = inp => `${fluidLabel(inp)} · ${inp.T_source} → ${inp.T_steam} °C · ${inp.compressor} · ${inp.compressor_Hz} Hz`;
function suggestName(inp) { return `${fluidLabel(inp)} · ${inp.T_source}/${inp.T_steam} °C · ${inp.compressor} · ${inp.compressor_Hz} Hz`; }
function summarize(r) {
  if (!r) return { ok: false };
  const o = { ok: true };
  for (const [k] of RESULT_ROWS) o[k] = k === "m_dot_kgh" ? r.m_dot_kgs * 3600 : k === "m_water_kgh" ? (r.m_water_kgs == null ? null : r.m_water_kgs * 3600)
    : k === "m_steam_kgh" ? (r.m_steam_kgs == null ? null : r.m_steam_kgs * 3600)
    : k === "P_discharge_bar" ? (r.P_discharge_bar ?? r.P_cond_bar) : r[k];   // max pressure: the compressor discharge
  o.in_envelope = !!r.in_envelope;
  for (const k in o) if (typeof o[k] === "number" && !isFinite(o[k])) o[k] = null;
  return o;
}

/* ---------------------------------------------------------------- solving a case away from the page */
function caseFluid(inp) {
  const k = fluidKey(inp);
  if (LIVE.fluid && LIVE.fluid.key === k) return LIVE.fluid;
  if (CS.fluids.has(k)) return CS.fluids.get(k);
  if (CS.fluids.size >= 6) { const [k0, f0] = CS.fluids.entries().next().value; if (f0 !== LIVE.fluid) f0.destroy(); CS.fluids.delete(k0); }
  const fl = new CoolPropFluid(CP, inp.comps.filter(c => c[1] > 0), inp.basis);
  fl.key = k; CS.fluids.set(k, fl); return fl;
}
const caseCfg = inp => ({
  T_source: inp.T_source, dTmin_source: inp.dTmin_source, T_steam: inp.T_steam, T_water_in: inp.T_water_in,
  T_water_out: inp.T_water_out, dTmin_steam: inp.dTmin_steam, dTmin_water: inp.dTmin_water, dT_superheat: inp.dT_superheat,
  T_ambient: inp.T_ambient, dTmin_ihx: inp.dTmin_ihx, compressor_Hz: inp.compressor_Hz, loss_dT_basis: "mean",
  loss_UA_cond_kWK: inp.UA_cond_W_K / 1000, loss_UA_subcool_kWK: inp.UA_subcool_W_K / 1000, loss_UA_evap_kWK: inp.UA_evap_W_K / 1000,
  preheat_subcool_share: inp.preheat_subcool_pct / 100,
  dP_K_cond: inp.dPK_cond || 0, dP_K_subcool: inp.dPK_subcool || 0, dP_K_evap: inp.dPK_evap || 0,
});
function solveInputs(inp, over = {}) {
  const fl = caseFluid(inp);
  return { fl, res: SGHP.solve(fl, Object.assign(caseCfg(inp), over), SGHP.compressorMap(inp.compressor), inp.UA_W_K / 1000) };
}

/* ---------------------------------------------------------------- storage */
async function casesBoot() {
  casesUI();
  const hash = (location.hash || "").match(/case=([A-Za-z0-9_-]+)/);
  if (hash) openCode(hash[1], "link");
  let db = null, user = null, downloads = null;
  try {
    if (window.claude && typeof window.claude.use === "function") {
      [db, user, downloads] = await Promise.all([claude.use("db"), claude.use("user"), claude.use("downloads")]);
    }
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
    }, err => { caseStatus(`The shared case library stopped updating (${err.code}). Reload the page to reconnect.`, "warn"); });
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
  if (code === "invalid_argument") { CS.canWrite = false; casesChanged(); return "You can view the shared cases but not change them (ask the page owner for Contributor access). Links still work."; }
  if (code === "quota_exceeded") return "The case library is full: delete some cases first.";
  return `Could not save (${code || (e && e.message) || e}).`;
}

/* ---------------------------------------------------------------- actions */
function currentResults() { return LIVE.result && LIVE.result.ok ? summarize(LIVE.result.r) : { ok: false }; }
async function saveCase(asNew) {
  const inp = caseInputs(), name = ($("case-name").value.trim() || suggestName(inp)).slice(0, 120);
  const note = $("case-note").value.trim().slice(0, 500);
  const now = new Date().toISOString();
  const res = currentResults();
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
    CS.compare = CS.compare.filter(x => x !== id);
    if (CS.loaded && CS.loaded.id === id) CS.loaded = null;
    caseStatus(`Deleted “${c.name}”.`); casesChanged();
  } catch (e) { caseStatus(writeError(e), "warn"); }
}
// delete every ticked case this viewer may delete (their own, or any as an editor); a second click confirms
async function deleteSelected() {
  const btn = $("case-del-sel"), ids = CS.compare.filter(id => CS.mayDelete && CS.mayDelete(id));
  if (!ids.length) { caseStatus("None of the ticked cases can be deleted by you (only your own, or any as an editor).", "warn"); return; }
  if (btn.dataset.armed !== "1") {
    btn.dataset.armed = "1"; btn.textContent = `Confirm: delete ${ids.length} case${ids.length > 1 ? "s" : ""}`;
    setTimeout(() => { if (btn.dataset.armed === "1") { btn.dataset.armed = ""; casesChanged(); } }, 5000);
    return;
  }
  btn.dataset.armed = ""; btn.disabled = true;
  let done = 0, err = null;
  for (const id of ids) {             // one write at a time
    try { await storeDelete(id); done++; CS.compare = CS.compare.filter(x => x !== id); delete CS.colors[id]; if (CS.loaded && CS.loaded.id === id) CS.loaded = null; }
    catch (e) { err = e; break; }
  }
  btn.disabled = false;
  const skipped = CS.compare.length;
  caseStatus(err ? `Deleted ${done} case${done === 1 ? "" : "s"}; then ${writeError(err)}`
    : `Deleted ${done} case${done === 1 ? "" : "s"}${skipped ? `; ${skipped} ticked case${skipped > 1 ? "s" : ""} not yours left in place` : ""}.`, err ? "warn" : "");
  casesChanged();
}
function toggleCompare(id, on) {
  if (on && !CS.compare.includes(id)) {
    if (CS.compare.length >= MAX_COMPARE) { caseStatus(`Compare up to ${MAX_COMPARE} cases: untick one first.`, "warn"); casesChanged(); return; }
    CS.compare.push(id);
    const used = new Set(Object.values(CS.colors)); CS.colors[id] = [0, 1, 2, 3, 4].find(k => !used.has(k));
  } else if (!on) { CS.compare = CS.compare.filter(x => x !== id); delete CS.colors[id]; }
  casesChanged();
}
// tick every case for comparison, up to the limit: those already ticked first, then the list order (newest first)
function selectAll() {
  const ids = CS.list.filter(c => { try { cleanInputs(c.inputs); return true; } catch (_) { return false; } }).map(c => c.id);
  for (const id of ids) { if (CS.compare.length >= MAX_COMPARE) break; if (!CS.compare.includes(id)) toggleCompare(id, true); }
  const left = ids.length - CS.compare.length;
  caseStatus(left > 0 ? `Comparing ${MAX_COMPARE} cases (the most); ${left} more not ticked.` : "");
  casesChanged();
}
// each compared case keeps its colour for as long as it is ticked, whatever its position
const caseColor = id => CASE_COLORS[CS.colors[id] != null ? CS.colors[id] : 0];
// re-order the comparison: every view (columns, tables, legend, diagrams, sweep) follows CS.compare
function moveCompare(id, to) {
  const from = CS.compare.indexOf(id); if (from < 0) return;
  to = Math.max(0, Math.min(CS.compare.length - 1, to)); if (to === from) return;
  CS.compare.splice(from, 1); CS.compare.splice(to, 0, id);
  setOrderBy("");                     // a hand re-order replaces the sort
  casesChanged();
}

/* "Order by": keep the comparison sorted by name or by a result, best first (higher first where neither is
   better). Cases not solved yet, or that failed, go last in their current order. */
const ORDER_KEY = "sghp-compare-order";
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
  const rest = CS.compare.filter(id => !order.includes(id)), next = order.concat(rest);
  if (next.join() === CS.compare.join()) return false;
  CS.compare = next; return true;
}
function startDrag(ev, col) {
  if (ev.button > 0) return;
  ev.preventDefault();
  const host = $("cmp-side"), id = col.dataset.id;
  col.classList.add("dragging");      // listen on the window: moving the column in the DOM would drop a pointer capture
  const move = e => {
    // place the column before the first sibling whose centre lies right of the pointer
    const others = [...host.children].filter(x => x !== col);
    const next = others.find(x => { const r = x.getBoundingClientRect(); return e.clientX < r.left + r.width / 2; });
    if (next !== col.nextSibling) host.insertBefore(col, next || null);
    const wrap = host.parentElement, r = wrap.getBoundingClientRect();           // scroll the strip near its edges
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
function caseCode(name, inp) { return b64url(JSON.stringify({ v: CASE_SCHEMA, n: name, i: inp })); }
function caseLink(name, inp) { return `${PAGE_URL}#case=${caseCode(name, inp)}`; }
function openCode(text, how) {
  const m = String(text).trim().match(/(?:case=)?([A-Za-z0-9_-]{20,})\s*$/);
  if (!m) { caseStatus("That doesn't look like a case link.", "warn"); return false; }
  let o, inp;
  try { o = JSON.parse(unb64url(m[1])); inp = cleanInputs(o.i); }
  catch (e) { caseStatus(`Could not read that case link${e && e.message ? `: ${e.message}` : ""}.`, "warn"); return false; }
  const name = typeof o.n === "string" ? o.n.slice(0, 120) : suggestName(inp);
  const go = () => {
    applyInputs(inp);
    CS.loaded = null; CS.fromLink = { name, json: inputsJSON(inp) };
    setCaseName(name, inp); $("case-note").value = "";
    caseStatus(`Opened “${name}” from a ${how}. Save it to add it to the cases.`);
    casesChanged();
  };
  if (CP) go(); else { const t = setInterval(() => { if (CP) { clearInterval(t); go(); } }, 200); }
  return true;
}
async function copyText(text, label) {
  const box = $("case-linkbox"), inp = $("case-linktext");
  inp.value = text; box.hidden = false; inp.select();
  try { await navigator.clipboard.writeText(text); caseStatus(`${label} copied.`); }
  catch (_) { caseStatus(`${label} is selected below: press Ctrl+C (⌘C) to copy it.`); }
}

/* ---------------------------------------------------------------- export */
function csvCell(v) { const s = v == null ? "" : String(v); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; }
async function exportCompare() {
  const cs = compareCases(); if (!cs.length || !CS.downloads) return;
  const rows = [["", ...cs.map(c => c.name)]];
  for (const [k, lab] of INPUT_ROWS) rows.push([lab, ...cs.map(c => inputValue(c.inputs, k))]);
  for (const [k, lab, d] of RESULT_ROWS) rows.push([lab, ...cs.map(c => { const r = solvedResults(c); return r && r.ok && r[k] != null ? (+r[k]).toFixed(d + 2) : ""; })]);
  rows.push(["Case link", ...cs.map(c => caseLink(c.name, c.inputs))]);
  const csv = rows.map(r => r.map(csvCell).join(",")).join("\n");
  const json = JSON.stringify({ schema: CASE_SCHEMA, model: MODEL_VERSION, coolprop: CP_VERSION, exported: new Date().toISOString(),
    cases: cs.map(c => ({ name: c.name, note: c.note || "", inputs: c.inputs, results: solvedResults(c) })) }, null, 2);
  try {
    await CS.downloads.save({ filename: "SGHP_case_comparison.csv", data: csv });
    await CS.downloads.save({ filename: "SGHP_cases.json", data: json });
  } catch (e) { if (e && e.code !== "declined") caseStatus(`Export failed (${e.code || e}).`, "warn"); }
}

/* ---------------------------------------------------------------- page */
/* The case name is automatic (built from the inputs, shown as the field's placeholder and kept up to date)
   until someone types one; clearing the field makes it automatic again. */
function setCaseName(name, inp) {
  CS.nameAuto = !name || name === suggestName(inp);
  $("case-name").value = CS.nameAuto ? "" : name;
}
function caseStatus(text, cls) { const s = $("case-status"); s.textContent = text || ""; s.className = "hint" + (cls ? " " + cls : ""); }
function casesUI() {
  // "Order by": a button and an option list (no native drop-down: some embedded viewers do not open them)
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
  $("case-copy").addEventListener("click", () => { const inp = caseInputs(); copyText(caseLink($("case-name").value.trim() || suggestName(inp), inp), "Link to the current inputs"); });
  $("case-open").addEventListener("click", () => { const b = $("case-openbox"); b.hidden = !b.hidden; if (!b.hidden) $("case-opentext").focus(); });
  $("case-opengo").addEventListener("click", () => { if (openCode($("case-opentext").value, "link")) { $("case-opentext").value = ""; $("case-openbox").hidden = true; } });
  $("case-export").addEventListener("click", exportCompare);
  $("case-all").addEventListener("click", selectAll);
  $("case-none").addEventListener("click", () => { CS.compare = []; CS.colors = {}; casesChanged(); });
  $("case-del-sel").addEventListener("click", deleteSelected);
  // sweep choices as toggle buttons (native drop-downs do not open in some embedded viewers)
  const group = (id, opts, cur, pick) => {
    const box = $(id); box.innerHTML = "";
    for (const [k, s] of Object.entries(opts)) {
      const b = document.createElement("button"); b.type = "button"; b.className = "tbtn"; b.dataset.v = k;
      b.textContent = s.label.replace(/, [^,]*$/, ""); b.setAttribute("aria-pressed", String(k === cur));
      b.addEventListener("click", () => { for (const o of box.children) o.setAttribute("aria-pressed", String(o === b)); pick(k); });
      box.append(b);
    }
  };
  group("sweep-x", SWEEPS, CS.sweep.x, k => { CS.sweep.x = k; runSweep(); });
  group("sweep-y", SWEEP_Y, CS.sweep.y, k => { CS.sweep.y = k; drawSweep(); });
}
let casesDirty = false;
function casesChanged() { if (casesDirty) return; casesDirty = true; requestAnimationFrame(() => { casesDirty = false; casesRender(); }); }
function casesOnRender() { casesHeader(); }

function casesHeader() {
  const inp = caseInputs(), j = inputsJSON(inp);
  const cur = $("case-current");
  if (CS.loaded) cur.textContent = j === CS.loaded.json ? `Showing “${CS.loaded.name}”` : `Showing “${CS.loaded.name}”, modified`;
  else if (CS.fromLink) cur.textContent = j === CS.fromLink.json ? `Showing “${CS.fromLink.name}” from a link (not saved)` : `Opened “${CS.fromLink.name}” from a link, modified`;
  else cur.textContent = "Current inputs are not saved as a case";
  const own = CS.loaded && CS.list.find(c => c.id === CS.loaded.id);
  $("case-save").textContent = own ? "Update case" : "Save case";
  $("case-saveas").hidden = !own;
  const ro = CS.store === "shared" && CS.canWrite === false;
  $("case-save").disabled = $("case-saveas").disabled = ro || !CP;
  $("case-name").placeholder = suggestName(inp);     // an automatic name follows the inputs; a typed one stays
  if (CS.nameAuto && document.activeElement !== $("case-name")) $("case-name").value = "";
  $("case-store").textContent = CS.store === "shared" ? (ro ? "Shared cases, view only" : "Shared with everyone who can open this page")
    : CS.store === "local" ? "Saved in this browser only" : "Loading cases…";
}

let renderSeq = 0;
async function casesRender() {
  const seq = ++renderSeq;
  casesHeader();
  // names of the people who saved cases (ids only are stored)
  const ids = [...new Set(CS.list.map(c => c.createdBy).filter(Boolean))];
  if (CS.user && ids.length) { try { const ps = await CS.user.profiles(ids); for (const id of ids) CS.names[id] = ps[id] ? (ps[id].isMe ? "you" : ps[id].name) : ""; } catch (_) {} }
  let editor = false;
  if (CS.user) { try { editor = await CS.user.canEdit(); } catch (_) {} }
  CS.mayDelete = id => { const c = CS.list.find(x => x.id === id); return !!c && (CS.store === "local" || (c.createdBy && c.createdBy === CS.me) || editor); };
  if (seq !== renderSeq) return;                       // a newer render is under way
  CS.compare = CS.compare.filter(id => CS.list.some(c => c.id === id)).slice(0, MAX_COMPARE);
  const box = $("case-list"); box.innerHTML = "";
  if (!CS.list.length) {
    const p = document.createElement("p"); p.className = "hint";
    p.textContent = CS.store ? "No cases yet. Set the inputs, then Save case." : "Loading cases…";
    box.append(p);
  }
  for (const c of CS.list) {
    let inp; try { inp = cleanInputs(c.inputs); } catch (_) { continue; }
    const ci = CS.compare.indexOf(c.id), r = c.results || {};
    const card = document.createElement("div"); card.className = "case-card" + (ci >= 0 ? " on" : "") + (CS.loaded && CS.loaded.id === c.id ? " loaded" : "");
    if (ci >= 0) card.style.setProperty("--case", `var(${caseColor(c.id)})`);
    const top = document.createElement("label"); top.className = "case-top";
    const cb = document.createElement("input"); cb.type = "checkbox"; cb.checked = ci >= 0; cb.setAttribute("aria-label", `Compare ${c.name}`);
    cb.addEventListener("change", () => toggleCompare(c.id, cb.checked));
    const nm = document.createElement("span"); nm.className = "case-name"; nm.textContent = c.name || "Untitled case";
    const cop = document.createElement("span"); cop.className = "case-cop"; cop.textContent = r.ok ? `COP ${fmt(r.COP, 2)}` : "no cycle";
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
    b("Copy link", () => copyText(caseLink(c.name, inp), `Link to “${c.name}”`));
    const mayDelete = CS.store === "local" || (c.createdBy && c.createdBy === CS.me) || editor;   // your own cases, or any as an editor
    if (mayDelete && !(CS.store === "shared" && CS.canWrite === false)) b("Delete", () => deleteCase(c.id), { "data-del": c.id });
    card.append(act);
    box.append(card);
  }
  const n = CS.list.length;
  $("case-all").hidden = n < 2 || CS.compare.length >= Math.min(n, MAX_COMPARE);
  $("case-all").textContent = n > MAX_COMPARE ? `Select ${MAX_COMPARE}` : "Select all";
  $("case-none").hidden = !CS.compare.length;
  const delBtn = $("case-del-sel");
  delBtn.hidden = !CS.compare.length || (CS.store === "shared" && CS.canWrite === false);
  if (delBtn.dataset.armed !== "1") delBtn.textContent = CS.compare.length > 1 ? `Delete selected (${CS.compare.length})` : "Delete selected";
  $("case-count").textContent = CS.compare.length ? `${CS.compare.length} of ${MAX_COMPARE} ticked to compare` : `Tick up to ${MAX_COMPARE} cases to compare them`;
  renderCompare();
}

/* ---------------------------------------------------------------- comparison */
const compareCases = () => CS.compare.map(id => CS.list.find(c => c.id === id)).filter(Boolean)
  .map(c => { try { return Object.assign({}, c, { inputs: cleanInputs(c.inputs) }); } catch (_) { return null; } }).filter(Boolean);
function solvedResults(c) { const s = CS.solved.get(c.id); return s && s.json === inputsJSON(c.inputs) ? summarize(s.case.ok ? s.case.r : null) : null; }
function inputValue(inp, k) {
  if (k === "fluid") return fluidLabel(inp);
  return inp[k];
}
let compareJob = 0;
function renderCompare() {
  const host = $("compare"), cs = compareCases();
  host.hidden = cs.length < 2;
  $("case-export").disabled = cs.length < 1;
  if (cs.length < 2) return;
  if (applyOrder(cs)) { casesChanged(); return; }       // sorted: render again in the new order
  // solve each compared case (states and paths for the diagrams), a case per frame
  const job = ++compareJob, todo = cs.filter(c => !solvedResults(c) && !(CS.solved.get(c.id) || {}).failed);
  compareTables(cs);
  const step = () => {
    if (job !== compareJob) return;
    if (!CP) { setTimeout(step, 300); return; }
    const c = todo.shift();
    if (c) {
      try {
        const { fl, res } = solveInputs(c.inputs);
        CS.solved.set(c.id, { json: inputsJSON(c.inputs), case: toCase(res, fl), dome: SGHP.dome(fl), label: fl.label });
      } catch (e) { CS.solved.set(c.id, { json: inputsJSON(c.inputs), case: { ok: false, msg: cpError(e) }, dome: null, failed: true }); }
      compareTables(cs);
      setTimeout(step, 0); return;
    }
    if (applyOrder(cs)) { casesChanged(); return; }     // every case solved: sort by its results
    compareTables(cs); sideBySide(cs); overlay("cmp-ts", cs, "ts"); overlay("cmp-ph", cs, "ph");
    runSweep();
  };
  $("cmp-ts").innerHTML = $("cmp-ph").innerHTML = `<div class="empty">Solving the cases…</div>`;
  sideBySide(cs);
  setTimeout(step, 0);
}

function compareTables(cs) {
  legend("lg-cmp", cs.map((c, i) => ({ color: caseColor(c.id), label: c.name })));
  const msgs = cs.map(c => { const s = CS.solved.get(c.id); return s && s.case && !s.case.ok ? `${c.name}: ${s.case.msg}` : null; }).filter(Boolean);
  const env = cs.filter(c => { const r = solvedResults(c); return r && r.ok && !r.in_envelope; }).map(c => c.name);
  const drift = cs.filter(c => { const r = solvedResults(c), s = c.results; return r && r.ok && s && s.ok && Math.abs(r.COP - s.COP) > 1e-4 * Math.abs(s.COP); }).map(c => c.name);
  const notes = [...msgs.map(m => `No cycle for ${m}`), ...(env.length ? [`Outside the compressor envelope (efficiencies extrapolated): ${env.join(", ")}`] : []),
    ...(drift.length ? [`Re-solved results differ from when these were saved (model or CoolProp changed): ${drift.join(", ")}`] : [])];
  $("cmp-notes").innerHTML = ""; for (const n of notes) { const s = document.createElement("span"); s.className = "chip warn"; s.textContent = n; $("cmp-notes").append(s); }
  // inputs that differ
  const ti = document.createElement("table"); ti.className = "cmp";
  const hi = ti.createTHead().insertRow(); const h0 = document.createElement("th"); h0.textContent = "Input"; hi.append(h0);
  for (const c of cs) { const th = document.createElement("th"); th.className = "ch"; const sw = document.createElement("i"); sw.className = "sw"; sw.style.background = `var(${caseColor(c.id)})`; const nm = document.createElement("span"); nm.textContent = c.name; th.append(sw, nm); hi.append(th); }
  const bi = ti.createTBody(); let ndiff = 0;
  for (const [k, lab] of INPUT_ROWS) {
    const vals = cs.map(c => inputValue(c.inputs, k));
    if (vals.every(v => String(v) === String(vals[0]))) continue;
    ndiff++;
    const tr = bi.insertRow(); tr.insertCell().textContent = lab;
    vals.forEach(v => { tr.insertCell().textContent = String(v); });
  }
  const same = INPUT_ROWS.filter(([k]) => { const v = cs.map(c => inputValue(c.inputs, k)); return v.every(x => String(x) === String(v[0])); }).map(([, l]) => l.replace(/, .*$/, ""));
  $("cmp-inputs").innerHTML = "";
  if (ndiff) $("cmp-inputs").append(ti); else $("cmp-inputs").innerHTML = `<p class="hint">These cases have identical inputs.</p>`;
  $("cmp-same").textContent = same.length && ndiff ? `Same in every case: ${same.join(", ")}.` : "";
}

/* side-by-side dashboard: one column per case, the headline gauges on scales shared by all the cases,
   key figures, the heat and power split, and the case's own T-s cycle */
function sideBySide(cs) {
  const host = $("cmp-side"); host.innerHTML = "";
  const rs = cs.map(c => { const s = CS.solved.get(c.id); return s && s.json === inputsJSON(c.inputs) ? s.case : null; });
  const ok = rs.filter(r => r && r.ok).map(r => r.r);
  const kwTop = ok.length ? dialTop(...ok.flatMap(r => [r.Q_cond_kW, r.Q_subcool_kW, r.W_electric_kW])) : 10;
  const pTop = ok.length ? dialTop(...ok.map(r => r.P_discharge_bar ?? r.P_cond_bar)) : 10;
  const powTop = ok.length ? Math.max(...ok.map(r => r.Q_heat_kW + r.W_electric_kW)) : 1;
  // one T-s frame for every column, so the cycles line up
  const pts = []; cs.forEach(c => { const s = CS.solved.get(c.id); if (s && s.case.ok) { s.case.paths.forEach(p => pts.push(...p.pts)); if (s.dome) pts.push(...s.dome.liq, ...s.dome.vap); } });
  const sx = pts.map(p => p[3]).filter(isFinite), sy = pts.map(p => p[0]).filter(isFinite);
  const frameTs = sx.length ? { xr: [Math.min(...sx) - 0.04 * (Math.max(...sx) - Math.min(...sx)), Math.max(...sx) + 0.04 * (Math.max(...sx) - Math.min(...sx))],
    yr: [Math.floor((Math.min(...sy) - 5) / 10) * 10, Math.ceil((Math.max(...sy) + 5) / 10) * 10] } : {};
  // the best value per metric, where higher or lower is better; only when the values differ
  const okRs = rs.map(cc => (cc && cc.ok ? cc.r : null));
  const mark = (f, better) => {
    const v = okRs.map(r => (r ? f(r) : null)), fin = v.filter(x => x != null && isFinite(x));
    if (!better || fin.length < 2) return () => "";
    const hi = Math.max(...fin), lo = Math.min(...fin), tol = 1e-9 * Math.max(1, Math.abs(hi));
    if (hi - lo <= tol) return () => "";
    const best = better > 0 ? hi : lo;
    return i => (v[i] != null && Math.abs(v[i] - best) <= tol ? "win" : "");
  };
  const M_ = {
    COP: mark(r => r.COP, 1), steam: mark(r => r.Q_cond_kW, 1), hw: mark(r => r.Q_subcool_kW, 1), el: mark(r => r.W_electric_kW, -1),
    pmax: mark(r => r.P_discharge_bar ?? r.P_cond_bar, -1), heat: mark(r => r.Q_heat_kW, 1),
  };
  const KV = [
    ["Heat delivered, kW", r => r.Q_heat_kW, 1, 1], ["Steam, kg/h", r => (r.m_steam_kgs == null ? NaN : r.m_steam_kgs * 3600), 1, 1],
    ["Hot water, kg/h", r => (r.m_water_kgs == null ? NaN : r.m_water_kgs * 3600), 0, 0],
    ["Low pressure, bar", r => r.P_evap_bar, 2, 0], ["Max pressure, bar", r => r.P_discharge_bar ?? r.P_cond_bar, 1, -1],
    ["Pressure ratio", r => r.PR, 2, -1], ["Discharge, °C", r => r.T_discharge, 1, -1], ["Isentropic efficiency", r => r.eta_is_overall, 3, 1],
    ["Volumetric efficiency", r => r.eta_vol, 3, 1], ["Mass flow, kg/h", r => r.m_dot_kgs * 3600, 0, 0],
    ["Suction density, kg/m³", r => r.rho_suction_kgm3, 2, 0],
    ["Suction flow, m³/h", r => r.V_suction_m3h, 1, -1], ["Discharge flow, m³/h", r => r.V_discharge_m3h, 1, -1],
    ["VHC, MJ/m³", r => r.VHC, 2, 1],
  ].map(([lab, f, d, better]) => ({ lab, f, d, cls: mark(f, better) }));
  // the left-most case is the baseline: the others show their change from it
  const base = okRs[0];
  const delta = (v, b, i) => {
    if (i === 0 || b == null || v == null || !isFinite(v) || !isFinite(b) || Math.abs(b) < 1e-12) return "";
    const pc = 100 * (v - b) / b;
    return Math.abs(pc) < 0.05 ? " (=)" : ` (${pc > 0 ? "+" : "−"}${Math.abs(pc).toFixed(1)} %)`;
  };
  cs.forEach((c, i) => {
    const ci = CS.compare.indexOf(c.id), col = document.createElement("div");
    col.className = "side-col"; col.style.setProperty("--case", `var(${caseColor(c.id)})`); col.dataset.id = c.id;
    // drag handle (pointer: mouse, pen and touch) plus left/right buttons for the keyboard
    const head = document.createElement("div"); head.className = "side-head";
    const grip = document.createElement("span"); grip.className = "grip"; grip.textContent = "⠿"; grip.title = "Drag to re-order"; grip.setAttribute("aria-hidden", "true");
    const h = document.createElement("h4"); h.textContent = c.name;
    const mv = (d, lab, sym) => { const b = document.createElement("button"); b.type = "button"; b.className = "icon mv"; b.textContent = sym; b.setAttribute("aria-label", `Move ${c.name} ${lab}`); b.disabled = (d < 0 && i === 0) || (d > 0 && i === cs.length - 1); b.addEventListener("click", () => moveCompare(c.id, i + d)); return b; };
    head.append(grip, h, mv(-1, "left", "‹"), mv(1, "right", "›"));
    if (i === 0) { const bt = document.createElement("span"); bt.className = "base-tag"; bt.textContent = "baseline"; head.append(bt); }
    grip.addEventListener("pointerdown", ev => startDrag(ev, col));
    const sub = document.createElement("div"); sub.className = "sub"; sub.textContent = caseSummary(c.inputs);
    col.append(head, sub);
    const cc = rs[i];
    if (!cc) { const e = document.createElement("div"); e.className = "empty"; e.textContent = "Solving…"; col.append(e); host.append(col); return; }
    if (!cc.ok) { const e = document.createElement("div"); e.className = "empty"; e.textContent = `No cycle: ${cc.msg}`; col.append(e); host.append(col); return; }
    const r = cc.r, n = copStars(r.COP);
    const cop = document.createElement("div"); cop.className = "side-cop";
    cop.className = "side-cop " + M_.COP(i);
    cop.innerHTML = `<span class="v">${fmt(r.COP, 2)}</span><span class="stars" role="img" aria-label="${n.toFixed(2)} of 5 stars">${starRow(n, `sst-${c.id}`)}</span><small>heating COP${delta(r.COP, base ? base.COP : null, i).replace(/^ \((.*)\)$/, " · $1 vs baseline")}${r.in_envelope ? "" : " · outside envelope"}</small>`;
    const dials = document.createElement("div"); dials.className = "side-dials";
    dials.innerHTML = `<div class="dial ${M_.steam(i)}">${dialSVG(r.Q_cond_kW, kwTop, "--s3", "kW steam")}<span class="dk">Steam, kW</span></div>
      <div class="dial ${M_.hw(i)}">${dialSVG(r.Q_subcool_kW, kwTop, "--s4", "kW hot water")}<span class="dk">Hot water, kW</span></div>
      <div class="dial ${M_.el(i)}">${dialSVG(r.W_electric_kW, kwTop, "--s2", "kW electrical")}<span class="dk">Electrical, kW</span></div>
      <div class="dial ${M_.pmax(i)}">${dialSVG(r.P_discharge_bar ?? r.P_cond_bar, pTop, "--s6", "bar maximum pressure")}<span class="dk">Max pressure, bar</span></div>`;
    // heat delivered and electrical input as one bar on a common length, so sizes compare across columns
    const bar = document.createElement("div"); bar.className = "side-bar"; bar.setAttribute("role", "img");
    bar.setAttribute("aria-label", `steam ${fmt(r.Q_cond_kW, 1)} kW, hot water ${fmt(r.Q_subcool_kW, 1)} kW, electrical ${fmt(r.W_electric_kW, 1)} kW`);
    for (const [v, cv] of [[r.Q_cond_kW, "--s3"], [r.Q_subcool_kW, "--s4"], [r.W_electric_kW, "--s2"]]) { const s = document.createElement("i"); s.style.width = `${(100 * v / powTop).toFixed(2)}%`; s.style.background = `var(${cv})`; bar.append(s); }
    const bl = document.createElement("div"); bl.className = "side-barlab"; bl.textContent = `Steam : hot water ${(100 * r.Q_cond_kW / r.Q_heat_kW).toFixed(0)} : ${(100 * r.Q_subcool_kW / r.Q_heat_kW).toFixed(0)} · ${fmt(r.Q_heat_kW, 1)} kW heat`;
    const kv = document.createElement("dl"); kv.className = "side-kv";
    for (const m of KV) {
      const dt = document.createElement("dt"), dd = document.createElement("dd"), cl = m.cls(i);
      dt.textContent = m.lab; dd.textContent = fmt(m.f(r), m.d);
      const dl = delta(m.f(r), base ? m.f(base) : null, i);
      if (dl) { const sm = document.createElement("small"); sm.className = "dlt"; sm.textContent = dl; dd.append(sm); }
      if (cl) dd.className = cl;
      kv.append(dt, dd);
    }
    const mini = document.createElement("div"); mini.className = "mini"; mini.id = `side-ts-${c.id}`;
    col.append(cop, dials, bar, bl, kv, mini);
    host.append(col);
    overlay(mini.id, [c], "ts", Object.assign({ W: 380, H: 300 }, frameTs));          // smaller frame: legible text in a narrow column
  });
}

function overlay(id, cs, mode, size = { W: 560, H: 340 }) {
  const host = $(id), X = mode === "ts" ? 3 : 2, Y = mode === "ts" ? 0 : 1;
  const items = cs.map(c => ({ c, s: CS.solved.get(c.id) })).filter(o => o.s && o.s.case.ok);
  if (!items.length) { host.innerHTML = `<div class="empty">No cycles to draw</div>`; return; }
  const all = [];
  for (const { s } of items) { s.case.paths.forEach(p => all.push(...p.pts)); if (s.dome) all.push(...s.dome.liq, ...s.dome.vap); }
  const xs = all.map(p => p[X]).filter(isFinite), ys = all.map(p => p[Y]).filter(v => isFinite(v) && (mode === "ts" || v > 0));
  let xr = [Math.min(...xs), Math.max(...xs)], yr = [Math.min(...ys), Math.max(...ys)];
  const xpad = (xr[1] - xr[0]) * 0.04; xr = [xr[0] - xpad, xr[1] + xpad];
  if (mode === "ts") yr = [Math.floor((yr[0] - 5) / 10) * 10, Math.ceil((yr[1] + 5) / 10) * 10]; else yr = [Math.max(0.8, yr[0] * 0.8), yr[1] * 1.25];
  if (size.xr) xr = size.xr;
  if (size.yr) yr = size.yr;
  const f = frame(host, { xr, yr, ylog: mode === "ph", xlab: mode === "ts" ? "Specific entropy, kJ/(kg·K)" : "Specific enthalpy, kJ/kg", ylab: mode === "ts" ? "Temperature, °C" : "Pressure, bar (log)", W: size.W, H: size.H });
  const seen = new Set(), hov = [], oneFluid = new Set(items.map(o => o.s.label)).size === 1;
  for (const { c, s } of items) {
    const col = css(caseColor(c.id));
    if (s.dome && !seen.has(s.label)) {          // each refrigerant's dome: grey for one fluid, else in its first case's colour
      seen.add(s.label);
      const dp = [...s.dome.liq, ...[...s.dome.vap].reverse()].map(p => [p[X], p[Y]]).filter(p => isFinite(p[0]) && isFinite(p[1]) && (mode === "ts" || p[1] > 0));
      el("path", { d: poly(dp, f.fx, f.fy), fill: "none", stroke: oneFluid ? css("--dome") : col, "stroke-opacity": oneFluid ? 1 : 0.45, "stroke-width": 1.2, "stroke-dasharray": "5 4" }, f.svg);
    }
    for (const p of s.case.paths) {
      el("path", { d: poly(p.pts.map(q => [q[X], q[Y]]), f.fx, f.fy), fill: "none", stroke: col, "stroke-width": 2.2, "stroke-linecap": "round", "stroke-linejoin": "round", "stroke-dasharray": p.c === "Compressor" ? "6 3" : "none" }, f.svg);
      for (const q of p.pts) hov.push({ x: f.fx(q[X]), y: f.fy(q[Y]), html: `<b>${escapeHTML(c.name)}</b><br>${p.c}<br>T ${fmt(q[0], 1)} °C · P ${fmt(q[1], 2)} bar<br>h ${fmt(q[2], 1)} kJ/kg · s ${fmt(q[3], 4)}` });
    }
  }
  hoverLayer(host, f, hov);
}
const escapeHTML = s => String(s).replace(/[&<>"']/g, ch => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[ch]));

/* sweep: one operating variable across its range for every compared case, a slice per frame */
function runSweep() {
  const cs = compareCases(); if (cs.length < 2 || !CP) return;
  const sx = SWEEPS[CS.sweep.x], key = CS.sweep.x + "|" + cs.map(c => c.id + inputsJSON(c.inputs)).sort().join("|");   // order-free: re-ordering does not re-solve
  if (CS.sweep.key === key && CS.sweep.data) { drawSweep(); return; }
  CS.sweep.key = key; CS.sweep.data = cs.map(c => ({ id: c.id, name: c.name, pts: [] }));
  const xs = []; for (let x = sx.lo; x <= sx.hi + 1e-9; x += sx.step) xs.push(x);
  const jobs = []; cs.forEach((c, i) => xs.forEach(x => jobs.push([i, c, x])));
  let n = 0; const total = jobs.length;
  const step = () => {
    if (CS.sweep.key !== key) return;
    const t0 = performance.now();
    while (jobs.length && performance.now() - t0 < 30) {
      const [i, c, x] = jobs.shift(); n++;
      let r = null;
      try { const { res } = solveInputs(c.inputs, { [CS.sweep.x]: x }); if (res.success) r = res.report; } catch (_) {}
      CS.sweep.data[i].pts.push([x, r]);
    }
    $("sweep-status").textContent = jobs.length ? `Solving ${n} of ${total} points…` : `${CS.sweep.data.reduce((a, d) => a + d.pts.filter(p => p[1]).length, 0)} of ${total} points solved`;
    drawSweep();
    if (jobs.length) requestAnimationFrame(step);
  };
  requestAnimationFrame(step);
}
function drawSweep() {
  const host = $("cmp-sweep"), d = CS.sweep.data; if (!d) return;
  const sx = SWEEPS[CS.sweep.x], sy = SWEEP_Y[CS.sweep.y];
  const ys = []; for (const s of d) for (const [, r] of s.pts) if (r) { const v = sy.f(r); if (isFinite(v)) ys.push(v); }
  if (!ys.length) { host.innerHTML = `<div class="empty">No solved points yet</div>`; return; }
  let lo = Math.min(...ys), hi = Math.max(...ys); const pad = (hi - lo) * 0.08 || Math.abs(hi) * 0.05 || 1; lo -= pad; hi += pad;
  const f = frame(host, { xr: [sx.lo, sx.hi], yr: [lo, hi], xlab: sx.label, ylab: sy.label, W: 1000, H: 360 });   // full-width panel
  const hov = [];
  for (const s of d) {
    if (!CS.compare.includes(s.id)) continue;
    const col = css(caseColor(s.id));
    let run = [];
    const flush = () => { if (run.length > 1) el("path", { d: poly(run, f.fx, f.fy), fill: "none", stroke: col, "stroke-width": 2.2, "stroke-linejoin": "round" }, f.svg); run = []; };
    for (const [x, r] of s.pts) {
      const v = r ? sy.f(r) : NaN;
      if (!isFinite(v)) { flush(); continue; }
      run.push([x, v]);
      el("circle", { cx: f.fx(x), cy: f.fy(v), r: 3, fill: col }, f.svg);
      hov.push({ x: f.fx(x), y: f.fy(v), html: `<b>${escapeHTML(s.name)}</b><br>${sx.label.replace(/, .*$/, "")} ${x}<br>${sy.label} ${fmt(v, sy.d)}` });
    }
    flush();
  }
  hoverLayer(host, f, hov);
}
