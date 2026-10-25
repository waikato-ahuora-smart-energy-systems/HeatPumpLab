
/* SGHP live page: CoolProp (WebAssembly) properties + the SGHP cycle model, solved in the browser.
   Assembled into index.html by build_sghp_live.py together with the dashboard renderers. */

/* ---------------------------------------------------------------- CoolProp loading */
const STATUS = msg => { const s = document.getElementById("status"); s.textContent = msg; };
let CP = null, CP_VERSION = "";

/* A C++ exception thrown out of CoolProp's WebAssembly into JavaScript leaves the WebAssembly stack
   pointer where the throw happened; after about 60 of them the stack overruns the heap and CoolProp
   aborts ("native code called abort()"). Every CoolProp call that can throw is therefore made through
   cpCall, which saves the stack pointer first and restores it on an exception. The build exposes the
   two stack helpers as sghpStackSave / sghpStackRestore (see build_sghp_live.py). */
function cpCall(fn) {
  const sp = CP && CP.sghpStackSave ? CP.sghpStackSave() : null;
  try { return fn(); }
  catch (e) { if (sp !== null) CP.sghpStackRestore(sp); throw new Error(cpError(e)); }
}
// an AbstractState whose methods all go through cpCall
function guarded(AS) {
  return new Proxy(AS, { get(t, k) { const v = t[k]; return typeof v === "function" ? (...a) => cpCall(() => v.apply(t, a)) : v; } });
}

function cpError(e) {
  if (typeof e === "number" && CP && typeof CP.getExceptionMessage === "function") {
    try { const m = CP.getExceptionMessage(e); return Array.isArray(m) ? m.join(": ") : String(m); } catch (_) {}
  }
  return (e && e.message) ? e.message : String(e);
}

async function loadCoolProp() {
  const locateFile = p => (p.endsWith(".wasm") ? "coolprop.wasm" : p);
  if (window.__cpFactoryPromise && !(typeof Module === "object" && Module)) {        // ES-module build (CoolProp 8)
    const factory = await Promise.race([window.__cpFactoryPromise,
      new Promise((_, rej) => setTimeout(() => rej(new Error("CoolProp loader did not start")), 20000))]);
    return await factory({ locateFile });
  }
  if (typeof Module === "function") return await Module({ locateFile });         // MODULARIZE build
  if (typeof createCoolPropModule === "function") return await createCoolPropModule({ locateFile });
  if (typeof Module === "object" && Module) {
    await window.__cpReadyPromise;
    return Module;
  }
  throw new Error("CoolProp module not found in the page");
}

/* ---------------------------------------------------------------- property adapter (MPa, C, kJ/kg) */
const FLUIDS = ["n-Butane", "IsoButane", "Propane", "Propylene", "n-Pentane", "Isopentane", "Neopentane", "n-Hexane",
  "Isohexane", "Cyclopentane", "Cyclohexane", "DimethylEther", "Ammonia", "Water", "CarbonDioxide", "R1233zd(E)",
  "R1234ze(Z)", "R1234ze(E)", "R1234yf", "R245fa", "R365MFC", "R134a", "R152A", "Toluene", "Acetone"];
const LABEL = n => ({ "IsoButane": "Isobutane", "DimethylEther": "Dimethyl ether", "CarbonDioxide": "CO₂", "R365MFC": "R365mfc", "R152A": "R152a" }[n] || n);

class CoolPropFluid {
  constructor(M, comps, basis = "mass") {
    this.M = M; this.comps = comps; this.basis = basis; this.notes = [];
    const names = comps.map(c => c[0]), x = comps.map(c => c[1] / 100);
    this.key = names.join("&") + "|" + x.join(";") + "|" + basis;
    this.label = comps.length === 1 ? LABEL(names[0]) : `${names.map(LABEL).join("/")} ${comps.map(c => +c[1].toFixed(2)).join("/")} (${basis})`;
    this.AS = this._make(names, x);
    this.SAT = this._make(names, x);
    this.isMix = names.length > 1;
    if (this.isMix) {   // PropsSI string with mole fractions, for fast phase-imposed single-phase flashes
      let z = this.AS.get_mole_fractions();
      if (z && typeof z.size === "function") { const a = []; for (let i = 0; i < z.size(); i++) a.push(z.get(i)); z.delete && z.delete(); z = a; }
      this.psFluid = "HEOS::" + names.map((n, i) => `${n}[${z[i]}]`).join("&");
    }
    this.satCache = new Map();
    this.pairs = M.input_pairs;
    try { this.T_crit = this.AS.T_critical() - 273.15; this.P_crit = this.AS.p_critical() / 1e6; }
    catch (e) { this.T_crit = NaN; this.P_crit = NaN; }
    if (!(this.T_crit > -100 && this.T_crit < 1000)) { this.T_crit = NaN; this.P_crit = NaN; }
    try { this.rho_crit = isFinite(this.T_crit) ? this.AS.rhomass_critical() : NaN; } catch (e) { this.rho_crit = NaN; }
    if (!(this.rho_crit > 0 && this.rho_crit < 5000)) this.rho_crit = NaN;
    this.cache = new Map(); this.dewCache = new Map(); this.fallbacks = 0;
    try { this.MW = this.AS.molar_mass(); } catch (e) { this.MW = NaN; }
    if (this.isMix) {
      this.env = this._envelope();
      if (!isFinite(this.T_crit) && this.env && this.env.crit) {
        this.T_crit = this.env.crit.T; this.P_crit = this.env.crit.P; this.rho_crit = this.env.crit.rho * this.MW;
        this.notes.push(`critical point estimated from the phase envelope: ${this.T_crit.toFixed(1)} C, ${(this.P_crit * 10).toFixed(2)} bar`);
      } else if (!isFinite(this.T_crit)) this.notes.push("critical point not available for this mixture");
      if (isFinite(this.T_crit) && !(this.rho_crit > 0) && this.env && this.env.crit) this.rho_crit = this.env.crit.rho * this.MW;
    }
  }

  /* ---- fallbacks for mixtures near the critical point (CoolProp's own routines may not converge there)
     1. Phase envelope: AbstractState.build_phase_envelope traces the dew and bubble lines through the
        critical point by continuation; get_phase_envelope_data gives T, p and the bulk ("vap") and
        incipient ("liq") phase densities along it. Bubble and dew points are interpolated from it
        (ln p against T, ln rho against T) when the pressure-quality or temperature-quality flash fails.
     2. Density-temperature evaluation: the Helmholtz equation of state is explicit in (rho, T), so a
        phase-imposed evaluation always succeeds; single-phase states and the interpolated saturation states
        are found by solving p(rho, T) = P for rho (and an outer search on T for given h or s).
     A temperature-quality point that fails is first re-solved exactly with pressure-quality flashes
     (which converge closer to the critical point), started from the envelope's pressure. The approximate
     routes (1 and 2) are counted in this.fallbacks and reported on the page. */
  _envelope() {
    let d;
    try { this.SAT.build_phase_envelope(""); d = this.SAT.get_phase_envelope_data(); } catch (e) { return null; }
    const arr = v => { if (!v) return []; if (Array.isArray(v)) return v.slice(); const a = []; for (let i = 0; i < v.size(); i++) a.push(v.get(i)); return a; };
    const T = arr(d.T).map(t => t - 273.15), P = arr(d.p).map(p => p / 1e6), rb = arr(d.rhomolar_vap), ri = arr(d.rhomolar_liq);
    if (T.length < 10) return null;
    // the trace starts on the dew line (bulk vapour, incipient liquid denser) and crosses the critical point,
    // where the bulk and incipient densities meet, onto the bubble line
    let ic = -1;
    for (let i = 1; i < T.length; i++) if ((ri[i - 1] - rb[i - 1]) > 0 && (ri[i] - rb[i]) <= 0) { ic = i; break; }
    if (ic < 0) return null;
    const f = (ri[ic - 1] - rb[ic - 1]) / ((ri[ic - 1] - rb[ic - 1]) - (ri[ic] - rb[ic]));
    const crit = { T: T[ic - 1] + f * (T[ic] - T[ic - 1]), P: P[ic - 1] + f * (P[ic] - P[ic - 1]), rho: rb[ic - 1] + f * (rb[ic] - rb[ic - 1]) };   // rho: mol/m3
    // each branch, ordered by rising pressure and cut where the pressure stops rising (retrograde region)
    // (points where T or p does not rise are skipped, so either can be the lookup key; this drops the
    // retrograde tip between the cricondentherm and the critical point, and any step-size jitter)
    const branch = (idx) => { const out = []; for (const i of idx) { const l = out[out.length - 1]; if (l && (P[i] <= l.P || T[i] <= l.T)) continue; if (rb[i] > 0) out.push({ T: T[i], P: P[i], lnP: Math.log(P[i]), lnr: Math.log(rb[i]) }); } return out; };
    const dew = branch([...Array(ic).keys()]);
    const bub = branch([...Array(T.length - ic).keys()].map(k => T.length - 1 - k));
    return { crit, dew, bub };
  }
  _interp(br, key, x, outs) {
    let lo = 0, hi = br.length - 1;
    if (!(x >= br[0][key] && x <= br[hi][key])) return null;
    while (hi - lo > 1) { const m = (lo + hi) >> 1; if (br[m][key] <= x) lo = m; else hi = m; }
    const t = (x - br[lo][key]) / (br[hi][key] - br[lo][key] || 1);
    const o = {}; for (const k of outs) o[k] = br[lo][k] + t * (br[hi][k] - br[lo][k]); return o;
  }
  /* Explicit Helmholtz evaluation at molar density and temperature. The phase is imposed (PropsSI
     "Dmolar|gas" / "Dmolar|liquid"): AbstractState.update(DmolarT) on a mixture also runs a phase-stability
     analysis, 15-400 ms a call in WebAssembly, against about 1.5 ms here. */
  _pDT(rhomol, TK, ph) {
    const p = cpCall(() => this.M.PropsSI("P", `Dmolar|${ph}`, rhomol, "T", TK, this.psFluid));
    if (!isFinite(p)) throw new Error("density-temperature evaluation failed");
    return p / 1e6;
  }
  _dt(rhomol, TK, ph) {
    const f = k => cpCall(() => this.M.PropsSI(k, `Dmolar|${ph}`, rhomol, "T", TK, this.psFluid));
    const st = { T: TK - 273.15, P: f("P") / 1e6, h: f("H") / 1000, s: f("S") / 1000, D: rhomol * this.MW, rho: rhomol };
    if (![st.P, st.h, st.s, st.D].every(isFinite)) throw new Error("density-temperature evaluation failed");
    return st;
  }
  _rhoAt(TK, P, rho0, ph) {  // p(rho, T) = P near rho0 (secant with a positive-density guard)
    let r0 = rho0, p0 = this._pDT(r0, TK, ph) - P, r1 = rho0 * 1.01, p1 = this._pDT(r1, TK, ph) - P;
    for (let k = 0; k < 60 && Math.abs(p1) > 1e-9 * P; k++) {
      let r2 = r1 - p1 * (r1 - r0) / (p1 - p0 || 1e-30);
      if (!(r2 > 0) || !isFinite(r2)) r2 = r1 * (p1 > 0 ? 0.7 : 1.3);
      r0 = r1; p0 = p1; r1 = r2; p1 = this._pDT(r1, TK, ph) - P;
    }
    if (!(Math.abs(p1) <= 1e-6 * P)) throw new Error("density solve did not converge");
    return r1;
  }
  _envSat(q, by, x) {        // bubble (q = 0) or dew (q = 1) point from the phase envelope at P (by "P") or T (by "T")
    const br = this.env && (q === 0 ? this.env.bub : this.env.dew);
    if (!br) throw new Error("no phase envelope for this mixture");
    const o = by === "P" ? this._interp(br, "lnP", Math.log(x), ["T", "lnr"]) : this._interp(br, "T", x, ["lnP", "lnr"]);
    if (!o) throw new Error(`outside the phase envelope (${by} ${x})`);
    const TK = (by === "P" ? o.T : x) + 273.15, P = by === "P" ? x : Math.exp(o.lnP);
    const ph = q === 0 ? "liquid" : "gas";
    const st = this._dt(this._rhoAt(TK, P, Math.exp(o.lnr), ph), TK, ph);
    st.P = P; st.q = q; st.phase = q === 0 ? "saturated_liquid" : "saturated_vapour"; st.fallback = true;
    this.fallbacks++; return st;
  }
  _tqByPressure(q, T) {     // temperature-quality point from pressure-quality flashes (these converge closer to
    const br = this.env && (q === 0 ? this.env.bub : this.env.dew);   // the critical point), started on the envelope
    const o = br && this._interp(br, "T", T, ["lnP"]);
    if (!o) return null;
    try {
      const f = lnP => this._flash(this.pairs.PQ_INPUTS, Math.exp(lnP) * 1e6, q, q).T - T;
      let x0 = o.lnP, f0 = f(x0), x1 = x0 + 0.002, f1 = f(x1);
      for (let k = 0; k < 30 && Math.abs(f1) > 1e-7; k++) { const x2 = x1 - f1 * (x1 - x0) / (f1 - f0 || 1e-30); x0 = x1; f0 = f1; x1 = x2; f1 = f(x1); }
      if (!(Math.abs(f1) <= 1e-5)) return null;
      const st = Object.assign({}, this._flash(this.pairs.PQ_INPUTS, Math.exp(x1) * 1e6, q, q)); st.T = T; st.fallback = "PQ";
      return st;
    } catch (e) { return null; }
  }
  _dtSingle(phase, P, key, val) {   // single-phase state at P from (rho, T) evaluations
    const sat = this._sat(P), edge = phase === "gas" ? sat.dew : sat.bub;
    // the root lies below the dew-point density (gas, hotter than the dew point) or above the bubble-point
    // density (liquid, colder), where p rises monotonically with density: bracket it there and use Brent
    const rhoE = edge.rho || edge.D / this.MW;
    const atT = TC => {
      const TK = TC + 273.15, g = lnr => this._pDT(Math.exp(lnr), TK, phase) - P;
      let lo = Math.log(rhoE) - 0.005, hi = lo + 0.01;     // (a little either side, for states right at the edge)
      if (phase === "gas") lo -= 12; else for (let k = 0; k < 20 && !(g(hi) > 0); k++) hi += 0.05;
      return this._dt(Math.exp(SGHP.brent(g, lo, hi, 1e-13)), TK, phase);
    };
    let st;
    if (key === "T") st = atT(val);
    else {
      const lo = phase === "gas" ? edge.T : edge.T - 250, hi = phase === "gas" ? edge.T + 300 : edge.T;
      const TC = SGHP.brent(t => atT(t)[key] - val, lo, hi, 1e-9);
      st = atT(TC);
    }
    st.P = P; st.q = -1; st.phase = phase === "gas" ? "superheated_vapour" : "liquid"; st.fallback = true;
    this.fallbacks++; return st;
  }
  _make(names, x) {
    const M = this.M, fluids = names.join("&");
    let AS;
    try { AS = guarded(cpCall(() => M.factory("HEOS", fluids))); }
    catch (e) {
      const msg = cpError(e);
      if (names.length > 1 && /binary pair|interaction/i.test(msg)) {
        for (let i = 0; i < names.length; i++) for (let j = i + 1; j < names.length; j++) {
          try {
            const a = cpCall(() => M.get_fluid_param_string(names[i], "CAS")), b = cpCall(() => M.get_fluid_param_string(names[j], "CAS"));
            cpCall(() => M.apply_simple_mixing_rule(a, b, "linear"));
          } catch (_) {}
        }
        if (!this.notes.includes("estimated interaction")) this.notes.push("no fitted interaction parameters: linear mixing rule");
        AS = guarded(cpCall(() => M.factory("HEOS", fluids)));
      } else throw new Error(msg);
    }
    if (names.length > 1) {
      const set = this.basis === "mass" ? "set_mass_fractions" : "set_mole_fractions";
      try { AS[set](x); }
      catch (e) {
        if (M.VectorDouble) { const v = new M.VectorDouble(); x.forEach(z => v.push_back(z)); AS[set](v); v.delete(); }
        else throw new Error(cpError(e));
      }
    }
    return AS;
  }
  destroy() { for (const a of [this.AS, this.SAT]) try { a.delete(); } catch (_) {} }
  _dewT(P) {
    const k = P.toPrecision(10);
    if (this.dewCache.has(k)) return this.dewCache.get(k);
    let T;
    try { this.SAT.update(this.pairs.PQ_INPUTS, P * 1e6, 1); T = this.SAT.T() - 273.15; }
    catch (_) { try { T = this._envSat(1, "P", P).T; } catch (__) { T = Infinity; } }
    this.dewCache.set(k, T); return T;
  }
  _flash(pair, v1, v2, qIn) {
    const k = (pair && pair.value !== undefined ? pair.value : String(pair)) + "|" + v1 + "|" + v2;
    const hit = this.cache.get(k); if (hit) return hit;
    const AS = this.AS;
    try { AS.update(pair, v1, v2); } catch (e) { throw new Error(cpError(e)); }
    const st = { T: AS.T() - 273.15, P: AS.p() / 1e6, h: AS.hmass() / 1000, s: AS.smass() / 1000, D: AS.rhomass() };
    let q = qIn != null ? qIn : AS.Q();
    if (!(isFinite(st.T) && isFinite(st.h) && isFinite(st.s) && isFinite(st.D))) throw new Error("nonfinite state");
    if (q != null && isFinite(q) && q >= -1e-12 && q <= 1 + 1e-12) {
      st.q = Math.min(Math.max(q, 0), 1);
      st.phase = st.q <= 1e-12 ? "saturated_liquid" : st.q >= 1 - 1e-12 ? "saturated_vapour" : "two_phase";
    } else {
      st.q = -1;
      st.phase = isFinite(this.P_crit) && st.P >= this.P_crit ? "supercritical" : (st.T >= this._dewT(st.P) - 1e-6 ? "superheated_vapour" : "liquid");
    }
    if (this.cache.size > 50000) this.cache.clear();
    this.cache.set(k, st); return st;
  }
  /* Mixtures: CoolProp's general flash runs a phase-stability analysis (60-400 ms per call in
     WebAssembly). Bubble and dew points at P decide the phase, two-phase states come from
     pressure-quality flashes, and single-phase states from phase-imposed PropsSI calls (~1 ms). */
  _sat(P) {
    const k = P.toPrecision(12);
    let v = this.satCache.get(k);
    if (!v) { v = { bub: this.PQ(P, 0), dew: this.PQ(P, 1) }; this.satCache.set(k, v); }
    return v;
  }
  _single(phase, P, key, val) {
    try { return this._singleFast(phase, P, key, val); }
    catch (e) { if (!this.isMix || !this.env) throw e; return this._dtSingle(phase, P, key, val); }
  }
  _singleFast(phase, P, key, val) {
    const M = this.M, fl = this.psFluid, Pin = `P|${phase}`, Pa = P * 1e6;
    const k = `${phase}|${key}|${P}|${val}`;
    const hit = this.cache.get(k); if (hit) return hit;
    let T;
    try {
      if (key === "T") T = val;
      else {
        try { T = M.PropsSI("T", Pin, Pa, key === "h" ? "H" : "S", val * 1000, fl) - 273.15; if (!isFinite(T)) throw new Error("nonfinite state"); }
        catch (e) {
          // CoolProp bounds this inversion at the mixture's mole-averaged maximum temperature (for example
          // 270 C for propane/isopentane 20/80); hot discharge gas is past it, so solve T from the forward call
          if (phase !== "gas") throw e;
          T = this._gasInvert(P, key, val);
        }
      }
      const TK = T + 273.15;
      const st = { T, P, h: M.PropsSI("H", Pin, Pa, "T", TK, fl) / 1000, s: M.PropsSI("S", Pin, Pa, "T", TK, fl) / 1000,
        D: M.PropsSI("D", Pin, Pa, "T", TK, fl), q: -1, phase: phase === "gas" ? "superheated_vapour" : "liquid" };
      if (![st.T, st.h, st.s, st.D].every(isFinite)) throw new Error("nonfinite state");
      this.cache.set(k, st); return st;
    } catch (e) { throw new Error(cpError(e)); }
  }
  _gasInvert(P, key, val) {
    const Pa = P * 1e6, K = key === "h" ? "H" : "S";
    const f = T => this.M.PropsSI(K, `P|gas`, Pa, "T", T + 273.15, this.psFluid) / 1000 - val;
    const g = T => { try { return f(T); } catch (_) { return NaN; } };
    let lo = isFinite(this.T_crit) ? this.T_crit - 20 : 100, flo = g(lo);
    while (!isFinite(flo) && lo < 600) { lo += 10; flo = g(lo); }     // no gas root yet (compressed liquid)
    if (!isFinite(flo) || flo > 0) throw new Error("gas state not bracketed");
    let hi = lo + 50, fhi = g(hi);
    while (isFinite(fhi) && fhi < 0 && hi < 800) { lo = hi; flo = fhi; hi += 50; fhi = g(hi); }
    if (!(isFinite(fhi) && fhi >= 0)) throw new Error("gas state not bracketed");
    return SGHP.brent(f, lo, hi, 1e-9);
  }
  _twoPhase(P, key, val, sat) {
    const f = q => this.PQ(P, q)[key] - val;
    const q = SGHP.brent(f, 0, 1, 1e-9);
    return this.PQ(P, Math.min(Math.max(q, 0), 1));
  }
  _mixFlash(P, key, val, pairFallback, v1, v2) {
    if (isFinite(this.P_crit) && P >= this.P_crit * 0.995) {
      // supercritical (gas cooler): one fluid phase; the phase-imposed evaluation skips the slow stability analysis
      if (P > this.P_crit) {
        // one fluid phase above the critical pressure, but a phase-imposed solve can still land on the other
        // (metastable) density root: liquid-like below the critical temperature, gas-like above it. The
        // order follows the state at the critical temperature on this isobar, and a root on the wrong side
        // of the critical temperature is used only if the other fails.
        const order = this._scOrder(P, key, val), got = [];
        for (const ph of order) {
          try {
            const st = Object.assign({}, this._singleFast(ph, P, key, val)); st.phase = "supercritical";
            if (ph === "liquid" ? st.T <= this.T_crit + 2 : st.T >= this.T_crit - 2) return st;
            got.push(st);
          } catch (_) {}
        }
        if (got.length) return got[0];
      }
      return this._flash(pairFallback, v1, v2);
    }
    let sat;
    try { sat = this._sat(P); } catch (e) { return this._flash(pairFallback, v1, v2); }
    if (val >= sat.dew[key]) return val === sat.dew[key] ? sat.dew : this._single("gas", P, key, val);
    if (val <= sat.bub[key]) return val === sat.bub[key] ? sat.bub : this._single("liquid", P, key, val);
    return this._twoPhase(P, key, val, sat);
  }
  _scOrder(P, key, val) {
    if (key === "T") return val < this.T_crit ? ["liquid", "gas"] : ["gas", "liquid"];
    const k = `scref|${P}`;
    let ref = this.cache.get(k);
    if (ref === undefined) {
      ref = null;
      for (const ph of ["gas", "liquid"]) { try { ref = this._singleFast(ph, P, "T", this.T_crit); break; } catch (_) {} }
      this.cache.set(k, ref);
    }
    return ref && val < ref[key] ? ["liquid", "gas"] : ["gas", "liquid"];
  }
  PT(P, T) { return this.isMix ? this._mixFlash(P, "T", T, this.pairs.PT_INPUTS, P * 1e6, T + 273.15) : this._flash(this.pairs.PT_INPUTS, P * 1e6, T + 273.15); }
  PH(P, h) { return this.isMix ? this._mixFlash(P, "h", h, this.pairs.HmassP_INPUTS, h * 1e3, P * 1e6) : this._flash(this.pairs.HmassP_INPUTS, h * 1e3, P * 1e6); }
  PS(P, s) { return this.isMix ? this._mixFlash(P, "s", s, this.pairs.PSmass_INPUTS, P * 1e6, s * 1e3) : this._flash(this.pairs.PSmass_INPUTS, P * 1e6, s * 1e3); }
  /* temperature (C) and mass density (kg/m3): the critical isochore of a supercritical high side */
  TD(T, D) {
    if (!this.isMix) return this._flash(this.pairs.DmassT_INPUTS, D, T + 273.15);
    const k = `TD|${T}|${D}`, hit = this.cache.get(k); if (hit) return hit;
    const st = this._dt(D / this.MW, T + 273.15, "gas");
    st.q = -1; st.phase = isFinite(this.P_crit) && st.P >= this.P_crit ? "supercritical" : "superheated_vapour";
    this.cache.set(k, st); return st;
  }
  PQ(P, q) { return this._satFlash(this.pairs.PQ_INPUTS, P * 1e6, q, "P", P, q); }
  TQ(T, q) { return this._satFlash(this.pairs.QT_INPUTS, q, T + 273.15, "T", T, q); }
  _satFlash(pair, v1, v2, by, x, q) {
    try { return this._flash(pair, v1, v2, q); }
    catch (e) {
      if (!this.isMix || !this.env) throw e;
      const k = `env|${by}|${x}|${q}`, hit = this.cache.get(k); if (hit) return hit;
      let st;
      if (q === 0 || q === 1) {
        st = by === "T" ? this._tqByPressure(q, x) : null;
        if (!st) st = this._envSat(q, by, x);
      }
      else {                 // inside the dome: interpolate between the envelope's bubble and dew states
        const b = this._satFlash(by === "P" ? this.pairs.PQ_INPUTS : this.pairs.QT_INPUTS, by === "P" ? v1 : 0, by === "P" ? 0 : v2, by, x, 0);
        const d = this._satFlash(by === "P" ? this.pairs.PQ_INPUTS : this.pairs.QT_INPUTS, by === "P" ? v1 : 1, by === "P" ? 1 : v2, by, x, 1);
        st = { T: b.T + q * (d.T - b.T), P: b.P + q * (d.P - b.P), h: b.h + q * (d.h - b.h), s: b.s + q * (d.s - b.s),
          D: 1 / (1 / b.D + q * (1 / d.D - 1 / b.D)), q, phase: "two_phase", fallback: true };
      }
      this.cache.set(k, st); return st;
    }
  }
}


