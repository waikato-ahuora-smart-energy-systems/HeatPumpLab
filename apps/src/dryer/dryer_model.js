
/* Advanced dryer heat pumps, live: one design of RHP-A (Design A, RHX1 in series) or SLTHP (RHX1 parallel to the
   low-temperature gas cooler), solved in the browser. Ports of HeatPumpLab RecuperativeHP (Design A: RHX1 only,
   two stages) and SplitHeatPump.evaluate, over a property object with PT/PH/PS/PQ/TQ(...) -> {T C, P MPa, h kJ/kg,
   s kJ/kg/K, D kg/m3, phase} (the steam-generating heat pump page's CoolPropFluid) and an air object with h(T), T(h)
   at 0.101 MPa. Units: MPa, C, kJ/kg, per kg of refrigerant. */
(function (root) {
  "use strict";

  // ---------------------------------------------------------------- compressor (components/compressor.py)
  // The six HeatPumpLab isentropic (shaft) efficiency models, clipped to 0.4-1.0 as Compressor.get_isentropic_efficiency:
  // 0 ideal, 1 constant 0.7/0.96, 2 polynomial (Wang & Cleland 2018), 3 Brendel et al. (needs the suction pressure, MPa),
  // 5 Pareto front (the default, and the study's model). 4, the first Pareto front, was removed; old cases using it load as 5
  const PF = [0.8695, 0.025392, 1.476, 0.29915, 2.9599, 0.2];
  const COMP_MODELS = ["Ideal (1.0)", "Constant", "Polynomial (Wang & Cleland)", "Brendel et al.", null, "Pareto front"];
  let compModel = 5, etaConst = 0.7 / 0.96;
  function setCompModel(i, eta = null) {
    compModel = [0, 1, 2, 3, 5].includes(+i) ? +i : 5;
    if (eta != null && isFinite(+eta)) etaConst = Math.min(Math.max(+eta, 0.4), 1);
  }
  function etaIs(pr, Pin = null) { return etaModel(compModel, etaConst, pr, Pin); }
  /* isentropic efficiency of model i (COMP_MODELS) at a stage pressure ratio; eConst for the Constant model. The MVR
     fans of the HP + MVR cycle use this with their own model, separate from the compressors' */
  function etaModel(model, eConst, pr, Pin = null) {
    pr = Math.max(pr, 1);
    const L = Math.log(pr);
    let e;
    switch (+model) {
      case 0: e = 1; break;
      case 1: e = Math.min(Math.max(+eConst, 0.4), 1); break;
      case 2: e = (0.65 + 0.015 * pr - 0.0015 * pr * pr) / 0.96; break;
      case 3: e = pr > 1 ? 0.66981 - 0.6 / Math.pow(pr - 0.01466, 8.38 * (Pin ?? 0.2)) - 0.00102 * Math.pow(pr, 1.8) : 0.4; break;
      default: { const [e1, a, b, d, L0, w] = PF; e = e1 - a * Math.pow(L, b) - d / (1 + Math.exp(-(L - L0) / w)); }
    }
    return Math.min(Math.max(e, 0.4), 1);
  }
  /* vapour at (P, T) where T may sit on the dew point: a pure fluid has no PT state exactly at saturation, so a
     superheater outlet at the dew point is the saturated vapour */
  function ptGas(pr, P, T) {
    try { return pr.PT(P, T); }
    catch (e) { const d = pr.PQ(P, 1); if (Math.abs(d.T - T) < 0.05) return d; throw e; }
  }
  // Brendel et al. volumetric efficiency
  const etaVol = pr => 1 - 0.08244 * Math.pow(Math.max(pr, 1) - 1, 0.72773);
  const GAS = new Set(["saturated_vapour", "superheated_vapour", "supercritical"]);
  const MAX_STAGE_PR = 20;

  // ---------------------------------------------------------------- numerics
  function brent(f, a, b, tol = 1e-10, maxit = 100) {
    let fa = f(a), fb = f(b);
    if (fa * fb > 0) throw new Error("root not bracketed");
    if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }
    let c = a, fc = fa, d = b - a, mflag = true;
    for (let it = 0; it < maxit; it++) {
      if (fb === 0 || Math.abs(b - a) < tol) return b;
      let s;
      if (fa !== fc && fb !== fc) s = a * fb * fc / ((fa - fb) * (fa - fc)) + b * fa * fc / ((fb - fa) * (fb - fc)) + c * fa * fb / ((fc - fa) * (fc - fb));
      else s = b - fb * (b - a) / (fb - fa);
      const q = (3 * a + b) / 4;
      if (!((s > Math.min(q, b) && s < Math.max(q, b))) || (mflag && Math.abs(s - b) >= Math.abs(b - c) / 2) ||
          (!mflag && Math.abs(s - b) >= Math.abs(c - d) / 2) || (mflag && Math.abs(b - c) < tol) || (!mflag && Math.abs(c - d) < tol)) {
        s = (a + b) / 2; mflag = true;
      } else mflag = false;
      const fs = f(s); d = c; c = b; fc = fb;
      if (fa * fs < 0) { b = s; fb = fs; } else { a = s; fa = fs; }
      if (Math.abs(fa) < Math.abs(fb)) { [a, b] = [b, a]; [fa, fb] = [fb, fa]; }
    }
    return b;
  }
  // minimum of f on [a, b] (golden section, then the best sample)
  function goldenMin(f, a, b, iters = 18) {
    const g = (Math.sqrt(5) - 1) / 2;
    let x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = f(x1), f2 = f(x2);
    for (let i = 0; i < iters; i++) {
      if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = f(x1); }
      else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = f(x2); }
    }
    return f1 < f2 ? { x: x1, f: f1 } : { x: x2, f: f2 };
  }
  // monotone table y(x), linear interpolation
  function table(xs, ys) {
    const o = xs.map((x, i) => [x, ys[i]]).sort((p, q) => p[0] - q[0]);
    const X = o.map(p => p[0]), Y = o.map(p => p[1]);
    return x => {
      if (x <= X[0]) return Y[0];
      if (x >= X[X.length - 1]) return Y[Y.length - 1];
      let lo = 0, hi = X.length - 1;
      while (hi - lo > 1) { const m = (lo + hi) >> 1; if (X[m] <= x) lo = m; else hi = m; }
      return Y[lo] + (Y[hi] - Y[lo]) * (x - X[lo]) / (X[hi] - X[lo] || 1);
    };
  }
  const linspace = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
  const lmtd = (dh, dc) => (dh > 0 && dc > 0 ? (Math.abs(dh - dc) < 1e-6 ? dh : (dh - dc) / Math.log(dh / dc)) : NaN);

  // ---------------------------------------------------------------- shared defaults (the dryer study)
  const DEFAULTS = {
    T_air_inlet: 15, T_air_tar: 200, T_source_in: 40, Tref_evap_inlet: 10, dT_GC_min: 20, dTmin_rhx: 5,
    motor_eff: 0.93, P_air: 0.101, T_source_out: null,   // source water outlet (C); null: isothermal source, as the study
    // RHP-A (Design A)
    lift: 80, lift2: 0, split: 0.5, gc_out_offset: 0, ihx_on: true, P_H_guess: null,   // lift2: RHX2 on the interstage gas (RHP-B, RHP-C)
    // SLTHP
    rhx1_cold_end: 25, P_H: 2.0, alpha: 0.3, n_profile: 41,
  };

  /* air at 0.101 MPa from a property function: hOfT(T C) -> kJ/kg. Tables both ways over 0-260 C. */
  function airModel(hOfT) {
    const T = linspace(-5, 280, 286), h = T.map(hOfT);
    return { h: table(T, h), T: table(h, T) };
  }
  /* primary + secondary air as one stream per kg of primary air: both are heated from Tin to Tmt (the secondary
     leaves there), so below Tmt the enthalpy rises (1 + phi) times as fast; above it only the primary is heated.
     Any model that takes an air object then sees the pair; its air flows are kg of primary air. phi <= 0, or Tmt
     outside (Tin, Tmax), returns the air unchanged. */
  function airSecondary(air, Tin, Tmt, phi, Tmax = Infinity) {
    if (!(phi > 0 && Tmt > Tin && Tmt < Tmax)) return air;
    const h0 = air.h(Tin), hm = air.h(Tmt), f = 1 + phi;
    const comp = airModel(T => T <= Tmt ? h0 + f * (air.h(T) - h0) : h0 + f * (hm - h0) + air.h(T) - hm);
    comp.W = air.W; comp.secondary = { phi, Tmt, Tin };
    return comp;
  }
  // share of the heat to the air that goes to the secondary stream (inlet to Tmt), for air heated Tin -> Tout
  function secondaryShare(air, Tin, Tmt, Tout, phi) {
    if (!(phi > 0 && Tmt > Tin && Tmt < Tout)) return 0;
    const h0 = air.h(Tin), qs = phi * (air.h(Tmt) - h0);
    return qs / (qs + air.h(Tout) - h0);
  }

  // ======================================================================== SLTHP (SplitHeatPump.evaluate)
  function slthpOnce(pr, air, cfg) {
    const c = Object.assign({}, DEFAULTS, cfg);
    const fail = (message, extra = {}) => Object.assign({ success: false, message, cycle: "S" }, extra);
    const T_superheat = c.T_source_in - (c.dT_src ?? c.dTmin_rhx), T_LT_out = c.T_air_inlet + c.dT_GC_min, T_RX_out = T_superheat + c.rhx1_cold_end;
    const P_H = c.P_H, lift = c.lift;
    let alpha = c.alpha;
    try {
      if (lift <= 0) alpha = 0;
      else if (!(alpha > 0 && alpha < 1)) return fail("RHX1 branch fraction outside (0, 1)", { P_H });
      c._scLow = scLow(pr, T_superheat);
      const lt = pr.PT(P_H, T_LT_out), rx = alpha > 0 ? pr.PT(P_H, T_RX_out) : lt;
      const h_mix = (1 - alpha) * lt.h + alpha * rx.h;
      const P_L = c.T_source_out == null ? evapPressureForInlet(pr, h_mix, P_H, c.Tref_evap_inlet, c.expander)
        : waterPressure(pr, c, P => expandTo(pr, P_H, h_mix, P, c.expander), c._scLow ? Math.min(P_H / SC_MIN_PR, c.P_L_max ?? Infinity) : Math.min(pr.TQ(T_superheat, 1).P, 0.98 * P_H), c.PL_guess).P_L;
      const h_ev = expandTo(pr, P_H, h_mix, P_L, c.expander), wx = h_mix - h_ev;
      const dew = c._scLow ? dewOrNull(pr, P_L) : pr.PQ(P_L, 1);   // supercritical low side: no dew point above P_crit
      // the study's optima can sit on this bound; within 0.05 K the superheater outlet is the dew point
      if (dew && dew.T > T_superheat + 0.05) return fail(`evaporator dew ${dew.T.toFixed(1)} C above superheater outlet ${T_superheat.toFixed(1)} C`, { P_H, P_L });
      const sup = dew && dew.T >= T_superheat ? dew : ptGas(pr, P_L, T_superheat), T_suc = T_superheat + lift;
      const suction = lift > 0 ? pr.PT(P_L, T_suc) : sup;
      if (!GAS.has(suction.phase)) return fail("suction not gas", { P_H, P_L });
      const comp = compressPowerShare(pr, suction, P_L, P_H, c.split);
      for (const k of ["out1", "out2"]) if (!GAS.has(comp[k].phase)) return fail(`wet compression (stage ${k === "out1" ? 1 : 2} outlet)`, { P_H, P_L });
      const h_dis = comp.out2.h, T_dis = comp.out2.T, q_rhx = suction.h - sup.h, w = comp.w1 + comp.w2;
      let h_split = null;
      if (alpha > 0) {
        h_split = rx.h + q_rhx / alpha;
        if (h_split > h_dis) return fail("RHX1 duty exceeds the branch (split above discharge)", { P_H, P_L, alpha_needed: q_rhx / (h_dis - rx.h) });
        if (h_split < lt.h) return fail("split below the low-temperature gas-cooler outlet", { P_H, P_L });
      }
      const q_air = h_dis - h_mix - q_rhx;
      const ha_in = air.h(c.T_air_inlet), ha_tar = air.h(c.T_air_tar), r = q_air / (ha_tar - ha_in);
      const q_LT = (1 - alpha) * ((h_split ?? h_dis) - lt.h), q_HT = q_air - q_LT, ha_mid = ha_in + q_LT / r;
      // exchanger profiles (tables along the isobars, as the Python model)
      const n = c.n_profile, h_lo = Math.min(lt.h, rx.h);
      const hTab = linspace(h_lo, h_dis, 90), hp = table(hTab, hTab.map(h => pr.PH(P_H, h).T));
      const fs = linspace(0, 1, n), minApp = (hot, cold) => Math.min(...fs.map(f => hot(f) - cold(f)));
      const h_top = h_split ?? h_dis, ap = {};
      ap.LT_gas_cooler = minApp(f => hp(lt.h + f * (h_top - lt.h)), f => air.T(ha_in + f * (ha_mid - ha_in)));
      let rlmtd = null;
      if (h_split != null) {
        ap.HT_gas_cooler = minApp(f => hp(h_split + f * (h_dis - h_split)), f => air.T(ha_mid + f * (ha_tar - ha_mid)));
        const hl = linspace(sup.h, suction.h, 30), lp = table(hl, hl.map(h => pr.PH(P_L, h).T));
        ap.RHX1 = minApp(f => hp(rx.h + f * (h_split - rx.h)), f => lp(sup.h + f * (suction.h - sup.h)));
        rlmtd = lmtd(hp(h_split) - T_suc, T_RX_out - T_superheat);
      }
      const margins = [ap.LT_gas_cooler - c.dT_GC_min];
      if (ap.HT_gas_cooler != null) margins.push(ap.HT_gas_cooler - c.dT_GC_min, ap.RHX1 - c.dTmin_rhx);
      const margin = Math.min(...margins);
      const T_split = h_split != null ? hp(h_split) : null;
      const states = {
        LT_gc_out: lt, RHX1_hot_out: alpha > 0 ? rx : null, valve_in: pr.PH(P_H, h_mix), evap_in: pr.PH(P_L, h_ev),
        evap_dew: dew, superheat_out: sup, suction, stage1_out: comp.out1, discharge: comp.out2, split: h_split != null ? pr.PH(P_H, h_split) : null,
      };
      // the REFPROP optimum sits on its approach constraints; with CoolProp a shortfall of a few hundredths of a
      // kelvin is a property difference, so it is reported, not failed
      const warnings = [];
      if (margin < -0.3) return fail(`approach violated by ${(-margin).toFixed(2)} K`, { P_H, P_L, approaches_K: ap, states, margin });
      if (margin < -1e-6) warnings.push(`approach short of the minimum by ${(-margin).toFixed(2)} K`);
      const src = sourceSide(pr, c, P_L, h_ev, sup.h, warnings);
      const xq = c.expander ? expWarn(pr, P_L, h_ev, warnings) : null;
      const v1 = 1 / suction.D, v2 = c.split > 0 ? 1 / comp.out1.D : 0;   // one compressor: no stage-2 suction
      const w_el = expEl(w, wx, c.motor_eff), cop = q_air / w_el;
      const pr1 = comp.P_I / P_L, pr2 = P_H / comp.P_I;
      return {
        success: true, message: "success", cycle: "S", warnings, COP: cop, VHC: q_air / (v1 + v2) / 1000, P_H, P_I: comp.P_I, P_L, states, expander: !!c.expander,
        report: {
          COP: cop, VHC: q_air / (v1 + v2) / 1000, Pmax_bar: P_H * 10, P_evap_bar: P_L * 10, P_int_bar: comp.P_I * 10,
          PR_1: pr1, PR_2: pr2, Overall_PR: P_H / P_L, eta_is_1: comp.eta1, eta_is_2: comp.eta2, n_vol_1: etaVol(pr1), n_vol_2: etaVol(pr2),
          T_evap_inlet: c.Tref_evap_inlet, T_dew_evap: dew ? dew.T : null, low_side_supercritical: !!c._scLow, T_superheat_out: T_superheat, T_suction: T_suc,
          Comp_1_T_out: comp.out1.T, T_comp_discharge: T_dis, T_split, T_LT_gc_out: T_LT_out, T_RHX1_hot_out: alpha > 0 ? T_RX_out : null,
          T_valve_in: states.valve_in.T, T_air_mid: air.T(ha_mid), alpha_RHX1: alpha,
          q_air_kJkg: q_air, q_rhx1_kJkg: q_rhx, w_shaft_kJkg: w, w_el_kJkg: w_el, q_evap_kJkg: sup.h - h_ev,
          expander: !!c.expander, eta_exp: c.expander ? c.expander.eta : null, w_exp_kJkg: wx, exp_out_quality: xq,
          q_HT_gc_kJkg: q_HT, q_LT_gc_kJkg: q_LT, air_per_ref_kgkg: r, rhx1_LMTD_K: rlmtd, min_approach_margin_K: margin,
          approaches_K: ap, v_suction_1: v1, v_suction_2: v2, rho_suction_1: suction.D, rho_suction_2: c.split > 0 ? comp.out1.D : null, power_frac_2: comp.w2 / (comp.w1 + comp.w2), ...src,
        },
      };
    } catch (e) { return fail(`property calculation failed: ${e && e.message ? e.message : e}`, { P_H }); }
  }

  /* Evaporator pressure from a finite source-water flow (water T_source_in -> T_source_out, counterflow, constant
     cp): the highest P_L whose evaporator and superheater stay at least dT_src (default dTmin_rhx) below the water line. hvOf(P) is
     the valve outlet enthalpy at P_L = P; the superheater outlet is T_source_in - dTmin_rhx, so the hot end sits on
     the minimum by construction and the cold end and the glide set the pressure. Two-phase points by quality. */
  /* source temperature at a fraction f of the evaporator duty from its cold end: water is linear; an air source
     passes source_curve(f) (humid air, condensing below its dew point) */
  function srcT(c, f) { return typeof c.source_curve === "function" ? c.source_curve(f) : c.T_source_out + f * (c.T_source_in - c.T_source_out); }
  /* Supercritical low side. When the heater outlet (source inlet less the source approach) is at or above the
     critical temperature there is no dew point to bound the low-side pressure: the "evaporator" is a gas heater along
     an isobar (above P_crit, or two-phase then gas below it), and its pressure is set by the source pinch alone, as
     the evaporating pressure is otherwise. scLow says whether a design takes this branch; subcritical designs keep
     the dew-point bound unchanged. */
  const SC_MIN_PR = 1.2;           // supercritical low side: the low-side pressure stays at least this ratio below P_H
  const SC_P_MAX = 30;             // and the gas-cooler pressure search stops at 300 bar (15 MPa / 150 bar subcritical)
  function scLow(pr, T_out) {
    if (isFinite(pr.T_crit) && T_out >= pr.T_crit) return true;
    try { pr.TQ(T_out, 1); return false; }
    catch (e) { if (isFinite(pr.T_crit) && T_out > pr.T_crit - 20) return true; throw e; }
  }
  // dew point at P, or null above the critical pressure (or where there is none)
  function dewOrNull(pr, P) {
    if (isFinite(pr.P_crit) && P >= pr.P_crit) return null;
    try { return pr.PQ(P, 1); } catch (e) { return null; }
  }
  // smallest source-minus-refrigerant difference along the low-side isobar from the valve outlet hv towards the
  // heater outlet, sampled in enthalpy (two-phase and supercritical alike) with a golden-section refinement at the
  // smallest. The outlet itself sits on the approach by construction, so (as the evaporator's dew-point samples) the
  // last 4 % of the duty is left out: otherwise the residual is flat at zero below the pressure the pinch sets.
  function isobarPinch(pr, c, P, hv, need = -Infinity) {
    const sup = ptGas(pr, P, c.T_source_in - (c.dT_src ?? c.dTmin_rhx)), span = sup.h - hv;
    const d = f => srcT(c, f) - pr.PH(P, hv + f * span).T;
    const fs = linspace(0, 0.96, 17), ds = fs.map(d);
    const i = ds.indexOf(Math.min(...ds));
    if (i === 0 || ds[i] > need + 2) return ds[i];
    return Math.min(ds[i], goldenMin(d, fs[i - 1], fs[Math.min(i + 1, fs.length - 1)], 9).f);
  }
  function waterPinch(pr, c, P, hv) {
    if (c._scLow) return isobarPinch(pr, c, P, hv, c.dT_src ?? c.dTmin_rhx);
    const Twi = c.T_source_in;
    const sup = ptGas(pr, P, Twi - (c.dT_src ?? c.dTmin_rhx)), bub = pr.PQ(P, 0), dew = pr.PQ(P, 1), span = sup.h - hv;
    const Tw = h => srcT(c, (h - hv) / span);
    let min = Tw(hv) - pr.PH(P, hv).T;
    const q0 = Math.min(Math.max((hv - bub.h) / (dew.h - bub.h), 0), 1);
    for (let i = 1; i <= 10; i++) { const st = pr.PQ(P, Math.min(1, q0 + (1 - q0) * i / 10)); min = Math.min(min, Tw(st.h) - st.T); }
    return min;
  }
  function waterPressure(pr, c, hvOf, P_top, PL_guess) {
    const resid = lnP => { const P = Math.exp(lnP); return waterPinch(pr, c, P, hvOf(P)) - (c.dT_src ?? c.dTmin_rhx); };
    const top = Math.log(P_top * (1 - 1e-7)), rTop = resid(top);
    if (rTop >= 0) return { P_L: Math.exp(top), capped: true };
    let b = top, rb = rTop, a = Math.log(PL_guess && PL_guess < P_top ? PL_guess : P_top * 0.5), ra = resid(a);
    if (ra < 0) { let st = 0.02; for (let k = 0; k < 40 && ra < 0; k++) { b = a; rb = ra; a -= st; st *= 2; ra = resid(a); } }
    else { let st = 0.02; for (let k = 0; k < 40; k++) { const x = Math.min(a + st, top); const rx = x === top ? rTop : resid(x); if (rx < 0) { b = x; rb = rx; break; } a = x; ra = rx; st *= 2; } }
    if (ra < 0 || rb > 0) throw new Error(`${c._scLow ? "gas-heater" : "evaporator"} pressure not bracketed against the ${c.source_name || "source water"}`);
    return { P_L: Math.exp(brent(resid, a, b, 1e-7)), capped: false };
  }

  /* Source water against the evaporator and superheater (counterflow, constant water cp): the smallest
     water-minus-refrigerant temperature difference, water per kg of refrigerant, and a warning below dTmin_rhx */
  function sourceSide(pr, c, P_L, h_in, h_out, warnings) {
    if (c.T_source_out == null) return { evap_pinch_K: null, water_per_ref_kgkg: null };
    const Twi = c.T_source_in, Two = c.T_source_out, n = 16, nm = c.source_name || "source water";
    let min = Infinity;
    for (let i = 0; i <= n; i++) {
      const f = i / n, h = h_in + f * (h_out - h_in);
      let T; try { T = pr.PH(P_L, h).T; } catch (e) { continue; }
      min = Math.min(min, srcT(c, f) - T);
    }
    if (min < (c.dT_src ?? c.dTmin_rhx) - 0.05) warnings.push(`evaporator pinch to the ${nm} ${min.toFixed(1)} K, below the ${(c.dT_src ?? c.dTmin_rhx)} K minimum: raise the ${nm} outlet temperature or lower the evaporator inlet`);
    return { evap_pinch_K: min, water_per_ref_kgkg: c.source_curve ? null : (h_out - h_in) / (4.18 * (Twi - Two)) };
  }

  /* SLTHP at a given RHX1 flow share; a share just too small for the RHX1 branch to carry its duty (a rounded
     study optimum on that bound) is raised to the smallest that does, up to 5 %, and reported */
  function slthp(pr, air, cfg) {
    let res = slthpOnce(pr, air, cfg);
    const a0 = cfg.alpha;
    if (res.success || !res.alpha_needed || !(a0 > 0) || res.alpha_needed > a0 * 1.05) return res;
    let a = a0;
    for (let k = 0; k < 8 && !res.success && res.alpha_needed; k++) {
      a = Math.max(res.alpha_needed * (1 + 1e-7), a * (1 + 1e-7));
      if (a > a0 * 1.05 || a >= 1) break;
      res = slthpOnce(pr, air, Object.assign({}, cfg, { alpha: a }));
    }
    if (res.success) res.warnings.unshift(`RHX1 flow share raised from ${(a0 * 100).toFixed(2)} % to ${(a * 100).toFixed(2)} %, the least at which the RHX1 branch carries its duty`);
    return res;
  }

  // P_L such that the valve outlet (isenthalpic, h_valve) sits at T_target (geometric bisection, as the Python model)
  function evapPressureForInlet(pr, h_valve, P_H, T_target, ex = null) {
    let lo = 0.02, hi = 0.98 * P_H;
    const T_at = P => pr.PH(P, expandTo(pr, P_H, h_valve, P, ex)).T;
    const tl = T_at(lo), th = T_at(hi);
    if (!(tl < T_target && T_target < th)) throw new Error(`evaporator inlet ${T_target} C not bracketed (${tl.toFixed(1)} .. ${th.toFixed(1)} C)`);
    for (let k = 0; k < 50; k++) {
      const mid = Math.sqrt(lo * hi);
      if (T_at(mid) < T_target) lo = mid; else hi = mid;
      if (hi / lo - 1 < 2e-6) break;
    }
    return Math.sqrt(lo * hi);
  }
  // two stages with the stage-2 shaft-power share (SLTHP)
  function compressPowerShare(pr, suction, P_L, P_H, share) {
    const stage = (P_in, s_in, h_in, P_out) => {
      const eta = etaIs(P_out / P_in, P_in), h_is = pr.PS(P_out, s_in).h, w = (h_is - h_in) / eta;
      return { w, eta, out: pr.PH(P_out, h_in + w) };
    };
    if (!(share > 0)) { const s = stage(P_L, suction.s, suction.h, P_H); return { P_I: P_H, w1: s.w, w2: 0, eta1: s.eta, eta2: null, out1: s.out, out2: s.out }; }
    let lo = Math.log(P_L * 1.0001), hi = Math.log(P_H / 1.0001);
    const evaluate = lnP => {
      const P_I = Math.exp(lnP), s1 = stage(P_L, suction.s, suction.h, P_I), s2 = stage(P_I, s1.out.s, s1.out.h, P_H);
      return [s2.w / (s1.w + s2.w), { P_I, w1: s1.w, w2: s2.w, eta1: s1.eta, eta2: s2.eta, out1: s1.out, out2: s2.out }];
    };
    const [flo] = evaluate(lo), [fhi] = evaluate(hi);
    if (!(fhi <= share && share <= flo)) throw new Error("compressor power split not reachable");
    let pack;
    for (let k = 0; k < 30; k++) {
      const mid = 0.5 * (lo + hi), [frac, p] = evaluate(mid); pack = p;
      if (frac > share) lo = mid; else hi = mid;
      if (Math.abs(frac - share) < 2e-4) break;
    }
    return pack;
  }

  /* Expansion from the high side: an isenthalpic valve, or with c.expander = {eta} an expander (turbine) of that
     isentropic efficiency, whose shaft work is returned as electricity at the motor efficiency (expEl). A flashing
     liquid leaves an expander two-phase: expWarn says so, as two-phase expanders are less proven than gas turbines. */
  function expandTo(pr, P_in, h_in, P_out, ex) {
    if (!ex || !(ex.eta > 0)) return h_in;
    pr = pr._base || pr;           // full property calls, not a fast temperature table
    const m = pr._expMemo || (pr._expMemo = new Map()), k = `${ex.eta}|${P_in}|${h_in}|${P_out}`;
    let h = m.get(k);
    if (h === undefined) { h = h_in - ex.eta * (h_in - pr.PS(P_out, pr.PH(P_in, h_in).s).h); if (m.size > 20000) m.clear(); m.set(k, h); }
    return h;
  }
  const expEl = (w, wx, motor) => w / motor - wx * motor;
  function expWarn(pr, P_out, h_out, warnings, name = "expander") {
    let o; try { o = (pr._base || pr).PH(P_out, h_out); } catch (e) { return null; }
    const wet = o.phase === "two_phase" || (o.q > 0 && o.q < 1);
    if (wet) warnings.push(`two-phase expansion: the ${name} outlet is ${(o.q * 100).toFixed(0)} % vapour (${o.T.toFixed(1)} C); two-phase expanders are less proven than gas turbines and their efficiency is uncertain`);
    return wet ? o.q : null;
  }

  // ======================================================================== RHP-A (Design A)
  /* Design A: evaporator -> superheater (to T_source - dTmin_rhx) -> [IHX cold side] -> RHX1 cold side (+lift)
     -> stage 1 -> stage 2 (no intercooling; ln(P_I/P_L) = (1 - split) ln(P_H/P_L)) -> gas cooler 2 (refrigerant to
     suction + dTmin_rhx) -> RHX1 hot side -> gas cooler 1 (to T_air_inlet + approach + offset) -> [IHX hot side]
     -> valve. The evaporator pressure puts the valve outlet at Tref_evap_inlet (the dew point may not exceed the
     superheater outlet); the gas-cooler pressure is the lowest at which the gas-cooler pinch reaches the approach
     (minus the 0.05 K search band of the Python solver). */
  const PINCH_BAND = 0.05 - 1e-5;
  function rhpaAtPressure(pr, air, c, P_H, PL_guess) {
    const T_sup = c.T_source_in - (c.dT_src ?? c.dTmin_rhx), T_gc_out = c.T_air_inlet + c.dT_GC_min + c.gc_out_offset;
    const ihx = !!c.ihx_on && c.T_source_in < T_gc_out;
    const lift = Math.max(c.lift, 0), lift2 = Math.max(c.lift2 || 0, 0), F2 = c.split;
    const gco = pr.PT(P_H, T_gc_out);
    // the valve inlet depends on P_L through the IHX duty
    const ihxAt = P_L => {
      const sup = ptGas(pr, P_L, T_sup);
      if (!ihx) return { sup, q: 0, hv: gco.h, he: expandTo(pr, P_H, gco.h, P_L, c.expander) };
      const hot = gco.h - pr.PT(P_H, sup.T + c.dTmin_rhx).h, cold = ptGas(pr, P_L, T_gc_out - c.dTmin_rhx).h - sup.h;
      const q = Math.max(Math.min(hot, cold), 0);
      return { sup, q, hv: gco.h - q, he: expandTo(pr, P_H, gco.h - q, P_L, c.expander) };
    };
    if (c._scLow == null) c._scLow = scLow(pr, T_sup);
    // dew point at the superheater outlet: upper bound on P_L; a supercritical low side is bounded below P_H instead
    const P_dew_max = c._scLow ? Math.min(P_H / SC_MIN_PR, c.P_L_max ?? Infinity) : pr.TQ(T_sup, 1).P;
    const resid = lnP => { const P = Math.exp(lnP), s = ihxAt(P); return pr.PH(P, s.he).T - c.Tref_evap_inlet; };
    let P_L, capped = false;
    const rHi = c.T_source_out == null ? resid(Math.log(P_dew_max * (1 - 1e-7))) : 0;
    if (c.T_source_out != null) ({ P_L, capped } = waterPressure(pr, c, P => ihxAt(P).he, P_dew_max, PL_guess));
    else if (rHi <= 0) { P_L = P_dew_max * (1 - 1e-7); capped = true; }
    else {
      // bracket from the previous P_L (small steps), else from half the dew-point pressure
      const top = Math.log(P_dew_max * (1 - 1e-7));
      let a = Math.log(PL_guess && PL_guess < P_dew_max ? PL_guess : P_dew_max * 0.5), ra = resid(a), b = top, rb = rHi;
      if (ra > 0) { let st = 0.01; for (let k = 0; k < 40 && ra > 0; k++) { b = a; rb = ra; a -= st; st *= 2; ra = resid(a); } }
      else { let st = 0.01; for (let k = 0; k < 40; k++) { const x = Math.min(a + st, top); const rx = x === top ? rHi : resid(x); if (rx > 0) { b = x; rb = rx; break; } a = x; ra = rx; st *= 2; } }
      if (ra > 0 || rb < 0) throw new Error("evaporator pressure not bracketed");
      P_L = Math.exp(brent(resid, a, b, 1e-8));
    }
    const ih = ihxAt(P_L), sup = ih.sup, dew = c._scLow ? dewOrNull(pr, P_L) : pr.PQ(P_L, 1);
    const ihxOut = ih.q > 0 ? pr.PH(P_L, sup.h + ih.q) : sup;           // IHX cold side outlet = RHX1 cold inlet
    const suction = lift > 0 ? pr.PT(P_L, ihxOut.T + lift) : ihxOut;
    const P_I = P_L * Math.exp((1 - F2) * Math.log(P_H / P_L));
    const stage = (P_in, s_in, h_in, P_out) => {
      const eta = etaIs(P_out / P_in, P_in), h_is = pr.PS(P_out, s_in).h, w = (h_is - h_in) / eta;
      return { w, eta, out: pr.PH(P_out, h_in + w), pr: P_out / P_in };
    };
    const s1 = stage(P_L, suction.s, suction.h, P_I);
    // RHX2 cold side heats the interstage gas by lift2 before stage 2 (RHP-B, RHP-C)
    const s2in = lift2 > 0 ? pr.PT(P_I, s1.out.T + lift2) : s1.out, q_rhx2 = s2in.h - s1.out.h;
    const s2 = stage(P_I, s2in.s, s2in.h, P_H);
    const dis = s2.out;
    // hot side, as hp.gas_coolers: gas cooler 3 to the stage-2 suction + dTmin_rhx, RHX2 hot side, gas cooler 2 to
    // the stage-1 suction + dTmin_rhx, RHX1 hot side, gas cooler 1 to T_gc_out
    let h = dis.h, T = dis.T, rx2In = null, rx2Out = null, rxIn = null, rxOut = null, q_rhx = 0;
    const gcs = [];
    if (lift2 > 0) {
      rx2In = pr.PT(P_H, s2in.T + c.dTmin_rhx); rx2Out = pr.PH(P_H, rx2In.h - q_rhx2);
      gcs.push({ name: "Gas cooler 3", hin: h, hout: rx2In.h, Tin: T, Tout: rx2In.T }); h = rx2Out.h; T = rx2Out.T;
    }
    if (lift > 0) {
      rxIn = pr.PT(P_H, suction.T + c.dTmin_rhx);
      q_rhx = suction.h - ihxOut.h;
      rxOut = pr.PH(P_H, rxIn.h - q_rhx);
      gcs.push({ name: "Gas cooler 2", hin: h, hout: rxIn.h, Tin: T, Tout: rxIn.T }); h = rxOut.h; T = rxOut.T;
    }
    gcs.push({ name: "Gas cooler 1", hin: h, hout: gco.h, Tin: T, Tout: gco.T });
    const q_air = gcs.reduce((s, g) => s + (g.hin - g.hout), 0);
    const ha_in = air.h(c.T_air_inlet), ha_out = air.h(c.T_air_tar), r = q_air / (ha_out - ha_in);
    let ha = ha_out;
    for (const g of gcs) { g.q = g.hin - g.hout; g.ha_out = ha; g.ha_in = ha - g.q / r; ha = g.ha_in; g.Ta_out = air.T(g.ha_out); g.Ta_in = air.T(g.ha_in); }
    return { P_H, P_L, P_I, capped, ihx, T_gc_out, gco, sup, dew, ihxOut, ihxQ: ih.q, hValve: ih.hv, hEv: ih.he, suction, s1, s2in, s2, dis, rxIn, rxOut, q_rhx, rx2In, rx2Out, q_rhx2, gcs, q_air, r, ha_in, ha_out };
  }
  // minimum refrigerant-air approach in each gas cooler: terminals and the internal minimum
  function pinches(pr, air, st) {
    return st.gcs.map(g => {
      if (!(g.q > 0)) return { hot: g.Tin - g.Ta_out, cold: g.Tout - g.Ta_in, internal: Infinity, at: 0, min: Math.min(g.Tin - g.Ta_out, g.Tout - g.Ta_in) };
      const d = f => pr.PH(st.P_H, g.hout + f * g.q).T - air.T(g.ha_in + f * g.q / st.r);
      // coarse samples then a golden-section refinement around the smallest
      const fs = linspace(0.001, 0.999, 11), ds = fs.map(d);
      let i = ds.indexOf(Math.min(...ds));
      const m = goldenMin(d, fs[Math.max(i - 1, 0)], fs[Math.min(i + 1, fs.length - 1)], 12);
      const hot = g.Tin - g.Ta_out, cold = g.Tout - g.Ta_in;
      return { hot, cold, internal: m.f, at: m.x, min: Math.min(hot, cold, m.f) };
    });
  }

  /* RHP cycles (RHP-A/B/C and the conventional cycle). A subcritical low side evaporates at the highest pressure
     the source pinch allows (rhpaSolve). A supercritical low side (heater outlet at or above the critical temperature)
     has no dew point to bound it, and the pinch alone can sit far above the best design (the gas-cooler pressure then
     climbs to its limit as P_L rises), so its pressure is chosen for COP: a cap on P_L, which the source pinch may
     lower further. From a previous design (P_L and P_H guesses) the search is local: P_L and 5 % either side, moved
     until the middle is best, then the parabola's vertex. Otherwise ln P_L is scanned over 0.4-2.5 P_crit (6 points)
     and refined by golden section. Designs over the gas-cooler pressure limit or the stage pressure ratio fail and so
     bound the search. The COP is flat near its optimum, so P_L is found to a few per cent. */
  function rhpa(pr, air, cfg) {
    let sc = false;
    if (cfg.P_L_max == null) {
      const T_out = (cfg.T_source_in ?? DEFAULTS.T_source_in) - (cfg.dT_src ?? cfg.dTmin_rhx ?? DEFAULTS.dTmin_rhx);
      try { sc = isFinite(pr.P_crit) && scLow(pr, T_out); } catch (e) { sc = false; }
    }
    if (!sc) return rhpaSolve(pr, air, cfg);
    const P_cap = cfg.P_max ?? SC_P_MAX, lo = Math.log(0.4 * pr.P_crit), hi = Math.log(Math.min(2.5 * pr.P_crit, P_cap / SC_MIN_PR));
    let PH = cfg.P_H_guess, PLg = cfg.PL_guess, n = 0, best = null, lastFail = null;
    const memo = new Map();
    let pinched = null;                   // a design whose P_L the source pinch set: every larger cap gives it again
    const f = lnP => {
      if (memo.has(lnP)) return memo.get(lnP);
      if (pinched && lnP >= pinched.ln) return pinched.COP;
      n++;
      const r = rhpaSolve(pr, air, Object.assign({}, cfg, { P_L_max: Math.exp(lnP), P_H_guess: PH, PL_guess: PLg }));
      if (r.success) {
        PH = r.P_H; PLg = r.P_L; if (!best || r.COP > best.COP) best = r;
        if (r.report.P_L_by === "source pinch" && (!pinched || Math.log(r.P_L) < pinched.ln)) pinched = { ln: Math.log(r.P_L), COP: r.COP };
      } else lastFail = r;
      const v = r.success ? r.COP : -Infinity; memo.set(lnP, v); return v;
    };
    const done = () => { best.report.P_L_by = best.report.P_L_by === "source pinch" ? "source pinch" : "COP"; best.report.P_L_search_solves = n; return best; };
    if (cfg.PL_guess > 0 && cfg.P_H_guess > 0 && Math.log(cfg.PL_guess) > lo && Math.log(cfg.PL_guess) < hi) {
      const h = 0.05;
      let x0 = Math.log(cfg.PL_guess), f0 = f(x0), fm = f(x0 - h), fp = f(x0 + h);
      for (let k = 0; k < 8 && isFinite(f0) && (fm > f0 || fp > f0); k++) {
        if (fp > f0 && fp >= fm) { fm = f0; x0 += h; f0 = fp; fp = f(x0 + h); }
        else { fp = f0; x0 -= h; f0 = fm; fm = f(x0 - h); }
      }
      if (isFinite(f0) && f0 >= fm && f0 >= fp) {
        if (isFinite(fm) && isFinite(fp)) {
          const den = fm - 2 * f0 + fp;
          if (den < 0) f(x0 + Math.min(Math.max(h * (fm - fp) / (2 * den), -h), h));
        } else {                          // a limit on one side: golden section across the cliff
          let a = x0 - h, b = x0 + h; const g = (Math.sqrt(5) - 1) / 2;
          let x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = f(x1), f2 = f(x2);
          for (let k = 0; k < 4; k++) {
            if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = f(x1); }
            else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = f(x2); }
          }
        }
        return done();
      }
    }
    const xs = linspace(lo, hi, 6), fs = xs.map(f);
    const i = fs.indexOf(Math.max(...fs));
    if (!best) return Object.assign(lastFail, { message: `${lastFail.message} (supercritical low side: no low-side pressure from ${(Math.exp(lo) * 10).toFixed(0)} to ${(Math.exp(hi) * 10).toFixed(0)} bar solves)` });
    // golden section for the largest COP between the scan neighbours of the best point
    let a = xs[Math.max(i - 1, 0)], b = xs[Math.min(i + 1, xs.length - 1)];
    const g = (Math.sqrt(5) - 1) / 2;
    let x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = f(x1), f2 = f(x2);
    for (let k = 0; k < 6; k++) {
      if (f1 > f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = f(x1); }
      else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = f(x2); }
    }
    return done();
  }
  function rhpaSolve(pr, air, cfg) {
    const c = Object.assign({}, DEFAULTS, cfg);
    const fail = (message, extra = {}) => Object.assign({ success: false, message, cycle: "A" }, extra);
    if ((c.T_air_tar - c.T_air_inlet) * 2 / 3 < c.lift + (c.lift2 || 0)) return fail(`RHX lifts ${c.lift + (c.lift2 || 0)} K are above two thirds of the air temperature rise (${((c.T_air_tar - c.T_air_inlet) * 2 / 3).toFixed(0)} K)`);
    if (!(c.gc_out_offset >= 0 && c.gc_out_offset <= 50)) return fail("gas-cooler outlet offset must be 0-50 K");
    try {
      c._scLow = scLow(pr, c.T_source_in - (c.dT_src ?? c.dTmin_rhx));
      const P_cap = c.P_max ?? (c._scLow ? SC_P_MAX : 15);
      let PL = c.PL_guess || null;
      const at = P => { const s = rhpaAtPressure(pr, air, c, P, PL); PL = s.P_L; s.pinch = pinches(pr, air, s); s.resid = Math.min(...s.pinch.map(p => p.min)) - c.dT_GC_min;
        // a discharge colder than the gas-cooler outlet (a gas cooler with negative duty) cannot heat the air: short of the approach
        if (s.gcs.some(g => g.q < -1e-9)) s.resid = Math.min(s.resid, -c.dT_GC_min - 50);
        return s; };
      const target = -PINCH_BAND;
      // bracket the band entry around the guess (pressure steps relative to the suction pressure, as the Python search)
      // a supercritical low side starts above the critical pressure, and steps at least 5 % in P_H: its P_L can sit
      // close below P_H, where steps relative to the suction pressure would crawl
      const sc = c._scLow && isFinite(pr.P_crit);
      // a guess from an earlier design may sit above the limit (a stale page guess, or a limit since lowered): start inside it
      let P0 = Math.min(c.P_H_guess || (sc ? 1.5 * pr.P_crit : 2.0), P_cap / 1.1), s0;
      try { s0 = at(P0); } catch (e) { P0 = 3.0; s0 = at(P0); }
      let lo = { P: P0, s: s0 }, hi = { P: P0, s: s0 };
      const step = 1.1;
      for (let k = 0; k < 60; k++) {
        if (lo.s.resid < target && hi.s.resid < target) {
          // still short of the approach: the pressure must rise, so a stage already above the PR limit cannot recover
          const f = Math.min(Math.max(c.split || 0, 0), 1);
          const prMax = c.max_stage_pr ?? MAX_STAGE_PR;
          if (Math.pow(hi.P / hi.s.P_L, Math.max(f, 1 - f)) > prMax) return fail(`stage pressure ratio above ${prMax}`);
          let Pn = hi.s.P_L + (hi.P - hi.s.P_L) * step;
          if (sc) Pn = Math.max(Pn, hi.P * 1.05);
          lo = hi; hi = { P: Pn, s: at(Pn) };
        } else if (lo.s.resid > target && hi.s.resid > target) {
          let Pn = lo.s.P_L + (lo.P - lo.s.P_L) / step;
          if (sc) Pn = Math.min(Pn, lo.P / 1.05);
          if (!(Pn > lo.s.P_L * 1.0001)) break;
          hi = lo; lo = { P: Pn, s: at(Pn) };
        } else break;
        if (hi.P > P_cap) break;
      }
      if (!(lo.s.resid <= target && target <= hi.s.resid)) return fail(`no gas-cooler pressure reaches the ${c.dT_GC_min} K approach (pinch ${(Math.max(lo.s.resid, hi.s.resid) + c.dT_GC_min).toFixed(1)} K at ${(hi.P * 10).toFixed(1)} bar)`);
      const lnP = brent(x => at(Math.exp(x)).resid - target, Math.log(lo.P), Math.log(hi.P), 1e-5, 60);
      const st = at(Math.exp(lnP));
      if (st.P_H > P_cap * (1 + 1e-9)) return fail(`gas-cooler pressure ${(st.P_H * 10).toFixed(0)} bar above the ${(P_cap * 10).toFixed(0)} bar limit`, { P_H: st.P_H, P_L: st.P_L });
      return finishRHPA(pr, air, c, st);
    } catch (e) { return fail(`property calculation failed: ${e && e.message ? e.message : e}`); }
  }

  function finishRHPA(pr, air, c, st) {
    const fail = (message, extra = {}) => Object.assign({ success: false, message, cycle: "A", P_H: st.P_H, P_L: st.P_L }, extra);
    const s1 = st.s1, s2 = st.s2;
    if (!GAS.has(st.suction.phase)) return fail("stage-1 suction is not gas: raise the RHX1 lift or the gas-cooler outlet offset");
    if (!GAS.has(s1.out.phase)) return fail("wet compression at the stage-1 outlet: raise the gas-cooler outlet offset (the study used 7-40 K for such designs)");
    if (!GAS.has(s2.out.phase)) return fail("wet compression at the stage-2 outlet");
    const prMax = c.max_stage_pr ?? MAX_STAGE_PR;
    if (s1.pr > prMax || s2.pr > prMax) return fail(`stage pressure ratio above ${prMax}`, { stage_pr: Math.max(s1.pr, s2.pr) });
    if (st.rxIn) {
      if (st.rxIn.T > st.dis.T + 1e-6) return fail("RHX1 hot inlet above the compressor discharge: gas cooler 2 would run backwards");
      if (!(st.rxOut.T > st.ihxOut.T)) return fail("RHX1 temperature cross");
    }
    if (st.rx2In) {
      if (st.rx2In.T > st.dis.T + 1e-6) return fail("RHX2 hot inlet above the compressor discharge: gas cooler 3 would run backwards");
      if (!(st.rx2Out.T > s1.out.T)) return fail("RHX2 temperature cross");
    }
    for (const g of st.gcs) if (g.q < -1e-6) return fail(`${g.name} duty is negative`);
    const wx = st.hValve - st.hEv, w = s1.w + s2.w, w_el = expEl(w, wx, c.motor_eff), cop = st.q_air / w_el;
    const v1 = 1 / st.suction.D, v2 = 1 / st.s2in.D;
    const minP = Math.min(...st.pinch.map(p => p.min));
    const warnings = [], src = sourceSide(pr, c, st.P_L, st.hEv, st.sup.h, warnings);
    const xq = c.expander ? expWarn(pr, st.P_L, st.hEv, warnings) : null;
    const lowSC = !!c._scLow, aboveCrit = lowSC && isFinite(pr.P_crit) && st.P_L >= pr.P_crit;
    if (lowSC && st.capped && c.P_L_max == null) warnings.push(`the source pinch does not bound the gas-heater pressure: held at P_H / ${SC_MIN_PR}`);
    const states = {
      evap_in: pr.PH(st.P_L, st.hEv), evap_dew: st.dew, superheat_out: st.sup, ihx_cold_out: st.ihx ? st.ihxOut : null,
      suction: st.suction, stage1_out: s1.out, stage2_in: st.s2in, discharge: st.dis, gc3_out: st.rx2In, rhx2_hot_out: st.rx2Out,
      gc2_out: st.rxIn, rhx1_hot_out: st.rxOut, gc1_out: st.gco,
      valve_in: pr.PH(st.P_H, st.hValve),
    };
    return {
      success: true, message: st.capped && !lowSC ? "success (evaporator outlet at the superheater temperature)" : "success", cycle: "A", warnings,
      COP: cop, VHC: st.q_air / (v1 + v2) / 1000, P_H: st.P_H, P_I: st.P_I, P_L: st.P_L, states, gcs: st.gcs, low_sc: lowSC, expander: !!c.expander,
      report: {
        COP: cop, VHC: st.q_air / (v1 + v2) / 1000, Pmax_bar: st.P_H * 10, P_evap_bar: st.P_L * 10, P_int_bar: st.P_I * 10,
        PR_1: s1.pr, PR_2: s2.pr, Overall_PR: st.P_H / st.P_L, eta_is_1: s1.eta, eta_is_2: s2.eta, n_vol_1: etaVol(s1.pr), n_vol_2: etaVol(s2.pr),
        // supercritical low side: no dew point above P_crit (null), and the "superheater outlet" is the gas-heater outlet
        low_side_supercritical: lowSC, low_side_above_Pcrit: aboveCrit, P_L_by: lowSC ? (st.capped ? "COP" : "source pinch") : "source pinch", P_crit_bar: isFinite(pr.P_crit) ? pr.P_crit * 10 : null, evap_in_quality: states.evap_in.q,
        T_evap_inlet: states.evap_in.T, T_dew_evap: st.dew ? st.dew.T : null, Tref_evap_outlet: st.dew ? st.dew.T : null, T_superheat_out: st.sup.T,
        T_ihx_out: st.ihx ? st.ihxOut.T : null, q_ihx_kJkg: st.ihxQ, ihx_on: st.ihx,
        T_suction: st.suction.T, Comp_1_T_out: s1.out.T, T_comp_discharge: st.dis.T, T_gc2_out: st.rxIn ? st.rxIn.T : null,
        T_rhx1_hot_out: st.rxOut ? st.rxOut.T : null, T_gc_out: st.T_gc_out, T_valve_in: states.valve_in.T,
        T_air_mid: st.gcs.length > 1 ? st.gcs[0].Ta_in : null, q_air_kJkg: st.q_air, q_rhx1_kJkg: st.q_rhx,
        w_shaft_kJkg: w, w_el_kJkg: w_el, q_evap_kJkg: st.sup.h - st.hEv, air_per_ref_kgkg: st.r,
        expander: !!c.expander, eta_exp: c.expander ? c.expander.eta : null, w_exp_kJkg: wx, exp_out_quality: xq,
        gc_minimum_approach: minP, gc_pinch: st.pinch, lift: c.lift, lift2: c.lift2 || 0, q_rhx2_kJkg: st.q_rhx2, split: c.split, gc_out_offset: c.gc_out_offset,
        v_suction_1: v1, v_suction_2: v2, rho_suction_1: st.suction.D, rho_suction_2: st.s2in.D,
        q_gc_kJkg: st.gcs.map(g => g.q), power_frac_2: s2.w / (s1.w + s2.w), ...src,
      },
    };
  }

  /* Diagram paths for a solved design: [{name, pts: [[T, P bar, h, s]]}], named as the study's cycle diagrams. */
  function paths(pr, res, n = 12) {
    const st = res.states, P_H = res.P_H, P_L = res.P_L, P_I = res.P_I;
    const pt = s => [s.T, s.P * 10, s.h, s.s];
    const iso = (P, a, b, k = n) => Array.from({ length: k + 1 }, (_, i) => { try { return pt(pr.PH(P, a + (b - a) * i / k)); } catch (e) { return null; } }).filter(Boolean);
    const comp = (P0, h0, P1, h1, k = 10) => Array.from({ length: k + 1 }, (_, i) => { const f = i / k; try { return pt(pr.PH(P0 * Math.pow(P1 / P0, f), h0 + f * (h1 - h0))); } catch (e) { return null; } }).filter(Boolean);
    const out = [];
    const hv = st.valve_in.h, he = st.evap_in.h;
    out.push({ name: res.expander ? "Expander" : "Expansion valve", pts: comp(P_H, hv, P_L, he) });
    out.push({ name: "Evaporator", pts: iso(P_L, he, st.superheat_out.h, 2 * n) });
    if (res.cycle === "A") {
      if (st.ihx_cold_out) out.push({ name: "Internal heat exchanger", pts: iso(P_L, st.superheat_out.h, st.ihx_cold_out.h) });
      const cin = st.ihx_cold_out || st.superheat_out;
      if (st.suction.h - cin.h > 1e-9) out.push({ name: "Recuperator", pts: iso(P_L, cin.h, st.suction.h) });
      out.push({ name: "Compressor", pts: comp(P_L, st.suction.h, P_I, st.stage1_out.h) });
      const s2in = st.stage2_in || st.stage1_out;
      if (s2in.h - st.stage1_out.h > 1e-9) out.push({ name: "Recuperator 2", pts: iso(P_I, st.stage1_out.h, s2in.h) });
      out.push({ name: "Compressor", pts: comp(P_I, s2in.h, P_H, st.discharge.h) });
      // hot side: [GC3, RHX2] [GC2, RHX1] GC1
      let h = st.discharge.h;
      if (st.gc3_out) {
        out.push({ name: "Gas cooler", pts: iso(P_H, h, st.gc3_out.h) });
        out.push({ name: "Recuperator 2", pts: iso(P_H, st.gc3_out.h, st.rhx2_hot_out.h) }); h = st.rhx2_hot_out.h;
      }
      if (st.gc2_out) {
        out.push({ name: "Gas cooler", pts: iso(P_H, h, st.gc2_out.h) });
        out.push({ name: "Recuperator", pts: iso(P_H, st.gc2_out.h, st.rhx1_hot_out.h) }); h = st.rhx1_hot_out.h;
      }
      out.push({ name: "Gas cooler", pts: iso(P_H, h, st.gc1_out.h) });
      if (st.ihx_cold_out) out.push({ name: "Internal heat exchanger", pts: iso(P_H, st.gc1_out.h, hv) });
    } else {
      out.push({ name: "RHX1 (cold)", pts: iso(P_L, st.superheat_out.h, st.suction.h) });
      out.push({ name: "Compressor 1", pts: comp(P_L, st.suction.h, P_I, st.stage1_out.h) });
      out.push({ name: "Compressor 2", pts: comp(P_I, st.stage1_out.h, P_H, st.discharge.h) });
      if (st.split) {
        out.push({ name: "HT gas cooler", pts: iso(P_H, st.discharge.h, st.split.h) });
        out.push({ name: "LT gas cooler", pts: iso(P_H, st.split.h, st.LT_gc_out.h) });
        out.push({ name: "RHX1 (hot)", pts: iso(P_H, st.split.h, st.RHX1_hot_out.h) });
        out.push({ name: "Mixing", pts: [pt(st.LT_gc_out), pt(st.valve_in), pt(st.RHX1_hot_out)] });
      } else out.push({ name: "LT gas cooler", pts: iso(P_H, st.discharge.h, st.LT_gc_out.h) });
    }
    return out.filter(p => p.pts.length > 1);
  }

  const api = { DEFAULTS, COMP_MODELS, setCompModel, getCompModel: () => compModel, getEtaConst: () => etaConst, etaIs, etaModel, etaVol, brent, goldenMin, table, airModel, airSecondary, secondaryShare, slthp, rhpa, rhpaAtPressure, pinches, paths, scLow, evapPressureForInlet, compressPowerShare, expandTo, expEl, expWarn };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.DRYER = api;
})(typeof self !== "undefined" ? self : this);

