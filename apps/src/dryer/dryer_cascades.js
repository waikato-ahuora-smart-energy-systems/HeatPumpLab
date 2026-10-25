
/* Advanced dryer heat pumps: multi-loop cycles, solved in the browser over CoolPropFluid property objects
   (PT/PH/PS/PQ/TQ -> {T C, P MPa, h kJ/kg, s kJ/kg/K, D kg/m3, phase}) and an air object {h(T), T(h)}.

   kong      Transcritical-transcritical CO2 cascade after Kong (2026, PhD thesis, Waikato), cycle TT-4 with the
             additional gas cooler (Fig. 58): the bottom cycle (IHX) rejects all its heat in the cascade gas cooler (CGC),
             superheating the top suction; the top cycle (IHX) heats the dryer air in the HT gas cooler and all the air
             (dryer + fluidised bed) in the MT gas cooler; both evaporators on the source water. The bottom cycle may run
             subcritical.
   cascade   n-level cascade after Hoang et al.: level 1 evaporates on the source water, level k+1 in level k's cascade
             condenser. Air passes the levels in series, 1 to n, each level's subcooler then its condenser. A lower
             level can subcool its liquid into the air and/or send a fraction of its desuperheating and condensing heat
             to the air (the rest drives the next level); the top level desuperheats, condenses and subcools into the
             air, at the lowest condensing temperature that keeps the air approach.
   booster   HP + MVR: a conventional heat pump bottom cycle (ammonia by default) whose condenser heats the air first
             and raises steam (sub-atmospheric when below 100 C); n mechanical-vapour-recompression stages at a fixed
             pressure ratio each (1.5) compress the steam, with water injection to just above saturation between
             stages, and the steam desuperheats and condenses into the air to the target; the condensate flashes back.

   Every result: {success, message, warnings, COP, VHC, report, loops: [{name, fluid, flow, paths}], hx: [{name, key,
   hot, cold, otherName}], q_air, w_shaft, q_source, ref: {rho, flow, PR}} with duties per kg of the reference flow
   (the top cycle for kong, level 1 for cascade, the bottom cycle for booster). */
(function (root) {
  "use strict";
  const D_ = (typeof module !== "undefined" && module.exports) ? require("./dryer_model.js") : root.DRYER;
  const { etaIs, etaModel, brent, expandTo, expEl, expWarn } = D_;
  const MOTOR = 0.93;

  // ------------------------------------------------------------------ helpers
  const lin = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
  function stage(pr, Pin, st, Pout, eta = null) {
    if (eta == null && stage.fixedEta) eta = stage.fixedEta;
    const e = eta ?? etaIs(Pout / Pin, Pin), his = pr.PS(Pout, st.s).h, w = (his - st.h) / e;
    return { w, eta: e, out: pr.PH(Pout, st.h + w), pr: Pout / Pin };
  }
  const isGas = st => st.phase === "superheated_vapour" || st.phase === "supercritical" || st.phase === "saturated_vapour";
  const pt = s => [s.T, s.P * 10, s.h, s.s];
  function iso(pr, P, h0, h1, n = 14) { const o = []; for (let i = 0; i <= n; i++) { try { o.push(pt(pr.PH(P, h0 + (h1 - h0) * i / n))); } catch (e) {} } return o; }
  function compPath(pr, P0, h0, P1, h1, n = 8) { const o = []; for (let i = 0; i <= n; i++) { const f = i / n; try { o.push(pt(pr.PH(P0 * Math.pow(P1 / P0, f), h0 + f * (h1 - h0)))); } catch (e) {} } return o; }
  // a stream cooled along an isobar from h_hot to h_cold, against air heated from Ta_in to Ta_out (counterflow);
  // per kg of that stream: [[q, T_ref, T_air]] from the cold end
  function airProfile(pr, air, P, hHot, hCold, TaIn, TaOut, n = 16) {
    const ha0 = air.h(TaIn), ha1 = air.h(TaOut), Q = hHot - hCold;
    return lin(0, 1, n).map(f => [f * Q, Tph(pr, P, hCold + f * Q), air.T(ha0 + f * (ha1 - ha0))]);
  }
  /* T(h) along an isobar from a table built once per (fluid, pressure): PT on a 0.5 K grid (fast, and dense in h
     where cp peaks near the pseudo-critical point), the saturation points added below the critical pressure, then
     linear in h. Used only where a model enables it (pr._fastT), for the many temperature lookups of pinch searches. */
  const ISO = new WeakMap();
  function isoTable(pr, P) {
    let m = ISO.get(pr); if (!m) { m = new Map(); ISO.set(pr, m); }
    const k = P.toFixed(7); let t = m.get(k);
    if (t) return t;
    const pts = [];
    for (let T = -40; T <= 450; T += 0.5) { try { const st = pr.PT(P, T); if (isFinite(st.h)) pts.push([st.h, T]); } catch (e) {} }
    if (isFinite(pr.P_crit) && P < pr.P_crit) for (const q of [0, 1]) { try { const st = pr.PQ(P, q); pts.push([st.h, st.T]); } catch (e) {} }
    pts.sort((a, b) => a[0] - b[0]);
    const H = [], T = [];
    for (const [h, tt] of pts) if (!H.length || h > H[H.length - 1] + 1e-9) { H.push(h); T.push(tt); }
    t = { H, T, ok: H.length > 20 };
    if (m.size > 400) m.clear();
    m.set(k, t); return t;
  }
  function fastT(pr, P, h) {
    const t = isoTable(pr._base || pr, P), H = t.H;
    if (!t.ok || h < H[0] || h > H[H.length - 1]) return pr.PH(P, h).T;
    let lo = 0, hi = H.length - 1;
    while (hi - lo > 1) { const md = (lo + hi) >> 1; if (H[md] <= h) lo = md; else hi = md; }
    return t.T[lo] + (t.T[hi] - t.T[lo]) * (h - H[lo]) / (H[hi] - H[lo]);
  }
  const Tph = (pr, P, h) => pr._fastT ? fastT(pr, P, h) : pr.PH(P, h).T;
  const fastFluid = pr => { if (pr._fastT) return pr; const f = Object.create(pr); f._fastT = true; f._base = pr; return f; };
  const minGap = prof => Math.min(...prof.map(p => p[1] - p[2]));
  /* outlet enthalpy of a gas cooler or condenser heating air from TaIn to TaOut: the coldest outlet that keeps the
     approach dT (the air flow follows from the duty) */
  function coolerOutlet(pr, air, P, hIn, TaIn, TaOut, dT) {
    const Tin = Tph(pr, P, hIn);
    if (Tin < TaOut + dT - 1e-9) throw new Error(`stream at ${Tin.toFixed(1)} C cannot heat air to ${TaOut.toFixed(1)} C with a ${dT} K approach`);
    const hLow = pr.PT(P, TaIn + dT).h;
    const g = h => minGap(airProfile(pr, air, P, hIn, h, TaIn, TaOut, 14)) - dT;
    if (g(hLow) >= -1e-6) return hLow;
    let lo = hLow, hi = hIn - 1e-6 * Math.abs(hIn);
    if (g(hi) < -0.02) throw new Error("no outlet keeps the air approach");
    if (g(hi) < 0) return hi;
    // near the pseudo-critical point the pinch moves inside and the feasible outlets need not be one interval:
    // scan up from the coldest outlet for the first feasible one, then refine below it
    const N = 40;
    for (let i = 1; i <= N; i++) { const h = hLow + (hi - hLow) * i / N; if (g(h) >= 0) { lo = hLow + (hi - hLow) * (i - 1) / N; hi = h; break; } }
    for (let k = 0; k < 40 && hi - lo > 0.01; k++) { const m = (lo + hi) / 2; if (g(m) >= 0) hi = m; else lo = m; }
    return hi;
  }
  /* counterflow heat exchanger, hot flow mh (kg per kg of the cold flow): the largest duty per kg cold that keeps dT */
  function maxDuty(prH, PH, hHin, mh, prC, PC, hCin, dT, n = 14) {
    const ThIn = Tph(prH, PH, hHin), TcIn = Tph(prC, PC, hCin);
    if (ThIn - TcIn < dT) return { q: 0, gap: ThIn - TcIn };
    const gap = q => { let m = Infinity; for (const f of lin(0, 1, n)) { const Tc = Tph(prC, PC, hCin + f * q), Th = Tph(prH, PH, hHin - q / mh + f * q / mh); m = Math.min(m, Th - Tc); } return m - dT; };
    let qmax = Math.min(mh * (hHin - prH.PT(PH, TcIn + dT).h), prC.PT(PC, ThIn - dT).h - hCin);
    if (!(qmax > 0)) return { q: 0, gap: ThIn - TcIn };
    if (gap(qmax) >= 0) return { q: qmax };
    // the pinch can sit inside (cp peaks near the pseudo-critical point), so the feasible duties need not be one
    // interval: scan down from qmax for the largest feasible duty, then refine above it
    let lo = 0, hi = qmax;
    for (let i = 39; i >= 0; i--) { const q = qmax * i / 40; if (q === 0 || gap(q) >= 0) { lo = q; hi = qmax * (i + 1) / 40; break; } }
    for (let k = 0; k < 40 && hi - lo > 1e-4 * qmax; k++) { const m = (lo + hi) / 2; if (gap(m) >= 0) lo = m; else hi = m; }
    return { q: lo };
  }
  function hxPairs(prH, PH, hHin, hHout, mh, prC, PC, hCin, hCout, mc, n = 14) {
    // [[q, T]] per kg basis for both sides from the cold end
    const Q = mc * (hCout - hCin);
    const hot = lin(0, 1, n).map(f => [f * Q, prH.PH(PH, hHout + f * (hHin - hHout)).T]);
    const cold = lin(0, 1, n).map(f => [f * Q, prC.PH(PC, hCin + f * (hCout - hCin)).T]);
    return { hot, cold };
  }
  // source temperature at a fraction f of the evaporator duty from its cold end (water linear, or an air source curve)
  const srcT = (c, f) => typeof c.source_curve === "function" ? c.source_curve(f) : (c.T_source_out ?? c.T_source_in) + f * (c.T_source_in - (c.T_source_out ?? c.T_source_in));
  function waterLine(Q, c, q0 = 0, f0 = 0, f1 = 1) {
    if (typeof c.source_curve !== "function") return [[q0, srcT(c, f0)], [q0 + Q, srcT(c, f1)]];
    return lin(0, 1, 16).map(t => [q0 + t * Q, srcT(c, f0 + t * (f1 - f0))]);
  }
  const fail = (message, extra = {}) => Object.assign({ success: false, message }, extra);

  // ================================================================== Kong et al. TT cascade
  /* Bottom cycle: evaporator -> IHX -> compressor -> cascade gas cooler (CGC, its only heat sink: it superheats the
     top suction) -> IHX -> valve. Top cycle: evaporator -> IHX -> CGC -> compressor -> HT gas cooler -> MT gas cooler
     -> IHX -> valve. Air: primary and secondary air pass the MT gas cooler together (inlet -> T_mt), the secondary
     leaves at T_mt and the primary goes on through the HT gas cooler to the target. The two gas coolers are one
     counterflow exchanger against that composite air curve; its outlet is the coldest that keeps the approach, and
     the air flow follows from the duty. air_ratio = secondary / primary air mass flow; r = bottom / top refrigerant
     mass flow. Kong's reference (Table 9, section 6.2.4): eta_s 0.70 x 0.96 mechanical, 150 bar, 20 K CO2-air and
     5 K CO2-CO2 approaches, evaporating 10 / 20 C, a third of the air to the fluidised bed (air_ratio 0.5): COP 2.55. */
  const KONG = { T_ev_b: 10, T_ev_t: 20, P_b: 14.62, P_t: 15.0, r: 0.437, T_mt: 102, dT_cgc: 5, dT_ihx_b: 5, dT_ihx_t: 5, dT_ev: 2, eta: 0.7, air_ratio: 0.5, dT_ht: null, gc_off: 0 };
  // composite air curve per kg primary air: [0, H1] both streams inlet -> T_mt, [H1, H1 + H2] primary only T_mt -> target
  function airComposite(air, Tin, Tmt, Tout, phi) {
    const ha = air.h(Tin), hm = air.h(Tmt), ho = air.h(Tout), H1 = (1 + phi) * (hm - ha), H2 = ho - hm, H = H1 + H2;
    const T = x => x <= H1 ? air.T(ha + x / (1 + phi)) : air.T(hm + (x - H1));
    return { H1, H2, H, T, Tf: f => T(f * H), f1: H1 / H };
  }
  // coldest gas-cooler outlet against a composite air curve that keeps the approach (scan, then bisect)
  // dT applies where both air streams are heated (the MT gas cooler), dT2 where only the primary is (HT gas cooler)
  function compositeOutlet(pr, P, hIn, comp, TaIn, dT, dT2 = dT, n = 32) {
    const gap = h => {
      let m = Infinity;
      for (let i = 0; i <= n; i++) { const f = i / n; m = Math.min(m, Tph(pr, P, h + f * (hIn - h)) - comp.Tf(f) - (f <= comp.f1 ? dT : dT2)); }
      if (comp.f1 < 1) m = Math.min(m, Tph(pr, P, h + comp.f1 * (hIn - h)) - comp.T(comp.H1) - Math.min(dT, dT2));
      return m;
    };
    let hLow; try { hLow = pr.PT(P, TaIn + dT).h; } catch (e) { hLow = hIn - 1; }
    if (gap(hLow) >= -1e-6) return { h: hLow, ok: true };
    let lo = hLow, hi = hIn - 1e-6 * Math.abs(hIn);
    if (gap(hi) < 0) return { h: hi, ok: false };
    const N = 40;
    for (let i = 1; i <= N; i++) { const h = hLow + (hi - hLow) * i / N; if (gap(h) >= 0) { lo = hLow + (hi - hLow) * (i - 1) / N; hi = h; break; } }
    for (let k = 0; k < 40 && hi - lo > 0.01; k++) { const m = (lo + hi) / 2; if (gap(m) >= 0) hi = m; else lo = m; }
    return { h: hi, ok: true };
  }
  function kong(fb, ft, air, cfg) {
    const c = Object.assign({ T_air_inlet: 15, T_air_tar: 200, dT_GC_min: 20, T_source_in: 40, T_source_out: 20 }, KONG, cfg);
    fb = fastFluid(fb); ft = fastFluid(ft);
    const warnings = [];
    stage.fixedEta = c.eta || null;
    try {
      const r = c.r, Tain = c.T_air_inlet, phi = Math.max(0, +c.air_ratio || 0);
      if (!(r > 0)) return fail("the bottom/top mass flow must be positive");
      if (!(c.T_mt > Tain && c.T_mt < c.T_air_tar)) return fail(`the secondary air target must lie between the air inlet (${Tain} C) and the dryer inlet target (${c.T_air_tar} C)`);
      const comp = airComposite(air, Tain, c.T_mt, c.T_air_tar, phi);
      const eb = fb.TQ(c.T_ev_b, 1), et = ft.TQ(c.T_ev_t, 1), Peb = eb.P, Pet = et.P;
      if (c.P_t <= Pet || c.P_b <= Peb) return fail("discharge pressure below the evaporating pressure");
      // IHX: cold side from the evaporator outlet, hot side from Thot (same flow); the largest duty both ends allow
      const ihx = (pr, Pc, cin, Ph, Thot, dT) => {
        if (Thot - dT <= cin.T) return 0;
        return Math.max(0, Math.min(pr.PT(Ph, Thot).h - pr.PT(Ph, cin.T + dT).h, pr.PT(Pc, Thot - dT).h - cin.h));
      };
      /* bottom cycle for a given top CGC inlet: its CGC outlet temperature Tb drives its own IHX, which raises its
         discharge and so the CGC outlet again. The physical state is the lowest self-consistent Tb (scan up from
         no IHX, then bisect); none means the bottom cycle brings more heat than the top suction can take. */
      const bottomAt = (ctIn, Tb) => {
        const qib = ihx(fb, Peb, eb, c.P_b, Tb, c.dT_ihx_b), cbIn = fb.PH(Peb, eb.h + qib), cb = stage(fb, Peb, cbIn, c.P_b);
        const cg = maxDuty(fb, c.P_b, cb.out.h, r, ft, Pet, ctIn.h, c.dT_cgc), hbo = cb.out.h - cg.q / r;
        return { qib, cbIn, cb, cg, hbo, f: Tph(fb, c.P_b, hbo) - Tb };
      };
      const solveBottom = ctIn => {
        let lo = eb.T + c.dT_ihx_b, A = bottomAt(ctIn, lo);
        if (A.f <= 0) return A;
        let hi = null;
        for (let T = lo + 5; T <= 400; T += 5) { const B = bottomAt(ctIn, T); if (B.f <= 0) { hi = T; break; } lo = T; A = B; }
        if (hi == null) return null;
        for (let k = 0; k < 30 && hi - lo > 0.005; k++) { const m = (lo + hi) / 2, M = bottomAt(ctIn, m); if (M.f > 0) { lo = m; A = M; } else hi = m; }
        return bottomAt(ctIn, hi);
      };
      /* top gas-cooler outlet: the coldest that keeps the approaches, raised by gc_off (K) when that pays through
         the top IHX (Kong: the outlet temperature is an optimisation variable, at least air inlet + approach) */
      const dTht = c.dT_ht ?? c.dT_GC_min;
      const gcOut = ctp => {
        const g = compositeOutlet(ft, c.P_t, ctp.out.h, comp, Tain, c.dT_GC_min, dTht);
        if (!g.ok || !(c.gc_off > 0)) return g;
        const T = Tph(ft, c.P_t, g.h) + c.gc_off;
        if (T >= ctp.out.T - 1) return { h: g.h, ok: false };
        return { h: ft.PT(c.P_t, T).h, ok: true };
      };
      // outer: the top gas-cooler outlet sets the top IHX duty
      let Tg = c._T0 ? c._T0[0] : Tain + c.dT_GC_min, S = null, conv = false;
      for (let it = 0; it < 40; it++) {
        const qit = ihx(ft, Pet, et, c.P_t, Tg, c.dT_ihx_t), ctIn = ft.PH(Pet, et.h + qit);
        const B = solveBottom(ctIn);
        if (!B) return fail("the bottom cycle brings more heat than the cascade gas cooler can pass to the top suction (lower the bottom/top mass flow or the bottom pressure)");
        const ts = ft.PH(Pet, ctIn.h + B.cg.q), ctp = stage(ft, Pet, ts, c.P_t);
        const go = gcOut(ctp);
        if (!go.ok) return fail(`the top discharge (${ctp.out.T.toFixed(0)} C) cannot heat the primary air to ${c.T_air_tar} C with a ${dTht} K approach (raise the top pressure or the bottom/top mass flow)`);
        const nTg = Tph(ft, c.P_t, go.h);
        S = { ctIn, cbIn: B.cbIn, cb: B.cb, cg: B.cg, ts, hbo: B.hbo, ctp, go };
        const done = Math.abs(nTg - Tg) < 0.01;
        Tg += (it < 3 ? 1 : 0.5) * (nTg - Tg);
        if (done) { conv = true; break; }
      }
      if (!conv) warnings.push("the top gas-cooler outlet did not settle to 0.01 K");
      const Tb = Tph(fb, c.P_b, S.hbo);
      const { ctIn, cbIn, cb, cg, ts, hbo, ctp, go } = S;
      if (!isGas(cbIn) || !isGas(ts)) return fail("suction not gas");
      // IHX with the converged outlets
      const hOut = go.h, Tout = Tph(ft, c.P_t, hOut);
      const qit = ihx(ft, Pet, et, c.P_t, Tout, c.dT_ihx_t), qib = ihx(fb, Peb, eb, c.P_b, Tph(fb, c.P_b, hbo), c.dT_ihx_b);
      const hvt = hOut - qit, hvb = hbo - qib;
      // valves, or expanders: the evaporator inlets
      const het = expandTo(ft, c.P_t, hvt, Pet, c.expander), heb = expandTo(fb, c.P_b, hvb, Peb, c.expander), wx = (hvt - het) + r * (hvb - heb);
      const qAir = ctp.out.h - hOut, mP = qAir / comp.H, hMid = hOut + mP * comp.H1;
      const qMT = mP * comp.H1, qHT = mP * comp.H2;
      const qPri = mP * (air.h(c.T_air_tar) - air.h(Tain)), qSec = phi * mP * (air.h(c.T_mt) - air.h(Tain));
      const w = ctp.w + r * cb.w, wEl = expEl(w, wx, c.motor_eff ?? MOTOR), cop = qAir / wEl;
      const qEt = et.h - het, qEb = r * (eb.h - heb);
      const xq = c.expander ? [expWarn(ft, Pet, het, warnings, "top expander"), expWarn(fb, Peb, heb, warnings, "bottom expander")] : null;
      if (cg.q <= 1e-6) warnings.push("the cascade gas cooler carries no heat: the bottom discharge is too cold for the top suction");
      if (Math.max(ctp.out.T, cb.out.T) > 320) warnings.push(`discharge at ${Math.max(ctp.out.T, cb.out.T).toFixed(0)} C: the bottom cycle brings more heat than the top suction can usefully take (lower the bottom/top mass flow)`);
      /* source water in series: it cools in the top evaporator first (the warmer one, 20 C), then in the bottom one,
         so the water between them sits at Tw_mid by the energy balance; each evaporator's cold end (the refrigerant
         evaporates at constant temperature) must keep the approach to the water leaving it */
      const Tw_in = c.T_source_in, Tw_out = c.T_source_out ?? c.T_source_in;
      const fMid = qEb / (qEt + qEb), Tw_mid = srcT(c, fMid);
      const gapEt = Tw_mid - c.T_ev_t, gapEb = Tw_out - c.T_ev_b;
      if (gapEt < c.dT_ev - 1e-6) warnings.push(`top evaporator: the ${c.source_name || 'source water'} leaves it at ${Tw_mid.toFixed(1)} C, ${gapEt.toFixed(1)} K above the ${c.T_ev_t} C evaporating temperature (approach ${c.dT_ev} K); lower the top evaporating temperature or raise the source outlet`);
      if (gapEb < c.dT_ev - 1e-6) warnings.push(`bottom evaporator: the ${c.source_name || 'source water'} leaves at ${Tw_out.toFixed(1)} C, ${gapEb.toFixed(1)} K above the ${c.T_ev_b} C evaporating temperature (approach ${c.dT_ev} K); lower the bottom evaporating temperature or raise the source outlet`);
      const ap = (hh, hc, Ti, To) => airProfile(ft, air, c.P_t, hh, hc, Ti, To, 16);
      const pHT = ap(ctp.out.h, hMid, c.T_mt, c.T_air_tar), pMT = ap(hMid, hOut, Tain, c.T_mt);
      const gaps = { HT: minGap(pHT), MT: minGap(pMT) };
      if (c._lite) return { success: true, COP: cop, warnings, report: { r, air_ratio: phi, gap_HT: gaps.HT, gap_MT: gaps.MT, T0: [Tg] } };
      const vt = 1 / ts.D, vb = 1 / cbIn.D;
      // state points for checks (per kg of top-cycle refrigerant; the bottom flow is r)
      const st = (pr, P, h) => { const x = pr.PH(P, h); return { T: x.T, P: x.P, h: x.h, s: x.s }; };
      const states = {
        top: { evap_out: st(ft, Pet, et.h), ihx_cold_out: st(ft, Pet, ctIn.h), suction: st(ft, Pet, ts.h), discharge: st(ft, c.P_t, ctp.out.h),
          htgc_out: st(ft, c.P_t, hMid), gc_out: st(ft, c.P_t, hOut), valve_in: st(ft, c.P_t, hvt), valve_out: st(ft, Pet, het) },
        bot: { evap_out: st(fb, Peb, eb.h), suction: st(fb, Peb, cbIn.h), discharge: st(fb, c.P_b, cb.out.h), cgc_out: st(fb, c.P_b, hbo),
          valve_in: st(fb, c.P_b, hvb), valve_out: st(fb, Peb, heb) },
        r, air_primary: mP, air_secondary: phi * mP,
      };
      const loops = [
        { name: "Top cycle", fluid: 1, flow: 1, paths: [
          { name: c.expander ? "Expander" : "Expansion valve", pts: compPath(ft, c.P_t, hvt, Pet, het) },
          { name: "Evaporator", pts: iso(ft, Pet, het, et.h) },
          ...(qit > 0 ? [{ name: "Internal heat exchanger", pts: iso(ft, Pet, et.h, ctIn.h) }] : []),
          { name: "Cascade gas cooler (cold)", pts: iso(ft, Pet, ctIn.h, ts.h) },
          { name: "Compressor", pts: compPath(ft, Pet, ts.h, c.P_t, ctp.out.h) },
          { name: "HT gas cooler", pts: iso(ft, c.P_t, ctp.out.h, hMid) },
          { name: "MT gas cooler", pts: iso(ft, c.P_t, hMid, hOut) },
          ...(qit > 0 ? [{ name: "Internal heat exchanger", pts: iso(ft, c.P_t, hOut, hvt) }] : []),
        ] },
        { name: "Bottom cycle", fluid: 0, flow: r, paths: [
          { name: c.expander ? "Expander" : "Expansion valve", pts: compPath(fb, c.P_b, hvb, Peb, heb) },
          { name: "Evaporator", pts: iso(fb, Peb, heb, eb.h) },
          ...(qib > 0 ? [{ name: "Internal heat exchanger", pts: iso(fb, Peb, eb.h, cbIn.h) }] : []),
          { name: "Compressor", pts: compPath(fb, Peb, cbIn.h, c.P_b, cb.out.h) },
          { name: "Cascade gas cooler (hot)", pts: iso(fb, c.P_b, cb.out.h, hbo) },
          ...(qib > 0 ? [{ name: "Internal heat exchanger", pts: iso(fb, c.P_b, hbo, hvb) }] : []),
        ] },
      ];
      const cgc = hxPairs(fb, c.P_b, cb.out.h, hbo, r, ft, Pet, ctIn.h, ts.h, 1);
      const hx = [
        { name: "HT gas cooler", key: "gc", hot: pHT.map(p => [qMT + p[0], p[1]]), cold: pHT.map(p => [qMT + p[0], p[2]]), otherName: "primary air" },
        { name: "MT gas cooler", key: "gc2", hot: pMT.map(p => [p[0], p[1]]), cold: pMT.map(p => [p[0], p[2]]), otherName: "primary + secondary air" },
        { name: "Cascade gas cooler", key: "rhx", hot: cgc.hot, cold: cgc.cold, otherName: "top-cycle suction" },
        // one water line through both evaporators: bottom from 0 to its duty, top stacked after it
        { name: "Top evaporator", key: "evap", hot: waterLine(qEt, c, qEb, fMid, 1), cold: iso(ft, Pet, het, et.h).map(p => [qEb + (p[2] - het), p[0]]), otherName: c.source_name || "source water", coldRef: true },
        { name: "Bottom evaporator", key: "evap2", hot: waterLine(qEb, c, 0, 0, fMid), cold: iso(fb, Peb, heb, eb.h).map(p => [r * (p[2] - heb), p[0]]), otherName: c.source_name || "source water", coldRef: true },
      ];
      return {
        success: true, message: "success", cycle: "K", warnings, states, COP: cop, VHC: qAir / (vt + r * vb) / 1000, loops, hx,
        q_air: qAir, w_shaft: w, q_source: qEt + qEb, ref: { rho: cbIn.D, flow: r, PR: c.P_b / Peb },
        Pmax_bar: Math.max(c.P_t, c.P_b) * 10, T_dis_max: Math.max(ctp.out.T, cb.out.T),
        air: [{ name: "Primary air", q: qPri, from: Tain, to: c.T_air_tar }, ...(phi > 0 ? [{ name: "Secondary air", q: qSec, from: Tain, to: c.T_mt }] : [])],
        report: {
          COP: cop, VHC: qAir / (vt + r * vb) / 1000, P_top_bar: c.P_t * 10, P_bot_bar: c.P_b * 10, P_evap_top_bar: Pet * 10, P_evap_bot_bar: Peb * 10,
          PR_top: c.P_t / Pet, PR_bot: c.P_b / Peb, eta_top: ctp.eta, eta_bot: cb.eta, T_dis_top: ctp.out.T, T_dis_bot: cb.out.T,
          T_suc_top: ts.T, q_cgc: cg.q, q_HT: qHT, q_MT: qMT, q_primary: qPri, q_secondary: qSec, share_sec: qSec / qAir,
          q_ihx_top: qit, q_ihx_bot: r * qib, gap_HT: gaps.HT, gap_MT: gaps.MT, T_gc_mid: Tph(ft, c.P_t, hMid), T_gc_out: Tout,
          T_bot_cgc_out: Tph(fb, c.P_b, hbo), T_water_mid: Tw_mid, gap_evap_top: gapEt, gap_evap_bot: gapEb, q_evap_top: qEt, q_evap_bot: qEb, gc_off: c.gc_off || 0, T_bot_valve: Tph(fb, c.P_b, hvb), primary_per_top: mP,
          bottom_supercritical: isFinite(fb.P_crit) ? c.P_b > fb.P_crit : null, w_el_kJkg: wEl, r, air_ratio: phi,
          expander: !!c.expander, w_exp_kJkg: wx, exp_out_quality: xq,
        },
      };
    } catch (e) { return fail(`property calculation failed: ${e && e.message ? e.message : e}`); } finally { stage.fixedEta = null; }
  }

  // ================================================================== n-level cascade (Hoang et al.)
  const CASC = { n: 2, T_cond: [100], cond_air: [true], sub_air: [true, true], dT_casc: 5, superheat: 5, ihx: [], dT_ihx: 5 };
  function cascade(fls, air, cfg) {
    const c = Object.assign({ T_air_inlet: 15, T_air_tar: 200, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, dTmin_rhx: 5 }, CASC, cfg);
    const n = c.n, warnings = [];
    if (fls.length < n) return fail("a refrigerant is needed for every level");
    try {
      const Tain = c.T_air_inlet;
      const ev1 = (c.T_source_out ?? c.T_source_in) - (c.dT_src ?? c.dTmin_rhx);          // level-1 evaporator
      // level k state builder for a condensing temperature Tc and evaporating Te
      // Tihx: with an IHX on the level, the suction is heated to the liquid entering the IHX less its approach
      const level = (pr, Te, Tc, sh, Tihx = null) => {
        const ev = pr.TQ(Te, 1), Pe = ev.P, bub = pr.TQ(Tc, 0), dew = pr.TQ(Tc, 1), Pc = Math.max(bub.P, dew.P);
        if (Pc <= Pe) throw new Error(`condensing ${Tc.toFixed(1)} C not above evaporating ${Te.toFixed(1)} C`);
        // suction superheat (an IHX against the liquid): at least sh, and enough for a dry discharge 5 K above dew
        const dry = s => { const o = stage(pr, Pe, pr.PT(Pe, Te + s), Pc).out; return o.T - dew.T - 5; };
        let shN = Math.max(sh, 0.5);
        if (dry(shN) < 0) {
          let hi = shN;
          while (dry(hi) < 0 && hi < Tc - Te - 5) hi += 10;
          if (dry(hi) < 0) throw new Error(`level with ${Te.toFixed(0)} -> ${Tc.toFixed(0)} C cannot compress dry (raise the superheat or change the fluid)`);
          let lo = shN; for (let k = 0; k < 30 && hi - lo > 0.05; k++) { const m = (lo + hi) / 2; if (dry(m) >= 0) hi = m; else lo = m; }
          shN = hi;
        }
        if (Tihx != null && Tihx - Te > shN) shN = Tihx - Te;
        const suc = pr.PT(Pe, Te + shN), qIhx = Math.max(0, suc.h - pr.PT(Pe, Te + Math.min(sh, shN)).h);
        const cp = stage(pr, Pe, suc, Pc);
        const lb = pr.PQ(Pc, 0);
        return { ev, Pe, Pc, suc, cp, hb: lb.h, Tb: lb.T, Td: dew.T, sh: shN, qIhx };
      };
      const Tcs = Array.from({ length: n - 1 }, (_, k) => c.T_cond[k] ?? 100);
      // a lower level that condenses into the air does so until the air leaves it at its condensing temperature less
      // the approach; the rest of its desuperheating and condensing heat drives the next level
      const fK = Array(n).fill(0.2);
      // a lower level's subcooling into the air is capped so the air leaving its subcooler stays at least the approach
      // below its condensing temperature (a big lower-level flow could otherwise overheat the air with liquid alone)
      const subMax = Array(n).fill(Infinity);
      const fA = k => (c.cond_air[k] ? fK[k] : 0);
      const sub = k => !!c.sub_air[k];
      const ihxOn = k => !!(c.ihx || [])[k];
      const Tliq = Array(n).fill(null);           // liquid entering each level's IHX (after any subcooling to the air)
      // for a top condensing temperature: fixed point on the air temperatures entering each level
      const solveAt = (Tcn, airIn0) => {
        let airIn = airIn0.slice(), L = null;
        for (let it = 0; it < 60; it++) {
          L = [];
          let Te = ev1;
          for (let k = 0; k < n; k++) {
            const Tc = k < n - 1 ? Tcs[k] : Tcn;
            L.push(level(fls[k], Te, Tc, c.superheat, ihxOn(k) && Tliq[k] != null ? Tliq[k] - c.dT_ihx : null));
            Te = Tc - c.dT_casc;
          }
          for (let k = 0; k < n; k++) {
            const l = L[k], pr = fls[k], top = k === n - 1;
            // subcooled into the air down to the air inlet + approach, but not below what the IHX needs
            const Tfloor = l.ev.T + l.sh + 5;
            l.hsub = sub(k) ? Math.min(l.hb, pr.PT(l.Pc, Math.min(Math.max(airIn[k] + c.dT_GC_min, Tfloor), l.Tb - 0.01)).h) : l.hb;
            if (!top && isFinite(subMax[k])) l.hsub = Math.max(l.hsub, l.hb - Math.max(subMax[k], 0));
            l.qa = top ? l.cp.out.h - l.hsub : fA(k) * (l.cp.out.h - l.hb) + (l.hb - l.hsub);
            l.qc = top ? 0 : (1 - fA(k)) * (l.cp.out.h - l.hb);
            l.hv = l.hsub - l.qIhx; l.he = expandTo(pr, l.Pc, l.hv, l.Pe, c.expander);
            l.ihxShort = pr.PH(l.Pc, l.hv).T < l.ev.T + 2;
            if (ihxOn(k)) { const t = pr.PH(l.Pc, l.hsub).T; l.dTliq = Tliq[k] == null ? Infinity : Math.abs(t - Tliq[k]); Tliq[k] = t; } else l.dTliq = 0;
          }
          const m = [1];
          for (let k = 0; k < n - 1; k++) m.push(m[k] * L[k].qc / (L[k + 1].suc.h - L[k + 1].qIhx - L[k + 1].he));
          const Q = L.map((l, k) => m[k] * l.qa), Qt = Q.reduce((x, y) => x + y, 0);
          const ha0 = air.h(Tain), mA = Qt / (air.h(c.T_air_tar) - ha0);
          const nIn = []; let ha = ha0;
          for (let k = 0; k < n; k++) { nIn.push(air.T(ha)); ha += Q[k] / mA; }
          let dfMax = 0;
          for (let k = 0; k < n - 1; k++) {
            if (!c.cond_air[k]) continue;
            const l = L[k], want = mA * (air.h(l.Tb - c.dT_GC_min) - air.h(nIn[k])) / m[k] - (l.hb - l.hsub);
            const f = Math.min(Math.max(want / (l.cp.out.h - l.hb), 0), 0.95);
            dfMax = Math.max(dfMax, Math.abs(f - fK[k])); fK[k] = fK[k] + 0.7 * (f - fK[k]);
          }
          let dsMax = 0;
          for (let k = 0; k < n - 1; k++) {
            if (!sub(k)) continue;
            const room = Math.max(0, mA * (air.h(L[k].Tb - c.dT_GC_min) - air.h(nIn[k])) / m[k]);
            const cur = L[k].hb - L[k].hsub;
            if (cur > room + 1e-9 || isFinite(subMax[k])) {
              const nv = isFinite(subMax[k]) ? subMax[k] + 0.7 * (room - subMax[k]) : room;
              dsMax = Math.max(dsMax, Math.abs(nv - (isFinite(subMax[k]) ? subMax[k] : cur)));
              subMax[k] = nv;
            }
          }
          const conv = Math.max(...nIn.map((v, k) => Math.abs(v - airIn[k]))) < 0.01 && dfMax < 1e-4 && dsMax < 1e-3 && Math.max(...L.map(l => l.dTliq)) < 0.02;
          airIn = nIn; Object.assign(L, { m, Q, mA, airIn, f: fK.slice() });
          if (conv) break;
        }
        const t = L[n - 1];
        L.gap = minGap(airProfile(fls[n - 1], air, t.Pc, t.cp.out.h, t.hsub, L.airIn[n - 1], c.T_air_tar, 18));
        return L;
      };
      const prTop = fls[n - 1], TcritTop = isFinite(prTop.T_crit) ? prTop.T_crit - 1 : c.T_air_tar + 60;
      let lo = (n > 1 ? Tcs[n - 2] : ev1) + c.dT_casc + 2, hi = Math.min(TcritTop, c.T_air_tar + 40), airIn = Array(n).fill(Tain), L;
      // scan up for the first condensing temperature that keeps the top air approach, then bisect
      let prev = lo, found = null, Lf = null, lastErr = '';
      for (let T = lo; T <= hi + 1e-9; T += 5) {
        let Lt = null; try { Lt = solveAt(T, airIn); } catch (e) { lastErr = e.message; }
        if (Lt) { airIn = Lt.airIn; if (Lt.gap >= c.dT_GC_min) { found = T; Lf = Lt; break; } }
        prev = T;
      }
      if (found == null) throw new Error(`no top condensing temperature up to ${hi.toFixed(0)} C heats the air to ${c.T_air_tar} C with a ${c.dT_GC_min} K approach${lastErr ? ` (${lastErr})` : ''}`);
      lo = prev; hi = found; L = Lf;
      for (let k = 0; k < 25 && hi - lo > 0.02; k++) {
        const mid = (lo + hi) / 2;
        let Lm = null; try { Lm = solveAt(mid, airIn); } catch (e) { Lm = null; }
        if (Lm && Lm.gap >= c.dT_GC_min) { hi = mid; L = Lm; airIn = Lm.airIn; } else lo = mid;
      }
      const Tcn = hi;
      const m = L.m, Q = L.Q;
      L.forEach((l, k) => { if (l.ihxShort) warnings.push(`level ${k + 1}: its IHX needs more liquid subcooling than is left`); });
      const w = L.reduce((s, l, k) => s + m[k] * l.cp.w, 0), wx = L.reduce((s, l, k) => s + m[k] * (l.hv - l.he), 0);
      const wEl = expEl(w, wx, c.motor_eff ?? MOTOR), qAir = Q.reduce((a, b) => a + b, 0), cop = qAir / wEl;
      const xq = c.expander ? L.map((l, k) => expWarn(fls[k], l.Pe, l.he, warnings, `level ${k + 1} expander`)) : null;
      // air approach of each level
      const loops = [], hx = [];
      let hamb = air.h(Tain);
      const mA = L.mA; airIn = L.airIn;
      L.forEach((l, k) => {
        const pr = fls[k], top = k === n - 1;
        const hv = l.hsub, hval = l.hv, hEvOut = l.suc.h - l.qIhx;
        const paths = [
          { name: c.expander ? "Expander" : "Expansion valve", pts: compPath(pr, l.Pc, hval, l.Pe, l.he) },
          { name: "Evaporator", pts: iso(pr, l.Pe, l.he, hEvOut) },
          ...(l.qIhx > 1e-6 ? [{ name: "Internal heat exchanger", pts: iso(pr, l.Pe, hEvOut, l.suc.h) }] : []),
          { name: "Compressor", pts: compPath(pr, l.Pe, l.suc.h, l.Pc, l.cp.out.h) },
        ];
        if (top) paths.push({ name: "Air condenser", pts: iso(pr, l.Pc, l.cp.out.h, hv) });
        else {
          paths.push({ name: "Cascade condenser", pts: iso(pr, l.Pc, l.cp.out.h, l.hb) });
          if (l.hb - hv > 1e-6) paths.push({ name: "Subcooler", pts: iso(pr, l.Pc, l.hb, hv) });
        }
        if (l.qIhx > 1e-6) paths.push({ name: "Internal heat exchanger", pts: iso(pr, l.Pc, hv, hval) });
        loops.push({ name: `Level ${k + 1}`, fluid: k, flow: m[k], paths });
        // air side of this level: subcooler then condenser share
        const Ta0 = air.T(hamb), Ta1 = air.T(hamb + Q[k] / mA); hamb += Q[k] / mA;
        if (Q[k] > 1e-9) {
          // the refrigerant duty to air, from the coldest point: subcooling, then f of desuperheat+condensing
          const pts = lin(0, 1, 16).map(f => {
            const q = f * l.qa; let T;
            const qs = l.hb - hv;
            if (q <= qs || top) T = pr.PH(l.Pc, hv + q).T;
            else { const fr = fA(k) > 1e-9 ? (q - qs) / fA(k) : 0; T = pr.PH(l.Pc, Math.min(l.hb + fr, l.cp.out.h)).T; }
            return [m[k] * q, T, air.T(air.h(Ta0) + f * Q[k] / mA)];
          });
          const gap = minGap(pts);
          if (!top && gap < c.dT_GC_min - 0.05) warnings.push(`level ${k + 1} air approach ${gap.toFixed(1)} K, below the ${c.dT_GC_min} K minimum: lower its condensing-to-air fraction or raise its condensing temperature`);
          if (gap < 0) throw new Error(`level ${k + 1} cannot heat the air from ${Ta0.toFixed(0)} to ${Ta1.toFixed(0)} C (temperature cross)`);
          hx.push({ name: `Level ${k + 1} to air`, key: k === n - 1 ? "gc" : `gc${k + 2}`, hot: pts.map(p => [p[0], p[1]]), cold: pts.map(p => [p[0], p[2]]), otherName: "dryer air", gap });
        }
        if (k < n - 1) {
          const nx = L[k + 1], prn = fls[k + 1], qc = m[k] * l.qc;
          // cascade condenser: level k's (1 - f) share of desuperheating+condensing against level k+1's evaporator
          const hot = lin(0, 1, 14).map(f => [f * qc, pr.PH(l.Pc, l.hb + f * (l.cp.out.h - l.hb)).T]);
          const hc0 = nx.he, hc1 = nx.suc.h - nx.qIhx;
          const cold = lin(0, 1, 14).map(f => [f * qc, prn.PH(nx.Pe, hc0 + f * (hc1 - hc0)).T]);
          hx.push({ name: `Cascade ${k + 1}→${k + 2}`, key: "rhx", hot, cold, otherName: `level ${k + 2} evaporator` });
        }
        if (l.qIhx > 1e-6) {
          // IHX: the liquid (after any subcooling to the air) heats the evaporator outlet vapour, counterflow
          const Q = m[k] * l.qIhx;
          const hot = lin(0, 1, 12).map(f => [f * Q, pr.PH(l.Pc, hval + f * l.qIhx).T]);
          const cold = lin(0, 1, 12).map(f => [f * Q, pr.PH(l.Pe, hEvOut + f * l.qIhx).T]);
          const gap = Math.min(...hot.map((p, i) => p[1] - cold[i][1]));
          if (ihxOn(k) && gap < c.dT_ihx - 0.05) warnings.push(`level ${k + 1} IHX approach ${gap.toFixed(1)} K, below its ${c.dT_ihx} K`);
          hx.push({ name: `Level ${k + 1} IHX`, key: "ihx", hot, cold, otherName: `level ${k + 1} suction`, gap });
        }
      });
      const l1 = L[0], qE = l1.suc.h - l1.qIhx - l1.he;
      hx.push({ name: "Level 1 evaporator", key: "evap", hot: waterLine(qE, c), cold: iso(fls[0], l1.Pe, l1.he, l1.suc.h - l1.qIhx).map(p => [p[2] - l1.he, p[0]]), otherName: c.source_name || "source water", coldRef: true });
      const vs = L.reduce((s, l, k) => s + m[k] / l.suc.D, 0);
      return {
        success: true, message: "success", cycle: "N", warnings, COP: cop, VHC: qAir / vs / 1000, loops, hx,
        q_air: qAir, w_shaft: w, q_source: qE, ref: { rho: l1.suc.D, flow: 1, PR: l1.Pc / l1.Pe },
        Pmax_bar: Math.max(...L.map(l => l.Pc)) * 10, T_dis_max: Math.max(...L.map(l => l.cp.out.T)),
        air: L.map((l, k) => ({ name: `Level ${k + 1}`, q: Q[k] })),
        report: {
          COP: cop, VHC: qAir / vs / 1000, levels: L.map((l, k) => ({ Te: l.ev.T, Tc: l.Tb, Pe_bar: l.Pe * 10, Pc_bar: l.Pc * 10, PR: l.Pc / l.Pe,
            eta: l.cp.eta, Tdis: l.cp.out.T, flow: m[k], q_air: Q[k], air_in: airIn[k], superheat: l.sh, q_ihx: l.qIhx, f_air: k < n - 1 ? fA(k) : 1 })), w_el_kJkg: wEl, T_cond_top: Tcn,
          expander: !!c.expander, w_exp_kJkg: wx, exp_out_quality: xq,
        },
      };
    } catch (e) { return fail(`${e && e.message ? e.message : e}`); }
  }

  // ================================================================== booster: bottom cycle + MVR steam
  const BOOST = { T_cond: 90, dT_sg: 5, n_mvr: 9, PR_mvr: 1.5, inject: true, sh_inject: 5, superheat: 5, sub_air: true, ihx: false, dT_ihx: 5 };
  /* HP + MVR. With split_fill (Auto), each fan's process share heats the air to its saturation less the approach;
     when that leaves the next stream short of its approach (the last fan's desuperheating needs the air a little
     cooler where it starts condensing), every cap is lowered by the smallest common margin that keeps it. */
  function booster(fb, water, air, cfg) {
    const r = boosterFill(fb, water, air, cfg);
    if (!r.success || !(r.report.sub_raise > 0.05) || (cfg && cfg.t1_drop != null)) return r;
    /* the subcooler pinched inside (a liquid whose cp climbs steeply, like propane near its critical point), and
       raising the liquid outlet fixed it; the other way is to keep the liquid at the air inlet plus the approach and
       let the bottom cycle heat the air less (the steam takes the rest). Keep whichever gives the higher COP. */
    const dT = cfg.dT_GC_min ?? 20, at = d => boosterFill(fb, water, air, { ...cfg, t1_drop: d, no_sub_raise: true });
    const ok = q => q.success && q.report.gap_bottom >= dT - 0.02;
    let hi = 2, rh = at(hi);
    while (!ok(rh) && hi < 64) { hi *= 2; rh = at(hi); }
    if (!ok(rh)) return r;
    let lo = 0;
    for (let k = 0; k < 12 && hi - lo > 0.05; k++) { const m = (lo + hi) / 2, rm = at(m); if (ok(rm)) { hi = m; rh = rm; } else lo = m; }
    return rh.COP > r.COP ? rh : r;
  }
  function boosterFill(fb, water, air, cfg) {
    if (cfg && cfg.split_fill && cfg.fill_margin == null) {
      const dT = cfg.dT_GC_min ?? 20, ok = r => r.success && r.report.gap_steam >= dT - 0.05;
      const r0 = boosterAt(fb, water, air, { ...cfg, fill_margin: 0 });
      if (!r0.success || ok(r0)) return r0;
      let hi = 2, rh = boosterAt(fb, water, air, { ...cfg, fill_margin: hi });
      while (!ok(rh) && hi < 64) { hi *= 2; rh = boosterAt(fb, water, air, { ...cfg, fill_margin: hi }); }
      if (!ok(rh)) return r0;
      let lo = 0;
      for (let k = 0; k < 18 && hi - lo > 0.01; k++) {
        const m = (lo + hi) / 2, rm = boosterAt(fb, water, air, { ...cfg, fill_margin: m });
        if (ok(rm)) { hi = m; rh = rm; } else lo = m;
      }
      rh.report.fill_margin = hi;
      return rh;
    }
    return boosterAt(fb, water, air, cfg);
  }
  function boosterAt(fb, water, air, cfg) {
    const c = Object.assign({ T_air_inlet: 15, T_air_tar: 200, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, dTmin_rhx: 5 }, BOOST, cfg);
    const warnings = [];
    // MVR fans: their own efficiency model (fan_model, fan_eta for Constant) and motor, else the compressors'
    const etaF = P => c.fan_model != null ? etaModel(c.fan_model, c.fan_eta ?? 0.8, c.PR_mvr, P) : etaIs(c.PR_mvr, P);
    const motF = c.fan_motor ?? c.motor_eff ?? MOTOR;
    try {
      const Tain = c.T_air_inlet;
      const Te = (c.T_source_out ?? c.T_source_in) - (c.dT_src ?? c.dTmin_rhx);
      const ev = fb.TQ(Te, 1), Pe = ev.P, cb = fb.TQ(c.T_cond, 0), Pc = Math.max(cb.P, fb.TQ(c.T_cond, 1).P);
      const evOut = c.superheat > 0 ? fb.PT(Pe, Math.min(Te + c.superheat, c.T_source_in - (c.dT_src ?? c.dTmin_rhx))) : ev;
      // bottom compressor suction: the evaporator outlet, or with the IHX on, heated by the liquid (after any
      // subcooling to the air) to its temperature less the IHX approach
      const compAt = hsub => {
        let suc = evOut;
        if (c.ihx) { const Ts = fb.PH(Pc, hsub).T - c.dT_ihx; if (Ts > evOut.T + 0.01) suc = fb.PT(Pe, Ts); }
        return { suc, cp: stage(fb, Pe, suc, Pc), qIhx: suc.h - evOut.h };
      };
      const hb = fb.PQ(Pc, 0).h;
      // steam: raised at T_cond - dT_sg, n fans at PR each. After fan i (all but the last) water is injected to
      // saturation + sh, then a fraction x_i of the exhaust goes to the process (condenses into the air at that
      // pressure) and the rest to fan i+1; the last fan's exhaust all goes to the process. The fan states do not
      // depend on the splits, only the masses do, so the chain is built once per kg of steam raised.
      const Ts = c.T_cond - c.dT_sg, sv = water.TQ(Ts, 1), sl = water.TQ(Ts, 0);
      const nS = Math.max(1, Math.round(c.n_mvr));
      const T1 = c.T_cond - c.dT_GC_min - (c.t1_drop || 0);     // air after the bottom cycle (approach at condensing, or lower)
      if (T1 <= Tain) return fail("condensing temperature too low to heat the air");
      const ha0 = air.h(Tain), ha1 = air.h(T1), ha2 = air.h(c.T_air_tar);
      let P = sv.P, st = sv;
      const fans = [], spaths = [];
      for (let i = 0; i < nS; i++) {
        const P2 = P * c.PR_mvr, e = etaF(P), s = stage(water, P, st, P2, e), last = i === nS - 1;
        spaths.push({ name: `MVR ${i + 1}`, pts: compPath(water, P, st.h, P2, s.out.h, 5) });
        let out = s.out, inj = 0;
        if (c.inject && !last) {
          const tgt = water.PT(P2, water.PQ(P2, 1).T + Math.max(c.sh_inject, 0.1));
          if (out.h > tgt.h) { inj = (out.h - tgt.h) / (tgt.h - sl.h); spaths.push({ name: "Water injection", pts: [pt(out), pt(tgt)] }); out = tgt; }
        }
        const lf = water.PQ(P2, 0), sat = water.PQ(P2, 1);
        fans.push({ P_in: P, P_out: P2, T_out: s.out.T, w: s.w, eta: e, inj, out, lf, Tsat: sat.T, q: out.h - lf.h });
        P = P2; st = out;
      }
      // masses per kg raised for process fractions x (last forced to 1)
      const chain = x => {
        let M = 1, J = 0, wS = 0;
        const E = [];
        fans.forEach((f, i) => {
          wS += M * f.w; J += M * f.inj;
          const Mo = M * (1 + f.inj), xi = i === nS - 1 ? 1 : Math.min(Math.max(+x[i] || 0, 0), 1);
          f.M = M; E.push(Mo * xi); M = Mo * (1 - xi);
        });
        const qS = E.reduce((a, e, i) => a + e * fans[i].q, 0);
        const g = sv.h + J * sl.h - E.reduce((a, e, i) => a + e * fans[i].lf.h, 0);   // generator, per kg raised
        return { E, J, wS, qS, g };
      };
      // split fractions: given, or each fan's exhaust sized to heat the air to its saturation less the approach
      // (split_fill); then the air per kg raised (a) is the root that brings the last fan's steam to the air target
      let x = Array.from({ length: nS }, (_, i) => i === nS - 1 ? 1 : Math.min(Math.max(+((c.split || [])[i]) || 0, 0), 1));
      if (c.split_fill && nS > 1) {
        const fill = a => {
          let M = 1, ha = ha1; const xs = [];
          fans.forEach((f, i) => {
            const Mo = M * (1 + f.inj);
            if (i === nS - 1) { xs.push(1); ha += Mo * f.q / a; return; }
            const cap = air.h(Math.min(f.Tsat - c.dT_GC_min - (c.fill_margin || 0), c.T_air_tar));
            const Ei = cap > ha ? Math.min(Mo, a * (cap - ha) / f.q) : 0;
            xs.push(Mo > 0 ? Ei / Mo : 0); ha += Ei * f.q / a; M = Mo - Ei;
          });
          return { xs, f: ha - ha2 };
        };
        let lo = Math.log(1e-3), hi = Math.log(1e5);
        if (fill(Math.exp(lo)).f > 0 && fill(Math.exp(hi)).f < 0) {
          for (let k = 0; k < 80; k++) { const m = (lo + hi) / 2; if (fill(Math.exp(m)).f > 0) lo = m; else hi = m; }
          x = fill(Math.exp(lo)).xs;
        }
      }
      const ch = chain(x), E = ch.E, wS = ch.wS, g = ch.g, qS = ch.qS;
      const lastF = fans[nS - 1], satF = { T: lastF.Tsat }, Pf = lastF.P_out, stF = lastF.out;
      if (satF.T < c.T_air_tar && stF.T < c.T_air_tar + c.dT_GC_min) warnings.push("final steam is cooler than the air target plus the approach: add MVR stages");
      // bottom condenser: subcooling to the air (air enters at Tain), then the air section up to T1, the rest raises steam
      // the subcooled liquid's T(h) bends (its cp rises with temperature), so with the liquid leaving at the air inlet
      // plus the approach the pinch can sit a little inside the subcooler: raise the liquid outlet until it holds
      const bottom = subOff => {
        const hsub = c.sub_air ? fb.PT(Pc, Math.min(Tain + c.dT_GC_min + subOff, cb.T - 0.01)).h : hb;
        const K = compAt(hsub), qCond = K.cp.out.h - hsub;
        // per kg bottom refrigerant: mA (ha1-ha0) + raised g = qCond ; mA (ha2-ha1) = raised qS
        const k2 = (ha2 - ha1) / qS;
        const mA = qCond / ((ha1 - ha0) + k2 * g), raised = mA * k2;
        const qA1 = mA * (ha1 - ha0);
        const p1 = lin(0, 1, 24).map(f => { const q = f * qA1; return [q, fb.PH(Pc, hsub + q).T, air.T(ha0 + q / mA)]; });
        return { hsub, qCond, mA, raised, qA1, p1, g1: minGap(p1), K };
      };
      let subOff = 0, B = bottom(0);
      for (let k = 0; k < 25 && c.sub_air && !c.no_sub_raise && B.g1 < c.dT_GC_min - 0.01 && B.mA > 0; k++) { subOff += (c.dT_GC_min - B.g1) * 1.5 + 0.005; B = bottom(subOff); }
      const { hsub, qCond, mA, raised, qA1, p1, g1 } = B, { suc, cp, qIhx } = B.K, hv = hsub - qIhx;
      const he = expandTo(fb, Pc, hv, Pe, c.expander), wx = hv - he;
      if (!isGas(cp.out)) return fail("wet compression in the bottom cycle");
      let gIhx = null;
      if (qIhx > 1e-6) { gIhx = fb.PH(Pc, hv).T - evOut.T; if (gIhx < c.dT_ihx - 0.05) warnings.push(`bottom IHX cold-end approach ${gIhx.toFixed(1)} K, below its ${c.dT_ihx} K`); }
      if (!(raised > 0) || !(mA > 0)) return fail("no steam split balances the condenser");
      const qA2 = mA * (ha2 - ha1), qAir = qA1 + qA2;
      const w = cp.w + raised * wS, wEl = expEl(cp.w, wx, c.motor_eff ?? MOTOR) + raised * wS / motF, cop = qAir / wEl;
      const xq = c.expander ? expWarn(fb, Pe, he, warnings) : null;
      // approach checks: the bottom cycle (above), then each process steam stream in turn (counterflow within each)
      const p2 = []; let qAcc = qA1, haAcc = ha1;
      fans.forEach((f, i) => {
        const Q = raised * E[i] * f.q; if (!(Q > 1e-9)) return;
        for (const fr of lin(0, 1, 12)) p2.push([qAcc + fr * Q, water.PH(f.P_out, f.lf.h + fr * f.q).T, air.T(haAcc + fr * Q / mA)]);
        qAcc += Q; haAcc += Q / mA;
      });
      const g2 = minGap(p2);
      if (g2 < 0) return fail(`the process steam cannot heat the air to ${c.T_air_tar} C (temperature cross): add MVR stages, send less exhaust to the process early, or raise the bottom condensing temperature`);
      if (g2 < c.dT_GC_min - 0.05) warnings.push(`steam-to-air approach ${g2.toFixed(1)} K, below the ${c.dT_GC_min} K minimum: add MVR stages or send less exhaust to the process early`);
      if (g1 < c.dT_GC_min - 0.05) warnings.push(`bottom-cycle air approach ${g1.toFixed(1)} K`);
      if (sv.P < 0.1013) warnings.push(`steam raised at ${(sv.P * 10).toFixed(2)} bar (vacuum)`);
      const qSG = raised * g, qE = evOut.h - he, mS = raised * E.reduce((a, b) => a + b, 0);
      const stages = fans.map((f, i) => ({ P_in: f.P_in, P_out: f.P_out, T_out: f.T_out, T_sat: f.Tsat, inj: f.inj, mass: f.M, eta: f.eta,
        split: x[i], to_process: raised * E[i], q_air: raised * E[i] * f.q }));
      const condPaths = [];
      fans.forEach((f, i) => {
        if (!(E[i] > 1e-9)) return;
        condPaths.push({ name: i === nS - 1 ? "Steam condenser" : `Process steam ${i + 1}`, pts: iso(water, f.P_out, f.out.h, f.lf.h) });
        condPaths.push({ name: "Flash", pts: [pt(f.lf), pt(water.PH(sv.P, f.lf.h))] });
      });
      const loops = [
        { name: "Bottom cycle", fluid: 0, flow: 1, paths: [
          { name: c.expander ? "Expander" : "Expansion valve", pts: compPath(fb, Pc, hv, Pe, he) },
          { name: "Evaporator", pts: iso(fb, Pe, he, evOut.h) },
          ...(qIhx > 1e-6 ? [{ name: "Internal heat exchanger", pts: iso(fb, Pe, evOut.h, suc.h) }, { name: "Internal heat exchanger", pts: iso(fb, Pc, hsub, hv) }] : []),
          { name: "Compressor", pts: compPath(fb, Pe, suc.h, Pc, cp.out.h) },
          { name: "Steam generator", pts: iso(fb, Pc, cp.out.h, hsub + qA1) },
          { name: "Air condenser", pts: iso(fb, Pc, hsub + qA1, hsub) },
        ] },
        { name: "Steam", fluid: "water", flow: mS, paths: [
          { name: "Steam generator", pts: [pt(sl), pt(sv)] },
          ...spaths,
          ...condPaths,
        ] },
      ];
      const sgHot = lin(0, 1, 14).map(f => [f * qSG, fb.PH(Pc, hsub + qA1 + f * (cp.out.h - hsub - qA1)).T]);
      const hx = [
        { name: "Bottom cycle to air", key: "gc2", hot: p1.map(p => [p[0], p[1]]), cold: p1.map(p => [p[0], p[2]]), otherName: "dryer air" },
        { name: "Steam to air", key: "gc", hot: p2.map(p => [p[0], p[1]]), cold: p2.map(p => [p[0], p[2]]), otherName: "dryer air" },
        { name: "Steam generator", key: "rhx", hot: sgHot, cold: [[0, Ts], [qSG, Ts]], otherName: `steam at ${Ts.toFixed(0)} C` },
        { name: "Evaporator", key: "evap", hot: waterLine(qE, c), cold: iso(fb, Pe, he, evOut.h).map(p => [p[2] - he, p[0]]), otherName: c.source_name || "source water", coldRef: true },
        ...(qIhx > 1e-6 ? [{ name: "Bottom IHX", key: "ihx", hot: lin(0, 1, 12).map(f => [f * qIhx, fb.PH(Pc, hv + f * qIhx).T]), cold: lin(0, 1, 12).map(f => [f * qIhx, fb.PH(Pe, evOut.h + f * qIhx).T]), otherName: "bottom suction" }] : []),
      ];
      return {
        success: true, message: "success", cycle: "M", warnings, COP: cop, VHC: qAir / (1 / suc.D) / 1000, loops, hx,
        q_air: qAir, w_shaft: w, q_source: qE, ref: { rho: suc.D, flow: 1, PR: Pc / Pe },
        Pmax_bar: Math.max(Pc, Pf) * 10, T_dis_max: Math.max(cp.out.T, ...stages.map(s => s.T_out)),
        air: [{ name: "Bottom cycle", q: qA1 }, { name: "Steam", q: qA2 }],
        report: {
          COP: cop, P_evap_bar: Pe * 10, P_cond_bar: Pc * 10, PR_bottom: Pc / Pe, eta_bottom: cp.eta, T_dis_bottom: cp.out.T,
          T_steam: Ts, P_steam_bar: sv.P * 10, P_final_bar: Pf * 10, T_sat_final: satF.T, T_steam_final: st.T, eta_mvr: etaF(sv.P), PR_mvr: c.PR_mvr, fan_motor: motF,
          steam_per_kg: raised, steam_final_per_kg: mS, w_bottom: cp.w, w_mvr: raised * wS, q_air_bottom: qA1, q_air_steam: qA2, q_sg: qSG,
          T_air_mid: T1, gap_bottom: g1, gap_steam: g2, stages, w_el_kJkg: wEl, splits: x.slice(0, -1), w_el_fans: raised * wS / motF, sub_raise: subOff, t1_drop: c.t1_drop || 0, q_ihx: qIhx, T_suction: suc.T, gap_ihx: gIhx,
          expander: !!c.expander, w_exp_kJkg: wx, exp_out_quality: xq,
        },
      };
    } catch (e) { return fail(`${e && e.message ? e.message : e}`); }
  }

  /* ------------------------------------------------------------------ TT cascade optimiser
     Maximises the COP over the bottom and top discharge pressures (MPa) and, when the air ratio is free
     (air_ratio 0), the bottom/top mass-flow ratio r; with an air ratio given, r is solved for it at every point.
     A compass search in scaled coordinates from the current design; infeasible points simply lose. Bounds: the
     evaporating pressure + 5 % to Pmax, r in [0.1, 30]. A generator: it yields {nEval, best} after each solve so a
     page can keep painting; the final value is the result. */
  function* kongOptimiseSteps(fb, ft, air, cfg, opt = {}) {
    const c = Object.assign({}, KONG, cfg), Pmax = opt.Pmax ?? 25, maxEval = opt.maxEval ?? 500;
    const freeR = true, dim = 4;
    const Peb = fb.TQ(c.T_ev_b, 1).P, Pet = ft.TQ(c.T_ev_t, 1).P;
    const lo = [Peb * 1.05, Pet * 1.05, Math.log(0.05), 0], hi = [Pmax, Pmax, Math.log(5), 80];
    const clip = x => x.map((v, i) => Math.min(Math.max(v, lo[i]), hi[i]));
    let rGuess = c.r, T0 = null;
    const cfgAt = x => ({ ...c, P_b: x[0], P_t: x[1], r: Math.exp(x[2]), gc_off: x[3] });
    let nEval = 0; const memo = new Map();
    const f = x => {
      const k = x.map(v => v.toFixed(6)).join();
      if (memo.has(k)) return memo.get(k);
      nEval++;
      let y = -Infinity, r = null;
      try { r = kong(fb, ft, air, { ...cfgAt(x), _lite: true }); if (r.success && isFinite(r.COP)) { y = r.COP; if (y > (memo.bestY ?? -Infinity)) { memo.bestY = y; T0 = r.report.T0; if (!freeR) rGuess = r.report.r; } } } catch (e) {}
      const o = { y, r }; memo.set(k, o); return o;
    };
    let x0 = clip([c.P_b, c.P_t, Math.log(c.r), c.gc_off || 0]), b0 = f(x0);
    const start = b0.y;
    yield { nEval, best: b0.y };
    // coarse grid: the COP surface has several local optima (bottom cycle sub- or transcritical, the HT approach
    // edge), so the compass search starts from the current design and from the best grid point
    let xg = null, bg = { y: -Infinity };
    for (const a of [0.3, 0.5, 0.7, 0.9]) for (const b of [0.6, 0.8, 1.0]) for (const r0 of [0.2, 0.35, 0.5, 0.8]) {
      const z = clip([lo[0] + a * (hi[0] - lo[0]), lo[1] + b * (hi[1] - lo[1]), Math.log(r0), 0]), o = f(z);
      if (o.y > bg.y) { bg = o; xg = z; }
      yield { nEval, best: Math.max(b0.y, bg.y) };
    }
    const tol = [0.002, 0.002, 0.002, 0.05];
    function* compass(x, best) {
      let step = [1.0, 1.0, 0.25, 8];          // MPa, MPa, ln r, K
      while (nEval < maxEval && step.some((s, i) => s > tol[i])) {
        let moved = false;
        for (let i = 0; i < dim && nEval < maxEval; i++) {
          if (step[i] <= tol[i]) continue;
          for (const d of [1, -1]) {
            const z = x.slice(); z[i] += d * step[i]; const zc = clip(z);
            if (zc[i] === x[i]) continue;
            const o = f(zc); yield best.y;
            if (o.y > best.y + 1e-7) {
              best = o; x = zc; moved = true;
              for (let k = 0; k < 20 && nEval < maxEval; k++) {
                const w = x.slice(); w[i] += d * step[i]; const wc = clip(w); if (wc[i] === x[i]) break;
                const ow = f(wc); yield best.y;
                if (ow.y > best.y + 1e-7) { best = ow; x = wc; } else break;
              }
              break;
            }
          }
        }
        if (!moved) step = step.map((s, i) => s > tol[i] ? s / 2 : s);
      }
      return { x, best, converged: !step.some((s, i) => s > tol[i]) };
    }
    let x = null, best = { y: -Infinity }, converged = false, overall = Math.max(b0.y, bg.y);
    for (const [xs, bs] of [[x0, b0], [xg, bg]]) {
      if (!xs || !isFinite(bs.y)) continue;
      const g = compass(xs, bs); let st;
      while (!(st = g.next()).done) { overall = Math.max(overall, st.value); yield { nEval, best: overall }; }
      if (st.value.best.y > best.y) ({ x, best, converged } = st.value);
    }
    if (!x) return fail("no feasible design found within the bounds", { nEval });
    // the optimum often sits on the edge of feasibility (the HT gas-cooler approach): round as the page shows the
    // values (0.1 bar, 0.001 in r) and step r up until the rounded design solves
    let P_b = Math.round(x[0] * 100) / 100, P_t = Math.round(x[1] * 100) / 100, rr = Math.round(Math.exp(x[2]) * 1000) / 1000, off = Math.round(x[3] * 10) / 10, chk = null;
    for (let k = 0; k < 60; k++) {
      chk = kong(fb, ft, air, { ...c, P_b, P_t, r: rr, gc_off: off });
      if (chk.success) break;
      rr = Math.round((rr + 0.001) * 1000) / 1000;
    }
    if (chk && chk.success) { x = [P_b, P_t, Math.log(rr), off]; best = { y: chk.COP, r: chk }; }
    const near = (v, b) => Math.abs(v - b) < 1e-3;
    return { success: true, COP: best.y, start, P_b: x[0], P_t: x[1], r: Math.exp(x[2]), gc_off: x[3], freeR, nEval,
      converged,
      atBound: { P_b: near(x[0], hi[0]) || near(x[0], lo[0]), P_t: near(x[1], hi[1]) || near(x[1], lo[1]), r: freeR && (near(x[2], lo[2]) || near(x[2], hi[2])) } };
  }
  function kongOptimise(fb, ft, air, cfg, opt = {}) { const g = kongOptimiseSteps(fb, ft, air, cfg, opt); let s; do { s = g.next(); } while (!s.done); return s.value; }

  /* ------------------------------------------------------------------ shared compass search
     f(x) -> {y, r} (y = -Infinity when infeasible); integer dimensions step by whole numbers and stop at a step of 1.
     Runs from each start in turn and keeps the best; yields the best COP so far after every solve. */
  function* compassSearch(f, starts, lo, hi, step0, tol, isInt, state, maxEval) {
    const dim = lo.length;
    const clip = x => x.map((v, i) => { const c = Math.min(Math.max(v, lo[i]), hi[i]); return isInt[i] ? Math.round(c) : c; });
    let X = null, B = { y: -Infinity }, conv = false;
    for (const [xs0, bs0] of starts) {
      if (!xs0 || !isFinite(bs0.y)) continue;
      let x = xs0, best = bs0, step = step0.slice();
      const live = () => step.some((s, i) => isInt[i] ? s >= 1 : s > tol[i]);
      while (state.nEval < maxEval && live()) {
        let moved = false;
        for (let i = 0; i < dim && state.nEval < maxEval; i++) {
          if (isInt[i] ? step[i] < 1 : step[i] <= tol[i]) continue;
          for (const d of [1, -1]) {
            const z = x.slice(); z[i] += d * step[i]; const zc = clip(z);
            if (zc[i] === x[i]) continue;
            const o = f(zc); yield Math.max(B.y, best.y);
            if (o.y > best.y + 1e-7) {
              best = o; x = zc; moved = true;
              for (let k = 0; k < 20 && state.nEval < maxEval; k++) {        // keep going while it pays
                const w = x.slice(); w[i] += d * step[i]; const wc = clip(w); if (wc[i] === x[i]) break;
                const ow = f(wc); yield Math.max(B.y, best.y);
                if (ow.y > best.y + 1e-7) { best = ow; x = wc; } else break;
              }
              break;
            }
          }
        }
        if (!moved) step = step.map((s, i) => isInt[i] ? (s >= 1 ? Math.floor(s / 2) : 0) : (s > tol[i] ? s / 2 : s));
      }
      if (best.y > B.y) { X = x; B = best; conv = !live(); }
    }
    return { x: X, best: B, converged: conv };
  }
  const memoF = (state, solve) => x => {
    const k = x.map(v => (+v).toFixed(4)).join();
    if (state.memo.has(k)) return state.memo.get(k);
    state.nEval++;
    let o = { y: -Infinity, r: null };
    try { const r = solve(x); if (r && r.success && isFinite(r.COP) && !r.bad) o = { y: r.COP, r }; else o = { y: -Infinity, r }; } catch (e) {}
    state.memo.set(k, o); return o;
  };

  /* ------------------------------------------------------------------ n-level cascade optimiser
     Maximises the COP over the lower levels' condensing temperatures (the top level's is solved for the air
     approach). A design counts only if it keeps every air approach and every IHX closes (no such warnings).
     Starts from the current design and from the best point of a grid of evenly spaced lifts. */
  function* cascadeOptimiseSteps(fls, air, cfg, opt = {}) {
    const c = Object.assign({ T_air_inlet: 15, T_air_tar: 200, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, dTmin_rhx: 5 }, CASC, cfg);
    const n = c.n, dim = n - 1, maxEval = opt.maxEval ?? 400, state = { nEval: 0, memo: new Map() };
    if (fls.length < n) return fail("a refrigerant is needed for every level");
    if (dim < 1) return fail("a cascade needs at least two levels");
    const ev1 = (c.T_source_out ?? c.T_source_in) - (c.dT_src ?? c.dTmin_rhx), top = c.T_air_tar - c.dT_GC_min + 30;
    const tc = k => isFinite(fls[k].T_crit) ? fls[k].T_crit - 2 : top;
    const lo = [], hi = [];
    for (let k = 0; k < dim; k++) {
      lo.push(ev1 + 10 + k * (c.dT_casc + 10));
      hi.push(Math.min(tc(k), top - (dim - 1 - k) * (c.dT_casc + 10)));
    }
    if (lo.some((v, k) => v >= hi[k])) return fail("no room between the source and the air target for the lower levels");
    const bad = r => r.success && r.warnings.some(w => /approach|IHX/.test(w));
    const f = memoF(state, x => { const r = cascade(fls, air, { ...c, T_cond: x.map(v => +v) }); if (r.success) r.bad = bad(r); return r; });
    const clip = x => x.map((v, i) => Math.min(Math.max(v, lo[i]), hi[i]));
    const x0 = clip(Array.from({ length: dim }, (_, k) => c.T_cond[k] ?? lo[k])), b0 = f(x0), start = b0.r && b0.r.success ? b0.r.COP : NaN, startOk = isFinite(b0.y);
    yield { nEval: state.nEval, best: b0.y };
    // grid: level k condensing at a fraction of the way from the source to the air target, lifts spread evenly
    // (with a bias towards the bottom or top)
    let xg = null, bg = { y: -Infinity };
    const span = top - ev1, fr = dim === 1 ? [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8] : [0.8, 0.9, 1, 1.1, 1.25];
    for (const a of fr) {
      const z = clip(Array.from({ length: dim }, (_, k) => dim === 1 ? ev1 + a * span : ev1 + span * Math.pow((k + 1) / n, a)));
      const o = f(z); if (o.y > bg.y) { bg = o; xg = z; }
      yield { nEval: state.nEval, best: Math.max(b0.y, bg.y) };
    }
    const g = compassSearch(f, [[x0, b0], [xg, bg]], lo, hi, Array(dim).fill(8), Array(dim).fill(0.05), Array(dim).fill(false), state, maxEval);
    let st, overall = Math.max(b0.y, bg.y);
    while (!(st = g.next()).done) { overall = Math.max(overall, st.value); yield { nEval: state.nEval, best: overall }; }
    let { x, best, converged } = st.value;
    if (!x) return fail("no design within the bounds keeps every approach", { nEval: state.nEval });
    // round as the page shows them (0.1 K), nudging up until the rounded design still keeps its approaches
    let xr = x.map(v => Math.round(v * 10) / 10);
    for (let k = 0; k < 20; k++) { const o = f(xr); if (isFinite(o.y)) { x = xr; best = o; break; } xr = xr.map(v => Math.round((v + 0.1) * 10) / 10); }
    const near = (v, b) => Math.abs(v - b) < 0.05;
    return { success: true, COP: best.y, start, startOk, T_cond: x.map(v => +v), T_cond_top: best.r.report.T_cond_top, nEval: state.nEval, converged,
      atBound: x.map((v, i) => near(v, lo[i]) || near(v, hi[i])), lo, hi };
  }
  function cascadeOptimise(fls, air, cfg, opt = {}) { const g = cascadeOptimiseSteps(fls, air, cfg, opt); let s; do { s = g.next(); } while (!s.done); return s.value; }

  /* ------------------------------------------------------------------ booster optimiser
     Maximises the COP over the bottom condensing temperature and the number of MVR stages. A design counts only
     if both the bottom-cycle and the steam air approaches hold. Scans every stage count on a coarse condensing
     grid, then refines from the best of that and from the current design. */
  function* boosterOptimiseSteps(fb, water, air, cfg, opt = {}) {
    const c = Object.assign({ T_air_inlet: 15, T_air_tar: 200, dT_GC_min: 20, T_source_in: 40, T_source_out: 30, dTmin_rhx: 5 }, BOOST, cfg);
    const maxEval = opt.maxEval ?? 1200, nMax = opt.nMax ?? 20, state = { nEval: 0, memo: new Map() };
    const Te = (c.T_source_out ?? c.T_source_in) - (c.dT_src ?? c.dTmin_rhx);
    const lo = [Math.max(c.T_air_inlet + c.dT_GC_min + 2, Te + 10), 1];
    const hi = [Math.min(isFinite(fb.T_crit) ? fb.T_crit - 2 : c.T_air_tar, c.T_air_tar + c.dT_GC_min), nMax];
    if (lo[0] >= hi[0]) return fail("no room for the bottom condensing temperature");
    const bad = r => r.report.gap_steam < c.dT_GC_min - 0.05 || r.report.gap_bottom < c.dT_GC_min - 0.05;
    const f = memoF(state, x => { const r = booster(fb, water, air, { ...c, T_cond: x[0], n_mvr: x[1] }); if (r.success) r.bad = bad(r); return r; });
    const clip = x => [Math.min(Math.max(x[0], lo[0]), hi[0]), Math.round(Math.min(Math.max(x[1], lo[1]), hi[1]))];
    const x0 = clip([c.T_cond, c.n_mvr]), b0 = f(x0), start = b0.r && b0.r.success ? b0.r.COP : NaN, startOk = isFinite(b0.y);
    yield { nEval: state.nEval, best: b0.y };
    // every stage count on a 5 K condensing grid, then a 1-D refinement of the condensing temperature for the best
    // three stage counts and for the current design's (the feasible region shifts with the count, so the
    // count is not stepped by the search itself)
    const bestByN = new Map();
    for (let m = lo[1]; m <= hi[1]; m++) for (let T = lo[0]; T <= hi[0] + 1e-9; T += 5) {
      const z = clip([T, m]), o = f(z), cur = bestByN.get(m);
      if (!cur || o.y > cur.o.y) bestByN.set(m, { z, o });
      yield { nEval: state.nEval, best: Math.max(b0.y, ...[...bestByN.values()].map(v => v.o.y)) };
    }
    const picks = [...bestByN.values()].filter(v => isFinite(v.o.y)).sort((p, q) => q.o.y - p.o.y).slice(0, 3);
    if (isFinite(b0.y) && !picks.some(v => v.z[1] === x0[1])) picks.push({ z: x0, o: b0 });
    let X = null, Bst = { y: -Infinity }, conv = true, overall = Math.max(b0.y, ...picks.map(v => v.o.y));
    for (const pk of picks) {
      const m = pk.z[1], f1 = z => f([z[0], m]);
      const starts = [[[pk.z[0]], pk.o]];
      if (pk.z[1] === x0[1] && pk.z !== x0 && isFinite(b0.y)) starts.push([[x0[0]], b0]);
      const g = compassSearch(f1, starts, [lo[0]], [hi[0]], [4], [0.05], [false], state, maxEval);
      let st;
      while (!(st = g.next()).done) { overall = Math.max(overall, st.value); yield { nEval: state.nEval, best: overall }; }
      if (st.value.best.y > Bst.y) { Bst = st.value.best; X = [st.value.x[0], m]; conv = st.value.converged; }
    }
    const st = { value: { x: X, best: Bst, converged: conv } };
    let { x, best, converged } = st.value;
    if (!x) return fail("no design within the bounds keeps both air approaches", { nEval: state.nEval });
    let xr = [Math.round(x[0] * 10) / 10, x[1]];
    for (let k = 0; k < 20; k++) { const o = f(xr); if (isFinite(o.y)) { x = xr; best = o; break; } xr = [Math.round((xr[0] + 0.1) * 10) / 10, xr[1]]; }
    const near = (v, b) => Math.abs(v - b) < 0.05;
    return { success: true, COP: best.y, start, startOk, T_cond: x[0], n_mvr: x[1], splits: best.r.report.splits, nEval: state.nEval, converged,
      atBound: { T_cond: near(x[0], lo[0]) || near(x[0], hi[0]), n_mvr: x[1] === hi[1] }, lo, hi };
  }
  function boosterOptimise(fb, water, air, cfg, opt = {}) { const g = boosterOptimiseSteps(fb, water, air, cfg, opt); let s; do { s = g.next(); } while (!s.done); return s.value; }

  // ================================================================== sCO2 reverse Brayton (Dalal, Klöppel et al.)
  /* Supercritical CO2 heat pump with a split compression and a preheater, after the architectures of Dalal,
     Klöppel, Rose and Walmsley (Euro sCO2 2027, "Evaluating sCO2 heat pump architectures for large-scale industrial
     air drying above 150 °C"). The whole cycle is above the critical pressure: the suction sits at P_L (80 bar).
       Low side: turbine (expander) -> LTHX (heated by the source) -> IHX cold side -> split. A fraction x (the
       secondary stream) goes to compressor 1; its hot discharge preheats the main stream (1 - x) in the preheater
       before compressor 2.
       Architecture 1: compressor 2 heats the air in the HTHX; the secondary stream, after the preheater, rejoins
       the main stream after the HTHX; together they give heat back in the IHX hot side, then expand in the turbine.
       Architecture 2 (with Held's split stream): the secondary stream, after the preheater, heats the cooler part of
       the air in HTHX 2 and the main stream the hotter part in HTHX 1; only the main stream passes the IHX hot side;
       the two rejoin before the turbine.
     The recuperators (IHX, preheater) are set by an effectiveness: the duty is that fraction of the largest the
     exchanger could pass (found along the profiles, so the pinch can sit inside near the pseudo-critical point). The
     LTHX (gas heater) heats the CO2 as far as the source allows with a minimum approach dT_src, as the gas cooler
     keeps the air approach. The HTHX(s) cool the CO2 to the
     coldest outlet that keeps the air approach. The discharge pressure P_H is a design variable (the COP usually
     rises above the lowest feasible pressure, as the hotter discharge gives more heat per kg); with P_H_min it is the
     lowest at which compressor 2's discharge reaches the dryer target plus the approach. Pressure drops: dP_hp per exchanger on the high side,
     dP_lp on the low side. Shaft work: compressors / eta_m less turbine x eta_m; electrical: shaft / motor.
     Duties per kg of total CO2 flow; the reference flow for the swept volume is the total suction volume of both
     compressors. */
  const BRAY = { arch: 1, P_L: 8.0, P_H: 17.0, P_H_min: false, x: 0.3, eps_ihx: 0.85, eps_pre: 0.85, eps_lthx: 0.85, eta_c: 0.85, eta_t: 0.85, eta_m: 0.99, dP_hp: 0.05, dP_lp: 0.01, Pmax: 30 };
  function braytonAt(pr, air, c, P_H) {
    const x = Math.min(Math.max(+c.x, 0.01), 0.95), arch = +c.arch === 2 ? 2 : 1;
    const Tain = c.T_air_inlet, Ttar = c.T_air_tar, a = c.dT_GC_min ?? 20;
    const dPh = Math.max(+c.dP_hp || 0, 0), dPl = Math.max(+c.dP_lp || 0, 0), PL = c.P_L;
    const Pto = PL + 2 * dPl, Pl1 = PL + dPl, Pc2 = PL - dPl, Ph1 = P_H - dPh, Pti = P_H - 2 * dPh;
    // LTHX (gas heater): the CO2 is heated as far as the source allows while keeping the source approach dT_src
    // everywhere along the exchanger (counterflow; the source profile is a fixed shape in the duty fraction)
    const dTs = c.dT_src ?? 5;
    const cap = c.q_lthx_max != null && isFinite(c.q_lthx_max) ? Math.max(+c.q_lthx_max, 0) : Infinity;
    const lthxOut = ht => Math.min(lthxApproach(ht), ht + cap);
    const lthxApproach = ht => {
      const hMax = pr.PT(Pl1, c.T_source_in - dTs).h;
      if (!(hMax > ht)) return ht;
      const ok = h1 => { for (const f of lin(0, 1, 16)) if (Tph(pr, Pl1, ht + f * (h1 - ht)) > srcT(c, f) - dTs + 1e-6) return false; return true; };
      if (ok(hMax)) return hMax;
      if (!ok(ht + 1e-6)) return ht;
      let lo = ht, hi = hMax;
      for (let k = 0; k < 40 && hi - lo > 0.01; k++) { const m = (lo + hi) / 2; if (ok(m)) lo = m; else hi = m; }
      return lo;
    };
    // start: the suction at the source inlet temperature less 20 K
    let hci = pr.PT(PL, c.T_source_in - 20).h, S = null, d = Infinity;
    for (let it = 0; it < 200; it++) {
      const ci = pr.PH(PL, hci);
      const s1 = stage(pr, PL, ci, P_H, c.eta_c);
      // preheater: secondary (x, hot, after compressor 1) against the main stream (1 - x, cold) before compressor 2
      const mh = x / (1 - x), qpMax = maxDuty(pr, P_H, s1.out.h, mh, pr, Pc2, hci, 0).q, qp = c.eps_pre * qpMax;
      const hc2 = hci + qp, hpho = s1.out.h - qp / mh;
      const s2 = stage(pr, Pc2, pr.PH(Pc2, hc2), P_H, c.eta_c), h2 = s2.out.h, T2 = s2.out.T;
      // HTHX(s): the coldest outlets that keep the air approach. Below the target plus the approach there is no
      // design; the air is then taken only to the discharge less the approach, so the discharge temperature stays
      // continuous in the pressure and the pressure search brackets the target
      const okHot = T2 >= Ttar + a - 1e-9, Tt = Math.max(Math.min(Ttar, T2 - a), Tain + 1);
      const cool = (P, hIn, TaIn, TaOut, dT) => { try { return coolerOutlet(pr, air, P, hIn, TaIn, TaOut, dT); } catch (e) { if (okHot) throw e; return hIn; } };
      let h3, h5 = null, Tmid = null, hHin, mHot;
      if (arch === 1) {
        h3 = cool(P_H, h2, Tain, Tt, Math.min(a, T2 - Tt));
        hHin = (1 - x) * h3 + x * hpho; mHot = 1;
      } else {
        const Tpho = Tph(pr, Ph1, hpho);
        const sec = Tm => { const hh5 = cool(Ph1, hpho, Tain, Tm, a); return { h5: hh5, mAir: x * (hpho - hh5) / (air.h(Tm) - air.h(Tain)) }; };
        const main = (Tm, mAir) => (1 - x) * h2 - mAir * (air.h(Tt) - air.h(Tm));          // (1 - x) h3 from the balance
        if (Tpho - a > Tain + 0.5 && Tt > Tain + 1) {
          // the lowest air temperature between the two HTHXs at which the main stream still keeps its approach:
          // both exchangers then sit at their pinch (the most heat to the air)
          const aa = Math.min(a, T2 - Tt);
          const G = Tm => { const s = sec(Tm), h3b = main(Tm, s.mAir) / (1 - x); return h3b - cool(P_H, h2, Tm, Tt, aa); };
          let lo = Tain + 0.5, hi = Math.min(Tpho - a, Tt - 0.5);
          if (G(hi) < 0) { Tmid = hi; } else if (G(lo) >= 0) { Tmid = lo; } else {
            for (let k = 0; k < 40 && hi - lo > 0.02; k++) { const m = (lo + hi) / 2; if (G(m) >= 0) hi = m; else lo = m; }
            Tmid = hi;
          }
          const s = sec(Tmid); h5 = s.h5; h3 = main(Tmid, s.mAir) / (1 - x);
        } else { Tmid = Tain; h5 = hpho; h3 = cool(P_H, h2, Tain, Tt, Math.min(a, T2 - Tt)); }
        hHin = h3; mHot = 1 - x;
      }
      // IHX: the hot side (all the flow in architecture 1, the main stream in 2) heats the low side after the LTHX
      const Phot = arch === 1 ? Ph1 : Ph1;
      const h1 = S ? S.h1 : pr.PT(Pl1, c.T_source_in - 10).h;
      const qiMax = maxDuty(pr, Phot, hHin, mHot, pr, Pl1, h1, 0).q, qi = c.eps_ihx * qiMax;
      const hci2 = h1 + qi, hho = hHin - qi / mHot;
      const hti = arch === 1 ? hho : (1 - x) * hho + x * h5;
      const ti = pr.PH(Pti, hti), tis = pr.PS(Pto, ti.s).h, ht = hti - c.eta_t * (hti - tis);
      const h1n = lthxOut(ht);
      S = { ci, s1, s2, hc2, hpho, h2, T2, okHot, h3, h5, Tmid, hHin, mHot, h1: h1n, qi, hho, hti, ht, qp, x, arch, Phot };
      d = Math.abs(hci2 - hci) + Math.abs(h1n - h1);
      if (d < 1e-3 && it > 2) break;
      hci = it < 3 ? hci2 : hci + (it < 40 ? 0.5 : 0.2) * (hci2 - hci);        // damped: the recuperators feed back
    }
    S.hci = hci; S.P = { PL, Pto, Pl1, Pc2, P_H, Ph1, Pti }; S.conv = d < 1e-2;
    return S;
  }
  // the lowest discharge pressure at which compressor 2 reaches the dryer target plus the approach (bisection; the
  // discharge temperature is continuous in the pressure, see braytonAt)
  function braytonMinP(pr, air, c) {
    const a = c.dT_GC_min ?? 20, g = P => { const S = braytonAt(pr, air, c, P); return { S, v: S.T2 - (c.T_air_tar + a) }; };
    let lo = c.P_L * 1.15, hi = Math.max(c.Pmax || 30, lo * 1.01);
    const gl = g(lo), gh = g(hi);
    if (gh.v < 0) return gh.S;
    if (gl.v >= 0) return gl.S;
    let S = gh.S;
    for (let k = 0; k < 30 && hi - lo > 1e-3; k++) { const m = (lo + hi) / 2, gm = g(m); if (gm.v >= 0) { hi = m; S = gm.S; } else lo = m; }
    return S;
  }
  function brayton(fl, air, cfg) {
    const c = Object.assign({ T_air_inlet: 15, T_air_tar: 200, dT_GC_min: 20, T_source_in: 40, T_source_out: 30 }, BRAY, cfg);
    const warnings = [], pr = fastFluid(fl), a = c.dT_GC_min ?? 20;
    try {
      if (!(c.P_L > 0)) return fail("set the suction pressure");
      if (c.P_L < pr.P_crit) warnings.push(`suction ${(c.P_L * 10).toFixed(0)} bar is below the critical pressure (${(pr.P_crit * 10).toFixed(1)} bar): the low side runs subcritical and the turbine may exhaust two-phase`);
      let S;
      const srcCheck = S => {
        const Tto = Tph(pr, S.P.Pl1, S.ht), Tlim = srcT(c, 0) - (c.dT_src ?? 5);
        return Tto > Tlim + 0.05 ? fail(`LTHX: the CO₂ leaves the turbine at ${Tto.toFixed(1)} °C, above the ${c.source_name || "source water"} outlet (${srcT(c, 0).toFixed(1)} °C) less the ${c.dT_src ?? 5} K source approach, so the source cannot heat it: raise the source outlet temperature, lower the source approach, or lower the suction pressure`) : null;
      };
      if (c.P_H_min) {
        S = braytonMinP(pr, air, c);
        const sc = srcCheck(S); if (sc) return sc;
        if (!S.okHot) return fail(`compressor 2 reaches only ${S.T2.toFixed(0)} °C at ${(S.P.P_H * 10).toFixed(0)} bar, short of the ${c.T_air_tar} °C target plus the ${a} K approach: raise the max pressure or the preheater/IHX effectiveness`);
      } else {
        if (!(c.P_H > c.P_L * 1.05)) return fail("the discharge pressure must be above the suction pressure");
        S = braytonAt(pr, air, c, c.P_H);
        const sc = srcCheck(S); if (sc) return sc;
        if (!S.okHot) {
          const m = c._quick ? { okHot: false } : braytonMinP(pr, air, c);          // the optimiser skips the search
          return fail(`compressor 2 discharges at ${S.T2.toFixed(1)} °C, below the ${c.T_air_tar} °C target plus the ${a} K approach: raise the discharge pressure${m.okHot ? ` to at least ${(m.P.P_H * 10).toFixed(1)} bar` : ""}`, { P_H: c.P_H });
        }
      }
      const { x, arch } = S, P = S.P, Tain = c.T_air_inlet, Ttar = c.T_air_tar;
      // duties and work per kg of total CO2
      const qA1 = (1 - x) * (S.h2 - S.h3), qA2 = arch === 2 ? x * (S.hpho - S.h5) : 0, qAir = qA1 + qA2;
      const wc = x * S.s1.w + (1 - x) * S.s2.w, wt = S.hti - S.ht, wShaft = wc / c.eta_m - wt * c.eta_m;
      const motor = c.motor_eff ?? MOTOR, wEl = wShaft / motor, cop = qAir / wEl, qE = S.h1 - S.ht;
      if (!(wShaft > 0)) return fail("the turbine gives back more than the compressors take");
      // the loop has settled when the cycle's energy balance closes (source + compressors = air + turbine)
      const bal = qE + wc - wt - qAir;
      if (Math.abs(bal) > 0.005 * qAir + 0.2) return fail(`the cycle does not settle (energy balance off by ${bal.toFixed(1)} kJ/kg): the IHX and preheater feed heat back to the suction without a steady state at this pressure and split; raise the discharge pressure or lower the split`);
      // LTHX against the source (water line or air curve): no temperature cross
      const coldT = f => Tph(pr, P.Pl1, S.ht + f * qE);
      let gapE = Infinity; for (const f of lin(0, 1, 16)) gapE = Math.min(gapE, srcT(c, f) - coldT(f));
      const dTs = c.dT_src ?? 5;
      if (gapE < dTs - 0.1) return fail(`LTHX approach ${gapE.toFixed(1)} K, below the ${dTs} K source approach`);
      // gaps on the air
      const pA1 = airProfile(pr, air, P.P_H, S.h2, S.h3, arch === 2 ? S.Tmid : Tain, Ttar, 16), gap1 = minGap(pA1);
      const pA2 = arch === 2 && qA2 > 1e-9 ? airProfile(pr, air, P.Ph1, S.hpho, S.h5, Tain, S.Tmid, 16) : null, gap2 = pA2 ? minGap(pA2) : null;
      if (gap1 < a - 0.1 || (gap2 != null && gap2 < a - 0.1)) warnings.push(`air approach ${Math.min(gap1, gap2 ?? Infinity).toFixed(1)} K, below the ${a} K set`);
      const vSuc = x / S.ci.D + (1 - x) / pr.PH(P.Pc2, S.hc2).D;
      const st = { suction: S.ci, c1_out: S.s1.out, c2_in: pr.PH(P.Pc2, S.hc2), c2_out: S.s2.out, hthx_out: pr.PH(P.P_H, S.h3), pre_hot_out: pr.PH(P.Ph1, S.hpho),
        ihx_hot_in: pr.PH(S.Phot, S.hHin), ihx_hot_out: pr.PH(S.Phot, S.hho), turb_in: pr.PH(P.Pti, S.hti), turb_out: pr.PH(P.Pto, S.ht), lthx_out: pr.PH(P.Pl1, S.h1) };
      if (arch === 2) st.hthx2_out = pr.PH(P.Ph1, S.h5);
      const paths = [
        { name: "Turbine (expander)", pts: compPath(pr, P.Pti, S.hti, P.Pto, S.ht) },
        { name: "LTHX evaporator (source)", pts: iso(pr, P.Pto, S.ht, S.h1) },
        { name: "Internal heat exchanger (IHX)", pts: iso(pr, P.Pl1, S.h1, S.hci) },
        { name: "Compressor 1 (secondary)", pts: compPath(pr, P.PL, S.hci, P.P_H, S.s1.out.h) },
        { name: "Preheater recuperator", pts: iso(pr, P.Pc2, S.hci, S.hc2) },
        { name: "Compressor 2 (main)", pts: compPath(pr, P.Pc2, S.hc2, P.P_H, S.h2) },
        { name: arch === 2 ? "HTHX 1 gas cooler" : "HTHX gas cooler", pts: iso(pr, P.P_H, S.h2, S.h3) },
        { name: "Preheater recuperator", pts: iso(pr, P.P_H, S.s1.out.h, S.hpho) },
        ...(arch === 2 ? [{ name: "HTHX 2 gas cooler", pts: iso(pr, P.Ph1, S.hpho, S.h5) }] : [{ name: "Mixing", pts: [pt(st.pre_hot_out), pt(st.ihx_hot_in)] }, { name: "Mixing", pts: [pt(st.hthx_out), pt(st.ihx_hot_in)] }]),
        { name: "Internal heat exchanger (IHX)", pts: iso(pr, S.Phot, S.hHin, S.hho) },
        ...(arch === 2 ? [{ name: "Mixing", pts: [pt(st.ihx_hot_out), pt(st.turb_in)] }, { name: "Mixing", pts: [pt(st.hthx2_out), pt(st.turb_in)] }] : []),
      ];
      const loops = [{ name: "sCO2 cycle", fluid: 0, flow: 1, paths }];
      const qI = S.qi, qP = (1 - x) * S.qp;
      const hx = [
        { name: arch === 2 ? "HTHX 1: main stream to air" : "HTHX: CO2 to air", key: "gc", hot: pA1.map(p => [p[0] * (1 - x), p[1]]), cold: pA1.map(p => [p[0] * (1 - x), p[2]]), otherName: "dryer air" },
        ...(pA2 ? [{ name: "HTHX 2: secondary stream to air", key: "gc2", hot: pA2.map(p => [p[0] * x, p[1]]), cold: pA2.map(p => [p[0] * x, p[2]]), otherName: "dryer air" }] : []),
        { name: "LTHX (source)", key: "evap", hot: waterLine(qE, c), cold: lin(0, 1, 16).map(f => [f * qE, coldT(f)]), otherName: c.source_name || "source water", coldRef: true },
        { name: "IHX", key: "ihx", hot: lin(0, 1, 14).map(f => [f * qI, Tph(pr, S.Phot, S.hho + f * qI / S.mHot)]), cold: lin(0, 1, 14).map(f => [f * qI, Tph(pr, P.Pl1, S.h1 + f * qI)]), otherName: "CO2 to the compressors" },
        { name: "Preheater", key: "rhx", hot: lin(0, 1, 14).map(f => [f * qP, Tph(pr, P.Ph1, S.hpho + f * qP / x)]), cold: lin(0, 1, 14).map(f => [f * qP, Tph(pr, P.Pc2, S.hci + f * S.qp)]), otherName: "main stream before compressor 2" },
      ];
      return {
        success: true, message: "success", cycle: "X", warnings, COP: cop, VHC: qAir / vSuc / 1000, loops, hx,
        q_air: qAir, w_shaft: wShaft, q_source: qE, ref: { rho: 1 / vSuc, flow: 1, PR: P.P_H / P.PL },
        Pmax_bar: P.P_H * 10, T_dis_max: Math.max(S.s1.out.T, S.T2),
        air: arch === 2 ? [{ name: "HTHX 2", q: qA2, to: S.Tmid }, { name: "HTHX 1", q: qA1 }] : [{ name: "HTHX", q: qA1 }],
        report: {
          COP: cop, COP_shaft: qAir / wShaft, arch, x, P_L_bar: P.PL * 10, P_H_bar: P.P_H * 10, PR: P.P_H / P.PL,
          T_suction: S.ci.T, T_c1_out: S.s1.out.T, T_c2_in: pr.PH(P.Pc2, S.hc2).T, T_c2_out: S.T2, T_hthx_out: st.hthx_out.T, T_hthx2_out: arch === 2 ? st.hthx2_out.T : null,
          T_air_mid: arch === 2 ? S.Tmid : null, T_turb_in: st.turb_in.T, T_turb_out: st.turb_out.T, T_lthx_out: st.lthx_out.T,
          balance_kJkg: bal, w_c1: x * S.s1.w, w_c2: (1 - x) * S.s2.w, w_turb: wt, w_shaft: wShaft, w_el_kJkg: wEl, q_air_1: qA1, q_air_2: qA2, q_ihx: qI, q_pre: qP, q_lthx: qE,
          gap_hthx: gap1, gap_hthx2: gap2, gap_lthx: gapE, eps: { ihx: c.eps_ihx, pre: c.eps_pre }, dT_src: c.dT_src ?? 5, q_lthx_cap: c.q_lthx_max != null && isFinite(c.q_lthx_max) ? +c.q_lthx_max : null, eta_c: c.eta_c, eta_t: c.eta_t, eta_m: c.eta_m,
        },
        states: st,
      };
    } catch (e) { return fail(`${e && e.message ? e.message : e}`); }
  }
  /* optimiser: the COP over the suction pressure P_L, the discharge pressure P_H and the secondary fraction x.
     Inner: a compass search in scaled coordinates over P_H (from just above the suction to Pmax) and x (0.02-0.9) at
     a fixed P_L; infeasible points lose. Outer: a compass search on P_L (steps 10, 5, 2.5, 1.25 bar), each trial
     re-optimising P_H and x from the best design so far (P_H scaled with P_L) to a coarser step, then a finish at
     the chosen P_L. A solve budget (600, or 300 for the slower architecture 2) caps the run. P_L has no lower bound tied to the
     critical pressure (the low side may run subcritical; the model warns of a two-phase turbine exhaust): it runs
     from 6 bar to Pmax / 1.15, and designs that fail (no discharge high enough under Pmax, a temperature cross)
     bound it in practice. opt_PL: false keeps P_L fixed. A generator that yields {nEval, best} after each solve,
     like the others. */
  function* braytonPHx(fl, air, cfg, PL, v0, step0, cnt, budget, tol = [0.01, 0.004]) {
    const Pmax = cfg.Pmax || BRAY.Pmax || 30, lo = [PL * 1.1, 0.02], hi = [Pmax, 0.9];
    const f = v => { const r = brayton(fl, air, { ...cfg, P_L: PL, P_H_min: false, P_H: v[0], x: v[1], _quick: true }); cnt.n++; return r.success ? r.COP : -Infinity; };
    if (!(lo[0] < hi[0])) return { COP: -Infinity };
    let v = [Math.min(Math.max(v0[0], lo[0]), hi[0]), Math.min(Math.max(v0[1], lo[1]), hi[1])];
    let best = f(v); yield { nEval: cnt.n, best: Math.max(best, cnt.best) };
    const start = best;
    // an infeasible start (pressure too low): walk the pressure up until it solves
    while (!isFinite(best) && v[0] < hi[0] && cnt.n < budget) { v = [Math.min(v[0] * 1.08, hi[0]), v[1]]; best = f(v); yield { nEval: cnt.n, best: Math.max(best, cnt.best) }; }
    if (!isFinite(best)) return { COP: -Infinity, start };
    let step = step0.slice();
    while ((step[0] > tol[0] || step[1] > tol[1]) && cnt.n < budget) {
      let moved = false;
      for (let d = 0; d < 2 && !moved; d++) for (const sg of [1, -1]) {
        const w = v.slice(); w[d] = Math.min(Math.max(w[d] + sg * step[d], lo[d]), hi[d]);
        if (w[d] === v[d]) continue;
        const fw = f(w); yield { nEval: cnt.n, best: Math.max(best, cnt.best) };
        if (fw > best + 1e-6) { best = fw; v = w; moved = true; break; }
      }
      if (!moved) step = step.map(t => t / 2);
    }
    return { COP: best, P_H: v[0], x: v[1], start, atBound: { P_H: v[0] >= hi[0] - 1e-6 || v[0] <= lo[0] + 1e-6, x: v[1] <= lo[1] + 1e-6 || v[1] >= hi[1] - 1e-6 } };
  }
  function* braytonOptimiseSteps(fl, air, cfg) {
    const c = Object.assign({}, BRAY, cfg), Pmax = c.Pmax || 30, optPL = c.opt_PL !== false;
    const cnt = { n: 0, best: -Infinity }, budget = c.opt_budget || (c.arch === 2 ? 300 : 600);
    const PLlo = 0.6, PLhi = Pmax / 1.15;
    let PL = Math.min(Math.max(c.P_L, PLlo), PLhi);
    const v0 = [c.P_H_min ? (c.P_H || 17) : c.P_H, c.x];
    // the design as set, then P_H and x at the set suction pressure
    const r0 = brayton(fl, air, { ...cfg, P_L: PL, P_H_min: false, P_H: v0[0], x: v0[1], _quick: true }); cnt.n++;
    const start = r0.success ? r0.COP : -Infinity;
    let cur = yield* braytonPHx(fl, air, cfg, PL, v0, [1.0, 0.08], cnt, budget);
    // no design at the set suction pressure (the max pressure too low for it): step the suction pressure down
    for (let k = 0; optPL && !isFinite(cur.COP) && PL > PLlo + 1e-9 && k < 12 && cnt.n < budget; k++) {
      PL = Math.max(PL - 1.0, PLlo);
      cur = yield* braytonPHx(fl, air, cfg, PL, [Math.min(v0[0], Pmax), v0[1]], [1.0, 0.08], cnt, budget);
    }
    if (!isFinite(cur.COP)) return { success: false, message: `no design up to the max pressure${optPL ? ` at a suction pressure from ${(PL * 10).toFixed(0)} to ${(Math.min(Math.max(c.P_L, PLlo), PLhi) * 10).toFixed(0)} bar` : ""}` };
    cur.P_L = PL; cnt.best = cur.COP;
    if (optPL) {
      let s = 1.0, dir = 1;
      while (s >= 0.125 - 1e-9 && cnt.n < budget) {
        let moved = false;
        for (const sg of [dir, -dir]) {
          const PLn = Math.min(Math.max(cur.P_L + sg * s, PLlo), PLhi);
          if (Math.abs(PLn - cur.P_L) < 1e-9) continue;
          const tr = yield* braytonPHx(fl, air, cfg, PLn, [cur.P_H * PLn / cur.P_L, cur.x], [0.25, 0.02], cnt, budget, [0.04, 0.008]);
          if (tr.COP > cur.COP + 1e-5) { cur = { ...tr, P_L: PLn }; cnt.best = cur.COP; dir = sg; moved = true; break; }
        }
        if (!moved) s /= 2;
      }
      // the P_L trials stop P_H and x at a coarser step: finish them at the chosen P_L
      const fin = yield* braytonPHx(fl, air, cfg, cur.P_L, [cur.P_H, cur.x], [0.04, 0.008], cnt, budget + 80);
      if (fin.COP > cur.COP) cur = { ...fin, P_L: cur.P_L };
    }
    return { success: true, COP: cur.COP, P_L: cur.P_L, P_H: cur.P_H, x: cur.x, start, startOk: isFinite(start), nEval: cnt.n,
      atBound: { P_L: optPL && (cur.P_L <= PLlo + 1e-6 || cur.P_L >= PLhi - 1e-6), P_H: cur.atBound.P_H, x: cur.atBound.x }, budgetHit: cnt.n >= budget };
  }

  const api = { KONG, CASC, BOOST, kong, kongOptimise, kongOptimiseSteps, cascade, booster,
    cascadeOptimise, cascadeOptimiseSteps, boosterOptimise, boosterOptimiseSteps, BRAY, brayton, braytonOptimiseSteps };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.DRYCASC = api;
})(typeof self !== "undefined" ? self : this);

