/* Screening range for the refrigerant charge of a solved design, without a heat exchanger design: the JavaScript port of
   HeatPumpLab/cycles/common/charge.py, which holds the method, its assumptions and its limits. Each exchanger is split
   into equal-duty segments, dA = dQ / (U dT) and dm = dA (V/A) rho, with U and V/A assumed by service. Specific charge in
   g per kW of heating, which does not depend on the capacity. */
(function (root) {
  "use strict";
  // case: service: [U W/m2/K, V/A m]. Air: U on the finned area. Recuperator: V/A of each refrigerant side.
  const GEOMETRY = {
    compact: { air: [60, 0.05e-3], liquid: [2500, 0.5e-3], recuperator: [800, 0.5e-3] },
    typical: { air: [40, 0.15e-3], liquid: [1500, 1.0e-3], recuperator: [500, 1.0e-3] },
    bulky: { air: [25, 0.30e-3], liquid: [800, 3.0e-3], recuperator: [250, 3.0e-3] },
  };
  const EXCHANGER_SHARE = { low: 2 / 3, central: 0.55, high: 0.435 };
  // bound: [geometry case, density index in a section's rho pairs: 0 homogeneous, 1 void fraction]
  const BOUNDS = { low: ["compact", 0], central: ["typical", 1], high: ["bulky", 1] };
  const N = 24;
  const mids = () => Array.from({ length: N }, (_, i) => (i + 0.5) / N);

  function voidDensity(rho, rl, rv) {
    const x = Math.min(Math.max((1 / rho - 1 / rl) / (1 / rv - 1 / rl), 0), 1);
    if (x === 0) return rl;
    if (x === 1) return rv;
    const a = 1 / (1 + (1 - x) / x * Math.pow(rv / rl, 2 / 3));
    return a * rv + (1 - a) * rl;
  }
  function pair(fl, P, h, memo) {
    const st = fl.PH(P, h), rho = st.D;
    if (st.phase !== "two_phase") return [rho, rho];
    const sat = q => { const k = `${P}:${q}`; if (!(k in memo)) memo[k] = fl.PQ(P, q).D; return memo[k]; };
    return [rho, voidDensity(rho, sat(0), sat(1))];
  }
  function interp(pts, q) {
    if (q <= pts[0][0]) return pts[0][1];
    for (let i = 1; i < pts.length; i++) {
      if (q <= pts[i][0]) { const [a, Ta] = pts[i - 1], [b, Tb] = pts[i]; return b > a ? Ta + (Tb - Ta) * (q - a) / (b - a) : Tb; }
    }
    return pts[pts.length - 1][1];
  }

  /* sections: [{name, service, duty, dT: [K], rho: [side: [[homogeneous, void fraction]]]}] */
  function estimate(sections, Qheat, mixture) {
    if (!(Qheat > 0)) throw new Error("charge estimate needs a positive heating duty");
    const totals = { low: 0, central: 0, high: 0 }, by = {};
    for (const s of sections) {
      if (!s.dT.length || Math.min(...s.dT) <= 0) throw new Error(`${s.name}: charge estimate needs a positive driving force`);
      const dQ = s.duty / s.dT.length, row = {};
      for (const b of Object.keys(BOUNDS)) {
        const [U, va] = GEOMETRY[BOUNDS[b][0]][s.service];
        let kg = 0;
        for (const side of s.rho) s.dT.forEach((d, i) => { kg += 1000 * dQ / (U * d) * va * side[i][BOUNDS[b][1]]; });
        row[b] = 1000 * kg / Qheat;
        totals[b] += row[b];
      }
      by[s.name] = row;
    }
    const sys = b => totals[b] / EXCHANGER_SHARE[b];
    return { exchanger_low: totals.low, exchanger_central: totals.central, exchanger_high: totals.high,
      system_low: sys("low"), system_central: sys("central"), system_high: sys("high"), by_exchanger: by, mixture: !!mixture };
  }

  /* a multi-loop result: each exchanger takes the path of the same name, as a gas cooler where supercritical. The
     economiser exchanger of the economiser arrangement is not counted (nor in the Python model). */
  function loopCharge(fl, r) {
    const paths = {};
    for (const L of r.loops) for (const p of L.paths) paths[p.name] = p;
    const memo = {}, secs = [];
    for (const ex of r.hx) {
      const name = ex.name.replace(/ to air$/, ""), path = paths[name] || paths[`${name} / gas cooler`];
      if (!path) continue;
      const hs = path.pts.map(p => p[2]), hC = Math.min(...hs), hH = Math.max(...hs), P = path.pts[0][1] / 10;
      const duty = ex.hot[ex.hot.length - 1][0], fs = mids();
      secs.push({ name: ex.name, service: /air/.test(ex.otherName) ? "air" : "liquid", duty,
        dT: fs.map(f => interp(ex.hot, f * duty) - interp(ex.cold, f * duty)),
        rho: [fs.map(f => pair(fl, P, hC + f * (hH - hC), memo))] });
    }
    return estimate(secs, r.q_air, fl.isMix);
  }

  /* an RHP result (RHP-A/B/C and the conventional cycle), per kg of refrigerant: the gas coolers against the air, the
     evaporator against the source (``src``: T_source_in, T_source_out, source_curve) and each recuperator's two sides */
  function rhpCharge(fl, air, r, src) {
    const st = r.states, p = r.report, memo = {}, fs = mids(), secs = [], R = p.air_per_ref_kgkg;
    const P_H = r.P_H, P_I = r.P_I, P_L = r.P_L;
    r.gcs.forEach((g, k) => {
      if (!(g.q > 1e-6)) return;
      secs.push({ name: g.name || `Gas cooler ${k + 1}`, service: "air", duty: g.q,
        dT: fs.map(f => fl.PH(P_H, g.hout + f * g.q).T - air.T(g.ha_in + f * g.q / R)),
        rho: [fs.map(f => pair(fl, P_H, g.hout + f * g.q, memo))] });
    });
    const To = src.T_source_out == null ? src.T_source_in : src.T_source_out;
    const srcT = f => src.source_curve ? src.source_curve(f) : To + f * (src.T_source_in - To);
    const hE = st.evap_in.h, hS = st.superheat_out.h;
    secs.push({ name: "Evaporator", service: src.source_curve ? "air" : "liquid", duty: hS - hE,
      dT: fs.map(f => srcT(f) - fl.PH(P_L, hE + f * (hS - hE)).T), rho: [fs.map(f => pair(fl, P_L, hE + f * (hS - hE), memo))] });
    // counterflow from the cold end: the hot side's outlet faces the cold side's inlet
    const recup = (name, Ph, hOut, hIn, Pc, cIn, cOut) => {
      const duty = cOut - cIn;
      if (!(duty > 1e-6)) return;
      secs.push({ name, service: "recuperator", duty,
        dT: fs.map(f => fl.PH(Ph, hOut + f * (hIn - hOut)).T - fl.PH(Pc, cIn + f * duty).T),
        rho: [fs.map(f => pair(fl, Ph, hOut + f * (hIn - hOut), memo)), fs.map(f => pair(fl, Pc, cIn + f * duty, memo))] });
    };
    const cIn = (st.ihx_cold_out || st.superheat_out).h;
    if (st.gc2_out) recup("RHX1", P_H, st.rhx1_hot_out.h, st.gc2_out.h, P_L, cIn, st.suction.h);
    if (st.gc3_out) recup("RHX2", P_H, st.rhx2_hot_out.h, st.gc3_out.h, P_I, st.stage1_out.h, st.stage2_in.h);
    if (st.ihx_cold_out) recup("IHX", P_H, st.valve_in.h, st.gc1_out.h, P_L, st.superheat_out.h, st.ihx_cold_out.h);
    return estimate(secs, p.q_air_kJkg, fl.isMix);
  }

  const api = { GEOMETRY, EXCHANGER_SHARE, voidDensity, estimate, loopCharge, rhpCharge };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.DRYCHARGE = api;
})(typeof globalThis !== "undefined" ? globalThis : this);
