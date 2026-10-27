/* Advanced dryer heat pumps: four more cycles, solved in the browser like dryer_cascades.js. Each is a port of the
   Python model in HeatPumpLab/cycles with the same equations and iterations, so the two agree to round-off.

   economiser    two-stage heat pump with an intermediate pressure: a flash tank, a closed economiser (vapour
                 injection) or parallel compression, or a single stage for comparison (EconomiserHeatPump). The high
                 pressure is the lowest that heats the air to the target with the approach; the intermediate pressure
                 is set by its saturation temperature T_int (null: the geometric mean of the evaporating and high
                 pressures).
   dualPressure  dual-pressure condensation after Dai et al. (2020): a low-pressure condenser heats the air from the
                 inlet to T_mid and a high-pressure one from T_mid to the target; the high-pressure liquid is
                 throttled into the low-pressure condenser (DualPressureHeatPump).
   stirling      reverse-Stirling units in series on the air, each a fraction of the Carnot COP between its hot side
                 (dT_hot above the air it delivers) and its cold side (dT_cold below the source outlet): a
                 semi-empirical block, not a physical Stirling model (StirlingHeatPump).
   hybrid        ammonia-water absorption-compression (Osenbrück) heat pump over the IAPWS G4-01 mixture model in
                 ammonia_water.js (HybridHeatPump).

   Results have the shape of DRYCASC's ({success, message, warnings, COP, VHC, report, loops, hx, q_air, w_shaft,
   q_source, ref, Pmax_bar, T_dis_max}); duties per kg of the basis flow named in each model (Stirling: per kg of
   primary air, with ref.air set because no compressor sets the size). Each optimiser is a generator that yields
   {nEval, best} after every solve and returns its outcome, like DRYCASC's. */
(function (root) {
  "use strict";
  const node = typeof module !== "undefined" && module.exports;
  const C_ = node ? require("./dryer_cascades.js") : root.DRYCASC;
  const AW = () => node ? require("./ammonia_water.js") : root.AMMWAT;
  const { stage, pt, iso, compPath, airProfile, minGap, coolerOutlet, lin, srcT, waterLine, fail, compassSearch } = C_._h;
  const MOTOR = 0.93;
  const GAS = new Set(["saturated_vapour", "superheated_vapour", "supercritical"]);
  const gasLike = (fl, st) => GAS.has(st.phase) || (isFinite(fl.P_crit) && st.P >= fl.P_crit && st.T >= fl.T_crit);
  const errMsg = e => `${e && e.message ? e.message : e}`;
  const jsRound = (x, d = 0) => Math.round(x * 10 ** d) / 10 ** d;
  const nameOf = fl => fl.label || fl.fluid || "the refrigerant";

  /* smallest suction superheat (K, at least shMin) whose compression from saturation at Pin (Tsat) to Pout
     discharges margin K above the dew point at Pout: overhanging fluids compress wet from near saturation */
  function drySuperheat(fl, Pin, Tsat, Pout, shMin, margin = 5) {
    const dew = isFinite(fl.P_crit) && Pout >= fl.P_crit ? -Infinity : fl.PQ(Pout, 1).T;
    const dry = sh => stage(fl, Pin, fl.PT(Pin, Tsat + sh), Pout).out.T - dew - margin;
    const sh = Math.max(shMin, 0.5);
    if (dry(sh) >= 0) return sh;
    let hi = sh;
    while (dry(hi) < 0 && hi < 120) hi += 10;
    if (dry(hi) < 0) throw new Error(`${nameOf(fl)} cannot compress dry from ${Tsat.toFixed(0)} C (raise the superheat or change the fluid)`);
    let lo = sh;
    for (let k = 0; k < 30 && hi - lo > 0.05; k++) { const m = (lo + hi) / 2; if (dry(m) >= 0) hi = m; else lo = m; }
    return hi;
  }
  /* the lowest high pressure (40 log steps from Plo to Pmax, then bisection in log P) whose design keeps the air
     approach: at(P) returns the design ({gap}) or throws */
  function lowestPressure(at, Plo, Pmax, dT, nScan = 40, tol = 1e-5, nBis = 30) {
    const tryAt = P => { try { return at(P); } catch (e) { return e instanceof Error ? e : new Error(errMsg(e)); } };
    let prev = Plo, found = null, last = null;
    for (let i = 0; i <= nScan; i++) {
      const P = Plo * Math.pow(Pmax / Plo, i / nScan), r = tryAt(P);
      if (!(r instanceof Error) && r.gap >= dT - 0.05) { found = [P, r]; break; }
      last = r; prev = P;
    }
    if (!found) return { last };
    let lo = prev, [hi, S] = found;
    for (let k = 0; k < nBis && hi / lo - 1 >= tol; k++) {
      const m = Math.sqrt(lo * hi), r = tryAt(m);
      if (!(r instanceof Error) && r.gap >= dT - 0.05) { hi = m; S = r; } else lo = m;
    }
    return { S };
  }
  const noPressure = (c, last) => `no high pressure up to ${c.P_max_bar} bar heats the air to ${c.T_air_tar} C with a ${c.dT_GC_min} K approach${last instanceof Error ? ` (${last.message})` : ""}`;
  // a one-variable optimiser: the start, a coarse grid and a compass search, rounded to 0.1 K
  function* optimise1(f, x0, lo, hi, grid, maxEval, state) {
    const b0 = x0 ? f(x0) : { y: -Infinity };
    yield { nEval: state.nEval, best: b0.y };
    let xg = null, bg = { y: -Infinity };
    for (const a of grid) {
      const z = [lo[0] + a * (hi[0] - lo[0])], o = f(z);
      if (o.y > bg.y) { bg = o; xg = z; }
      yield { nEval: state.nEval, best: Math.max(b0.y, bg.y) };
    }
    const g = compassSearch(f, [[x0, b0], [xg, bg]], lo, hi, [8], [0.05], [false], state, maxEval);
    let st, overall = Math.max(b0.y, bg.y);
    while (!(st = g.next()).done) { overall = Math.max(overall, st.value); yield { nEval: state.nEval, best: overall }; }
    let { x, best, converged } = st.value;
    if (!x) return { x: null, b0 };
    const xr = [jsRound(x[0], 1)], o = f(xr);
    if (isFinite(o.y)) { x = xr; best = o; }
    return { x, best, converged, b0 };
  }
  const memo1 = (state, solve) => x => {
    const k = (+x[0]).toFixed(4);
    if (state.memo.has(k)) return state.memo.get(k);
    state.nEval++;
    let r; try { r = solve(+x[0]); } catch (e) { r = fail(errMsg(e)); }
    const o = { y: r.success ? r.COP : -Infinity, r };
    state.memo.set(k, o); return o;
  };
  const air0 = { T_air_inlet: 15, T_air_tar: 170, dT_GC_min: 20, T_source_in: 60, T_source_out: 50, dT_src: 5 };
  const evapT = c => (c.T_source_out ?? c.T_source_in) - c.dT_src;
  const prWarn = (prs, dis, warnings) => {
    if (Math.max(...prs) > 20) warnings.push(`stage pressure ratio ${Math.max(...prs).toFixed(0)}, beyond the compressor efficiency correlations (fitted up to 20): treat the efficiency, and the COP, with caution`);
    if (Math.max(...dis) > 250) warnings.push(`discharge at ${Math.max(...dis).toFixed(0)} C: check the refrigerant's and the oil's stability limits`);
  };
  const evapHx = (fl, c, Pe, he, hs, m, qE) => ({ name: "Evaporator", key: "evap", hot: waterLine(qE, c), coldRef: true,
    cold: iso(fl, Pe, he, hs).map(p => [m * (p[2] - he), p[0]]), otherName: c.source_name || "source water" });

  // ================================================================== economiser / flash tank / parallel compression
  /* Per kg of the high-pressure flow. flash: the condenser liquid flashes at P_i; the flash vapour (x) joins the
     low-stage discharge at the high-stage suction and the saturated liquid expands to the evaporator. economiser: a
     fraction y of the condenser liquid expands to P_i and evaporates (superheat_inj above saturation) against the
     rest, which is subcooled to dT_econ above the saturation temperature at P_i; the vapour joins the low-stage
     discharge. parallel: as flash, but the flash vapour (heated by the liquid if a dry discharge needs it) has its
     own compressor to the high pressure, and the low stage compresses straight to it. single: one stage. */
  const ECON = { arrangement: "flash", T_int: null, dT_econ: 5, superheat_inj: 5, superheat: 5, P_max_bar: 60 };
  const ARR = ["flash", "economiser", "parallel", "single"];
  function econCycle(fl, air, c, P_H, Pe, Te, P_i) {
    const arr = c.arrangement, two = arr !== "single";
    if (two && !(Pe < P_i && P_i < P_H)) throw new Error("the intermediate pressure must lie between the evaporating and high pressures");
    const loOut = arr === "flash" || arr === "economiser" ? P_i : P_H;
    const sh = drySuperheat(fl, Pe, Te, loOut, Math.max(c.superheat, 0.5));
    const suc = fl.PT(Pe, Te + sh), low = stage(fl, Pe, fl.PT(Pe, Te + sh), loOut);
    let lq = null, vq = null, hInj = null, par = null;
    if (two) {
      lq = fl.PQ(P_i, 0); vq = fl.PQ(P_i, 1);
      // the parallel compressor's suction is heated (superheat_inj, more if a dry discharge needs it) by the liquid
      const shInj = arr === "parallel" ? drySuperheat(fl, P_i, vq.T, P_H, c.superheat_inj) : arr === "economiser" ? Math.max(c.superheat_inj, 0) : 0;
      hInj = shInj > 0 ? fl.PT(P_i, vq.T + shInj).h : vq.h;
      if (arr === "parallel") par = stage(fl, P_i, fl.PH(P_i, hInj), P_H);
    }
    // the condenser inlet, its outlet and the injected fraction that outlet implies, at an injected fraction v
    const evaluate = v => {
      let hi = null, hIn;
      if (arr === "single") hIn = low.out.h;
      else if (arr === "flash" || arr === "economiser") {
        const mix = fl.PH(P_i, (1 - v) * low.out.h + v * hInj);
        if (!gasLike(fl, mix)) throw new Error("the high-stage suction is wet (raise the superheat or the intermediate pressure)");
        hi = stage(fl, P_i, mix, P_H); hIn = hi.out.h;
      } else hIn = (1 - v) * low.out.h + v * par.out.h;
      const h3 = coolerOutlet(fl, air, P_H, hIn, c.T_air_inlet, c.T_air_tar, c.dT_GC_min);
      let vNew;
      if (arr === "single") vNew = 0;
      else if (arr === "economiser") { const hSub = Math.min(h3, fl.PT(P_H, vq.T + c.dT_econ).h); vNew = h3 > hSub ? (h3 - hSub) / (hInj - hSub) : 0; }
      else { const dh = hInj - vq.h; vNew = Math.min(Math.max((h3 - lq.h) / (vq.h - lq.h + dh), 0), 1); }
      return { hIn, h3, hi, vNew };
    };
    // iterate directly; where that does not settle, the smallest self-consistent fraction by a scan and bisection
    let v = 0, E = evaluate(0), settled = !two || E.vNew < 1e-12;
    for (let k = 0; k < 25 && !settled; k++) {
      const vn = E.vNew, En = evaluate(vn);
      if (Math.abs(En.vNew - vn) < 1e-8) settled = true;
      v = vn; E = En;
    }
    if (!settled) {
      const g = u => evaluate(u).vNew - u;
      let lo = 0, hi = null;
      for (let k = 1; k <= 40; k++) { const u = k / 40; if (g(u) <= 0) { hi = u; break; } lo = u; }
      if (hi == null) throw new Error("no self-consistent injected fraction");
      for (let k = 0; k < 40 && hi - lo >= 1e-8; k++) { const m = (lo + hi) / 2; if (g(m) > 0) lo = m; else hi = m; }
      v = hi; E = evaluate(v);
      if (Math.abs(E.vNew - v) < 1e-5) v = E.vNew;
      E = evaluate(v);
    }
    const { hIn, h3, hi } = E;
    const x = arr === "flash" || arr === "parallel" ? v : 0, y = arr === "economiser" ? v : 0;
    let mLow, hValve;
    if (arr === "single") { mLow = 1; hValve = h3; }
    else if (arr === "flash" || arr === "parallel") { mLow = 1 - x; hValve = x > 0 ? lq.h : h3; }
    else { mLow = 1 - y; hValve = y > 0 ? h3 - y / (1 - y) * (hInj - h3) : h3; }
    const he = hValve, qAir = hIn - h3;
    const w = mLow * low.w + (hi ? hi.w : 0) + (par ? x * par.w : 0), qSrc = mLow * (suc.h - he);
    const gap = minGap(airProfile(fl, air, P_H, hIn, h3, c.T_air_inlet, c.T_air_tar, 18));
    return { P_H, P_i, Pe, Te, sh, suc, low, hi, par, x, y, hInj, hIn, h3, hValve, he, mLow, qAir, w, qSrc, gap };
  }
  function economiser(fl, air, cfg) {
    const c = Object.assign({}, air0, ECON, cfg), warnings = [];
    if (!ARR.includes(c.arrangement)) return fail(`arrangement must be one of ${ARR.join(", ")}`);
    try {
      const Te = evapT(c), Pe = fl.TQ(Te, 1).P, Pmax = c.P_max_bar / 10;
      const Tcc = isFinite(fl.T_crit) ? fl.T_crit : Infinity;
      // floor of the pressure scan: condensing well below the air target, and below the critical point
      const T0 = Math.min(Math.max(Te + 10, c.T_air_tar - 60), Tcc - 1);
      const Plo = T0 > Te + 5 ? fl.TQ(T0, 1).P : Pe * 2;
      if (Plo >= Pmax) return fail(`${nameOf(fl)} needs more than ${c.P_max_bar} bar to reach ${c.T_air_tar} C`);
      const Pint = P_H => c.arrangement === "single" ? null : c.T_int != null ? fl.TQ(c.T_int, 1).P : Math.sqrt(Pe * P_H);
      const { S, last } = lowestPressure(P_H => econCycle(fl, air, c, P_H, Pe, Te, Pint(P_H)), Plo, Pmax, c.dT_GC_min);
      if (!S) return fail(noPressure(c, last));
      return econResult(fl, air, c, S, warnings);
    } catch (e) { return fail(errMsg(e), { warnings }); }
  }
  function econResult(fl, air, c, S, warnings) {
    const wEl = S.w / (c.motor_eff ?? MOTOR), cop = S.qAir / wEl, { P_H, P_i, Pe } = S, arr = c.arrangement;
    const supercritical = isFinite(fl.P_crit) && P_H >= fl.P_crit;
    const hiIn = S.hi ? S.hi.out.h - S.hi.w : null;
    let vs = S.mLow / S.suc.D;
    if (S.hi) vs += 1 / fl.PH(P_i, hiIn).D;
    if (S.par) vs += S.x / fl.PH(P_i, S.hInj).D;
    const VHC = S.qAir / vs / 1000;
    const dis = [S.low.out.T, ...(S.hi ? [S.hi.out.T] : []), ...(S.par ? [S.par.out.T] : [])];
    const prs = [S.low.pr, ...(S.hi ? [S.hi.pr] : []), ...(S.par ? [S.par.pr] : [])];
    prWarn(prs, dis, warnings);
    const gcName = "Condenser" + (supercritical ? " / gas cooler" : "");
    const fromHigh = arr === "single" || arr === "economiser";   // the main liquid expands from the high pressure
    const paths = [
      ...(arr === "flash" || arr === "parallel" ? [{ name: "Flash valve", pts: compPath(fl, P_H, S.h3, P_i, S.h3) }] : []),
      { name: "Expansion valve", pts: compPath(fl, fromHigh ? P_H : P_i, S.hValve, Pe, S.he) },
      { name: "Evaporator", pts: iso(fl, Pe, S.he, S.suc.h) },
      { name: "Low-stage compressor", pts: compPath(fl, Pe, S.suc.h, S.low.out.P, S.low.out.h) },
      ...(S.hi ? [{ name: "High-stage compressor", pts: compPath(fl, P_i, hiIn, P_H, S.hi.out.h) }] : []),
      ...(S.par ? [{ name: "Parallel compressor", pts: compPath(fl, P_i, S.hInj, P_H, S.par.out.h) }] : []),
      ...(arr === "economiser" && S.y > 0 ? [{ name: "Injection valve", pts: compPath(fl, P_H, S.h3, P_i, S.h3) },
        { name: "Economiser (injected side)", pts: iso(fl, P_i, S.h3, S.hInj) }, { name: "Economiser (main liquid)", pts: iso(fl, P_H, S.h3, S.hValve) }] : []),
      { name: gcName, pts: iso(fl, P_H, S.hIn, S.h3) },
    ];
    const prof = airProfile(fl, air, P_H, S.hIn, S.h3, c.T_air_inlet, c.T_air_tar, 18);
    const hx = [{ name: `${gcName} to air`, key: "gc", hot: prof.map(p => [p[0], p[1]]), cold: prof.map(p => [p[0], p[2]]), otherName: "dryer air" },
      evapHx(fl, c, Pe, S.he, S.suc.h, S.mLow, S.qSrc)];
    if (arr === "economiser" && S.y > 0) {
      const q = S.y * (S.hInj - S.h3);           // per kg high-pressure flow
      hx.push({ name: "Economiser", key: "ihx", hot: lin(0, 1, 14).map(f => [f * q, fl.PH(P_H, S.hValve + f * (S.h3 - S.hValve)).T]),
        cold: lin(0, 1, 14).map(f => [f * q, fl.PH(P_i, S.h3 + f * (S.hInj - S.h3)).T]), otherName: "injected vapour" });
    }
    const report = { COP: cop, VHC, arrangement: arr, fluid: nameOf(fl), P_high_bar: P_H * 10, P_int_bar: P_i != null ? P_i * 10 : null, P_evap_bar: Pe * 10,
      T_int: P_i != null ? fl.PQ(P_i, 1).T : null, T_evap: S.Te, superheat: S.sh, supercritical,
      flows: { high: 1, low: S.mLow, injected: arr === "flash" || arr === "parallel" ? S.x : S.y }, x_flash: S.x, y_econ: S.y, w_el_kJkg: wEl,
      q_air: S.qAir, q_evap: S.qSrc, gap_air: S.gap, T_dis_low: S.low.out.T, T_dis_high: S.hi ? S.hi.out.T : null, T_dis_parallel: S.par ? S.par.out.T : null,
      PR_low: S.low.pr, PR_high: S.hi ? S.hi.pr : S.par ? S.par.pr : null,
      air_per_ref: S.qAir / (air.h(c.T_air_tar) - air.h(c.T_air_inlet)), air_W: air.W, balance: S.qAir - S.w - S.qSrc };
    return { success: true, message: "success", cycle: "E", warnings, COP: cop, VHC, loops: [{ name: "Cycle", fluid: 0, flow: 1, paths }], hx,
      q_air: S.qAir, w_shaft: S.w, q_source: S.qSrc, ref: { rho: 1 / vs, flow: 1, PR: Math.max(...prs) },
      Pmax_bar: P_H * 10, T_dis_max: Math.max(...dis), report };
  }
  // best-COP intermediate saturation temperature
  function* economiserOptimiseSteps(fl, air, cfg, opt = {}) {
    const c = Object.assign({}, air0, ECON, cfg), maxEval = opt.maxEval ?? 120, state = { nEval: 0, memo: new Map() };
    if (c.arrangement === "single") { const r = economiser(fl, air, c); return r.success ? { success: true, COP: r.COP, start: r.COP, startOk: true, T_int: null, nEval: 1, atBound: false } : r; }
    const Te = evapT(c), top = Math.min(c.T_air_tar, isFinite(fl.T_crit) ? fl.T_crit - 2 : c.T_air_tar);
    const lo = [Te + 5], hi = [top - 5];
    if (lo[0] >= hi[0]) return fail("no room for an intermediate pressure between the evaporator and the condenser");
    const f = memo1(state, T => economiser(fl, air, { ...c, T_int: T }));
    const base = economiser(fl, air, { ...c, T_int: null });
    const x0 = base.success ? [Math.min(Math.max(base.report.T_int, lo[0]), hi[0])] : null;
    const o = yield* optimise1(f, x0, lo, hi, [0.2, 0.35, 0.5, 0.65, 0.8], maxEval, state);
    if (!o.x) return fail("no intermediate pressure gives a feasible design", { nEval: state.nEval });
    const near = v => Math.abs(v - lo[0]) < 0.05 || Math.abs(v - hi[0]) < 0.05;
    return { success: true, COP: o.best.y, start: base.success ? base.COP : NaN, startOk: base.success, T_int: o.x[0], nEval: state.nEval, converged: o.converged, atBound: near(o.x[0]), lo, hi };
  }

  // ================================================================== dual-pressure condensation (Dai et al.)
  /* Per kg of evaporator flow: evaporator -> low-stage compressor (Pe -> P_L) -> split: 1 - a to the low-pressure
     condenser, a through the high-stage compressor to the high-pressure condenser and a valve back to P_L, where it
     joins the low-pressure condenser inlet; all the flow leaves that condenser for the evaporator valve. The air
     passes the low-pressure condenser (inlet -> T_mid), then the high-pressure one (T_mid -> target). a follows
     from the air balance, T_mid is the highest the low-pressure condenser allows and P_H the lowest that heats the
     air to the target. */
  const DUAL = { T_cond_low: 110, superheat: 5, P_max_bar: 60 };
  function dualAt(fl, air, c, P_H, Pe, Te, P_L, low) {
    const { sh, suc, lo } = low, h1 = lo.out.h, hi = stage(fl, P_L, lo.out, P_H), h2 = hi.out.h;
    const Tain = c.T_air_inlet, Ttar = c.T_air_tar, haIn = air.h(Tain), haTar = air.h(Ttar), TdewL = fl.PQ(P_L, 1).T;
    const sections = Tmid => {
      const h3H = coolerOutlet(fl, air, P_H, h2, Tmid, Ttar, c.dT_GC_min), qH = h2 - h3H, hm = air.h(Tmid);
      let a = 0.5;
      for (let k = 0; k < 60; k++) {
        const hInL = (1 - a) * h1 + a * h3H, h3L = coolerOutlet(fl, air, P_L, hInL, Tain, Tmid, c.dT_GC_min), qL = hInL - h3L;
        const aNew = qL * (haTar - hm) / (qH * (hm - haIn));          // the same air flow through both condensers
        if (!(aNew > 0 && aNew < 1)) throw new Error("no split of the flow balances the two condensers");
        if (Math.abs(aNew - a) < 1e-9) { a = aNew; break; }
        a = a + 0.8 * (aNew - a);
      }
      const hInL = (1 - a) * h1 + a * h3H, h3L = coolerOutlet(fl, air, P_L, hInL, Tain, Tmid, c.dT_GC_min);
      return { a, h3H, qH, hInL, h3L, qL: hInL - h3L };
    };
    const top = TdewL - c.dT_GC_min;
    if (top <= Tain + 1) throw new Error(`the low condensing temperature (${TdewL.toFixed(0)} C) is too close to the air inlet`);
    let best = null, loT = Tain + 1, hiT = null;
    for (let k = 0; k <= 20; k++) {
      const T = top - (top - Tain - 1) * k / 20;
      try { best = [T, sections(T)]; loT = T; break; } catch (e) { hiT = T; }
    }
    if (!best) throw new Error("the low-pressure condenser cannot take any of the air heating");
    if (hiT != null) for (let k = 0; k < 25 && hiT - loT >= 0.02; k++) {
      const m = (loT + hiT) / 2;
      try { best = [m, sections(m)]; loT = m; } catch (e) { hiT = m; }
    }
    const [Tmid, S] = best, a = S.a, he = S.h3L;
    const gapH = minGap(airProfile(fl, air, P_H, h2, S.h3H, Tmid, Ttar, 18)), gapL = minGap(airProfile(fl, air, P_L, S.hInL, S.h3L, Tain, Tmid, 18));
    return { P_H, P_L, Pe, Te, sh, suc, lo, hi, Tmid, h1, h2, he, qAir: a * S.qH + S.qL, w: lo.w + a * hi.w, qSrc: suc.h - he,
      gap: Math.min(gapH, gapL), gapH, gapL, ...S };
  }
  function dualPressure(fl, air, cfg) {
    const c = Object.assign({}, air0, DUAL, cfg), warnings = [];
    try {
      const Te = evapT(c), Pe = fl.TQ(Te, 1).P;
      if (isFinite(fl.T_crit) && c.T_cond_low >= fl.T_crit - 1) return fail(`the low condensing temperature must be below ${nameOf(fl)}'s critical temperature`);
      if (c.T_cond_low <= Te + 5) return fail("the low condensing temperature must be above the evaporating temperature");
      const P_L = fl.TQ(c.T_cond_low, 1).P;
      const sh = drySuperheat(fl, Pe, Te, P_L, c.superheat), suc = fl.PT(Pe, Te + sh);
      const low = { sh, suc, lo: stage(fl, Pe, suc, P_L) };
      const Pmax = c.P_max_bar / 10, Plo = P_L * 1.05;
      if (Plo >= Pmax) return fail(`${nameOf(fl)} needs more than ${c.P_max_bar} bar`);
      const { S, last } = lowestPressure(P_H => dualAt(fl, air, c, P_H, Pe, Te, P_L, low), Plo, Pmax, c.dT_GC_min);
      if (!S) return fail(noPressure(c, last));
      return dualResult(fl, air, c, S, warnings);
    } catch (e) { return fail(errMsg(e), { warnings }); }
  }
  function dualResult(fl, air, c, S, warnings) {
    const wEl = S.w / (c.motor_eff ?? MOTOR), cop = S.qAir / wEl, { P_H, P_L, Pe, a } = S;
    const supercritical = isFinite(fl.P_crit) && P_H >= fl.P_crit;
    const prs = [S.lo.pr, S.hi.pr], dis = [S.lo.out.T, S.hi.out.T];
    prWarn(prs, dis, warnings);
    const vs = 1 / S.suc.D + a / S.lo.out.D, VHC = S.qAir / vs / 1000;
    const hpName = "High-pressure condenser" + (supercritical ? " / gas cooler" : "");
    const paths = [
      { name: "Expansion valve", pts: compPath(fl, P_L, S.h3L, Pe, S.he) },
      { name: "Evaporator", pts: iso(fl, Pe, S.he, S.suc.h) },
      { name: "Low-stage compressor", pts: compPath(fl, Pe, S.suc.h, P_L, S.h1) },
      { name: "High-stage compressor", pts: compPath(fl, P_L, S.h1, P_H, S.h2) },
      { name: hpName, pts: iso(fl, P_H, S.h2, S.h3H) },
      { name: "High-pressure valve", pts: compPath(fl, P_H, S.h3H, P_L, S.h3H) },
      { name: "Low-pressure condenser", pts: iso(fl, P_L, S.hInL, S.h3L) },
    ];
    const pH = airProfile(fl, air, P_H, S.h2, S.h3H, S.Tmid, c.T_air_tar, 18), pL = airProfile(fl, air, P_L, S.hInL, S.h3L, c.T_air_inlet, S.Tmid, 18);
    const hx = [
      { name: "Low-pressure condenser to air", key: "gc2", hot: pL.map(p => [p[0], p[1]]), cold: pL.map(p => [p[0], p[2]]), otherName: "dryer air" },
      { name: `${hpName} to air`, key: "gc", hot: pH.map(p => [a * p[0], p[1]]), cold: pH.map(p => [a * p[0], p[2]]), otherName: "dryer air" },
      evapHx(fl, c, Pe, S.he, S.suc.h, 1, S.qSrc),
    ];
    const report = { COP: cop, VHC, fluid: nameOf(fl), P_high_bar: P_H * 10, P_low_cond_bar: P_L * 10, P_evap_bar: Pe * 10, T_cond_low: c.T_cond_low,
      T_evap: S.Te, T_mid: S.Tmid, split_high: a, superheat: S.sh, supercritical, w_el_kJkg: wEl, q_air: S.qAir, q_air_low: S.qL, q_air_high: a * S.qH,
      q_evap: S.qSrc, gap_air: S.gap, gap_low: S.gapL, gap_high: S.gapH, T_dis_low: dis[0], T_dis_high: dis[1], PR_low: prs[0], PR_high: prs[1],
      air_per_ref: S.qAir / (air.h(c.T_air_tar) - air.h(c.T_air_inlet)), air_W: air.W, balance: S.qAir - S.w - S.qSrc };
    return { success: true, message: "success", cycle: "D", warnings, COP: cop, VHC, loops: [{ name: "Cycle", fluid: 0, flow: 1, paths }], hx,
      q_air: S.qAir, w_shaft: S.w, q_source: S.qSrc, ref: { rho: 1 / vs, flow: 1, PR: Math.max(...prs) }, Pmax_bar: P_H * 10, T_dis_max: Math.max(...dis),
      air: [{ name: "low-pressure condenser", q: S.qL, to: S.Tmid }, { name: "high-pressure", q: a * S.qH }], report };
  }
  // best-COP low condensing temperature
  function* dualPressureOptimiseSteps(fl, air, cfg, opt = {}) {
    const c = Object.assign({}, air0, DUAL, cfg), maxEval = opt.maxEval ?? 120, state = { nEval: 0, memo: new Map() };
    const Te = evapT(c), top = Math.min(c.T_air_tar, isFinite(fl.T_crit) ? fl.T_crit - 2 : c.T_air_tar);
    const lo = [Math.max(Te + 10, c.T_air_inlet + c.dT_GC_min + 5)], hi = [top - 5];
    if (lo[0] >= hi[0]) return fail("no room for a low condensing temperature between the source and the air target");
    const f = memo1(state, T => dualPressure(fl, air, { ...c, T_cond_low: T }));
    const x0 = [Math.min(Math.max(c.T_cond_low, lo[0]), hi[0])];
    const o = yield* optimise1(f, x0, lo, hi, [0.2, 0.35, 0.5, 0.65, 0.8], maxEval, state);
    if (!o.x) return fail("no low condensing temperature gives a feasible design", { nEval: state.nEval });
    const near = v => Math.abs(v - lo[0]) < 0.05 || Math.abs(v - hi[0]) < 0.05;
    return { success: true, COP: o.best.y, start: o.b0.y, startOk: isFinite(o.b0.y), T_cond_low: o.x[0], nEval: state.nEval, converged: o.converged, atBound: near(o.x[0]), lo, hi };
  }

  // ================================================================== reverse-Stirling units
  /* COP_unit = eta_carnot * T_hot / (T_hot - T_cold) (K), eta_carnot a constant (0.37 from the RFC Maasdam case:
     COP 2 at about 428 K / 348 K) or eta_of_ratio(T_hot / T_cold). n_units in series on the air heat bands between
     the boundaries (equal steps when null); duties per kg of primary air. */
  const STIR = { n_units: 1, boundaries: null, eta_carnot: 0.37, eta_of_ratio: null, dT_hot: 20, dT_cold: 5 };
  const DEMONSTRATED_T_HOT = 200;
  const stirBounds = c => c.boundaries != null ? c.boundaries.map(Number) : Array.from({ length: c.n_units - 1 }, (_, k) => c.T_air_inlet + (c.T_air_tar - c.T_air_inlet) * (k + 1) / c.n_units);
  function stirling(air, cfg) {
    const c = Object.assign({}, air0, STIR, cfg), warnings = [], n = c.n_units;
    if (!(n >= 1 && n <= 4)) return fail("1 to 4 Stirling units in series");
    const bounds = stirBounds(c);
    if (bounds.length !== n - 1) return fail(`${n} units need ${n - 1} boundary temperatures`);
    const T = [c.T_air_inlet, ...bounds, c.T_air_tar];
    if (T.some((t, i) => i > 0 && t <= T[i - 1])) return fail("the unit boundaries must rise from the air inlet to the target");
    const Tcold = (c.T_source_out ?? c.T_source_in) - c.dT_cold;
    const units = []; let qAir = 0, W = 0;
    for (let k = 0; k < n; k++) {
      const Thot = T[k + 1] + c.dT_hot, Th = Thot + 273.15, Tc = Tcold + 273.15;
      if (Th <= Tc) return fail(`unit ${k + 1}: its hot side (${Thot.toFixed(0)} C) is not above its cold side (${Tcold.toFixed(0)} C)`);
      const eta = typeof c.eta_of_ratio === "function" ? c.eta_of_ratio(Th / Tc) : c.eta_carnot;
      if (!(eta > 0 && eta <= 1)) return fail(`unit ${k + 1}: Carnot fraction ${eta} outside (0, 1]`);
      const carnot = Th / (Th - Tc), cop = eta * carnot, q = air.h(T[k + 1]) - air.h(T[k]), w = q / cop;
      qAir += q; W += w;
      units.push({ T_air_in: T[k], T_air_out: T[k + 1], T_hot: Thot, T_cold: Tcold, carnot_COP: carnot, eta_carnot: eta, COP: cop, q_air: q, w_el: w, q_source: q - w, lift: Thot - Tcold });
      if (Thot > DEMONSTRATED_T_HOT) warnings.push(`unit ${k + 1} hot side at ${Thot.toFixed(0)} C, above the about ${DEMONSTRATED_T_HOT} C of the published field data`);
    }
    const cop = qAir / W;
    const report = { COP: cop, units, n_units: n, boundaries: bounds, T_cold: Tcold, w_el_kJkg: W, q_air: qAir, q_source: qAir - W, air_W: air.W,
      eta_model: typeof c.eta_of_ratio === "function" ? "eta_of_ratio" : `constant ${c.eta_carnot}`, balance: 0 };
    const hx = units.map((u, k) => ({ name: `Unit ${k + 1} to air`, key: ["gc", "gc2", "gc3", "gc4"][k], hot: [[0, u.T_hot], [u.q_air, u.T_hot]],
      cold: [[0, u.T_air_in], [u.q_air, u.T_air_out]], otherName: "dryer air" }));
    hx.push({ name: "Cold side (source)", key: "evap", hot: waterLine(qAir - W, c), cold: [[0, Tcold], [qAir - W, Tcold]], otherName: c.source_name || "source water", coldRef: true });
    return { success: true, message: "success", cycle: "R", warnings, COP: cop, VHC: NaN, loops: [], hx, q_air: qAir, w_shaft: W, q_source: qAir - W,
      ref: { air: true, rho: NaN, flow: 1, PR: 1 }, Pmax_bar: NaN, T_dis_max: Math.max(...units.map(u => u.T_hot)),
      air: units.map((u, k) => ({ name: `unit ${k + 1}`, q: u.q_air, to: u.T_air_out })), report };
  }
  // best-COP air temperatures between the units
  function* stirlingOptimiseSteps(air, cfg, opt = {}) {
    const c = Object.assign({}, air0, STIR, cfg), n = c.n_units, maxEval = opt.maxEval ?? 200, state = { nEval: 0 };
    if (n <= 1) { const r = stirling(air, c); return r.success ? { success: true, COP: r.COP, start: r.COP, startOk: true, boundaries: [], nEval: 1, atBound: false } : r; }
    const lo = Array(n - 1).fill(c.T_air_inlet + 1), hi = Array(n - 1).fill(c.T_air_tar - 1);
    const f = x => {
      state.nEval++;
      if (x.some((b, i) => i > 0 && b <= x[i - 1])) return { y: -Infinity, r: null };
      const r = stirling(air, { ...c, boundaries: x.map(Number) });
      return { y: r.success ? r.COP : -Infinity, r };
    };
    const x0 = stirBounds(c), b0 = f(x0);
    yield { nEval: state.nEval, best: b0.y };
    const g = compassSearch(f, [[x0, b0]], lo, hi, Array(n - 1).fill(10), Array(n - 1).fill(0.05), Array(n - 1).fill(false), state, maxEval);
    let st; while (!(st = g.next()).done) yield { nEval: state.nEval, best: st.value };
    let { x, best, converged } = st.value;
    if (!x) return fail("no set of unit boundaries gives a feasible design");
    const xr = x.map(v => jsRound(v, 1)), o = f(xr);
    if (isFinite(o.y)) { x = xr; best = o; }
    return { success: true, COP: best.y, start: b0.y, startOk: isFinite(b0.y), boundaries: x.map(Number), nEval: state.nEval, converged, atBound: false };
  }

  // ================================================================== ammonia-water hybrid (Osenbrück)
  /* Per kg of solution through the absorber. Desorber (P_L): the throttled solution is heated by the source to
     T_des_out = T_source_in - dT_src, where q_des has boiled off; the separator sends the ammonia-rich vapour to the
     compressor and the lean liquid to the pump, less any part injected at the compressor suction (wet compression,
     to hold the discharge at T_dis_max); both are absorbed in the absorber (P_H), which heats the air; the absorber
     outlet is the coldest that keeps the approach, P_H the lowest that heats the air to the target. The solution keeps
     one composition through the absorber and the desorber, which behave as a condenser and an evaporator with a
     large glide (Jensen et al. 2015). Properties: ammonia_water.js. */
  const HYB = { w_rich: 0.5, q_des: 0.4, T_air_tar: 120, dT_GC_min: 10, T_source_out: 40, P_max_bar: 50, T_dis_max: 170,
    wet_compression: true, limit_discharge: false, eta_pump: 0.6 };
  const STREAMS = new Map();
  const awStream = w => { const k = w.toFixed(12); let s = STREAMS.get(k); if (!s) { if (STREAMS.size > 400) STREAMS.clear(); s = AW().stream(w); STREAMS.set(k, s); } return s; };
  // T(h) from the stream's isobar table, as the Python model's fast lookups (airProfile and coolerOutlet use PH(..).T)
  const fastAW = R => ({ PT: (P, T) => R.PT(P, T), PH: (P, h) => ({ T: R.T_ph_fast(P, h) }) });
  function hybrid(air, cfg) {
    const c = Object.assign({}, air0, HYB, cfg), warnings = [];
    if (!(c.w_rich > 0 && c.w_rich < 1) || !(c.q_des > 0 && c.q_des < 1)) return fail("w_rich and q_des must lie between 0 and 1");
    try {
      const R = awStream(c.w_rich), Rf = fastAW(R), Tdo = c.T_source_in - c.dT_src;
      const desOut = R.TQ(Tdo, c.q_des), P_L = desOut.P;            // low pressure: q_des boiled off at Tdo
      const V = awStream(desOut.w_vap), L = awStream(desOut.w_liq);
      const vap = V.PQ(P_L, 1), liq = L.PQ(P_L, 0);
      if (Math.abs(vap.T - Tdo) > 0.05 || Math.abs(liq.T - Tdo) > 0.05) throw new Error("the separator phases are not at the desorber outlet temperature");
      const Pmax = c.P_max_bar / 10;
      if (P_L * 1.2 >= Pmax) return fail(`the low pressure (${(P_L * 10).toFixed(1)} bar) leaves no room below ${c.P_max_bar} bar`);
      const wv = desOut.w_vap, wl = desOut.w_liq;
      // the separator vapour (q_des kg) compressed with m kg of lean liquid injected at the suction
      const compression = (P_H, m) => {
        if (m <= 0) return [stage(V, P_L, vap, P_H), 0];
        const n = c.q_des + m, Cm = awStream((c.q_des * wv + m * wl) / n), suc = Cm.PH(P_L, (c.q_des * vap.h + m * liq.h) / n);
        return [stage(Cm, P_L, suc, P_H), m];
      };
      const at = P_H => {
        let [comp, mInj] = compression(P_H, 0);
        if (c.wet_compression && comp.out.T > c.T_dis_max) {
          const mHi = 1 - c.q_des - 1e-6;
          if (compression(P_H, mHi)[0].out.T > c.T_dis_max) throw new Error(`even injecting all the lean liquid cannot hold the discharge at ${c.T_dis_max} C`);
          // regula falsi (Illinois) on the discharge temperature, which falls smoothly with m
          let a = 0, b = mHi, ga = comp.out.T - c.T_dis_max, gb = compression(P_H, b)[0].out.T - c.T_dis_max, side = 0, best = b, m = b, gm = gb;
          for (let k = 0; k < 40; k++) {
            m = b - gb * (b - a) / (gb - ga); gm = compression(P_H, m)[0].out.T - c.T_dis_max;
            if (gm <= 0) best = m;
            if (Math.abs(gm) < 0.02) break;
            if ((gm > 0) === (ga > 0)) { a = m; ga = gm; if (side === -1) gb /= 2; side = -1; }
            else { b = m; gb = gm; if (side === 1) ga /= 2; side = 1; }
          }
          [comp, mInj] = compression(P_H, Math.abs(gm) < 0.02 ? m : best);
        }
        const mc = c.q_des + mInj, mp = 1 - mc, wPump = (P_H - P_L) * 1e3 / liq.D / c.eta_pump, hPump = liq.h + wPump;
        const hIn = mc * comp.out.h + mp * hPump;
        const hOut = coolerOutlet(Rf, air, P_H, hIn, c.T_air_inlet, c.T_air_tar, c.dT_GC_min);
        const gap = minGap(airProfile(Rf, air, P_H, hIn, hOut, c.T_air_inlet, c.T_air_tar, 18));
        return { P_H, comp, mInj, mc, mp, wPump, hPump, hIn, hOut, gap };
      };
      const { S, last } = lowestPressure(at, P_L * 1.2, Pmax, c.dT_GC_min, 24, 2e-4, 20);
      if (!S) return fail(noPressure(c, last));
      const P_H = S.P_H, hDi = S.hOut, qSrc = desOut.h - hDi;
      const pts = lin(0, 1, 16).map(f => [f * qSrc, R.T_ph_fast(P_L, hDi + f * qSrc), srcT(c, f)]);
      const gapSrc = Math.min(...pts.map(p => p[2] - p[1]));
      if (gapSrc < c.dT_src - 0.05) return fail(`the desorber comes within ${gapSrc.toFixed(1)} K of the source water (needs ${c.dT_src} K): raise the degassing or the ammonia fraction, or lower the source outlet`);
      const w = S.mc * S.comp.w + S.mp * S.wPump, qAir = S.hIn - S.hOut, wEl = w / (c.motor_eff ?? MOTOR), cop = qAir / wEl, Tdis = S.comp.out.T;
      if (Tdis > c.T_dis_max + 0.05) {
        const msg = `compressor discharge at ${Tdis.toFixed(0)} C, above the ${c.T_dis_max} C limit`;
        if (c.limit_discharge) return fail(msg);
        warnings.push(msg);
      }
      return hybridResult(air, c, S, desOut, vap, liq, P_L, qAir, w, wEl, cop, qSrc, gapSrc, pts, warnings);
    } catch (e) { return fail(errMsg(e), { warnings }); }
  }
  function hybridResult(air, c, S, desOut, vap, liq, P_L, qAir, w, wEl, cop, qSrc, gapSrc, srcPts, warnings) {
    const R = awStream(c.w_rich), P_H = S.P_H, comp = S.comp;
    const absIn = R.PH(P_H, S.hIn), absOut = R.PH(P_H, S.hOut), desIn = R.PH(P_L, S.hOut);
    const [Tb, Td] = R.bubble_dew(P_H);
    const prof = airProfile(fastAW(R), air, P_H, S.hIn, S.hOut, c.T_air_inlet, c.T_air_tar, 18);
    const vs = S.mInj > 0 ? S.mc / R.PH(P_L, comp.out.h - comp.w).D : c.q_des / vap.D;
    const VHC = qAir / vs / 1000;
    // the diagrams follow the solution (rich composition) through the absorber, valve and desorber
    const sucSt = S.mInj > 0 ? awStream((c.q_des * desOut.w_vap + S.mInj * desOut.w_liq) / S.mc) : awStream(desOut.w_vap);
    const paths = [
      { name: "Absorber (to air)", pts: iso(R, P_H, S.hIn, S.hOut, 10) },
      { name: "Expansion valve", pts: [pt(absOut), pt(desIn)] },
      { name: "Desorber (source)", pts: iso(R, P_L, S.hOut, desOut.h, 10) },
      { name: "Compressor", pts: [pt(sucSt.PH(P_L, comp.out.h - comp.w)), pt(comp.out)] },
      { name: "Solution pump", pts: [pt(liq), pt(Object.assign({}, liq, { P: P_H, h: S.hPump }))] },
    ];
    const hx = [
      { name: "Absorber to air", key: "gc", hot: prof.map(p => [p[0], p[1]]), cold: prof.map(p => [p[0], p[2]]), otherName: "dryer air" },
      { name: "Desorber", key: "evap", hot: srcPts.map(p => [p[0], p[2]]), cold: srcPts.map(p => [p[0], p[1]]), otherName: c.source_name || "source water", coldRef: true },
    ];
    const report = { COP: cop, VHC, w_rich: c.w_rich, q_des: c.q_des, w_vapour: desOut.w_vap, w_lean: desOut.w_liq, P_low_bar: P_L * 10, P_high_bar: P_H * 10,
      PR: P_H / P_L, T_des_out: desOut.T, T_des_in: desIn.T, T_abs_in: absIn.T, T_abs_out: absOut.T, T_bubble_high: Tb - 273.15, T_dew_high: Td - 273.15,
      T_dis: comp.out.T, eta_comp: comp.eta, w_comp: S.mc * comp.w, w_pump: S.mp * S.wPump, m_injected: S.mInj, wet_compression: S.mInj > 0, w_el_kJkg: wEl,
      q_air: qAir, q_source: qSrc, gap_air: S.gap, gap_source: gapSrc, air_per_solution: qAir / (air.h(c.T_air_tar) - air.h(c.T_air_inlet)), air_W: air.W,
      balance: qAir - w - qSrc, properties: "IAPWS G4-01 (Tillner-Roth and Friend), ammonia_water.js" };
    return { success: true, message: "success", cycle: "Y", warnings, COP: cop, VHC, loops: [{ name: "Solution", fluid: 0, flow: 1, paths }], hx,
      q_air: qAir, w_shaft: w, q_source: qSrc, ref: { rho: 1 / vs, flow: 1, PR: P_H / P_L }, Pmax_bar: P_H * 10, T_dis_max: comp.out.T, report };
  }
  // best-COP solution composition and degassing
  function* hybridOptimiseSteps(air, cfg, opt = {}) {
    const c = Object.assign({}, air0, HYB, cfg), maxEval = opt.maxEval ?? 80, state = { nEval: 0 }, memo = new Map();
    const lo = [0.2, 0.05], hi = [0.9, 0.9];
    const f = x => {
      const k = `${(+x[0]).toFixed(4)},${(+x[1]).toFixed(4)}`;
      if (memo.has(k)) return memo.get(k);
      state.nEval++;
      const r = hybrid(air, { ...c, w_rich: +x[0], q_des: +x[1] }), ok = r.success && !(c.limit_discharge && r.T_dis_max > c.T_dis_max);
      const o = { y: ok ? r.COP : -Infinity, r }; memo.set(k, o); return o;
    };
    const x0 = [Math.min(Math.max(c.w_rich, lo[0]), hi[0]), Math.min(Math.max(c.q_des, lo[1]), hi[1])], b0 = f(x0), starts = [[x0, b0]];
    yield { nEval: state.nEval, best: b0.y };
    if (!isFinite(b0.y)) for (const w of [0.35, 0.5, 0.65]) for (const q of [0.2, 0.4, 0.6]) {
      const o = f([w, q]); if (isFinite(o.y)) starts.push([[w, q], o]);
      yield { nEval: state.nEval, best: Math.max(...starts.map(s => s[1].y)) };
    }
    const g = compassSearch(f, starts, lo, hi, [0.08, 0.08], [0.002, 0.002], [false, false], state, maxEval);
    let st; while (!(st = g.next()).done) yield { nEval: state.nEval, best: st.value };
    let { x, best, converged } = st.value;
    if (!x) return fail("no solution composition and degassing gives a feasible design", { nEval: state.nEval });
    const xr = [jsRound(x[0], 3), jsRound(x[1], 3)], o = f(xr);
    if (isFinite(o.y)) { x = xr; best = o; }
    const near = (v, i) => Math.abs(v - lo[i]) < 1e-3 || Math.abs(v - hi[i]) < 1e-3;
    return { success: true, COP: best.y, start: b0.y, startOk: isFinite(b0.y), w_rich: +x[0], q_des: +x[1], nEval: state.nEval, converged,
      atBound: { w_rich: near(x[0], 0), q_des: near(x[1], 1) } };
  }

  const api = { ECON, DUAL, STIR, DEMONSTRATED_T_HOT, drySuperheat, economiser, economiserOptimiseSteps, dualPressure, dualPressureOptimiseSteps,
    stirling, stirlingOptimiseSteps, HYB, hybrid, hybridOptimiseSteps };
  if (node) module.exports = api; else root.DRYMORE = api;
})(typeof self !== "undefined" ? self : this);
