/* Ammonia-water mixture properties for the hybrid absorption-compression heat pump, in plain JavaScript: a port of
   HeatPumpLab/utils/ammonia_water.py (AmmoniaWater, AmmoniaWaterStream), which evaluates the same model through teqp.

   Model     IAPWS Guideline G4-01 (2001), Tillner-Roth and Friend, J. Phys. Chem. Ref. Data 27 (1998) 63, exactly as in
             teqp 0.23.2 (AmmoniaWaterTillnerRoth): reduced Helmholtz energy at tau = Tred(x)/T, delta = rho/rhored(x)
             of x * ammonia (Tillner-Roth 1993) + (1 - x) * water (IAPWS-95) + the departure function; x is the ammonia
             mole fraction, R = 8.314471 J/mol/K. Analytic tau, delta and x derivatives (chain rule through Tred(x),
             rhored(x)), and the second delta derivative for dp/drho.
   Ideal gas CoolProp's pure-fluid ideal parts as teqp converts them (IdealHelmholtz), so h and s carry CoolProp's
             reference states, as in Python: alpha_ig = sum x_i (alpha_ig,i(T, rho) + ln x_i).
   Algorithms as in Python: densities by Newton on p(rho) with the same root checks; equilibrium at (T, p) by Newton
             on the ln K equations (Python: teqp's mix_VLE_Tp, same equations and solution), damped successive
             substitution as the fallback; bubble/dew points by secant on T or by pressure substitution; Brent's method
             (scipy's brentq, ported) where Python uses brentq. Starting values use the pure-fluid vapour-pressure
             ancillaries (Python: CoolProp saturation tables) - converged results do not depend on them.

   Units of the model object: T K, p Pa, rho mol/m3, h J/mol, s J/mol/K. Streams (fixed ammonia mass fraction w):
   T C, P MPa, h kJ/kg, s kJ/kg/K, D kg/m3; a state is {T, P, h, s, D, q, phase, w_liq, w_vap[, x_mol, y_mol]} with
   q = -1 outside the two-phase region. Methods take (P, x) like the cycles' fluids and throw Error on failure. */
(function (root) {
  "use strict";
  const M_NH3 = 17.03026e-3, M_H2O = 18.015268e-3;         // kg/mol (IAPWS G4-01, Table 1)
  const R = 8.314471, X_MIN = 1e-10, T_K = 273.15;

  const mole_fraction = w => { const a = w / M_NH3; return a / (a + (1 - w) / M_H2O); };
  const mass_fraction = x => x * M_NH3 / (x * M_NH3 + (1 - x) * M_H2O);
  const molarMass = x => x * M_NH3 + (1 - x) * M_H2O;
  const clampX = x => Math.min(Math.max(x, X_MIN), 1 - X_MIN);

  // ------------------------------------------------------------------ coefficients (teqp's hard-coded JSON)
  const F64 = a => Float64Array.from(a);
  const power = (n, t, d, l) => ({ n: F64(n), t: F64(t), d: F64(d), l: Int32Array.from(l), N: n.length });
  const NH3 = power(   // Tillner-Roth et al. (1993): n delta^d tau^t exp(-delta^l), no exponential where l = 0
    [0.04554431, 0.7238548, 0.0122947, -1.858814, 2.141882e-11, -0.0143002, 0.3441324, -0.2873571, 2.352589e-05,
      -0.03497111, 0.02397852, 0.001831117, -0.04085375, 0.2379275, -0.03548972, -0.1823729, 0.02281556, -0.006663444,
      -0.008847486, 0.002272635, -0.0005588655],
    [-0.5, 0.5, 1, 1.5, 3, 0, 3, 4, 4, 5, 3, 5, 6, 8, 8, 10, 10, 5, 7.5, 15, 30],
    [2, 1, 4, 1, 15, 3, 3, 1, 8, 2, 1, 8, 1, 2, 3, 2, 4, 3, 1, 2, 4],
    [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3]);
  const H2O = power(   // IAPWS-95 (Wagner and Pruss 2002), polynomial and exponential terms
    [0.012533547935523, 7.8957634722828, -8.7803203303561, 0.31802509345418, -0.26145533859358, -0.0078199751687981,
      0.0088089493102134, -0.66856572307965, 0.20433810950965, -6.6212605039687e-05, -0.19232721156002,
      -0.25709043003438, 0.16074868486251, -0.040092828925807, 3.9343422603254e-07, -7.5941377088144e-06,
      0.00056250979351888, -1.5608652257135e-05, 1.1537996422951e-09, 3.6582165144204e-07, -1.3251180074668e-12,
      -6.2639586912454e-10, -0.10793600908932, 0.017611491008752, 0.22132295167546, -0.40247669763528,
      0.58083399985759, 0.0049969146990806, -0.031358700712549, -0.74315929710341, 0.4780732991548,
      0.020527940895948, -0.13636435110343, 0.014180634400617, 0.0083326504880713, -0.029052336009585,
      0.038615085574206, -0.020393486513704, -0.0016554050063734, 0.0019955571979541, 0.00015870308324157,
      -1.638856834253e-05, 0.043613615723811, 0.034994005463765, -0.076788197844621, 0.022446277332006,
      -6.2689710414685e-05, -5.5711118565645e-10, -0.19905718354408, 0.31777497330738, -0.11841182425981],
    [-0.5, 0.875, 1, 0.5, 0.75, 0.375, 1, 4, 6, 12, 1, 5, 4, 2, 13, 9, 3, 4, 11, 4, 13, 1, 7, 1, 9, 10, 10, 3, 7, 10,
      10, 6, 10, 10, 1, 2, 3, 4, 8, 6, 9, 8, 16, 22, 23, 23, 10, 50, 44, 46, 50],
    [1, 1, 1, 2, 2, 3, 4, 1, 1, 1, 2, 2, 3, 4, 4, 5, 7, 9, 10, 11, 13, 15, 1, 2, 2, 2, 3, 4, 4, 4, 5, 6, 6, 7, 9, 9, 9,
      9, 9, 10, 10, 12, 3, 4, 4, 5, 14, 3, 6, 6, 6],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
      2, 2, 2, 2, 3, 3, 3, 3, 4, 6, 6, 6, 6]);
  // IAPWS-95 Gaussian terms: n delta^d tau^t exp(-eta (delta - eps)^2 - beta (tau - gamma)^2)
  const GA = { n: [-31.306260323435, 31.546140237781, -2521.3154341695], t: [0, 1, 4], d: [3, 3, 3],
    eta: [20, 20, 20], eps: [1, 1, 1], beta: [150, 150, 250], gam: [1.21, 1.21, 1.25] };
  // IAPWS-95 non-analytic terms: n Delta^b delta Psi
  const NA = { n: [-0.14874640856724, 0.31806110878444], A: [0.32, 0.32], B: [0.2, 0.2], C: [28, 32], D: [700, 800],
    a: [3.5, 3.5], b: [0.85, 0.95], beta: [0.3, 0.3] };
  // departure function (terms 1..14 of aw.hpp; the first has no exponential)
  const DEP = power(
    [-1.855822E-02, 5.258010E-02, 3.552874E-10, 5.451379E-06, -5.998546E-13, -3.687808E-06, 0.2586192, -1.368072E-08,
      1.226146E-02, -7.181443E-02, 9.970849E-02, 1.0584086E-03, -0.1963687, -0.7777897],
    [1.5, 0.5, 6.5, 1.75, 15, 6, -1, 4, 3.5, 0, -1, 8, 7.5, 4],
    [4, 5, 15, 12, 12, 15, 4, 15, 4, 5, 6, 10, 6, 2],
    [0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 2, 2, 2]);
  const DEP_GAMMA = 0.5248379;
  const TcN = 405.40, TcW = 647.096, kT = 0.9648407, ALPHA = 1.125455;
  const vcN = 0.01703026 / 225, vcW = 0.018015268 / 322, kV = 1.2395117, BETA = 0.8978069;
  // ideal gas, teqp-converted from CoolProp: ln rho + a1 + a2/T + c + aT ln T + sum n ln(1 - exp(-theta/T))
  const IG = [
    { a1: -16.11892003782038, a2: 2271.5462320021797, c: 18.015806484253137, aT: -3.0,
      n: [2.224, 3.148, 0.9579], th: [1646.0000000000002, 3965.0, 7231.0] },
    { a1: -18.11153368769503, a2: 4324.6787995634495, c: 19.458390148113008, aT: -3.00632,
      n: [0.012436, 0.97315, 1.2795, 0.96956, 0.24873],
      th: [832.9999962983201, 2289.00000119312, 5009.00000151968, 5982.000000404159, 17800.000014507998] }];
  // pure vapour-pressure ancillaries (starting values only): p = pc exp(Tr/T sum n (1 - T/Tr)^t)
  const ANC = [
    { Tt: 195.495, Tr: 405.4, pc: 11333000.0, n: [0.0016131703802769548, 0.01339903644956835, -6.4517151211338994,
      -4.349105787320862, 1.8295554336003688, 7.011533373126273], t: [0.108, 0.435, 0.971, 4.023, 4.944, 17.494] },
    { Tt: 273.16, Tr: 647.096, pc: 22064000.0, n: [-9.75639641045262, 3.3357600887120102, -1.10029278432831,
      0.02037617155190105, -2.6668589845604367, 6.676721087238668], t: [1.018, 1.206, 2.327, 5.753, 4.215, 14.951] }];

  // ------------------------------------------------------------------ residual Helmholtz energy
  // Accumulators: value, delta*d/ddelta, tau*d/dtau, delta^2*d2/ddelta2.
  const acc = new Float64Array(4);
  const dpow = new Float64Array(7);
  function powerTerms(P, lnt, lnd, out) {
    let A = 0, Ad = 0, At = 0, Add = 0;
    const n = P.n, t = P.t, d = P.d, l = P.l;
    for (let i = 0; i < P.N; i++) {
      const L = l[i], dl = L ? dpow[L] : 0;
      const v = n[i] * Math.exp(t[i] * lnt + d[i] * lnd - dl);
      const k = d[i] - L * dl;
      A += v; Ad += v * k; At += v * t[i]; Add += v * (k * (k - 1) - L * L * dl);
    }
    out[0] = A; out[1] = Ad; out[2] = At; out[3] = Add;
  }
  function waterExtra(tau, delta, out) {         // Gaussian and non-analytic terms, added to out
    let A = 0, Ad = 0, At = 0, Add = 0;
    for (let i = 0; i < 3; i++) {
      const de = delta - GA.eps[i], dt = tau - GA.gam[i];
      const v = GA.n[i] * Math.pow(delta, GA.d[i]) * Math.pow(tau, GA.t[i]) *
        Math.exp(-GA.eta[i] * de * de - GA.beta[i] * dt * dt);
      const k = GA.d[i] - 2 * GA.eta[i] * delta * de;
      A += v; Ad += v * k; At += v * (GA.t[i] - 2 * GA.beta[i] * tau * dt);
      Add += v * (k * k - GA.d[i] - 2 * GA.eta[i] * delta * delta);
    }
    const dm1 = delta - 1, S = dm1 * dm1, tm1 = tau - 1;
    for (let i = 0; i < 2; i++) {
      const a = NA.a[i], b = NA.b[i], B = NA.B[i], Aa = NA.A[i], C = NA.C[i], D = NA.D[i], be = NA.beta[i];
      const k = 1 / (2 * be);
      const Sk1 = Math.pow(S, k - 1), Sa1 = Math.pow(S, a - 1);
      const theta = (1 - tau) + Aa * Math.pow(S, k);
      const Dl = theta * theta + B * Math.pow(S, a);
      const Dl_d = dm1 * (Aa * theta * (2 / be) * Sk1 + 2 * B * a * Sa1);
      const Dl_dd = Aa * theta * (2 / be) * Sk1 + 2 * B * a * Sa1 + 4 * B * a * (a - 1) * Sa1 +
        2 * Aa * Aa / (be * be) * Math.pow(S, 2 * k - 1) + Aa * theta * (4 / be) * (k - 1) * Sk1;
      const Db = Math.pow(Dl, b), Db1 = Math.pow(Dl, b - 1), Db2 = Math.pow(Dl, b - 2);
      const Db_d = b * Db1 * Dl_d, Db_dd = b * (Db1 * Dl_dd + (b - 1) * Db2 * Dl_d * Dl_d);
      const Db_t = -2 * theta * b * Db1;
      const Psi = Math.exp(-C * S - D * tm1 * tm1);
      const Psi_d = -2 * C * dm1 * Psi, Psi_dd = (2 * C * S - 1) * 2 * C * Psi, Psi_t = -2 * D * tm1 * Psi;
      const n = NA.n[i];
      const v = n * Db * delta * Psi;
      const vd = n * (Db * (Psi + delta * Psi_d) + Db_d * delta * Psi);
      const vdd = n * (Db * (2 * Psi_d + delta * Psi_dd) + 2 * Db_d * (Psi + delta * Psi_d) + Db_dd * delta * Psi);
      const vt = n * delta * (Db_t * Psi + Db * Psi_t);
      if (!(Number.isFinite(v) && Number.isFinite(vd) && Number.isFinite(vdd) && Number.isFinite(vt))) continue;
      A += v; Ad += delta * vd; At += tau * vt; Add += delta * delta * vdd;
    }
    out[0] += A; out[1] += Ad; out[2] += At; out[3] += Add;
  }

  // the last evaluation: alphar, delta*ar_delta, tau*ar_tau, delta^2*ar_deltadelta, d(alphar)/dx at (T, rho)
  const RES = { T: NaN, rho: NaN, x: NaN, ar: 0, d: 0, t: 0, dd: 0, ax: 0 };
  const aN = new Float64Array(4), aW = new Float64Array(4), s0 = new Float64Array(4), s1 = new Float64Array(4),
    s2 = new Float64Array(4);
  const DEP0 = power([...DEP.n.slice(0, 6)], [...DEP.t.slice(0, 6)], [...DEP.d.slice(0, 6)], [...DEP.l.slice(0, 6)]);
  const DEP1 = power([...DEP.n.slice(6, 13)], [...DEP.t.slice(6, 13)], [...DEP.d.slice(6, 13)], [...DEP.l.slice(6, 13)]);
  const DEP2 = power([...DEP.n.slice(13)], [...DEP.t.slice(13)], [...DEP.d.slice(13)], [...DEP.l.slice(13)]);

  function evalRes(T, rho, x) {
    x = clampX(x);
    if (T === RES.T && rho === RES.rho && x === RES.x) return RES;
    const x1 = 1 - x, xa = Math.pow(x, ALPHA), xb = Math.pow(x, BETA);
    const Tred = TcN * x * x + TcW * x1 * x1 + 2.0 * x * (1.0 - xa) * kT / 2.0 * (TcN + TcW);
    const dTred = 2 * TcN * x - 2 * TcW * x1 + (1 - (1 + ALPHA) * xa) * kT * (TcN + TcW);
    const vred = vcN * x * x + vcW * x1 * x1 + 2.0 * x * (1.0 - xb) * kV / 2.0 * (vcN + vcW);
    const dvred = 2 * vcN * x - 2 * vcW * x1 + (1 - (1 + BETA) * xb) * kV * (vcN + vcW);
    const tau = Tred / T, delta = rho / (1 / vred);
    const lnt = Math.log(tau), lnd = Math.log(delta);
    dpow[1] = delta; dpow[2] = delta * delta; dpow[3] = dpow[2] * delta; dpow[4] = dpow[2] * dpow[2];
    dpow[5] = dpow[4] * delta; dpow[6] = dpow[3] * dpow[3];
    powerTerms(NH3, lnt, lnd, aN);
    powerTerms(H2O, lnt, lnd, aW);
    waterExtra(tau, delta, aW);
    powerTerms(DEP0, lnt, lnd, s0);
    powerTerms(DEP1, lnt, lnd, s1);
    powerTerms(DEP2, lnt, lnd, s2);
    const xg = Math.pow(x, DEP_GAMMA), F = x * (1 - xg), dF = 1 - (1 + DEP_GAMMA) * xg;
    let ar = 0, d = 0, t = 0, dd = 0;
    for (let j = 0; j < 4; j++) {
      const sum = s0[j] + x * s1[j] + x * x * s2[j];
      acc[j] = x * aN[j] + x1 * aW[j] + F * sum;
    }
    ar = acc[0]; d = acc[1]; t = acc[2]; dd = acc[3];
    const depx = dF * (s0[0] + x * s1[0] + x * x * s2[0]) + F * (s1[0] + 2 * x * s2[0]);
    RES.T = T; RES.rho = rho; RES.x = x;
    RES.ar = ar; RES.d = d; RES.t = t; RES.dd = dd;
    // d(alphar)/dx at T, rho: through tau (Tred) and delta (vred), and explicitly
    RES.ax = t * dTred / Tred + d * dvred / vred + (aN[0] - aW[0] + depx);
    return RES;
  }

  // ideal gas: alpha_ig and -T d(alpha_ig)/dT at constant rho
  const IGR = { a0: 0, a10: 0 };
  function evalIdeal(T, rho, x) {
    x = clampX(x);
    let a0 = 0, a10 = 0;
    for (let i = 0; i < 2; i++) {
      const g = IG[i], zi = i === 0 ? x : 1 - x;
      let ai = Math.log(rho) + g.a1 + g.a2 / T + g.c + g.aT * Math.log(T), bi = g.a2 / T - g.aT;
      for (let j = 0; j < g.n.length; j++) {
        const u = g.th[j] / T;
        ai += g.n[j] * Math.log(1.0 - Math.exp(-u));
        bi += g.n[j] * u / (Math.exp(u) - 1);
      }
      a0 += zi * (ai + Math.log(zi));
      a10 += zi * bi;
    }
    IGR.a0 = a0; IGR.a10 = a10;
    return IGR;
  }

  // ------------------------------------------------------------------ Brent (scipy's brentq)
  function brentq(f, xa, xb, xtol = 2e-12, rtol = 4 * 2.220446049250313e-16, maxiter = 100) {
    let xpre = xa, xcur = xb, xblk = 0, fblk = 0, spre = 0, scur = 0;
    let fpre = f(xpre), fcur = f(xcur);
    if (fpre === 0) return xpre;
    if (fcur === 0) return xcur;
    if ((fpre > 0) === (fcur > 0)) throw new Error("brentq: f(a) and f(b) must have different signs");
    for (let i = 0; i < maxiter; i++) {
      if (fpre !== 0 && fcur !== 0 && (fpre > 0) !== (fcur > 0)) {
        xblk = xpre; fblk = fpre; spre = scur = xcur - xpre;
      }
      if (Math.abs(fblk) < Math.abs(fcur)) {
        xpre = xcur; xcur = xblk; xblk = xpre;
        fpre = fcur; fcur = fblk; fblk = fpre;
      }
      const delta = (xtol + rtol * Math.abs(xcur)) / 2, sbis = (xblk - xcur) / 2;
      if (fcur === 0 || Math.abs(sbis) < delta) return xcur;
      if (Math.abs(spre) > delta && Math.abs(fcur) < Math.abs(fpre)) {
        let stry;
        if (xpre === xblk) stry = -fcur * (xcur - xpre) / (fcur - fpre);
        else {
          const dpre = (fpre - fcur) / (xpre - xcur), dblk = (fblk - fcur) / (xblk - xcur);
          stry = -fcur * (fblk * dblk - fpre * dpre) / (dblk * dpre * (fblk - fpre));
        }
        if (2 * Math.abs(stry) < Math.min(Math.abs(spre), 3 * Math.abs(sbis) - delta)) { spre = scur; scur = stry; }
        else { spre = sbis; scur = sbis; }
      } else { spre = sbis; scur = sbis; }
      xpre = xcur; fpre = fcur;
      xcur += Math.abs(scur) > delta ? scur : (sbis > 0 ? delta : -delta);
      fcur = f(xcur);
    }
    throw new Error("brentq: failed to converge");
  }

  // ------------------------------------------------------------------ the model
  // pure vapour pressures for starting values: ln p linear in 1/T between 120 ancillary points (as Python tabulates
  // CoolProp's saturation pressures), extended along the end segments' slopes
  const PSAT = ANC.map(a => {
    const anc = T => { const th = 1 - T / a.Tr; let s = 0; for (let i = 0; i < a.n.length; i++) s += a.n[i] * Math.pow(th, a.t[i]); return a.pc * Math.exp(a.Tr / T * s); };
    const Ts = []; for (let i = 0; i < 120; i++) Ts.push(a.Tt + 0.5 + (a.Tr - 1.0 - a.Tt - 0.5) * i / 119);
    Ts.reverse();
    return { invT: Float64Array.from(Ts, t => 1 / t), lnp: Float64Array.from(Ts, t => Math.log(anc(t))) };
  });
  function psatPure(i, T) {
    const { invT, lnp } = PSAT[i], u = 1 / T, N = invT.length;
    let k;
    if (u <= invT[0]) k = 0;
    else if (u >= invT[N - 1]) k = N - 2;
    else { let lo = 0, hi = N; while (lo < hi) { const m = (lo + hi) >> 1; if (invT[m] < u) lo = m + 1; else hi = m; } k = Math.min(lo - 1, N - 2); }
    return Math.exp(lnp[k] + (lnp[k + 1] - lnp[k]) * (u - invT[k]) / (invT[k + 1] - invT[k]));
  }

  function pressure(T, rho, x) { return rho * R * T * (1.0 + evalRes(T, rho, x).d); }
  function dpdrho(T, rho, x) { const r = evalRes(T, rho, x); return R * T * (1.0 + 2.0 * r.d + r.dd); }

  /* Molar density of the liquid or vapour root at T (K), p (Pa), or null; rho0 overrides the starting value (used to
     warm-start the equilibrium Newton). */
  function density(T, p, x, phase, rho0) {
    const liq = phase === "liquid";
    let rho = rho0 !== undefined ? rho0 : (liq ? 50000.0 - 15000.0 * x : p / (R * T)), ok = false;
    for (let it = 0; it < 100; it++) {
      const r = evalRes(T, rho, x);
      const pr = rho * R * T * (1.0 + r.d), d = R * T * (1.0 + 2.0 * r.d + r.dd);
      if (!(d > 0 && Number.isFinite(pr))) { rho *= liq ? 1.05 : 0.9; continue; }
      let nw = rho - (pr - p) / d;
      if (nw <= 0) nw = rho / 2;
      if (liq && nw < 0.25 * rho) nw = 0.5 * rho;
      if (Math.abs(nw - rho) <= 1e-11 * rho) { rho = nw; ok = true; break; }
      rho = nw;
    }
    if (!ok) return null;
    if (!(dpdrho(T, rho, x) > 0) || !(Math.abs(pressure(T, rho, x) - p) <= 1e-6 * p + 1e-3)) return null;
    const Zc = p / (rho * R * T);
    if (liq && Zc > 0.3) return null;
    if (!liq && Zc < 0.3) return null;
    return rho;
  }

  function h_s(T, rho, x) {
    const r = evalRes(T, rho, x), ar = r.ar, ad = r.d, at = r.t, ig = evalIdeal(T, rho, x);
    const a10 = at + ig.a10;
    return [R * T * (a10 + 1.0 + ad), R * (a10 - ar - ig.a0)];
  }

  const LP = new Float64Array(2);
  function lnPhiInto(T, rho, x, out) {
    const r = evalRes(T, rho, x), xc = RES.x, base = r.ar + r.d - Math.log(1 + r.d);
    out[0] = base + (1 - xc) * r.ax;
    out[1] = base - xc * r.ax;
    return out;
  }
  function ln_phi(T, rho, x) { const o = lnPhiInto(T, rho, x, LP); return [o[0], o[1]]; }

  /* Newton on F = ln(x_i phi_i^L) - ln(y_i phi_i^V) = 0 in (x, y), each phase's density solved at (T, p); the
     composition derivatives of ln phi by finite differences (they only steer). */
  const LL = new Float64Array(2), LV = new Float64Array(2), LH = new Float64Array(2);
  function vleNewton(T, p, x, y, rL, rV) {
    for (let it = 0; it < 60; it++) {
      rL = density(T, p, x, "liquid", rL) || density(T, p, x, "liquid");
      rV = density(T, p, y, "vapour", rV) || density(T, p, y, "vapour");
      if (rL === null || rV === null) return null;
      lnPhiInto(T, rL, x, LL); lnPhiInto(T, rV, y, LV);
      const F1 = LL[0] + Math.log(x) - LV[0] - Math.log(y);
      const F2 = LL[1] + Math.log(1 - x) - LV[1] - Math.log(1 - y);
      if (!(Number.isFinite(F1) && Number.isFinite(F2))) return null;
      if (Math.abs(F1) < 1e-14 && Math.abs(F2) < 1e-14) return [x, y, rL, rV];   // (the noise floor is ~1e-13)
      const hx = 1e-7 * Math.min(x, 1 - x), hy = 1e-7 * Math.min(y, 1 - y);
      const rLh = density(T, p, x + hx, "liquid", rL), rVh = density(T, p, y + hy, "vapour", rV);
      if (rLh === null || rVh === null) return null;
      lnPhiInto(T, rLh, x + hx, LH);
      const J11 = (LH[0] - LL[0]) / hx + 1 / x, J21 = (LH[1] - LL[1]) / hx - 1 / (1 - x);
      lnPhiInto(T, rVh, y + hy, LH);
      const J12 = -(LH[0] - LV[0]) / hy - 1 / y, J22 = -(LH[1] - LV[1]) / hy + 1 / (1 - y);
      const det = J11 * J22 - J12 * J21;
      if (!(det !== 0 && Number.isFinite(det))) return null;
      let dx = -(J22 * F1 - J12 * F2) / det, dy = -(-J21 * F1 + J11 * F2) / det;
      let lam = 1;
      while (lam > 1e-3 && !(x + lam * dx > 0 && x + lam * dx < 1 && y + lam * dy > 0 && y + lam * dy < 1)) lam /= 2;
      dx *= lam; dy *= lam;
      x += dx; y += dy;
      if (Math.abs(dx) < 1e-11 && Math.abs(dy) < 1e-11) {         // quadratic: this last step lands at the noise
        rL = density(T, p, x, "liquid", rL); rV = density(T, p, y, "vapour", rV);
        return rL === null || rV === null ? null : [x, y, rL, rV];
      }
    }
    return null;
  }

  /* Coexisting liquid and vapour ammonia mole fractions at T (K) and p (Pa): [x, y, rhoL, rhoV], or null. */
  function equilibrium_TP(T, p, x0, y0) {
    if (x0 == null || y0 == null) {
      const Ka = psatPure(0, T) / p, Kw = psatPure(1, T) / p;
      if (!(Ka > 1 && 1 > Kw)) return null;
      x0 = (1 - Kw) / (Ka - Kw);
      y0 = Ka * x0;
    }
    let x = x0, y = y0;
    let rL = density(T, p, x, "liquid"), rV = density(T, p, y, "vapour");
    if (rL === null || rV === null) return null;
    if (x > 0 && x < 1 && y > 0 && y < 1) {
      const r = vleNewton(T, p, x, y, rL, rV);
      if (r) {
        const [xs, ys, dL, dV] = r;
        if (Math.abs(ys - xs) > 1e-6 && ys > xs && Math.abs(pressure(T, dL, xs) / p - 1) < 1e-7 &&
            Math.abs(pressure(T, dV, ys) / p - 1) < 1e-7 && dL > dV) return r;
      }
    }
    let conv = false;
    for (let it = 0; it < 400; it++) {
      rL = density(T, p, x, "liquid");
      rV = density(T, p, y, "vapour");
      if (rL === null || rV === null) return null;
      lnPhiInto(T, rL, x, LL); lnPhiInto(T, rV, y, LV);
      const Ka = Math.exp(LL[0] - LV[0]), Kw = Math.exp(LL[1] - LV[1]);
      if (!(Ka > 1 && 1 > Kw)) return null;
      const xn = (1 - Kw) / (Ka - Kw), yn = Ka * xn;
      if (!(xn > 0 && xn < 1 && yn > 0 && yn < 1)) return null;
      if (Math.abs(xn - x) < 1e-10 && Math.abs(yn - y) < 1e-10) { x = xn; y = yn; conv = true; break; }
      x += 0.5 * (xn - x); y += 0.5 * (yn - y);
    }
    if (!conv || Math.abs(y - x) < 1e-6) return null;
    return [x, y, density(T, p, x, "liquid"), density(T, p, y, "vapour")];
  }

  function raoultT(p, x, dew) {
    return brentq(T => {
      const pa = psatPure(0, T), pw = psatPure(1, T);
      return dew ? (x / pa + (1 - x) / pw) * p - 1 : x * pa + (1 - x) * pw - p;
    }, 180.0, 900.0);
  }

  const P_ = (p) => (p / 1e6).toFixed(3), X_ = (x) => x.toFixed(4);
  /* Bubble (dew=false: x is the liquid) or dew (dew=true: x is the vapour) temperature (K) at p (Pa), and the
     incipient phase's mole fraction: successive substitution on the K-values, secant on temperature. */
  const T_RETRY = [0.99, 1.01, 0.98, 1.02, 0.97, 1.03, 0.95, 1.05];   // saturation starts tried when the estimate has no root
  function saturationT(p, x, dew, T0) {
    let T = T0 != null ? T0 : raoultT(p, x, dew);
    const base = dew ? "vapour" : "liquid", other = dew ? "liquid" : "vapour";
    const LB = new Float64Array(2), LI = new Float64Array(2);
    function resid(T, inc) {
      const rB = density(T, p, x, base);
      if (rB === null) return [null, inc];
      lnPhiInto(T, rB, x, LB);
      if (inc === null) {
        const K0 = psatPure(0, T) / p, K1 = psatPure(1, T) / p;
        inc = dew ? (x / K0) / (x / K0 + (1 - x) / K1) : K0 * x / (K0 * x + K1 * (1 - x));
      }
      let S = NaN;
      for (let it = 0; it < 60; it++) {
        const rI = density(T, p, inc, other);
        if (rI === null) return [null, inc];
        lnPhiInto(T, rI, inc, LI);
        const t0 = x * Math.exp(LB[0] - LI[0]), t1 = (1 - x) * Math.exp(LB[1] - LI[1]);
        S = t0 + t1;
        const nw = t0 / S;
        if (Math.abs(nw - inc) < 1e-11) { inc = nw; break; }
        inc = nw;
      }
      return [Math.log(S), inc];
    }
    let [g0, inc] = resid(T, null);
    // a phase root is missing at the estimate (Raoult's law is far off for this mixture): try nearby temperatures
    for (const f of T_RETRY) { if (g0 !== null) break; [g0, inc] = resid(T * f, null); if (g0 !== null) T *= f; }
    if (g0 === null) throw new Error(`ammonia-water ${dew ? "dew" : "bubble"} point at ${P_(p)} MPa, x=${X_(x)}: ` +
      "no phase root at the starting temperature");
    let T1 = T * ((g0 > 0) === dew ? 1.002 : 0.998);
    let [g1, inc1] = resid(T1, inc);
    for (let it = 0; it < 60; it++) {
      if (g1 === null) { T1 = (T + T1) / 2; [g1, inc1] = resid(T1, inc); continue; }
      if (Math.abs(g1) < 1e-11) return [T1, inc1];
      if (g1 === g0) break;
      let T2 = T1 - g1 * (T1 - T) / (g1 - g0);
      T2 = Math.min(Math.max(T2, T1 - 30), T1 + 30);
      T = T1; g0 = g1; inc = inc1;
      T1 = T2;
      [g1, inc1] = resid(T1, inc);
    }
    if (g1 !== null && Math.abs(g1) < 1e-8) return [T1, inc1];
    throw new Error(`ammonia-water ${dew ? "dew" : "bubble"} point at ${P_(p)} MPa, x=${X_(x)} did not converge`);
  }

  /* Bubble or dew pressure (Pa) at T (K), and the incipient phase's mole fraction. */
  function saturationP(T, x, dew) {
    const pa = psatPure(0, T), pw = psatPure(1, T);
    let p, inc;
    if (dew) { p = 1.0 / (x / pa + (1 - x) / pw); inc = (x / pa) / (x / pa + (1 - x) / pw); }
    else { p = x * pa + (1 - x) * pw; inc = x * pa / p; }
    const base = dew ? "vapour" : "liquid", other = dew ? "liquid" : "vapour";
    const LB = new Float64Array(2), LI = new Float64Array(2);
    for (let it = 0; it < 500; it++) {
      const rB = density(T, p, x, base), rI = density(T, p, inc, other);
      if (rB === null || rI === null) { p *= ((rI === null) !== dew) ? 0.7 : 1.3; continue; }
      lnPhiInto(T, rB, x, LB); lnPhiInto(T, rI, inc, LI);
      const t0 = x * Math.exp(LB[0] - LI[0]), t1 = (1 - x) * Math.exp(LB[1] - LI[1]), S = t0 + t1;
      const ni = t0 / S, pn = dew ? p / S : p * S;
      if (Math.abs(pn / p - 1) < 1e-10 && Math.abs(ni - inc) < 1e-9) return [pn, ni];
      p = pn; inc = ni;
    }
    throw new Error(`ammonia-water ${dew ? "dew" : "bubble"} pressure at ${T.toFixed(2)} K, x=${X_(x)} did not converge`);
  }

  const model = {
    pressure, density: (T, p, x, phase) => density(T, p, x, phase), h_s, ln_phi, equilibrium_TP,
    bubble_T: (p, x) => saturationT(p, x, false), dew_T: (p, y) => saturationT(p, y, true),
    bubble_p: (T, x) => saturationP(T, x, false), dew_p: (T, y) => saturationP(T, y, true),
    psat_pure: (name, T) => psatPure(name === "Ammonia" ? 0 : 1, T),
  };

  // ------------------------------------------------------------------ streams of fixed composition
  class AmmoniaWaterStream {
    constructor(w, m) {
      if (!(w > 0 && w < 1)) throw new Error("the ammonia mass fraction must lie between 0 and 1");
      this.model = m || model;
      this.w = w;
      this.x = mole_fraction(w);
      this.name = this.fluid = `ammonia-water (${(w * 100).toFixed(1)} % NH3 by mass)`;
      this.backend = "IAPWS G4-01 (JavaScript)";
      this.P_crit = NaN; this.T_crit = NaN; this.fast = false;
      this._memo = new Map(); this._sat = new Map(); this._iso = new Map();
    }
    fast_view() { return this; }

    // T(h) along an isobar: liquid from 120 K below the bubble point, the two-phase region in 1 K steps (warm-started),
    // vapour to 160 K above the dew point; linear between the points
    _iso_table(P) {
      const key = Math.round(P * 1e9);
      const got = this._iso.get(key);
      if (got) return got;
      const [Tb, Td] = this.bubble_dew(P);
      const pts = [];
      for (let T = Tb - 120.0; T < Tb; T += 4.0) {
        try { const st = this._single(T, P, "liquid"); pts.push([st.h, st.T]); } catch (e) { /* skip */ }
      }
      pts.push([this._single(Tb, P, "liquid").h, Tb - T_K]);
      const n = Math.max(Math.trunc((Td - Tb) / 1.0), 8);
      let guess = null;
      for (let i = 1; i < n; i++) {
        const Tt = Tb + (Td - Tb) * i / n;
        try { const st = this._two_phase(Tt, P, guess); guess = [st.x_mol, st.y_mol]; pts.push([st.h, st.T]); }
        catch (e) { guess = null; }
      }
      pts.push([this._single(Td, P, "vapour").h, Td - T_K]);
      for (let T = Td + 4.0; T <= Td + 160.0; T += 4.0) {
        try { const st = this._single(T, P, "vapour"); pts.push([st.h, st.T]); } catch (e) { /* skip */ }
      }
      pts.sort((a, b) => a[0] - b[0] || a[1] - b[1]);
      const H = [], Ts = [];
      for (const [h, t] of pts) if (!H.length || h > H[H.length - 1] + 1e-9) { H.push(h); Ts.push(t); }
      const tab = { H: Float64Array.from(H), T: Float64Array.from(Ts) };
      if (this._iso.size > 200) this._iso.clear();
      this._iso.set(key, tab);
      return tab;
    }

    /** Bubble and dew temperatures (K) at P (MPa), with the incipient vapour and liquid mole fractions. */
    bubble_dew(P) {
      const key = Math.round(P * 1e10);
      let r = this._sat.get(key);
      if (!r) {
        const m = this.model, p = P * 1e6;
        const [Tb, yb] = m.bubble_T(p, this.x), [Td, xd] = m.dew_T(p, this.x);
        r = [Tb, Td, yb, xd];
        this._sat.set(key, r);
      }
      return r;
    }

    _single(T, P, phase) {
      const m = this.model, p = P * 1e6;
      let rho = m.density(T, p, this.x, phase);
      if (rho === null) {
        const other = phase === "liquid" ? "vapour" : "liquid";
        rho = m.density(T, p, this.x, other);
        if (rho === null) throw new Error(`${this.name}: no single-phase state at ${(T - T_K).toFixed(2)} C, ${P.toFixed(4)} MPa`);
        phase = other;
      }
      const [h, s] = m.h_s(T, rho, this.x), M = molarMass(this.x), liq = phase === "liquid";
      return { T: T - T_K, P, h: h / M / 1000, s: s / M / 1000, D: rho * M, q: -1.0,
        phase: liq ? "liquid" : "superheated_vapour", w_liq: liq ? this.w : null, w_vap: liq ? null : this.w };
    }

    _two_phase(T, P, guess) {
      const m = this.model, p = P * 1e6;
      if (!guess) {
        const [Tb, Td, yb, xd] = this.bubble_dew(P);
        const f = Td > Tb ? Math.min(Math.max((T - Tb) / (Td - Tb), 0.0), 1.0) : 0.5;
        guess = [this.x + f * (xd - this.x), yb + f * (this.x - yb)];
      }
      let eq = m.equilibrium_TP(T, p, guess[0], guess[1]);
      if (eq === null) eq = m.equilibrium_TP(T, p);
      if (eq === null) throw new Error(`${this.name}: no two-phase equilibrium at ${(T - T_K).toFixed(2)} C, ${P.toFixed(4)} MPa`);
      const [x, y, rL, rV] = eq;
      const V = Math.min(Math.max((this.x - x) / (y - x), 0.0), 1.0);
      const [hL, sL] = m.h_s(T, rL, x), [hV, sV] = m.h_s(T, rV, y);
      const ML = molarMass(x), MV = molarMass(y), nL = 1.0 - V, nV = V, mass = nL * ML + nV * MV;
      return { T: T - T_K, P, h: (nL * hL + nV * hV) / mass / 1000, s: (nL * sL + nV * sV) / mass / 1000,
        D: mass / (nL / rL + nV / rV), q: nV * MV / mass, phase: "two_phase",
        w_liq: mass_fraction(x), w_vap: mass_fraction(y), x_mol: x, y_mol: y };
    }

    _memoised(key, fn) {
      let st = this._memo.get(key);
      if (st === undefined) {
        st = fn();
        if (this._memo.size > 100000) this._memo.clear();
        this._memo.set(key, st);
      }
      return st;
    }

    PT(P, T) {
      return this._memoised("PT" + P + "," + T, () => {
        const Tk = T + T_K, [Tb, Td] = this.bubble_dew(P);
        if (Tk <= Tb) return this._single(Tk, P, "liquid");
        if (Tk >= Td) return this._single(Tk, P, "vapour");
        return this._two_phase(Tk, P);
      });
    }

    _solve_T(P, prop, value) {
      const [Tb, Td] = this.bubble_dew(P);
      const vb = this._single(Tb, P, "liquid")[prop], vd = this._single(Td, P, "vapour")[prop];
      let lo, hi, region;
      if (value <= vb) {
        lo = Tb - 1.0; hi = Tb; region = "liquid";
        while (this._single(lo, P, "liquid")[prop] > value) {
          hi = lo; lo -= 25.0;
          if (lo < 180.0) throw new Error(`${this.name}: ${prop}=${value} below the liquid range at ${P} MPa`);
        }
      } else if (value >= vd) {
        lo = Td; hi = Td + 1.0; region = "vapour";
        while (this._single(hi, P, "vapour")[prop] < value) {
          lo = hi; hi += 25.0;
          if (hi > 1200.0) throw new Error(`${this.name}: ${prop}=${value} above the vapour range at ${P} MPa`);
        }
      } else { lo = Tb; hi = Td; region = "two"; }
      const f = T => {
        if (region === "two") {
          if (T <= Tb + 1e-9) return vb - value;
          if (T >= Td - 1e-9) return vd - value;
          return this._two_phase(T, P)[prop] - value;
        }
        return this._single(T, P, region)[prop] - value;
      };
      const T = brentq(f, lo, hi, 1e-10, 1e-13);
      if (region === "two")
        return Tb < T && T < Td ? this._two_phase(T, P) : this._single(T, P, T <= Tb ? "liquid" : "vapour");
      return this._single(T, P, region);
    }

    PH(P, h) { return this._memoised("PH" + P + "," + h, () => this._solve_T(P, "h", h)); }
    PS(P, s) { return this._memoised("PS" + P + "," + s, () => this._solve_T(P, "s", s)); }

    /** Saturated liquid (q=0), saturated vapour (q=1) or a vapour mass fraction q in between, at P (MPa). */
    PQ(P, q) {
      return this._memoised("PQ" + P + "," + q, () => {
        const [Tb, Td] = this.bubble_dew(P);
        if (q <= 0) return Object.assign(this._single(Tb, P, "liquid"), { q: 0.0, phase: "saturated_liquid" });
        if (q >= 1) return Object.assign(this._single(Td, P, "vapour"), { q: 1.0, phase: "saturated_vapour" });
        const T = brentq(T => this._two_phase(T, P).q - q, Tb + 1e-7, Td - 1e-7, 1e-10);
        return this._two_phase(T, P);
      });
    }

    /** The saturation state at T (C) and vapour mass fraction q (0: bubble, 1: dew). */
    TQ(T, q) {
      return this._memoised("TQ" + T + "," + q, () => {
        const m = this.model, Tk = T + T_K;
        let p = q <= 0 ? m.bubble_p(Tk, this.x)[0] : q >= 1 ? m.dew_p(Tk, this.x)[0] : null;
        if (p === null) {
          const lo = m.dew_p(Tk, this.x)[0], hi = m.bubble_p(Tk, this.x)[0];
          p = Math.exp(brentq(lp => this.PT(Math.exp(lp) / 1e6, T).q - q, Math.log(lo) + 1e-9, Math.log(hi) - 1e-9, 1e-12));
        }
        return this.PQ(p / 1e6, q);
      });
    }

    is_gas(st) { return st.phase === "saturated_vapour" || st.phase === "superheated_vapour"; }
    Tph(P, h) { return this.PH(P, h).T; }

    /** Temperature (C) at (P, h) from the isobar table (a PH flash outside it). */
    T_ph_fast(P, h) {
      const { H, T } = this._iso_table(P), N = H.length;
      if (h < H[0] || h > H[N - 1]) return this.PH(P, h).T;
      let lo = 0, hi = N;
      while (lo < hi) { const m = (lo + hi) >> 1; if (h < H[m]) hi = m; else lo = m + 1; }
      const i = Math.min(lo - 1, N - 2);
      return T[i] + (T[i + 1] - T[i]) * (h - H[i]) / (H[i + 1] - H[i]);
    }
  }

  const api = { M_NH3, M_H2O, R_MIX: R, mole_fraction, mass_fraction, model, AmmoniaWaterStream,
    stream: (w, m) => new AmmoniaWaterStream(w, m), brentq, backend: "IAPWS G4-01 (JavaScript)" };
  if (typeof module !== "undefined" && module.exports) module.exports = api; else root.AMMWAT = api;
})(typeof self !== "undefined" ? self : this);
