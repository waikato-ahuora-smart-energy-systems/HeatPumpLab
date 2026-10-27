"""Ammonia-water mixture properties for the hybrid absorption-compression heat pump.

The mixture model is the IAPWS Guideline G4-01 (2001) formulation of Tillner-Roth and Friend (J. Phys. Chem. Ref. Data
27 (1998) 63), evaluated through NIST's teqp (``teqp.AmmoniaWaterTillnerRoth``: the residual Helmholtz energy and its
derivatives). CoolProp has no ammonia-water mixture, so this module adds what a cycle model needs on top of teqp:

* the ideal-gas parts, from CoolProp's pure-fluid data (IAPWS-95 water and CoolProp's ammonia) through teqp's
  ``IdealHelmholtz`` model. Enthalpy and entropy therefore use CoolProp's reference states for each pure component,
  which cancel in any balance that conserves each component;
* densities from temperature and pressure (Newton on p(rho) for the liquid and the vapour root);
* bubble and dew points, and the equilibrium liquid and vapour of a two-phase state, from the fugacity coefficients
  (for a binary the coexisting compositions at given T and p do not depend on the overall composition);
* PT, PH and PS flashes at a fixed overall composition, and saturation states at a given vapour fraction.

Mass basis throughout: ``w`` is the ammonia mass fraction. Units follow HeatPumpLab: T in C, P in MPa, h in kJ/kg,
s in kJ/kg/K, D in kg/m3. A state is a dict like the cycles' fluid states, plus the phase compositions
(``w_liq``, ``w_vap``) and the vapour mass fraction ``q`` (-1 outside the two-phase region).

teqp is an optional dependency: ``pip install "heatpumplab[ammonia-water]"``.
"""
from __future__ import annotations

import json
import math

import numpy as np
from scipy.optimize import brentq

M_NH3 = 17.03026e-3          # kg/mol (IAPWS G4-01, Table 1)
M_H2O = 18.015268e-3
R_MIX = 8.314471             # J/mol/K, the gas constant of the mixture model
_X_MIN = 1e-10               # the Tillner-Roth model rejects mole fractions of exactly zero
T_K = 273.15


def mole_fraction(w: float) -> float:
    """Ammonia mole fraction from its mass fraction."""
    a = w / M_NH3
    return a / (a + (1.0 - w) / M_H2O)


def mass_fraction(x: float) -> float:
    """Ammonia mass fraction from its mole fraction."""
    return x * M_NH3 / (x * M_NH3 + (1.0 - x) * M_H2O)


def _molar_mass(x: float) -> float:
    return x * M_NH3 + (1.0 - x) * M_H2O


class AmmoniaWater:
    """The ammonia-water property model (see the module docstring). One instance serves any composition."""

    def __init__(self):
        try:
            import teqp
            import CoolProp.CoolProp as CP
        except ImportError as exc:   # pragma: no cover - depends on the optional extra
            raise ImportError('ammonia-water properties need teqp: pip install "heatpumplab[ammonia-water]"') from exc
        self.res = teqp.AmmoniaWaterTillnerRoth()

        def ideal(name):
            d = json.loads(CP.get_fluid_param_string(name, "JSON"))
            return teqp.convert_CoolProp_idealgas(json.dumps(d[0] if isinstance(d, list) else d), 0)

        self.ig = teqp.make_model({"kind": "IdealHelmholtz", "model": [ideal("Ammonia"), ideal("Water")]})
        self._CP = CP
        self._psat_tabs: dict = {}

    # ---- single phase at (T, rho) -----------------------------------------------------

    @staticmethod
    def _z(x: float) -> np.ndarray:
        x = min(max(x, _X_MIN), 1.0 - _X_MIN)
        return np.array([x, 1.0 - x])

    def pressure(self, T: float, rho: float, x: float) -> float:
        """Pressure (Pa) at T (K), molar density rho (mol/m3) and ammonia mole fraction x."""
        return rho * R_MIX * T * (1.0 + self.res.get_Ar01(T, rho, self._z(x)))

    def _dpdrho(self, T, rho, z):
        return R_MIX * T * (1.0 + 2.0 * self.res.get_Ar01(T, rho, z) + self.res.get_Ar02(T, rho, z))

    def density(self, T: float, p: float, x: float, phase: str) -> float | None:
        """Molar density (mol/m3) of the liquid or vapour root at T (K) and p (Pa), or None when that root does not
        exist (no mechanically stable solution reached from the phase's starting guess)."""
        z = self._z(x)
        rho = 50000.0 - 15000.0 * x if phase == "liquid" else p / (R_MIX * T)
        for _ in range(100):
            try:
                pr = rho * R_MIX * T * (1.0 + self.res.get_Ar01(T, rho, z))
                d = self._dpdrho(T, rho, z)
            except RuntimeError:
                return None
            if not (d > 0 and math.isfinite(pr)):
                # outside the stable branch: step towards the phase's side
                rho = rho * (1.05 if phase == "liquid" else 0.9)
                continue
            step = (pr - p) / d
            new = rho - step
            if new <= 0:
                new = rho / 2
            if phase == "liquid" and new < 0.25 * rho:
                new = 0.5 * rho
            if abs(new - rho) <= 1e-11 * rho:
                rho = new
                break
            rho = new
        else:
            return None
        try:
            if self._dpdrho(T, rho, z) <= 0 or abs(self.pressure(T, rho, x) - p) > 1e-6 * p + 1e-3:
                return None
        except RuntimeError:
            return None
        # a liquid guess that slid to a gas-like density (or the reverse) is not that phase's root
        Zc = p / (rho * R_MIX * T)
        if phase == "liquid" and Zc > 0.3:
            return None
        if phase == "vapour" and Zc < 0.3:
            return None
        return rho

    def h_s(self, T: float, rho: float, x: float) -> tuple[float, float]:
        """Molar enthalpy (J/mol) and entropy (J/mol/K) at T (K), rho (mol/m3), x."""
        z = self._z(x)
        r, i = self.res, self.ig
        a10 = r.get_Ar10(T, rho, z) + i.get_Ar10(T, rho, z)
        h = R_MIX * T * (a10 + 1.0 + r.get_Ar01(T, rho, z))
        s = R_MIX * (a10 - r.get_Ar00(T, rho, z) - i.get_Ar00(T, rho, z))
        return h, s

    def ln_phi(self, T: float, rho: float, x: float) -> np.ndarray:
        z = self._z(x)
        return np.log(np.asarray(self.res.get_fugacity_coefficients(T, rho * z)))

    # ---- phase equilibrium ------------------------------------------------------------

    def _psat_pure(self, name: str, T: float) -> float:
        """A starting estimate of a pure component's vapour pressure (Pa): ln p linear in 1/T between CoolProp
        saturation points tabulated once (triple to critical point), extended along the end segments' slopes."""
        tab = self._psat_tabs.get(name)
        if tab is None:
            CP = self._CP
            Tc, Tt = CP.PropsSI("Tcrit", name), CP.PropsSI("Ttriple", name)
            Ts = [Tt + 0.5 + (Tc - 1.0 - Tt - 0.5) * i / 119 for i in range(120)]
            invT = [1.0 / t for t in reversed(Ts)]
            lnp = [math.log(CP.PropsSI("P", "T", t, "Q", 0, name)) for t in reversed(Ts)]
            tab = self._psat_tabs[name] = (invT, lnp)
        invT, lnp = tab
        u = 1.0 / T
        if u <= invT[0]:
            i = 0
        elif u >= invT[-1]:
            i = len(invT) - 2
        else:
            i = min(int(np.searchsorted(invT, u)) - 1, len(invT) - 2)
        return math.exp(lnp[i] + (lnp[i + 1] - lnp[i]) * (u - invT[i]) / (invT[i + 1] - invT[i]))

    def _raoult_T(self, p: float, x: float, dew: bool) -> float:
        def f(T):
            pa, pw = self._psat_pure("Ammonia", T), self._psat_pure("Water", T)
            return (x / pa + (1 - x) / pw) * p - 1 if dew else x * pa + (1 - x) * pw - p
        return brentq(f, 180.0, 900.0)

    def equilibrium_TP(self, T: float, p: float, x0: float | None = None, y0: float | None = None):
        """Coexisting liquid and vapour ammonia mole fractions at T (K) and p (Pa): (x, y, rhoL, rhoV), or None when
        no two phases exist there. teqp's Newton solver (``mix_VLE_Tp``) from the given guesses (or Raoult's law),
        checked for a non-trivial solution with both phases at p; damped successive substitution on the K-values
        (x = (1 - K_w) / (K_a - K_w), y = K_a x for a binary) if it fails."""
        if x0 is None or y0 is None:
            pa, pw = self._psat_pure("Ammonia", T), self._psat_pure("Water", T)
            Ka, Kw = pa / p, pw / p
            if not Ka > 1 > Kw:
                return None
            x0 = (1 - Kw) / (Ka - Kw)
            y0 = Ka * x0
        x, y = x0, y0
        rL = self.density(T, p, x, "liquid")
        rV = self.density(T, p, y, "vapour")
        if rL is None or rV is None:
            return None
        try:
            r = self.res.mix_VLE_Tp(T, p, rL * self._z(x), rV * self._z(y))
            vL, vV = np.asarray(r.rhovecL, dtype=float), np.asarray(r.rhovecV, dtype=float)
            if np.all(np.isfinite(vL)) and np.all(np.isfinite(vV)) and np.all(vL > 0) and np.all(vV > 0):
                xs, ys, dL, dV = vL[0] / vL.sum(), vV[0] / vV.sum(), vL.sum(), vV.sum()
                if (abs(ys - xs) > 1e-6 and ys > xs and abs(self.pressure(T, dL, xs) / p - 1) < 1e-7
                        and abs(self.pressure(T, dV, ys) / p - 1) < 1e-7 and dL > dV):
                    return float(xs), float(ys), float(dL), float(dV)
        except RuntimeError:
            pass
        for _ in range(400):
            rL = self.density(T, p, x, "liquid")
            rV = self.density(T, p, y, "vapour")
            if rL is None or rV is None:
                return None
            K = np.exp(self.ln_phi(T, rL, x) - self.ln_phi(T, rV, y))
            Ka, Kw = K
            if not Ka > 1 > Kw:
                return None
            xn = (1 - Kw) / (Ka - Kw)
            yn = Ka * xn
            if not (0 < xn < 1 and 0 < yn < 1):
                return None
            if abs(xn - x) < 1e-10 and abs(yn - y) < 1e-10:
                x, y = xn, yn
                break
            x, y = x + 0.5 * (xn - x), y + 0.5 * (yn - y)
        else:
            return None
        if abs(y - x) < 1e-6:
            return None       # trivial solution
        return x, y, self.density(T, p, x, "liquid"), self.density(T, p, y, "vapour")

    def _saturation(self, p: float, x: float, dew: bool, T0: float | None = None) -> tuple[float, float]:
        """Bubble (dew=False: x is the liquid) or dew (dew=True: x is the vapour) temperature (K) at p (Pa), and the
        incipient phase's ammonia mole fraction. Successive substitution on the K-values, secant on temperature."""
        T = T0 if T0 is not None else self._raoult_T(p, x, dew)
        inc = None

        def resid(T, inc):
            """ln of the sum the K-values must make 1, and the incipient composition."""
            base = "vapour" if dew else "liquid"
            rB = self.density(T, p, x, base)
            if rB is None:
                return None, inc
            lpB = self.ln_phi(T, rB, x)
            if inc is None:
                pa, pw = self._psat_pure("Ammonia", T), self._psat_pure("Water", T)
                K = np.array([pa / p, pw / p])
                inc = (x / K[0]) / (x / K[0] + (1 - x) / K[1]) if dew else K[0] * x / (K[0] * x + K[1] * (1 - x))
            for _ in range(60):
                rI = self.density(T, p, inc, "liquid" if dew else "vapour")
                if rI is None:
                    return None, inc
                lpI = self.ln_phi(T, rI, inc)
                zv = np.array([x, 1 - x])
                if dew:
                    t = zv * np.exp(lpB - lpI)            # x_i = y_i phi_V,i / phi_L,i
                else:
                    t = zv * np.exp(lpB - lpI)            # y_i = x_i phi_L,i / phi_V,i
                S = t.sum()
                new = t[0] / S
                if abs(new - inc) < 1e-11:
                    inc = new
                    break
                inc = new
            return math.log(S), inc

        g0, inc = resid(T, inc)
        if g0 is None:
            raise ValueError(f"ammonia-water {'dew' if dew else 'bubble'} point at {p / 1e6:.3f} MPa, x={x:.4f}: "
                             "no phase root at the starting temperature")
        T1 = T * (1.002 if (g0 > 0) == dew else 0.998)
        g1, inc1 = resid(T1, inc)
        for _ in range(60):
            if g1 is None:
                T1 = (T + T1) / 2
                g1, inc1 = resid(T1, inc)
                continue
            if abs(g1) < 1e-11:
                return T1, inc1
            if g1 == g0:
                break
            T2 = T1 - g1 * (T1 - T) / (g1 - g0)
            T2 = min(max(T2, T1 - 30), T1 + 30)
            T, g0, inc = T1, g1, inc1
            T1 = T2
            g1, inc1 = resid(T1, inc)
        if g1 is not None and abs(g1) < 1e-8:
            return T1, inc1
        raise ValueError(f"ammonia-water {'dew' if dew else 'bubble'} point at {p / 1e6:.3f} MPa, x={x:.4f} did not "
                         "converge")

    def bubble_T(self, p: float, x: float) -> tuple[float, float]:
        """Bubble temperature (K) at p (Pa) for a liquid of ammonia mole fraction x, and the vapour's mole fraction."""
        return self._saturation(p, x, dew=False)

    def dew_T(self, p: float, y: float) -> tuple[float, float]:
        """Dew temperature (K) at p (Pa) for a vapour of ammonia mole fraction y, and the liquid's mole fraction."""
        return self._saturation(p, y, dew=True)

    def _saturation_p(self, T: float, x: float, dew: bool) -> tuple[float, float]:
        """Bubble (dew=False: x is the liquid) or dew (dew=True: x is the vapour) pressure (Pa) at T (K), and the
        incipient phase's mole fraction: p <- p * sum(K x) for a bubble point, p <- p / sum(x / K) for a dew point,
        with the incipient composition updated at each step."""
        pa, pw = self._psat_pure("Ammonia", T), self._psat_pure("Water", T)
        zv = np.array([x, 1 - x])
        if dew:
            p = 1.0 / (x / pa + (1 - x) / pw)
            inc = (x / pa) / (x / pa + (1 - x) / pw)
        else:
            p = x * pa + (1 - x) * pw
            inc = x * pa / p
        base, other = ("vapour", "liquid") if dew else ("liquid", "vapour")
        for _ in range(500):
            rB = self.density(T, p, x, base)
            rI = self.density(T, p, inc, other)
            if rB is None or rI is None:
                # a root is missing at this pressure: move towards the phase that is missing
                p *= 0.7 if (rI is None) != dew else 1.3
                continue
            t = zv * np.exp(self.ln_phi(T, rB, x) - self.ln_phi(T, rI, inc))
            S = t.sum()
            new_inc = t[0] / S
            p_new = p * S if not dew else p / S
            if abs(p_new / p - 1) < 1e-10 and abs(new_inc - inc) < 1e-9:
                return p_new, new_inc
            p, inc = p_new, new_inc
        raise ValueError(f"ammonia-water {'dew' if dew else 'bubble'} pressure at {T:.2f} K, x={x:.4f} did not converge")

    def bubble_p(self, T: float, x: float) -> tuple[float, float]:
        """Bubble pressure (Pa) at T (K) for a liquid of mole fraction x, and the vapour's mole fraction."""
        return self._saturation_p(T, x, dew=False)

    def dew_p(self, T: float, y: float) -> tuple[float, float]:
        """Dew pressure (Pa) at T (K) for a vapour of mole fraction y, and the liquid's mole fraction."""
        return self._saturation_p(T, y, dew=True)


class AmmoniaWaterStream:
    """An ammonia-water stream of fixed overall ammonia mass fraction ``w``, with the property interface of the cycles'
    ``Fluid`` (PT, PH, PS, PQ, TQ, Tph, T_ph_fast; C, MPa, kJ/kg). Results are memoised."""

    backend = "teqp (IAPWS G4-01)"

    def __init__(self, w: float, model: AmmoniaWater | None = None):
        if not 0.0 < w < 1.0:
            raise ValueError("the ammonia mass fraction must lie between 0 and 1")
        self.model = model if model is not None else AmmoniaWater()
        self.w = w
        self.x = mole_fraction(w)
        self.name = self.fluid = f"ammonia-water ({w * 100:.1f} % NH3 by mass)"
        self.P_crit = float("nan")
        self.T_crit = float("nan")
        self.fast = False
        self._memo: dict = {}
        self._sat: dict = {}
        self._iso: dict = {}

    def fast_view(self):
        return self

    # ---- fast T(h) along an isobar ------------------------------------------------------

    def _iso_table(self, P: float):
        """T(h) along the isobar at this composition: liquid from 120 K below the bubble point, the two-phase region
        in 1 K steps (each equilibrium warm-started from the last) and vapour to 160 K above the dew point; linear
        between the points. Built once per pressure, for the many lookups of a pinch search."""
        key = round(P, 9)
        tab = self._iso.get(key)
        if tab is not None:
            return tab
        Tb, Td, _, _ = self.bubble_dew(P)
        pts = []
        T = Tb - 120.0
        while T < Tb:
            try:
                st = self._single(T, P, "liquid")
                pts.append((st["h"], st["T"]))
            except ValueError:
                pass
            T += 4.0
        pts.append((self._single(Tb, P, "liquid")["h"], Tb - T_K))
        n = max(int((Td - Tb) / 1.0), 8)
        guess = None
        for i in range(1, n):
            Tt = Tb + (Td - Tb) * i / n
            try:
                st = self._two_phase(Tt, P, guess)
                guess = (st["x_mol"], st["y_mol"])
                pts.append((st["h"], st["T"]))
            except ValueError:
                guess = None
        pts.append((self._single(Td, P, "vapour")["h"], Td - T_K))
        T = Td + 4.0
        while T <= Td + 160.0:
            try:
                st = self._single(T, P, "vapour")
                pts.append((st["h"], st["T"]))
            except ValueError:
                pass
            T += 4.0
        pts.sort()
        H, Ts = [], []
        for h, t in pts:
            if not H or h > H[-1] + 1e-9:
                H.append(h)
                Ts.append(t)
        tab = (H, Ts)
        if len(self._iso) > 200:
            self._iso.clear()
        self._iso[key] = tab
        return tab

    # ---- saturation -------------------------------------------------------------------

    def bubble_dew(self, P: float) -> tuple[float, float, float, float]:
        """Bubble and dew temperatures (K) at P (MPa), with the incipient vapour and liquid mole fractions."""
        key = round(P, 10)
        if key not in self._sat:
            m, p = self.model, P * 1e6
            Tb, yb = m.bubble_T(p, self.x)
            Td, xd = m.dew_T(p, self.x)
            self._sat[key] = (Tb, Td, yb, xd)
        return self._sat[key]

    # ---- states -----------------------------------------------------------------------

    def _single(self, T: float, P: float, phase: str) -> dict:
        m, p = self.model, P * 1e6
        rho = m.density(T, p, self.x, phase)
        if rho is None:
            other = "vapour" if phase == "liquid" else "liquid"
            rho = m.density(T, p, self.x, other)
            if rho is None:
                raise ValueError(f"{self.name}: no single-phase state at {T - T_K:.2f} C, {P:.4f} MPa")
            phase = other
        h, s = m.h_s(T, rho, self.x)
        M = _molar_mass(self.x)
        return dict(T=float(T - T_K), P=P, h=float(h / M / 1000), s=float(s / M / 1000), D=float(rho * M), q=-1.0,
                    phase="liquid" if phase == "liquid" else "superheated_vapour", w_liq=self.w if phase == "liquid"
                    else None, w_vap=self.w if phase != "liquid" else None)

    def _two_phase(self, T: float, P: float, guess: tuple | None = None) -> dict:
        m, p = self.model, P * 1e6
        if guess is None:
            # start between this stream's bubble point (liquid at x, incipient vapour yb) and dew point (vapour at x,
            # incipient liquid xd): far closer than Raoult's law for this strongly non-ideal mixture
            Tb, Td, yb, xd = self.bubble_dew(P)
            f = min(max((T - Tb) / (Td - Tb), 0.0), 1.0) if Td > Tb else 0.5
            guess = (self.x + f * (xd - self.x), yb + f * (self.x - yb))
        eq = m.equilibrium_TP(T, p, *guess)
        if eq is None:
            eq = m.equilibrium_TP(T, p)
        if eq is None:
            raise ValueError(f"{self.name}: no two-phase equilibrium at {T - T_K:.2f} C, {P:.4f} MPa")
        x, y, rL, rV = eq
        V = (self.x - x) / (y - x)                 # vapour mole fraction of the overall stream (lever rule)
        V = min(max(V, 0.0), 1.0)
        hL, sL = m.h_s(T, rL, x)
        hV, sV = m.h_s(T, rV, y)
        ML, MV = _molar_mass(x), _molar_mass(y)
        nL, nV = 1.0 - V, V
        mass = nL * ML + nV * MV
        q = nV * MV / mass
        h = (nL * hL + nV * hV) / mass / 1000
        s = (nL * sL + nV * sV) / mass / 1000
        vol = nL / rL + nV / rV                        # m3 per mol of overall stream
        return dict(T=float(T - T_K), P=P, h=float(h), s=float(s), D=float(mass / vol), q=float(q), phase="two_phase",
                    w_liq=float(mass_fraction(x)), w_vap=float(mass_fraction(y)), x_mol=float(x), y_mol=float(y))

    def _memoised(self, key, fn):
        st = self._memo.get(key)
        if st is None:
            st = fn()
            if len(self._memo) > 100000:
                self._memo.clear()
            self._memo[key] = st
        return st

    def PT(self, P: float, T: float) -> dict:
        def calc():
            Tk = T + T_K
            Tb, Td, _, _ = self.bubble_dew(P)
            if Tk <= Tb:
                return self._single(Tk, P, "liquid")
            if Tk >= Td:
                return self._single(Tk, P, "vapour")
            return self._two_phase(Tk, P)
        return self._memoised(("PT", P, T), calc)

    def _solve_T(self, P: float, prop: str, value: float) -> dict:
        """The state at P whose ``prop`` (h or s) equals ``value``: liquid, two-phase or vapour by the saturated
        values, then a bracketed root in temperature (h and s rise with T at fixed P)."""
        Tb, Td, _, _ = self.bubble_dew(P)
        vb = self._single(Tb, P, "liquid")[prop]
        vd = self._single(Td, P, "vapour")[prop]
        if value <= vb:
            lo, hi, region = Tb - 1.0, Tb, "liquid"
            while self._single(lo, P, "liquid")[prop] > value:
                hi, lo = lo, lo - 25.0
                if lo < 180.0:
                    raise ValueError(f"{self.name}: {prop}={value} below the liquid range at {P} MPa")
        elif value >= vd:
            lo, hi, region = Td, Td + 1.0, "vapour"
            while self._single(hi, P, "vapour")[prop] < value:
                lo, hi = hi, hi + 25.0
                if hi > 1200.0:
                    raise ValueError(f"{self.name}: {prop}={value} above the vapour range at {P} MPa")
        else:
            lo, hi, region = Tb, Td, "two"

        def f(T):
            if region == "two":
                if T <= Tb + 1e-9:
                    return vb - value
                if T >= Td - 1e-9:
                    return vd - value
                return self._two_phase(T, P)[prop] - value
            return self._single(T, P, region)[prop] - value

        T = brentq(f, lo, hi, xtol=1e-10, rtol=1e-13)
        if region == "two":
            return self._two_phase(T, P) if Tb < T < Td else (self._single(T, P, "liquid" if T <= Tb else "vapour"))
        st = self._single(T, P, region)
        return st

    def PH(self, P: float, h: float) -> dict:
        return self._memoised(("PH", P, h), lambda: self._solve_T(P, "h", h))

    def PS(self, P: float, s: float) -> dict:
        return self._memoised(("PS", P, s), lambda: self._solve_T(P, "s", s))

    def PQ(self, P: float, q: float) -> dict:
        """Saturated liquid (q=0), saturated vapour (q=1) or a vapour mass fraction q in between, at P (MPa)."""
        def calc():
            Tb, Td, _, _ = self.bubble_dew(P)
            if q <= 0:
                st = self._single(Tb, P, "liquid")
                st.update(q=0.0, phase="saturated_liquid")
                return st
            if q >= 1:
                st = self._single(Td, P, "vapour")
                st.update(q=1.0, phase="saturated_vapour")
                return st
            T = brentq(lambda T: self._two_phase(T, P)["q"] - q, Tb + 1e-7, Td - 1e-7, xtol=1e-10)
            return self._two_phase(T, P)
        return self._memoised(("PQ", P, q), calc)

    def TQ(self, T: float, q: float) -> dict:
        """The saturation state at T (C) and vapour mass fraction q (0: bubble, 1: dew)."""
        def calc():
            m, Tk = self.model, T + T_K
            p = m.bubble_p(Tk, self.x)[0] if q <= 0 else m.dew_p(Tk, self.x)[0] if q >= 1 else None
            if p is None:
                lo, hi = m.dew_p(Tk, self.x)[0], m.bubble_p(Tk, self.x)[0]
                p = math.exp(brentq(lambda lp: self.PT(math.exp(lp) / 1e6, T)["q"] - q, math.log(lo) + 1e-9,
                                    math.log(hi) - 1e-9, xtol=1e-12))
            return self.PQ(p / 1e6, q)
        return self._memoised(("TQ", T, q), calc)

    def is_gas(self, st: dict) -> bool:
        return st.get("phase") in ("saturated_vapour", "superheated_vapour")

    def Tph(self, P: float, h: float) -> float:
        return self.PH(P, h)["T"]

    def T_ph_fast(self, P: float, h: float) -> float:
        """Temperature (C) at (P, h) from the isobar table (a PH flash outside it)."""
        import bisect
        H, Ts = self._iso_table(P)
        if h < H[0] or h > H[-1]:
            return self.PH(P, h)["T"]
        i = min(bisect.bisect_right(H, h) - 1, len(H) - 2)
        return Ts[i] + (Ts[i + 1] - Ts[i]) * (h - H[i]) / (H[i + 1] - H[i])
