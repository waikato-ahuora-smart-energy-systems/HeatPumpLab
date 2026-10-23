"""CoolProp (v8) property backend that can stand in for REFPROP's PropertyPackage.

``CoolPropPackage`` takes the same fluid strings as ``PropertyPackage`` ("Butane",
"Propane;Isopentane|0.2;0.8|mass") and returns the same dictionaries from ``get_std_properties``:
P (MPa), T (degC), h (kJ/kg), s (kJ/kg/K), D (kg/m3), Cp (kJ/kg/K), q, q_refprop and phase, with
REFPROP's phase labels. Failures raise ``RefpropError`` so existing error handling applies.

Backend: CoolProp HEOS (the same reference equations of state as REFPROP for the common
refrigerants). Reference state: IIR for each component (h = 200 kJ/kg, s = 1 kJ/kg/K, saturated
liquid at 0 C), REFPROP's default for the hydrocarbon refrigerants; water and air keep CoolProp's
default. Absolute h and s may still differ from REFPROP where its fluid file uses another default;
differences of state (and every cycle result) do not depend on the reference.

Mixtures use CoolProp's HEOS mixture model with its own binary interaction parameters, which can
differ from REFPROP's (and are estimated for pairs without fitted data); two-phase flashes for
mixtures are less robust than REFPROP's. Compare against REFPROP before relying on mixture results.

CoolProp is the default backend. Select REFPROP with ``make_property_package(fluid, backend="refprop")``
or the environment variable ``HEATPUMPLAB_PROPERTY_BACKEND=refprop`` (see ``property_backend``).
"""
from __future__ import annotations

import functools
import math
import os
import re

from .properties_package import RefpropError
from .property_cache import active_property_call_cache

try:
    import CoolProp.CoolProp as _CP
except ImportError:  # pragma: no cover - optional dependency
    _CP = None

DEFAULT_VAL = -9999990.0
ENV_VAR = "HEATPUMPLAB_PROPERTY_BACKEND"

# REFPROP fluid-file names used in HeatPumpLab -> CoolProp names (others are passed through,
# CoolProp also accepts many aliases and R-numbers).
NAME_MAP = {
    "butane": "n-Butane", "nbutane": "n-Butane", "n-butane": "n-Butane",
    "isobutane": "IsoButane", "ibutane": "IsoButane",
    "pentane": "n-Pentane", "npentane": "n-Pentane", "n-pentane": "n-Pentane",
    "isopentane": "Isopentane", "ipentane": "Isopentane",
    "neopentane": "Neopentane", "hexane": "n-Hexane", "nhexane": "n-Hexane", "isohexane": "Isohexane",
    "cyclopentane": "Cyclopentane", "cyclohexane": "Cyclohexane", "propane": "Propane",
    "propylen": "Propylene", "propylene": "Propylene", "water": "Water", "air": "Air",
    "co2": "CarbonDioxide", "ammonia": "Ammonia", "nh3": "Ammonia", "dme": "DimethylEther",
    "r1234zee": "R1234ze(E)", "r1234ze(e)": "R1234ze(E)", "r1233zde": "R1233zd(E)",
    "r1233zd(e)": "R1233zd(E)", "r1234yf": "R1234yf", "r1336mzzz": "R1336mzz(Z)",
    "r245fa": "R245fa", "r134a": "R134a", "r600": "n-Butane", "r600a": "IsoButane",
    "r290": "Propane", "r601": "n-Pentane", "r601a": "Isopentane", "r718": "Water", "r744": "CarbonDioxide",
}
_KEEP_DEFAULT_REFERENCE = {"Water", "Air", "CarbonDioxide", "Nitrogen", "Oxygen", "Argon"}

GAS_LIKE = frozenset({"saturated_vapour", "superheated_vapour", "supercritical"})


def coolprop_available() -> bool:
    return _CP is not None


def coolprop_version() -> str | None:
    if _CP is None:
        return None
    try:
        return _CP.get_global_param_string("version")
    except Exception:  # pragma: no cover
        return getattr(__import__("CoolProp"), "__version__", None)


def parse_fluid(fluid: str) -> tuple[list[str], list[float], str]:
    """HeatPumpLab fluid string -> (CoolProp names, fractions, basis 'mol' | 'mass')."""
    parts = fluid.split("|")
    names = [n.strip() for n in parts[0].split(";")]
    if len(parts) == 1:
        if len(names) != 1:
            raise ValueError("a mixture needs fractions: 'A;B|0.3;0.7|mass'")
        fractions = [1.0]
    else:
        fractions = [float(v) for v in parts[1].split(";")]
    basis = parts[2].strip().lower() if len(parts) == 3 else "mol"
    if basis not in ("mol", "mass"):
        raise ValueError("Fluid composition basis must be 'mol' or 'mass'.")
    if len(fractions) != len(names):
        raise ValueError("Fluid composition length does not match component count.")
    if any(v < 0 for v in fractions) or sum(fractions) <= 0:
        raise ValueError("Fluid composition must be non-negative with a positive sum.")
    total = sum(fractions)
    return [coolprop_name(n) for n in names], [v / total for v in fractions], basis


def coolprop_name(name: str) -> str:
    key = name.strip()
    if key.upper().endswith(".FLD"):
        key = key[:-4]
    return NAME_MAP.get(key.lower().replace(" ", ""), key)


def _namespaces():
    return _namespaces_for(id(_CP))


@functools.lru_cache(maxsize=None)
def _namespaces_for(_cp_id: int):
    """Where CoolProp may define its input-pair and phase constants (v6/7: CoolProp.CoolProp and
    CoolProp.constants; v8 nanobind: possibly the package or enum classes on these modules)."""
    spaces = [_CP]
    import importlib
    for mod in ("CoolProp", "CoolProp.constants", "CoolProp._constants"):
        try:
            spaces.append(importlib.import_module(mod))
        except Exception:
            pass
    out = []
    for sp in spaces:
        if sp is None or sp in out:
            continue
        out.append(sp)
        for attr in dir(sp):                     # enum classes, e.g. input_pairs / phases
            if attr.startswith("_"):
                continue
            obj = getattr(sp, attr, None)
            if isinstance(obj, type) and obj not in out:
                out.append(obj)
    return tuple(out)


def _const(name: str, index_function: str | None = None, index_name: str | None = None):
    """A CoolProp constant by name, searched across the namespaces above; then the index function."""
    return _const_for(id(_CP), name, index_function, index_name)


@functools.lru_cache(maxsize=None)
def _const_for(_cp_id: int, name: str, index_function: str | None, index_name: str | None):
    """``_const``, cached per CoolProp module (tests swap in a fake one)."""
    spaces = _namespaces()
    for sp in spaces:
        value = getattr(sp, name, None)
        if value is not None and not callable(value):
            return value
    if index_function:
        for sp in spaces:
            fn = getattr(sp, index_function, None)
            if callable(fn):
                try:
                    return fn(index_name or name)
                except Exception:
                    pass
    raise RefpropError(f"CoolProp constant '{name}' not found (searched {len(spaces)} namespaces)")


def _as_int(value):
    try:
        return int(value)
    except Exception:
        return value


def _brent(f, a: float, b: float, xtol: float = 1e-10, maxiter: int = 100) -> float:
    from scipy.optimize import brentq
    return brentq(f, a, b, xtol=xtol, maxiter=maxiter)


class _MixtureFlash:
    """Robust, fast flashes for a mixture (a port of the HeatPumpLab web apps' CoolProp adapter).

    CoolProp's general mixture flash runs a phase-stability analysis on every call (tens of ms) and can
    misclassify the phase near the critical point. Instead:

    - bubble and dew points at P come from pressure-quality flashes, cached per pressure;
    - the phase follows from comparing the input h, s or T with them; two-phase states are found on the
      quality between the bubble and dew points;
    - single-phase states use phase-imposed flashes, which skip the stability analysis;
    - above the critical pressure the liquid-like and gas-like roots are tried in the order set by the state
      at the critical temperature on that isobar, so the root on the right side of T_crit is used;
    - where the saturation flashes fail (close to the critical point), bubble and dew points are interpolated
      from the traced phase envelope and evaluated from the explicit (density, temperature) equation of state.

    Units: P in MPa, T in degC, h in kJ/kg, s in kJ/kg/K, D in kg/m3 (the package's units)."""

    def __init__(self, backend: str, names: list[str], fractions: list[float], basis: str,
                 T_crit: float, P_crit: float, molar_mass_kg: float):
        def make():
            AS = _CP.AbstractState(backend, "&".join(names))
            (AS.set_mass_fractions if basis == "mass" else AS.set_mole_fractions)(fractions)
            return AS
        self.AS, self.ph, self.SAT = make(), make(), make()
        self.MW = molar_mass_kg
        self.T_crit, self.P_crit = T_crit, P_crit
        self.PQ_I = _const("PQ_INPUTS", "get_input_pair_index")
        self.QT_I = _const("QT_INPUTS", "get_input_pair_index")
        self.PT_I = _const("PT_INPUTS", "get_input_pair_index")
        self.HP_I = _const("HmassP_INPUTS", "get_input_pair_index")
        self.PS_I = _const("PSmass_INPUTS", "get_input_pair_index")
        self.DT_I = _const("DmolarT_INPUTS", "get_input_pair_index")
        self.GAS = _as_int(_const("iphase_gas", "get_phase_index", "phase_gas"))
        self.LIQ = _as_int(_const("iphase_liquid", "get_phase_index", "phase_liquid"))
        self.sat_cache: dict = {}
        self.sc_ref: dict = {}
        self._env = None
        self._env_built = False
        self.fallbacks = 0

    # ---- states --------------------------------------------------------------------

    @staticmethod
    def _pack(AS, phase: str, q: float) -> dict:
        try:
            cp = AS.cpmass() / 1000.0
        except Exception:
            cp = DEFAULT_VAL
        st = {"T": AS.T() - 273.15, "P": AS.p() / 1e6, "h": AS.hmass() / 1000.0, "s": AS.smass() / 1000.0,
              "D": AS.rhomass(), "Cp": cp, "phase": phase, "q": q}
        if not all(math.isfinite(st[k]) for k in ("T", "P", "h", "s", "D")):
            raise RefpropError("nonfinite mixture state")
        return st

    @staticmethod
    def _sat_label(q: float) -> str:
        return "saturated_liquid" if q <= 1e-12 else "saturated_vapour" if q >= 1 - 1e-12 else "two_phase"

    def _imposed(self, phase: str, pair, v1: float, v2: float):
        self.ph.specify_phase(self.GAS if phase == "gas" else self.LIQ)
        try:
            self.ph.update(pair, v1, v2)
        finally:
            self.ph.unspecify_phase()
        return self.ph

    # ---- saturation ------------------------------------------------------------------

    def pq(self, P: float, q: float) -> dict:
        try:
            self.AS.update(self.PQ_I, P * 1e6, q)
            st = self._pack(self.AS, self._sat_label(q), q)
            st["P"] = P
            return st
        except Exception:
            return self._sat_fallback("P", P, q)

    def tq(self, T: float, q: float) -> dict:
        try:
            self.AS.update(self.QT_I, q, T + 273.15)
            st = self._pack(self.AS, self._sat_label(q), q)
            st["T"] = T
            return st
        except Exception:
            if q in (0.0, 1.0):
                st = self._tq_by_pressure(q, T)
                if st is not None:
                    return st
            return self._sat_fallback("T", T, q)

    def sat(self, P: float) -> tuple[dict, dict]:
        key = round(P, 12)
        hit = self.sat_cache.get(key)
        if hit is None:
            hit = (self.pq(P, 0.0), self.pq(P, 1.0))
            if len(self.sat_cache) > 20000:
                self.sat_cache.clear()
            self.sat_cache[key] = hit
        return hit

    def _sat_fallback(self, by: str, x: float, q: float) -> dict:
        if q in (0.0, 1.0):
            return self._env_sat(q, by, x)
        b, d = self._sat_fallback(by, x, 0.0), self._sat_fallback(by, x, 1.0)
        mix = {k: b[k] + q * (d[k] - b[k]) for k in ("T", "P", "h", "s")}
        mix.update(D=1.0 / (1.0 / b["D"] + q * (1.0 / d["D"] - 1.0 / b["D"])), Cp=DEFAULT_VAL,
                   phase="two_phase", q=q)
        return mix

    def _tq_by_pressure(self, q: float, T: float):
        """A temperature-quality point from pressure-quality flashes (these converge closer to the critical
        point), started from the envelope's pressure."""
        env = self.envelope()
        br = env and (env["bub"] if q == 0.0 else env["dew"])
        o = br and self._interp(br, "T", T, ("lnP",))
        if not o:
            return None
        try:
            def f(lnP):
                self.AS.update(self.PQ_I, math.exp(lnP) * 1e6, q)
                return self.AS.T() - 273.15 - T
            x0 = o["lnP"]
            f0 = f(x0)
            x1 = x0 + 0.002
            f1 = f(x1)
            for _ in range(30):
                if abs(f1) <= 1e-7:
                    break
                x0, f0, x1 = x1, f1, x1 - f1 * (x1 - x0) / ((f1 - f0) or 1e-30)
                f1 = f(x1)
            if abs(f1) > 1e-5:
                return None
            self.AS.update(self.PQ_I, math.exp(x1) * 1e6, q)
            st = self._pack(self.AS, self._sat_label(q), q)
            st["T"] = T
            return st
        except Exception:
            return None

    # ---- phase envelope fallback -----------------------------------------------------------

    def envelope(self):
        if self._env_built:
            return self._env
        self._env_built = True
        try:
            self.SAT.build_phase_envelope("")
            d = self.SAT.get_phase_envelope_data()
            T = [t - 273.15 for t in d.T]
            P = [p / 1e6 for p in d.p]
            rb, ri = list(d.rhomolar_vap), list(d.rhomolar_liq)
        except Exception:
            return None
        if len(T) < 10:
            return None
        # the trace starts on the dew line and crosses the critical point, where the bulk and incipient
        # densities meet, onto the bubble line
        ic = next((i for i in range(1, len(T)) if ri[i - 1] - rb[i - 1] > 0 >= ri[i] - rb[i]), -1)
        if ic < 0:
            return None
        f = (ri[ic - 1] - rb[ic - 1]) / ((ri[ic - 1] - rb[ic - 1]) - (ri[ic] - rb[ic]))
        crit = {"T": T[ic - 1] + f * (T[ic] - T[ic - 1]), "P": P[ic - 1] + f * (P[ic] - P[ic - 1]),
                "rho": rb[ic - 1] + f * (rb[ic] - rb[ic - 1])}

        def branch(idx):
            out = []
            for i in idx:      # rising T and P only: drops the retrograde tip and step jitter
                if out and (P[i] <= out[-1]["P"] or T[i] <= out[-1]["T"]):
                    continue
                if rb[i] > 0:
                    out.append({"T": T[i], "P": P[i], "lnP": math.log(P[i]), "lnr": math.log(rb[i])})
            return out

        self._env = {"crit": crit, "dew": branch(range(ic)), "bub": branch(range(len(T) - 1, ic - 1, -1))}
        return self._env

    @staticmethod
    def _interp(br, key: str, x: float, outs):
        if not br or not (br[0][key] <= x <= br[-1][key]):
            return None
        lo, hi = 0, len(br) - 1
        while hi - lo > 1:
            m = (lo + hi) // 2
            if br[m][key] <= x:
                lo = m
            else:
                hi = m
        t = (x - br[lo][key]) / ((br[hi][key] - br[lo][key]) or 1.0)
        return {k: br[lo][k] + t * (br[hi][k] - br[lo][k]) for k in outs}

    def _p_dt(self, rhomol: float, TK: float, phase: str) -> float:
        return self._imposed(phase, self.DT_I, rhomol, TK).p() / 1e6

    def _dt_state(self, rhomol: float, TK: float, phase: str) -> dict:
        return self._pack(self._imposed(phase, self.DT_I, rhomol, TK), phase, DEFAULT_VAL)

    def _rho_at(self, TK: float, P: float, rho0: float, phase: str) -> float:
        r0, r1 = rho0, rho0 * 1.01
        p0, p1 = self._p_dt(r0, TK, phase) - P, self._p_dt(r1, TK, phase) - P
        for _ in range(60):
            if abs(p1) <= 1e-9 * P:
                break
            r2 = r1 - p1 * (r1 - r0) / ((p1 - p0) or 1e-30)
            if not (r2 > 0 and math.isfinite(r2)):
                r2 = r1 * (0.7 if p1 > 0 else 1.3)
            r0, p0, r1 = r1, p1, r2
            p1 = self._p_dt(r1, TK, phase) - P
        if abs(p1) > 1e-6 * P:
            raise RefpropError("density solve did not converge")
        return r1

    def _env_sat(self, q: float, by: str, x: float) -> dict:
        env = self.envelope()
        if not env:
            raise RefpropError("saturation flash failed and no phase envelope is available")
        br = env["bub"] if q == 0.0 else env["dew"]
        o = self._interp(br, "lnP", math.log(x), ("T", "lnr")) if by == "P" else self._interp(br, "T", x, ("lnP", "lnr"))
        if not o:
            raise RefpropError(f"outside the phase envelope ({by} = {x})")
        TK = (o["T"] if by == "P" else x) + 273.15
        P = x if by == "P" else math.exp(o["lnP"])
        phase = "liquid" if q == 0.0 else "gas"
        st = self._dt_state(self._rho_at(TK, P, math.exp(o["lnr"]), phase), TK, phase)
        st.update(P=P, q=q, phase=self._sat_label(q))
        self.fallbacks += 1
        return st

    # ---- single-phase and general flashes --------------------------------------------

    def single(self, phase: str, P: float, key: str, val: float) -> dict:
        label = "superheated_vapour" if phase == "gas" else "liquid"
        try:
            if key == "T":
                AS = self._imposed(phase, self.PT_I, P * 1e6, val + 273.15)
            elif key == "h":
                AS = self._imposed(phase, self.HP_I, val * 1e3, P * 1e6)
            else:
                AS = self._imposed(phase, self.PS_I, P * 1e6, val * 1e3)
            st = self._pack(AS, label, DEFAULT_VAL)
        except Exception:
            if phase == "gas" and key != "T":
                st = self.single("gas", P, "T", self._gas_invert(P, key, val))
            else:
                st = self._dt_single(phase, P, key, val)
        st["P"] = P
        if key == "T":
            st["T"] = val
        return st

    def _gas_invert(self, P: float, key: str, val: float) -> float:
        """T of a gas state at P with the given h or s, from phase-imposed P-T flashes (CoolProp bounds its own
        inversion at the mixture's mole-averaged maximum temperature, which hot discharge gas can pass)."""
        def g(T):
            try:
                AS = self._imposed("gas", self.PT_I, P * 1e6, T + 273.15)
                return (AS.hmass() if key == "h" else AS.smass()) / 1000.0 - val
            except Exception:
                return float("nan")
        lo = self.T_crit - 20 if math.isfinite(self.T_crit) else 100.0
        flo = g(lo)
        while not math.isfinite(flo) and lo < 600:
            lo += 10
            flo = g(lo)
        if not math.isfinite(flo) or flo > 0:
            raise RefpropError("gas state not bracketed")
        hi = lo + 50
        fhi = g(hi)
        while math.isfinite(fhi) and fhi < 0 and hi < 800:
            lo, hi = hi, hi + 50
            fhi = g(hi)
        if not (math.isfinite(fhi) and fhi >= 0):
            raise RefpropError("gas state not bracketed")
        return _brent(g, lo, hi, 1e-9)

    def _dt_single(self, phase: str, P: float, key: str, val: float) -> dict:
        """Single-phase state at P from (density, temperature) evaluations, bracketed beyond the saturation
        density on the requested side, where p rises monotonically with density."""
        bub, dew = self.sat(P)
        edge = dew if phase == "gas" else bub
        rho_e = edge["D"] / self.MW

        def at_T(TC):
            TK = TC + 273.15

            def gfun(lnr):
                return self._p_dt(math.exp(lnr), TK, phase) - P
            lo = math.log(rho_e) - 0.005
            hi = lo + 0.01
            if phase == "gas":
                lo -= 12
            else:
                for _ in range(20):
                    if gfun(hi) > 0:
                        break
                    hi += 0.05
            return self._dt_state(math.exp(_brent(gfun, lo, hi, 1e-13)), TK, phase)

        if key == "T":
            st = at_T(val)
        else:
            lo, hi = (edge["T"], edge["T"] + 300) if phase == "gas" else (edge["T"] - 250, edge["T"])
            st = at_T(_brent(lambda t: at_T(t)[key] - val, lo, hi, 1e-9))
        st.update(P=P, q=DEFAULT_VAL, phase="superheated_vapour" if phase == "gas" else "liquid")
        self.fallbacks += 1
        return st

    def _sc_order(self, P: float, key: str, val: float):
        if key == "T":
            return ("liquid", "gas") if val < self.T_crit else ("gas", "liquid")
        ref = self.sc_ref.get(P, ...)
        if ref is ...:
            ref = None
            for ph in ("gas", "liquid"):
                try:
                    ref = self.single(ph, P, "T", self.T_crit)
                    break
                except Exception:
                    pass
            self.sc_ref[P] = ref
        return ("liquid", "gas") if ref is not None and val < ref[key] else ("gas", "liquid")

    def flash(self, P: float, key: str, val: float):
        """State at P with T, h or s = val; None to ask the caller for CoolProp's general flash."""
        if math.isfinite(self.P_crit) and P >= 0.995 * self.P_crit:
            if P > self.P_crit:
                got = []
                for ph in self._sc_order(P, key, val):
                    try:
                        st = self.single(ph, P, key, val)
                    except Exception:
                        continue
                    st["phase"] = "supercritical"
                    if (st["T"] <= self.T_crit + 2) if ph == "liquid" else (st["T"] >= self.T_crit - 2):
                        return st
                    got.append(st)
                if got:
                    return got[0]
            return None
        try:
            bub, dew = self.sat(P)
        except Exception:
            # no bubble or dew point here (the retrograde region just below the critical pressure, where the
            # traced envelope stops): take the liquid-like or gas-like root by the side of T_crit, as above P_crit
            if not math.isfinite(self.T_crit):
                return None
            for ph in self._sc_order(P, key, val):
                try:
                    st = self.single(ph, P, key, val)
                except Exception:
                    continue
                if (st["T"] <= self.T_crit + 2) if ph == "liquid" else (st["T"] >= self.T_crit - 2):
                    return st
            return None
        if val >= dew[key]:
            return dict(dew) if val == dew[key] else self.single("gas", P, key, val)
        if val <= bub[key]:
            return dict(bub) if val == bub[key] else self.single("liquid", P, key, val)
        q = _brent(lambda x: self.pq(P, x)[key] - val, 0.0, 1.0, 1e-10)
        return self.pq(P, min(max(q, 0.0), 1.0))


_MIXTURE_FLASHES: dict = {}


def _apply_linear_mixing_rules(names) -> None:
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            try:
                cas_a, cas_b = (_CP.get_fluid_param_string(n, "CAS") for n in (a, b))
                _CP.apply_simple_mixing_rule(cas_a, cas_b, "linear")
            except Exception:
                pass          # a pair that already has parameters keeps them
_CRITICAL_POINTS: dict = {}


class CoolPropPackage:
    """CoolProp stand-in for ``PropertyPackage`` (same units, keys and phase labels)."""

    DEFAULT_VAL = DEFAULT_VAL
    TC_to_K = 273.15

    def __init__(self, fld_name: str, pp_instance: str = "", backend: str = "HEOS",
                 reference_state: str | None = "IIR"):
        if _CP is None:
            raise RefpropError("CoolProp is not installed (pip install 'coolprop>=8')")
        self.fld_name = fld_name
        self.pp_instance = pp_instance
        self.names, fractions, self.composition_basis = parse_fluid(fld_name)
        self.n_components = len(self.names)
        self.reference_notes: list[str] = []
        if reference_state:
            for n in self.names:
                if n in _KEEP_DEFAULT_REFERENCE:
                    continue
                try:
                    setter = getattr(_CP, "set_reference_state", None) or getattr(_CP, "set_reference_stateS")
                    setter(n, reference_state)
                except Exception as exc:  # e.g. IIR undefined for a fluid with Tc < 0 C
                    self.reference_notes.append(f"{n}: kept default reference ({exc})")
        try:
            try:
                self.AS = _CP.AbstractState(backend, "&".join(self.names))
            except Exception as exc:
                # no fitted binary interaction parameters for a pair (e.g. isopentane/neopentane): estimate
                # them with CoolProp's linear mixing rule, as the web apps do, and note it
                if self.n_components < 2 or not re.search(r"binary pair|interaction", str(exc), re.I):
                    raise
                _apply_linear_mixing_rules(self.names)
                self.reference_notes.append("no fitted binary interaction parameters: linear mixing rule")
                self.AS = _CP.AbstractState(backend, "&".join(self.names))
            if self.n_components > 1:
                if self.composition_basis == "mass":
                    self.AS.set_mass_fractions(fractions)
                else:
                    self.AS.set_mole_fractions(fractions)
        except Exception as exc:
            raise RefpropError(f"CoolProp could not set up '{fld_name}': {exc}") from exc
        self.fld_composition = list(self.AS.get_mole_fractions()) if self.n_components > 1 else [1.0]
        self.molar_mass = float(self.AS.molar_mass()) * 1000.0      # g/mol
        crit_key = (id(_CP), backend, fld_name)   # a mixture's critical-point search costs tens of ms
        if crit_key not in _CRITICAL_POINTS:
            _CRITICAL_POINTS[crit_key] = self._critical_point()
        self.P_crit, self.T_crit, self.D_crit = _CRITICAL_POINTS[crit_key]
        try:
            self.P_max = float(self.AS.pmax()) / 1e6
        except Exception:
            self.P_max = float("nan")
        self._mix = None
        if self.n_components > 1:          # one flash engine (and phase envelope) per fluid, shared by its streams
            key = (id(_CP), backend, fld_name)
            self._mix = _MIXTURE_FLASHES.get(key)
            if self._mix is None:
                self._mix = _MixtureFlash(backend, self.names, fractions, self.composition_basis,
                                          self.T_crit, self.P_crit, self.molar_mass / 1000.0)
                if not math.isfinite(self.T_crit):
                    env = self._mix.envelope()
                    if env:
                        c = env["crit"]
                        self._mix.T_crit, self._mix.P_crit = c["T"], c["P"]
                _MIXTURE_FLASHES[key] = self._mix
            if not math.isfinite(self.T_crit) and math.isfinite(self._mix.T_crit):
                self.T_crit, self.P_crit = self._mix.T_crit, self._mix.P_crit
        self._pairs = {
            "PT": (_const("PT_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e6, b + 273.15)),
            "TP": (_const("PT_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e6, a + 273.15)),
            "PH": (_const("HmassP_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e3, a * 1e6)),
            "HP": (_const("HmassP_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e3, b * 1e6)),
            "PS": (_const("PSmass_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e6, b * 1e3)),
            "SP": (_const("PSmass_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e6, a * 1e3)),
            "PQ": (_const("PQ_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e6, b)),
            "QP": (_const("PQ_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e6, a)),
            "TQ": (_const("QT_INPUTS", "get_input_pair_index"), lambda a, b: (b, a + 273.15)),
            "QT": (_const("QT_INPUTS", "get_input_pair_index"), lambda a, b: (a, b + 273.15)),
            "TD": (_const("DmassT_INPUTS", "get_input_pair_index"), lambda a, b: (b, a + 273.15)),
            "DT": (_const("DmassT_INPUTS", "get_input_pair_index"), lambda a, b: (a, b + 273.15)),
            "PD": (_const("DmassP_INPUTS", "get_input_pair_index"), lambda a, b: (b, a * 1e6)),
            "DP": (_const("DmassP_INPUTS", "get_input_pair_index"), lambda a, b: (a, b * 1e6)),
            "TS": (_const("SmassT_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e3, a + 273.15)),
            "ST": (_const("SmassT_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e3, b + 273.15)),
        }
        self._phase = {}
        for cname, label in (("iphase_liquid", "liquid"), ("iphase_gas", "superheated_vapour"),
                             ("iphase_supercritical", "supercritical"),
                             ("iphase_supercritical_gas", "superheated_vapour"),
                             ("iphase_supercritical_liquid", "liquid"), ("iphase_twophase", "two_phase")):
            try:
                self._phase[_as_int(_const(cname, "get_phase_index", cname[1:]))] = label
            except RefpropError:
                pass                          # fall back to quality/temperature below

    def _sat_T(self, P_Pa: float) -> float:
        """Dew temperature (degC) at P, for phase inference when CoolProp's phase index is unavailable."""
        if P_Pa >= self.P_crit * 1e6:
            return -273.15
        st = _CP.AbstractState(self.AS.backend_name(), "&".join(self.names)) if self.n_components == 1 else None
        try:
            probe = st or self.AS
            probe.update(self._pairs["PQ"][0], P_Pa, 1.0)
            return probe.T() - 273.15
        except Exception:
            return float("inf")

    # ---- critical point ------------------------------------------------------------

    def _critical_point(self):
        try:
            return (float(self.AS.p_critical()) / 1e6, float(self.AS.T_critical()) - 273.15,
                    float(self.AS.rhomass_critical()))
        except Exception:
            pass
        try:   # mixtures: CoolProp's critical point search
            pts = self.AS.all_critical_points()
            if pts:
                pt = max(pts, key=lambda c: c.T)
                return pt.p / 1e6, pt.T - 273.15, pt.rhomolar * self.molar_mass / 1000.0
        except Exception:
            pass
        return float("nan"), float("nan"), float("nan")

    # ---- properties ----------------------------------------------------------------

    def get_std_properties(self, propsIn_str: str, prop1: float, prop2: float, comp=None) -> dict:
        """PropertyPackage-compatible state dict (MPa, degC, kJ/kg, kJ/kg/K, kg/m3)."""
        cache = active_property_call_cache()
        if cache is None:
            return self._std(propsIn_str, prop1, prop2)
        key = ("coolprop", id(self), self.fld_name, propsIn_str, float(prop1), float(prop2))
        hit, result = cache.lookup(key)
        if hit:
            return dict(result)
        result = self._std(propsIn_str, prop1, prop2)
        cache.store(key, dict(result))
        return dict(result)

    _MIX_KEYS = {"PT": ("P", "T"), "TP": ("T", "P"), "PH": ("P", "h"), "HP": ("h", "P"),
                 "PS": ("P", "s"), "SP": ("s", "P")}

    def _mix_std(self, pair: str, a: float, b: float):
        """Mixture state through the shared flash engine; None when it defers to the general flash."""
        a, b = float(a), float(b)
        if pair in ("PQ", "QP"):
            P, q = (a, b) if pair == "PQ" else (b, a)
            st = self._mix.pq(P, q)
        elif pair in ("TQ", "QT"):
            T, q = (a, b) if pair == "TQ" else (b, a)
            st = self._mix.tq(T, q)
        else:
            first, second = self._MIX_KEYS[pair]
            P, val = (a, b) if first == "P" else (b, a)
            st = self._mix.flash(P, second if first == "P" else first, val)
            if st is None:
                return None
            st = dict(st)
            st["P"] = P
            st[{"T": "T", "H": "h", "S": "s"}[pair.replace("P", "")]] = val      # the input comes back exactly
        phase, q = st["phase"], st["q"]
        if phase not in ("two_phase", "saturated_liquid", "saturated_vapour"):
            q = DEFAULT_VAL
        q_refprop = {"liquid": -1.0, "superheated_vapour": 998.0, "supercritical": 999.0}.get(phase, q)
        cp = st.get("Cp", DEFAULT_VAL)
        return {"P": st["P"], "h": st["h"], "q": q, "q_refprop": q_refprop, "phase": phase, "T": st["T"],
                "s": st["s"], "Cp": cp if cp > 0 else DEFAULT_VAL, "D": st["D"]}

    def _std(self, pair: str, a: float, b: float) -> dict:
        if pair not in self._pairs:
            raise RefpropError(f"Unsupported CoolProp input pair '{pair}'")
        if self._mix is not None and (pair in self._MIX_KEYS or pair in ("PQ", "QP", "TQ", "QT")):
            try:
                st = self._mix_std(pair, a, b)
            except Exception:
                st = None              # the general flash below decides, and reports any failure
            if st is not None:
                return st
        index, convert = self._pairs[pair]
        v1, v2 = convert(float(a), float(b))
        try:
            try:
                self.AS.update(index, v1, v2)
            except Exception as exc:
                P_in = (a if pair[0] == "P" else b) if "P" in pair else None
                if pair in ("PT", "TP") and "aturation" in str(exc):
                    # a P-T state exactly on the saturation line (a pure fluid) is ambiguous to CoolProp;
                    # REFPROP returns the subcooled-liquid root there, so impose the liquid phase and flash again
                    self.AS.specify_phase(_as_int(_const("iphase_liquid", "get_phase_index", "phase_liquid")))
                    try:
                        self.AS.update(index, v1, v2)
                    finally:
                        self.AS.unspecify_phase()
                elif P_in is not None and math.isfinite(self.P_crit) and abs(P_in / self.P_crit - 1) < 1e-9:
                    # CoolProp's single-phase solvers can fail exactly at the critical pressure; a pressure
                    # 1e-9 lower gives the same state (the input pressure is still reported)
                    w1, w2 = convert(*((P_in * (1 - 1e-9), b) if pair[0] == "P" else (a, P_in * (1 - 1e-9))))
                    self.AS.update(index, w1, w2)
                else:
                    raise
            T, P = self.AS.T(), self.AS.p()
            h, s, D = self.AS.hmass(), self.AS.smass(), self.AS.rhomass()
            ph = self.AS.phase()
            try:
                Q = self.AS.Q()
            except Exception:
                Q = float("nan")
            try:
                cp = self.AS.cpmass() / 1000.0
            except Exception:
                cp = DEFAULT_VAL
        except Exception as exc:
            err = RefpropError(f"CoolProp {pair} flash failed for {self.fld_name} at ({a}, {b}): {exc}")
            err.property_context = {"fluid": self.fld_name, "input_pair": pair, "inputs": [a, b],
                                    "backend": "CoolProp"}
            raise err from exc
        if not all(math.isfinite(x) for x in (T, P, h, s, D)):
            raise RefpropError(f"CoolProp {pair} returned a nonfinite state for {self.fld_name}")
        # an input pressure or temperature comes back exactly, as from REFPROP: CoolProp's flashes recompute
        # them and differ in the last digits, which breaks exact-equality checks (a zero pressure drop
        # becoming 5e-10 MPa, say) in the cycle solvers
        if "P" in pair:
            P = (a if pair[0] == "P" else b) * 1e6
        if "T" in pair:
            T = (a if pair[0] == "T" else b) + 273.15
        if "H" in pair:          # likewise an input enthalpy or entropy (J/kg, J/kg/K here)
            h = (a if pair[0] == "H" else b) * 1e3
        if "S" in pair:
            s = (a if pair[0] == "S" else b) * 1e3
        phase = self._phase.get(_as_int(ph), "unknown")
        if phase == "unknown":               # phase constants unavailable: infer from quality
            phase = "two_phase" if math.isfinite(Q) and 0.0 <= Q <= 1.0 else (
                "superheated_vapour" if math.isfinite(self.T_crit) and T - 273.15 >= self._sat_T(P) else "liquid")
        if pair in ("PQ", "QP", "TQ", "QT"):
            q_in = b if pair in ("PQ", "TQ") else a
            Q = q_in
        if phase == "two_phase" or (math.isfinite(Q) and 0.0 <= Q <= 1.0 and pair in ("PQ", "QP", "TQ", "QT")):
            if Q <= 1e-12:
                phase, Q = "saturated_liquid", 0.0
            elif Q >= 1 - 1e-12:
                phase, Q = "saturated_vapour", 1.0
            else:
                phase = "two_phase"
        q = Q if phase in ("two_phase", "saturated_liquid", "saturated_vapour") else DEFAULT_VAL
        q_refprop = {"liquid": -1.0, "superheated_vapour": 998.0, "supercritical": 999.0}.get(phase, q)
        return {"P": P / 1e6, "h": h / 1000.0, "q": q, "q_refprop": q_refprop, "phase": phase,
                "T": T - 273.15, "s": s / 1000.0, "Cp": cp if cp > 0 else DEFAULT_VAL, "D": D}

    # REFPROPdll-style outputs this backend supports (mass basis, the package's units)
    _SPECIFIED = {
        "T": lambda AS: AS.T() - 273.15, "P": lambda AS: AS.p() / 1e6, "H": lambda AS: AS.hmass() / 1000.0,
        "S": lambda AS: AS.smass() / 1000.0, "D": lambda AS: AS.rhomass(), "CP": lambda AS: AS.cpmass() / 1000.0,
        "CV": lambda AS: AS.cvmass() / 1000.0, "VIS": lambda AS: AS.viscosity() * 1e6,
        "TCX": lambda AS: AS.conductivity(), "PRANDTL": lambda AS: AS.Prandtl(),
        "W": lambda AS: AS.speed_sound(), "Q": lambda AS: AS.Q(),
    }

    def get_specified_properties(self, propsOut: str, propsIn: str, prop1: float, prop2: float):
        """REFPROPdll-style selected outputs (``"PRANDTL;VIS;TCX"``) for one state, as a list.

        Viscosity is in uPa s and conductivity in W/(m K), as REFPROP's MASS SI units. Only the outputs in
        ``_SPECIFIED`` are supported; REFPROP-only outputs (ODP, GWP, ...) raise ``RefpropError``."""
        st = self._std(propsIn, prop1, prop2)     # sets the state, with the same input-pair handling
        imposed = False
        if self._mix is not None:                 # the mixture engine does not leave self.AS at the state
            if st["phase"] in ("two_phase", "saturated_liquid", "saturated_vapour"):
                self.AS.update(self._pairs["PQ"][0], st["P"] * 1e6, st["q"])
            else:
                gas = st["phase"] != "liquid" and not (st["phase"] == "supercritical" and st["T"] < self.T_crit)
                self.AS.specify_phase(_as_int(_const("iphase_gas" if gas else "iphase_liquid", "get_phase_index",
                                                     "phase_gas" if gas else "phase_liquid")))
                imposed = True
                self.AS.update(self._pairs["PT"][0], st["P"] * 1e6, st["T"] + 273.15)
        out = []
        try:
            for name in propsOut.split(";"):
                fn = self._SPECIFIED.get(name.strip().upper())
                if fn is None:
                    raise RefpropError(f"CoolProp backend does not provide '{name}'; use the REFPROP backend")
                try:
                    out.append(float(fn(self.AS)))
                except Exception as exc:
                    raise RefpropError(f"CoolProp could not evaluate {name} for {self.fld_name}: {exc}") from exc
        finally:
            if imposed:
                self.AS.unspecify_phase()
        return out

    def _set_RP_lib(self, pp_instance: str) -> None:
        """REFPROP keeps separate library copies per instance; one CoolProp state serves them all."""
        self.pp_instance = pp_instance

    def get_component_molar_masses(self) -> tuple[float, ...]:
        out = []
        for n in self.names:
            out.append(float(_CP.PropsSI("molar_mass", n)) * 1000.0)
        return tuple(out)


# ----------------------------------------------------------------------------------------
# Backend selection
# ----------------------------------------------------------------------------------------


def property_backend(backend: str | None = None) -> str:
    """'refprop' or 'coolprop': the argument, else $HEATPUMPLAB_PROPERTY_BACKEND, else 'coolprop'."""
    value = (backend or os.environ.get(ENV_VAR) or "coolprop").strip().lower()
    if value not in ("refprop", "coolprop"):
        raise ValueError(f"property backend must be 'refprop' or 'coolprop', not '{value}'")
    return value


def make_property_package(fluid: str, pp_instance: str = "", backend: str | None = None):
    """PropertyPackage (REFPROP) or CoolPropPackage for ``fluid``."""
    if property_backend(backend) == "coolprop":
        return CoolPropPackage(fluid, pp_instance)
    from .properties_package import PropertyPackage
    return PropertyPackage(fluid, pp_instance)
