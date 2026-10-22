"""Shared building blocks of the multi-loop dryer heat pumps (TT CO2 cascade, n-level cascade).

Ports of the helpers in the "Advanced dryer heat pumps" web app (``DRYER`` and ``DRYCASC`` scripts), so the Python
models give the same numbers as the app over the same property library:

- ``LoopFluid``: one refrigerant through ``make_property_package`` (CoolProp or REFPROP); PT/PH/PS/PQ/TQ states as
  dicts in MPa, C, kJ/kg, kJ/kg/K, kg/m3, plus a fast T(h) table per isobar for pinch searches.
- ``AirModel``: air enthalpy per kg of dry air as a 1 K table both ways (dry air from the property package, or humid
  air at a fixed humidity ratio from CoolProp's HAPropsSI).
- compressor isentropic efficiency models (HeatPumpLab's six, model 5 = Pareto front v2), ``stage``;
- the optional expander (``expand_to``, ``exp_el``, ``exp_warn``): it replaces an expansion valve, its shaft work is
  credited at the motor efficiency, and a two-phase outlet is flagged;
- exchanger helpers (``air_profile``, ``max_duty``, ``min_gap``) and path builders for T-s / P-h diagrams.

Units: MPa, C, kJ/kg, kJ/kg/K, kg/m3; duties per kg of refrigerant.
"""
from __future__ import annotations

import bisect
import math
from dataclasses import asdict, dataclass, field

from ...components.compressor import PARETO_FRONT_V2_COEFFS
from ...utils.coolprop_properties import make_property_package, property_backend

GAS_PHASES = frozenset({"superheated_vapour", "supercritical", "saturated_vapour"})
DEFAULT_MOTOR_EFF = 0.93
P_ATM_KPA = 101.325


def js_round(x: float, decimals: int = 0) -> float:
    """Math.round(x * 10**d) / 10**d (half up, as the app rounds; Python's round() is half-even)."""
    f = 10.0 ** decimals
    return math.floor(x * f + 0.5) / f


def js_num(x) -> str:
    """A number as JavaScript prints it in a template string (250, not 250.0), for messages that match the app's."""
    x = float(x)
    return str(int(x)) if x.is_integer() and abs(x) < 1e15 else repr(x)


def lin(a: float, b: float, n: int) -> list[float]:
    return [a + (b - a) * i / (n - 1) for i in range(n)]


# ----------------------------------------------------------------------------------------
# Refrigerant
# ----------------------------------------------------------------------------------------


class LoopFluid:
    """One refrigerant of a loop: property states as dicts {T, P, h, s, D, q, phase}.

    ``fast_view()`` returns a view whose ``Tph`` reads a T(h) table built once per isobar (PT on a 0.5 K grid from
    -40 to 450 C, plus the saturation points below the critical pressure, linear in h), as the app does for the many
    temperature lookups of pinch searches. Outside the table, or with the table off, ``Tph`` is a PH flash.
    """

    def __init__(self, fluid: str, backend: str | None = None, package=None):
        self.fluid = fluid
        self.backend = property_backend(backend) if package is None else getattr(package, "backend", "custom")
        self.pp = package if package is not None else make_property_package(fluid, backend=self.backend)
        self.P_crit = float(getattr(self.pp, "P_crit", float("nan")))
        self.T_crit = float(getattr(self.pp, "T_crit", float("nan")))
        self.fast = False
        self._cache: dict = {}
        self._iso: dict = {}

    def fast_view(self) -> "LoopFluid":
        if self.fast:
            return self
        view = object.__new__(LoopFluid)
        view.__dict__.update(self.__dict__)      # shares the property package and both caches
        view.fast = True
        return view

    # ---- states ------------------------------------------------------------------------

    def _state(self, pair: str, a: float, b: float) -> dict:
        key = (pair, a, b)
        hit = self._cache.get(key)
        if hit is not None:
            return hit
        d = self.pp.get_std_properties(pair, a, b, getattr(self.pp, "fld_composition", None))
        st = {k: d.get(k) for k in ("T", "P", "h", "s", "D", "q", "phase")}
        if len(self._cache) > 200000:
            self._cache.clear()
        self._cache[key] = st
        return st

    def PT(self, P, T):
        return self._state("PT", P, T)

    def PH(self, P, h):
        return self._state("PH", P, h)

    def PS(self, P, s):
        return self._state("PS", P, s)

    def PQ(self, P, q):
        return self._state("PQ", P, q)

    def TQ(self, T, q):
        return self._state("TQ", T, q)

    def is_gas(self, st: dict) -> bool:
        """Saturated or superheated vapour, or supercritical (CoolProp labels a state above the critical pressure but
        below the critical temperature 'liquid'; such a state is not compressor suction)."""
        if st.get("phase") in GAS_PHASES:
            return True
        return math.isfinite(self.P_crit) and st["P"] >= self.P_crit and st["T"] >= self.T_crit

    # ---- temperature along an isobar ---------------------------------------------------

    def _iso_table(self, P: float):
        key = f"{P:.7f}"
        t = self._iso.get(key)
        if t is not None:
            return t
        pts = []
        T = -40.0
        while T <= 450.0:
            try:
                st = self.PT(P, T)
                if math.isfinite(st["h"]):
                    pts.append((st["h"], T))
            except Exception:
                pass
            T += 0.5
        if math.isfinite(self.P_crit) and P < self.P_crit:
            for q in (0.0, 1.0):
                try:
                    st = self.PQ(P, q)
                    pts.append((st["h"], st["T"]))
                except Exception:
                    pass
        pts.sort(key=lambda p: p[0])
        H, TT = [], []
        for h, tt in pts:
            if not H or h > H[-1] + 1e-9:
                H.append(h)
                TT.append(tt)
        t = (H, TT, len(H) > 20)
        if len(self._iso) > 400:
            self._iso.clear()
        self._iso[key] = t
        return t

    def Tph(self, P: float, h: float) -> float:
        if not self.fast:
            return self.PH(P, h)["T"]
        H, TT, ok = self._iso_table(P)
        if not ok or h < H[0] or h > H[-1]:
            return self.PH(P, h)["T"]
        hi = bisect.bisect_right(H, h)
        hi = min(max(hi, 1), len(H) - 1)
        lo = hi - 1
        return TT[lo] + (TT[hi] - TT[lo]) * (h - H[lo]) / (H[hi] - H[lo])


# ----------------------------------------------------------------------------------------
# Air
# ----------------------------------------------------------------------------------------


class AirModel:
    """Air enthalpy (kJ per kg of dry air) against temperature, tabled at 1 K from -5 to 280 C, linear both ways
    (clamped at the ends), as the app's ``DRYER.airModel``."""

    def __init__(self, h_of_T, W: float = 0.0, label: str = "air"):
        self.T_grid = lin(-5.0, 280.0, 286)
        self.h_grid = [h_of_T(T) for T in self.T_grid]
        self.W = W
        self.label = label
        order = sorted(range(len(self.h_grid)), key=lambda i: self.h_grid[i])
        self._hT = (self.T_grid, self.h_grid)
        self._Th = ([self.h_grid[i] for i in order], [self.T_grid[i] for i in order])

    @staticmethod
    def _interp(X, Y, x):
        if x <= X[0]:
            return Y[0]
        if x >= X[-1]:
            return Y[-1]
        hi = bisect.bisect_right(X, x)
        hi = min(max(hi, 1), len(X) - 1)
        lo = hi - 1
        return Y[lo] + (Y[hi] - Y[lo]) * (x - X[lo]) / ((X[hi] - X[lo]) or 1.0)

    def h(self, T: float) -> float:
        return self._interp(*self._hT, T)

    def T(self, h: float) -> float:
        return self._interp(*self._Th, h)


def dry_air(P_kPa: float = P_ATM_KPA, backend: str | None = None) -> AirModel:
    """Dry air through the property package (CoolProp or REFPROP)."""
    pp = make_property_package("Air", "snk", backend=backend)
    comp = getattr(pp, "fld_composition", None)
    P = P_kPa / 1000.0
    return AirModel(lambda T: pp.get_std_properties("PT", P, T, comp)["h"], 0.0, "dry air")


def humid_air(T_in: float, RH: float, P_kPa: float = P_ATM_KPA) -> AirModel:
    """Humid air with the humidity ratio of air at T_in and RH (%), per kg of dry air (CoolProp HAPropsSI; the
    moisture above saturation at T_in is capped there). The app's dryer air (inlet 15 C, 60 % RH by default)."""
    try:
        import CoolProp.CoolProp as CP
    except ImportError as exc:  # pragma: no cover - CoolProp is a dependency
        raise RuntimeError("humid air needs CoolProp (HAPropsSI)") from exc
    P = P_kPa * 1000.0
    W = CP.HAPropsSI("W", "T", T_in + 273.15, "P", P, "R", min(max(RH, 0.0), 100.0) / 100.0)
    return AirModel(lambda T: CP.HAPropsSI("H", "T", T + 273.15, "P", P, "W", W) / 1000.0, W, f"humid air ({RH:g} % RH)")


def air_secondary(air: AirModel, T_in: float, T_mt: float, phi: float, T_max: float = math.inf) -> AirModel:
    """Primary + secondary air as one stream per kg of primary air: both are heated from T_in to T_mt (where the
    secondary leaves), so below T_mt the enthalpy rises (1 + phi) times as fast. phi <= 0, or T_mt outside
    (T_in, T_max), returns the air unchanged (the app's ``DRYER.airSecondary``)."""
    if not (phi > 0 and T_in < T_mt < T_max):
        return air
    h0, hm, f = air.h(T_in), air.h(T_mt), 1.0 + phi
    comp = AirModel(lambda T: h0 + f * (air.h(T) - h0) if T <= T_mt else h0 + f * (hm - h0) + air.h(T) - hm,
                    air.W, f"{air.label} + secondary")
    comp.secondary = dict(phi=phi, T_mt=T_mt, T_in=T_in)
    return comp


# ----------------------------------------------------------------------------------------
# Compressor and expander
# ----------------------------------------------------------------------------------------

PF_V1 = (0.87096, 0.01991, 2.63312)


def eta_is(pr: float, model: int = 5, eta_const: float = 0.7 / 0.96, P_in: float | None = None) -> float:
    """Isentropic (shaft) efficiency at a stage pressure ratio, clipped to 0.4-1.0: 0 ideal, 1 constant eta_const,
    2 polynomial (Wang & Cleland 2018), 3 Brendel et al. (needs the suction pressure, MPa), 4 Pareto front v1,
    5 Pareto front v2 (HeatPumpLab's and the app's default)."""
    pr = max(pr, 1.0)
    L = math.log(pr)
    if model == 0:
        e = 1.0
    elif model == 1:
        e = min(max(eta_const, 0.4), 1.0)
    elif model == 2:
        e = (0.65 + 0.015 * pr - 0.0015 * pr * pr) / 0.96
    elif model == 3:
        e = 0.66981 - 0.6 / (pr - 0.01466) ** (8.38 * (P_in if P_in is not None else 0.2)) - 0.00102 * pr ** 1.8 \
            if pr > 1 else 0.4
    elif model == 4:
        e1, a, b = PF_V1
        e = e1 - a * L ** b
    else:
        e1, a, b, d, L0, w = PARETO_FRONT_V2_COEFFS["isentropic"]
        e = e1 - a * L ** b - d / (1 + math.exp(-(L - L0) / w))
    return min(max(e, 0.4), 1.0)


def stage(fl: LoopFluid, P_in: float, st: dict, P_out: float, eta: float | None = None,
          model: int = 5, eta_const: float = 0.7 / 0.96) -> dict:
    """One compression from state st at P_in to P_out: {w, eta, out, pr}."""
    e = eta if eta is not None else eta_is(P_out / P_in, model, eta_const, P_in)
    h_is = fl.PS(P_out, st["s"])["h"]
    w = (h_is - st["h"]) / e
    return dict(w=w, eta=e, out=fl.PH(P_out, st["h"] + w), pr=P_out / P_in)


def expand_to(fl: LoopFluid, P_in: float, h_in: float, P_out: float, eta: float | None) -> float:
    """Outlet enthalpy: an isenthalpic valve (eta None), or an expander of isentropic efficiency eta."""
    if eta is None or not eta > 0:
        return h_in
    return h_in - eta * (h_in - fl.PS(P_out, fl.PH(P_in, h_in)["s"])["h"])


def exp_el(w: float, wx: float, motor: float) -> float:
    """Electrical input: compressor shaft work over the motor efficiency, less the expander work at it."""
    return w / motor - wx * motor


def exp_warn(fl: LoopFluid, P_out: float, h_out: float, warnings: list, name: str = "expander"):
    """Expander outlet quality if two-phase (and a warning), else None."""
    try:
        o = fl.PH(P_out, h_out)
    except Exception:
        return None
    q = o.get("q")
    wet = o.get("phase") == "two_phase" or (q is not None and 0 < q < 1)
    if wet:
        warnings.append(f"two-phase expansion: the {name} outlet is {q * 100:.0f} % vapour ({o['T']:.1f} C); "
                        "two-phase expanders are less proven than gas turbines and their efficiency is uncertain")
    return q if wet else None


# ----------------------------------------------------------------------------------------
# Exchangers and paths
# ----------------------------------------------------------------------------------------


def air_profile(fl: LoopFluid, air: AirModel, P: float, h_hot: float, h_cold: float, Ta_in: float, Ta_out: float,
                n: int = 16) -> list:
    """A stream cooled along an isobar from h_hot to h_cold against air heated Ta_in -> Ta_out (counterflow), per kg
    of the stream: [[q, T_ref, T_air]] from the cold end."""
    ha0, ha1, Q = air.h(Ta_in), air.h(Ta_out), h_hot - h_cold
    return [[f * Q, fl.Tph(P, h_cold + f * Q), air.T(ha0 + f * (ha1 - ha0))] for f in lin(0.0, 1.0, n)]


def min_gap(profile: list) -> float:
    return min(p[1] - p[2] for p in profile)


def max_duty(prH: LoopFluid, PH: float, hHin: float, mh: float, prC: LoopFluid, PC: float, hCin: float, dT: float,
             n: int = 14) -> dict:
    """Counterflow exchanger, hot flow mh per kg of the cold flow: the largest duty per kg cold that keeps dT
    (scanned down from the end-limited duty, as the pinch can sit inside near the pseudo-critical point)."""
    ThIn, TcIn = prH.Tph(PH, hHin), prC.Tph(PC, hCin)
    if ThIn - TcIn < dT:
        return dict(q=0.0, gap=ThIn - TcIn)
    fs = lin(0.0, 1.0, n)

    def gap(q):
        m = math.inf
        for f in fs:
            Tc = prC.Tph(PC, hCin + f * q)
            Th = prH.Tph(PH, hHin - q / mh + f * q / mh)
            m = min(m, Th - Tc)
        return m - dT

    qmax = min(mh * (hHin - prH.PT(PH, TcIn + dT)["h"]), prC.PT(PC, ThIn - dT)["h"] - hCin)
    if not qmax > 0:
        return dict(q=0.0, gap=ThIn - TcIn)
    if gap(qmax) >= 0:
        return dict(q=qmax)
    lo, hi = 0.0, qmax
    for i in range(39, -1, -1):
        q = qmax * i / 40
        if q == 0 or gap(q) >= 0:
            lo, hi = q, qmax * (i + 1) / 40
            break
    k = 0
    while k < 40 and hi - lo > 1e-4 * qmax:
        m = (lo + hi) / 2
        if gap(m) >= 0:
            lo = m
        else:
            hi = m
        k += 1
    return dict(q=lo)


def _pt(st: dict) -> list:
    return [st["T"], st["P"] * 10, st["h"], st["s"]]


def iso_path(fl: LoopFluid, P: float, h0: float, h1: float, n: int = 14) -> list:
    """[[T C, P bar, h, s]] along an isobar (points that fail are skipped)."""
    out = []
    for i in range(n + 1):
        try:
            out.append(_pt(fl.PH(P, h0 + (h1 - h0) * i / n)))
        except Exception:
            pass
    return out


def comp_path(fl: LoopFluid, P0: float, h0: float, P1: float, h1: float, n: int = 8) -> list:
    """[[T C, P bar, h, s]] of a compression or expansion (log pressure, linear enthalpy)."""
    out = []
    for i in range(n + 1):
        f = i / n
        try:
            out.append(_pt(fl.PH(P0 * (P1 / P0) ** f, h0 + f * (h1 - h0))))
        except Exception:
            pass
    return out


def plain_state(fl: LoopFluid, P: float, h: float) -> dict:
    x = fl.PH(P, h)
    return dict(T=x["T"], P=x["P"], h=x["h"], s=x["s"])


def source_T(T_in: float, T_out: float | None, f: float) -> float:
    """Source water temperature at a fraction f of the evaporator duty from its cold end (linear)."""
    T_out = T_in if T_out is None else T_out
    return T_out + f * (T_in - T_out)


@dataclass
class LoopResult:
    """Result of a multi-loop dryer heat pump. Duties per kg of the cycle's basis flow (see each model).

    paths: {loop name: [{"name": component, "pts": [[T C, P bar, h, s], ...]}]} for T-s / P-h diagrams;
    exchangers: [{"name", "hot": [[q, T]], "cold": [[q, T]], ...}] temperature-duty profiles.
    """

    success: bool
    message: str
    COP: float | None = None          # heat to air / electrical input (compressors less any expander credit)
    VHC: float | None = None          # MJ/m3, heat to air per m3 of suction gas, all compressors
    Pmax_bar: float | None = None
    T_dis_max: float | None = None
    q_air: float | None = None        # kJ per kg of basis flow
    w_shaft: float | None = None
    q_source: float | None = None
    report: dict = field(default_factory=dict)
    states: dict = field(default_factory=dict)
    paths: dict = field(default_factory=dict)
    exchangers: list = field(default_factory=list)
    warnings: list = field(default_factory=list)

    def to_dict(self):
        return asdict(self)


def compass_search(f, starts, lo, hi, step0, tol, state: dict, max_eval: int):
    """The app's shared compass search: f(x) -> {y, r} (y = -inf when infeasible). Runs from each start in turn and
    keeps the best: each coordinate stepped both ways, a successful step repeated while it pays, steps halved when
    nothing moves. Returns (x, best, converged)."""
    dim = len(lo)

    def clip(x):
        return [min(max(v, lo[i]), hi[i]) for i, v in enumerate(x)]

    X, B, conv = None, dict(y=-math.inf), False
    for xs0, bs0 in starts:
        if xs0 is None or not math.isfinite(bs0["y"]):
            continue
        x, best, step = list(xs0), bs0, list(step0)

        def live():
            return any(s > tol[i] for i, s in enumerate(step))

        while state["nEval"] < max_eval and live():
            moved = False
            for i in range(dim):
                if state["nEval"] >= max_eval:
                    break
                if step[i] <= tol[i]:
                    continue
                for d in (1, -1):
                    z = list(x)
                    z[i] += d * step[i]
                    zc = clip(z)
                    if zc[i] == x[i]:
                        continue
                    o = f(zc)
                    if o["y"] > best["y"] + 1e-7:
                        best, x, moved = o, zc, True
                        for _ in range(20):
                            if state["nEval"] >= max_eval:
                                break
                            w = list(x)
                            w[i] += d * step[i]
                            wc = clip(w)
                            if wc[i] == x[i]:
                                break
                            ow = f(wc)
                            if ow["y"] > best["y"] + 1e-7:
                                best, x = ow, wc
                            else:
                                break
                        break
            if not moved:
                step = [s / 2 if s > tol[i] else s for i, s in enumerate(step)]
        if best["y"] > B["y"]:
            X, B, conv = x, best, not live()
    return X, B, conv
