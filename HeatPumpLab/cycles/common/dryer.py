"""Shared building blocks of the dryer heat pump cycles ported from the "Advanced dryer heat pumps" web app
(the page's DRYER and DRYCASC scripts): property access, the dryer-air model, compressors and expanders, the
exchanger helpers and the page's compass search.

The cycle modules (``TTCascadeHeatPump``, ``CascadeHeatPump``, ``MVRHeatPump``, ``BraytonHeatPump``) follow the
page's JavaScript closely, so that both give the same results on the same property data (CoolProp 8 in the page;
CoolProp or REFPROP here, through ``make_property_package``). Units follow HeatPumpLab: MPa, C, kJ/kg, kJ/kg/K,
kg/m3, per kg of the reference flow.
"""
from __future__ import annotations

import bisect
import math
from dataclasses import asdict, dataclass, field
from typing import Callable, Sequence

from ...components.compressor import pareto_front_efficiency
from ...utils.coolprop_properties import make_property_package, property_backend
from ...utils.properties_package import RefpropError

GAS_PHASES = frozenset({"saturated_vapour", "superheated_vapour", "supercritical"})
TWO_PHASE = frozenset({"two_phase", "saturated_liquid", "saturated_vapour"})
MOTOR = 0.93                       # default motor efficiency of the page
DEFAULT_MOTOR_EFF = MOTOR
P_ATM = 101325.0                   # Pa
P_ATM_KPA = P_ATM / 1000.0

# errors a cycle solve turns into a failed result (as the page's try/catch does)
SOLVE_ERRORS = (RefpropError, ValueError, ZeroDivisionError, ArithmeticError)


def lin(a: float, b: float, n: int) -> list[float]:
    """n evenly spaced values from a to b (inclusive), evaluated as the page does."""
    return [a + (b - a) * i / (n - 1) for i in range(n)]


def js_round(x: float, decimals: int = 0):
    """JavaScript's Math.round(x * 10**d) / 10**d (halves up; Python's round() is half-even). An int for d = 0."""
    if decimals == 0:
        return math.floor(x + 0.5)
    f = 10.0 ** decimals
    return math.floor(x * f + 0.5) / f


def js_num(x) -> str:
    """A number as JavaScript prints it in a template string (250, not 250.0), for messages that match the app's."""
    x = float(x)
    return str(int(x)) if x.is_integer() and abs(x) < 1e15 else repr(x)


class Table:
    """Monotone table y(x), linear between the points, clamped at the ends (the page's ``table``)."""

    def __init__(self, xs: Sequence[float], ys: Sequence[float]):
        pts = sorted(zip(xs, ys), key=lambda p: p[0])
        self.X = [p[0] for p in pts]
        self.Y = [p[1] for p in pts]

    def __call__(self, x: float) -> float:
        X, Y = self.X, self.Y
        if x <= X[0]:
            return Y[0]
        if x >= X[-1]:
            return Y[-1]
        lo = bisect.bisect_right(X, x) - 1
        hi = lo + 1
        return Y[lo] + (Y[hi] - Y[lo]) * (x - X[lo]) / ((X[hi] - X[lo]) or 1)


# ----------------------------------------------------------------------------------------
# Dryer air
# ----------------------------------------------------------------------------------------


class AirModel:
    """Dryer air as h(T) and T(h) tables (kJ/kg dry air, C), 1 K steps from -5 to 280 C, linear between them and
    clamped at the ends (the page's ``DRYER.airModel``). ``air.h(T)`` and ``air.T(h)`` look them up.

    ``AirModel.humid`` is the page's air: humid air (CoolProp HAPropsSI) at a fixed moisture content, set from the
    inlet temperature and relative humidity (default 15 C, 60 %, 101.325 kPa); enthalpy per kg of dry air, so air
    flows are kg of dry air. ``AirModel.dry`` uses W = 0, and ``dry_air`` dry air from the property package. Any
    h(T) function can be given instead.
    """

    def __init__(self, h_of_T: Callable[[float], float], W: float = 0.0, label: str = "air"):
        self.T_grid = lin(-5.0, 280.0, 286)
        self.h_grid = [h_of_T(t) for t in self.T_grid]
        self.h = Table(self.T_grid, self.h_grid)
        self.T = Table(self.h_grid, self.T_grid)
        self.W = W
        self.label = label

    @classmethod
    def humid(cls, T_in: float = 15.0, RH: float = 60.0, P: float = P_ATM, W: float | None = None,
              label: str | None = None) -> "AirModel":
        """Humid air with the moisture of air at T_in and RH % (or the moisture content W kg/kg dry air); P in Pa."""
        try:
            from CoolProp.CoolProp import HAPropsSI
        except ImportError as exc:   # pragma: no cover - CoolProp is a dependency
            raise RefpropError("the humid-air model needs CoolProp (pip install 'coolprop>=8')") from exc
        if W is None:
            rh = min(max(RH, 0.0), 100.0) / 100.0
            W = HAPropsSI("W", "T", T_in + 273.15, "P", P, "R", rh) if rh > 0 else 0.0
            try:
                W = min(W, HAPropsSI("W", "T", T_in + 273.15, "P", P, "R", 1.0))   # above saturation: saturated
            except ValueError:
                pass
        return cls(lambda t: HAPropsSI("H", "T", t + 273.15, "P", P, "W", W) / 1000.0, W,
                   label or (f"humid air ({RH:g} % RH)" if W else "dry air"))

    @classmethod
    def dry(cls, P: float = P_ATM) -> "AirModel":
        return cls.humid(P=P, W=0.0)


def humid_air(T_in: float, RH: float, P_kPa: float = P_ATM_KPA) -> AirModel:
    """Humid air with the humidity ratio of air at T_in and RH (%), per kg of dry air (CoolProp HAPropsSI): the
    app's dryer air (inlet 15 C, 60 % RH by default)."""
    return AirModel.humid(T_in, RH, P_kPa * 1000.0, label=f"humid air ({RH:g} % RH)")


def dry_air(P_kPa: float = P_ATM_KPA, backend: str | None = None) -> AirModel:
    """Dry air through the property package (CoolProp or REFPROP)."""
    pp = make_property_package("Air", "snk", backend=backend)
    comp = getattr(pp, "fld_composition", None)
    P = P_kPa / 1000.0
    return AirModel(lambda T: pp.get_std_properties("PT", P, T, comp)["h"], 0.0, "dry air")


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
# Working fluids
# ----------------------------------------------------------------------------------------


def _finite(x) -> bool:
    return x is not None and isinstance(x, (int, float)) and math.isfinite(x)


class Fluid:
    """A working fluid with the page's property interface: PT, PH, PS, PQ and TQ return dicts with T (C), P (MPa),
    h (kJ/kg), s (kJ/kg/K), D (kg/m3), q (quality, or a negative number outside the dome) and phase.

    Properties come from ``make_property_package`` (CoolProp or REFPROP: ``backend``), or from a ready
    PropertyPackage-like ``package``. Results are memoised per fluid.

    The page's fast temperature lookup along an isobar (``fastT``) is a T(h) table built once per pressure from PT
    states on a 0.5 K grid from -40 to 450 C, plus the saturation points below the critical pressure, linear in h,
    for the many temperature lookups of the pinch searches: ``T_ph_fast`` always uses it, and ``Tph`` uses it on a
    ``fast_view()`` (a view sharing the package and caches) and is a PH flash otherwise.
    """

    def __init__(self, fluid: str, backend: str | None = None, package=None):
        self.name = self.fluid = fluid
        self.backend = property_backend(backend) if package is None else getattr(package, "backend", "custom")
        self.pp = package if package is not None else make_property_package(fluid, backend=self.backend)
        self.P_crit = float(getattr(self.pp, "P_crit", float("nan")))
        self.T_crit = float(getattr(self.pp, "T_crit", float("nan")))
        self.fast = False
        self._memo: dict = {}
        self._iso: dict = {}

    def fast_view(self) -> "Fluid":
        if self.fast:
            return self
        view = object.__new__(type(self))
        view.__dict__.update(self.__dict__)      # shares the property package and both caches
        view.fast = True
        return view

    def _flash(self, pair: str, a: float, b: float) -> dict:
        key = (pair, a, b)
        st = self._memo.get(key)
        if st is None:
            r = self.pp.get_std_properties(pair, a, b, getattr(self.pp, "fld_composition", None))
            q = r.get("q")
            st = {"T": r["T"], "P": r["P"], "h": r["h"], "s": r["s"], "D": r["D"], "phase": r.get("phase"),
                  "q": q if r.get("phase") in TWO_PHASE and _finite(q) else -1.0}
            if not all(_finite(st[k]) for k in ("T", "h", "s", "D")):
                raise RefpropError(f"{self.name}: nonfinite {pair} state at ({a}, {b})")
            if len(self._memo) > 200000:
                self._memo.clear()
            self._memo[key] = st
        return st

    def PT(self, P: float, T: float) -> dict:
        return self._flash("PT", P, T)

    def PH(self, P: float, h: float) -> dict:
        return self._flash("PH", P, h)

    def PS(self, P: float, s: float) -> dict:
        return self._flash("PS", P, s)

    def PQ(self, P: float, q: float) -> dict:
        return self._flash("PQ", P, q)

    def TQ(self, T: float, q: float) -> dict:
        return self._flash("TQ", T, q)

    def is_gas(self, st: dict) -> bool:
        """Saturated or superheated vapour, or supercritical, or above both the critical pressure and temperature
        (a state CoolProp may label otherwise); suitable compressor suction."""
        if st.get("phase") in GAS_PHASES:
            return True
        return math.isfinite(self.P_crit) and st["P"] >= self.P_crit and st["T"] >= self.T_crit

    # ---- fast T(h) along an isobar --------------------------------------------------

    def _iso_table(self, P: float):
        key = round(P, 7)
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
        H, Ts = [], []
        for h, tt in pts:
            if not H or h > H[-1] + 1e-9:
                H.append(h)
                Ts.append(tt)
        t = (H, Ts, len(H) > 20)
        if len(self._iso) > 400:
            self._iso.clear()
        self._iso[key] = t
        return t

    def T_ph_fast(self, P: float, h: float) -> float:
        H, Ts, ok = self._iso_table(P)
        if not ok or h < H[0] or h > H[-1]:
            return self.PH(P, h)["T"]
        lo = min(bisect.bisect_right(H, h) - 1, len(H) - 2)
        hi = lo + 1
        return Ts[lo] + (Ts[hi] - Ts[lo]) * (h - H[lo]) / (H[hi] - H[lo])

    def Tph(self, P: float, h: float) -> float:
        return self.T_ph_fast(P, h) if self.fast else self.PH(P, h)["T"]


def T_ph(fl: Fluid, P: float, h: float, fast: bool = False) -> float:
    """Temperature at (P, h): the fast isobar table when ``fast`` or on a fast view (the page's ``Tph``)."""
    return fl.T_ph_fast(P, h) if fast else fl.Tph(P, h)


def is_gas(st: dict) -> bool:
    return st.get("phase") in GAS_PHASES


def point(st: dict) -> list[float]:
    """A path point [T C, P bar, h kJ/kg, s kJ/kg/K]."""
    return [st["T"], st["P"] * 10, st["h"], st["s"]]


def plain_state(fl: Fluid, P: float, h: float) -> dict:
    x = fl.PH(P, h)
    return dict(T=x["T"], P=x["P"], h=x["h"], s=x["s"])


def iso_path(fl: Fluid, P: float, h0: float, h1: float, n: int = 14) -> list:
    """[[T C, P bar, h, s]] along an isobar (points that fail are skipped)."""
    out = []
    for i in range(n + 1):
        try:
            out.append(point(fl.PH(P, h0 + (h1 - h0) * i / n)))
        except Exception:
            pass
    return out


def comp_path(fl: Fluid, P0: float, h0: float, P1: float, h1: float, n: int = 8) -> list:
    """[[T C, P bar, h, s]] of a compression or expansion (log pressure, linear enthalpy)."""
    out = []
    for i in range(n + 1):
        f = i / n
        try:
            out.append(point(fl.PH(P0 * (P1 / P0) ** f, h0 + f * (h1 - h0))))
        except Exception:
            pass
    return out


# ----------------------------------------------------------------------------------------
# Compressors and expanders
# ----------------------------------------------------------------------------------------

COMP_MODELS = ("Ideal (1.0)", "Constant", "Polynomial (Wang & Cleland)", "Brendel et al.", None,   # 4 removed
               "Pareto front")


def eta_model(model: int, eta_const: float, pressure_ratio: float, P_in: float | None = None) -> float:
    """Isentropic (shaft) efficiency of HeatPumpLab compressor model ``model`` (COMP_MODELS) at a stage pressure
    ratio, clipped to 0.4-1.0: 0 ideal, 1 constant ``eta_const``, 2 polynomial (Wang & Cleland 2018), 3 Brendel
    et al. (suction pressure P_in, MPa), 5 Pareto front (the default; 4, the first Pareto front, was removed)."""
    pr = max(pressure_ratio, 1.0)
    m = int(model)
    if m == 0:
        e = 1.0
    elif m == 1:
        e = min(max(float(eta_const), 0.4), 1.0)
    elif m == 2:
        e = (0.65 + 0.015 * pr - 0.0015 * pr * pr) / 0.96
    elif m == 3:
        P = 0.2 if P_in is None else P_in
        e = 0.66981 - 0.6 / (pr - 0.01466) ** (8.38 * P) - 0.00102 * pr ** 1.8 if pr > 1 else 0.4
    elif m == 5:
        e = pareto_front_efficiency(pr, "isentropic")
    elif m == 4:
        raise ValueError("compressor model 4 (first Pareto front) was removed; use 5 (Pareto front)")
    else:
        raise ValueError(f"unknown compressor model {model}")
    return min(max(e, 0.4), 1.0)


def eta_is(pr: float, model: int = 5, eta_const: float = 0.7 / 0.96, P_in: float | None = None) -> float:
    """``eta_model`` with the pressure ratio first (model 5, the Pareto front, by default)."""
    return eta_model(model, eta_const, pr, P_in)


def stage(fl: Fluid, P_in: float, st: dict, P_out: float, eta: float | None = None,
          model: int = 5, eta_const: float = 0.7 / 0.96) -> dict:
    """One compression from state ``st`` at P_in to P_out at isentropic efficiency ``eta`` (or the compressor
    model's at this pressure ratio): {w, eta, out, pr}."""
    e = eta if eta is not None else eta_is(P_out / P_in, model, eta_const, P_in)
    h_is = fl.PS(P_out, st["s"])["h"]
    w = (h_is - st["h"]) / e
    return {"w": w, "eta": e, "out": fl.PH(P_out, st["h"] + w), "pr": P_out / P_in}


def dry_superheat(fl: Fluid, P_in: float, T_sat: float, P_out: float, sh_min: float, model: int = 5,
                  eta_const: float = 0.7 / 0.96, margin: float = 5.0) -> float:
    """Smallest suction superheat (K, at least sh_min) whose compression from saturation at P_in (T_sat) to P_out
    discharges ``margin`` K above the dew point at P_out: overhanging fluids compress wet from near saturation."""
    supercritical = math.isfinite(fl.P_crit) and P_out >= fl.P_crit
    dew_out = -math.inf if supercritical else fl.PQ(P_out, 1.0)["T"]

    def dry(sh):
        return stage(fl, P_in, fl.PT(P_in, T_sat + sh), P_out, None, model, eta_const)["out"]["T"] - dew_out - margin

    sh = max(sh_min, 0.5)
    if dry(sh) >= 0:
        return sh
    hi = sh
    while dry(hi) < 0 and hi < 120:
        hi += 10
    if dry(hi) < 0:
        raise ValueError(f"{fl.fluid} cannot compress dry from {T_sat:.0f} C (raise the superheat or change the fluid)")
    lo = sh
    for _ in range(30):
        if hi - lo <= 0.05:
            break
        m = (lo + hi) / 2
        if dry(m) >= 0:
            hi = m
        else:
            lo = m
    return hi


def expand_to(fl: Fluid, P_in: float, h_in: float, P_out: float, eta: float | None) -> float:
    """Outlet enthalpy from the high side: an isenthalpic valve (eta None or 0) or an expander of isentropic
    efficiency ``eta``."""
    if not eta or not eta > 0:
        return h_in
    return h_in - eta * (h_in - fl.PS(P_out, fl.PH(P_in, h_in)["s"])["h"])


def expander_electric(w: float, wx: float, motor: float) -> float:
    """Electrical input: compressor shaft work over the motor less the expander work returned through it."""
    return w / motor - wx * motor


def expander_warning(fl: Fluid, P_out: float, h_out: float, warnings: list, name: str = "expander"):
    """Quality at a two-phase expander outlet (and a warning), else None."""
    try:
        o = fl.PH(P_out, h_out)
    except Exception:
        return None
    q = o.get("q", -1.0)
    wet = o.get("phase") == "two_phase" or (q is not None and 0 < q < 1)
    if wet:
        warnings.append(f"two-phase expansion: the {name} outlet is {q * 100:.0f} % vapour ({o['T']:.1f} C); "
                        "two-phase expanders are less proven than gas turbines and their efficiency is uncertain")
    return q if wet else None


# ----------------------------------------------------------------------------------------
# Exchangers
# ----------------------------------------------------------------------------------------


def air_profile(fl: Fluid, air: AirModel, P: float, h_hot: float, h_cold: float, Ta_in: float, Ta_out: float,
                n: int = 16, fast: bool = False) -> list:
    """A stream cooled along an isobar from h_hot to h_cold against air heated from Ta_in to Ta_out (counterflow):
    [[q, T_stream, T_air]] per kg of the stream from the cold end."""
    ha0, ha1 = air.h(Ta_in), air.h(Ta_out)
    Q = h_hot - h_cold
    return [[f * Q, T_ph(fl, P, h_cold + f * Q, fast), air.T(ha0 + f * (ha1 - ha0))] for f in lin(0.0, 1.0, n)]


def min_gap(profile: list) -> float:
    return min(p[1] - p[2] for p in profile)


def cooler_outlet(fl: Fluid, air: AirModel, P: float, h_in: float, Ta_in: float, Ta_out: float, dT: float,
                  fast: bool = False) -> float:
    """Outlet enthalpy of a gas cooler or condenser heating air from Ta_in to Ta_out: the coldest outlet that keeps
    the approach dT (the air flow follows from the duty)."""
    T_in = T_ph(fl, P, h_in, fast)
    if T_in < Ta_out + dT - 1e-9:
        raise ValueError(f"stream at {T_in:.1f} C cannot heat air to {Ta_out:.1f} C with a {dT} K approach")
    h_low = fl.PT(P, Ta_in + dT)["h"]

    def g(h):
        return min_gap(air_profile(fl, air, P, h_in, h, Ta_in, Ta_out, 14, fast)) - dT

    if g(h_low) >= -1e-6:
        return h_low
    lo, hi = h_low, h_in - 1e-6 * abs(h_in)
    if g(hi) < -0.02:
        raise ValueError("no outlet keeps the air approach")
    if g(hi) < 0:
        return hi
    # near the pseudo-critical point the pinch moves inside and the feasible outlets need not be one interval:
    # scan up from the coldest outlet for the first feasible one, then refine below it
    N = 40
    for i in range(1, N + 1):
        h = h_low + (hi - h_low) * i / N
        if g(h) >= 0:
            lo = h_low + (hi - h_low) * (i - 1) / N
            hi = h
            break
    k = 0
    while k < 40 and hi - lo > 0.01:
        m = (lo + hi) / 2
        if g(m) >= 0:
            hi = m
        else:
            lo = m
        k += 1
    return hi


def max_duty(fh: Fluid, P_h: float, h_h_in: float, m_h: float, fc: Fluid, P_c: float, h_c_in: float, dT: float,
             n: int = 14, fast: bool = False) -> float:
    """Counterflow exchanger, hot flow m_h (kg per kg of the cold flow): the largest duty per kg of the cold flow
    that keeps the approach dT everywhere."""
    Th_in, Tc_in = T_ph(fh, P_h, h_h_in, fast), T_ph(fc, P_c, h_c_in, fast)
    if Th_in - Tc_in < dT:
        return 0.0
    fs = lin(0.0, 1.0, n)

    def gap(q):
        m = math.inf
        for f in fs:
            Tc = T_ph(fc, P_c, h_c_in + f * q, fast)
            Th = T_ph(fh, P_h, h_h_in - q / m_h + f * q / m_h, fast)
            m = min(m, Th - Tc)
        return m - dT

    qmax = min(m_h * (h_h_in - fh.PT(P_h, Tc_in + dT)["h"]), fc.PT(P_c, Th_in - dT)["h"] - h_c_in)
    if not qmax > 0:
        return 0.0
    if gap(qmax) >= 0:
        return qmax
    # the pinch can sit inside (cp peaks near the pseudo-critical point), so the feasible duties need not be one
    # interval: scan down from qmax for the largest feasible duty, then refine above it
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
    return lo


def source_T(T_in: float, T_out: float | None, f: float, curve: Callable[[float], float] | None = None) -> float:
    """Heat-source temperature at a fraction f of the source duty from its cold end: a water line from T_out to
    T_in, or ``curve(f)`` (e.g. an air source)."""
    if curve is not None:
        return curve(f)
    To = T_in if T_out is None else T_out
    return To + f * (T_in - To)


# ----------------------------------------------------------------------------------------
# Results of the multi-loop cycles
# ----------------------------------------------------------------------------------------


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


# ----------------------------------------------------------------------------------------
# Compass search (the page's optimisers)
# ----------------------------------------------------------------------------------------


def compass_search(f, starts, lo, hi, step0, tol, is_int, state, max_eval, score=None, empty=None):
    """The page's compassSearch. f(x) returns an outcome scored by ``score`` (default ``o[0]``; -inf when
    infeasible); integer dimensions step by whole numbers and stop at a step of 1. Each coordinate is stepped both
    ways, a successful step is repeated while it pays, and the steps halve when nothing moves. Runs from each start
    (x, outcome) in turn and keeps the best. Returns (x, best outcome, converged); ``empty`` is the outcome when no
    start is feasible (default (-inf, None))."""
    if score is None:
        def score(o):
            return o[0]
    dim = len(lo)

    def clip(x):
        out = []
        for i, v in enumerate(x):
            c = min(max(v, lo[i]), hi[i])
            out.append(js_round(c) if is_int[i] else c)
        return out

    X, B, conv = None, (empty if empty is not None else (-math.inf, None)), False
    for xs0, bs0 in starts:
        if xs0 is None or not math.isfinite(score(bs0)):
            continue
        x, best, step = list(xs0), bs0, list(step0)

        def live():
            return any((s >= 1) if is_int[i] else (s > tol[i]) for i, s in enumerate(step))

        while state["nEval"] < max_eval and live():
            moved = False
            for i in range(dim):
                if state["nEval"] >= max_eval:
                    break
                if (step[i] < 1) if is_int[i] else (step[i] <= tol[i]):
                    continue
                for d in (1, -1):
                    z = list(x)
                    z[i] += d * step[i]
                    zc = clip(z)
                    if zc[i] == x[i]:
                        continue
                    o = f(zc)
                    if score(o) > score(best) + 1e-7:
                        best, x, moved = o, zc, True
                        for _ in range(20):                       # keep going while it pays
                            if state["nEval"] >= max_eval:
                                break
                            w = list(x)
                            w[i] += d * step[i]
                            wc = clip(w)
                            if wc[i] == x[i]:
                                break
                            ow = f(wc)
                            if score(ow) > score(best) + 1e-7:
                                best, x = ow, wc
                            else:
                                break
                        break
            if not moved:
                step = [(math.floor(s / 2) if s >= 1 else 0) if is_int[i] else (s / 2 if s > tol[i] else s)
                        for i, s in enumerate(step)]
        if score(best) > score(B):
            X, B, conv = x, best, not live()
    return X, B, conv


def memo_objective(state: dict, solve):
    """The page's memoF: f(x) -> (COP, result), -inf for a failed or flagged ('bad') design; memoised on x."""
    def f(x):
        k = ",".join(f"{float(v):.4f}" for v in x)
        if k in state["memo"]:
            return state["memo"][k]
        state["nEval"] += 1
        o = (-math.inf, None)
        try:
            r = solve(x)
            if r is not None and r.success and r.COP is not None and math.isfinite(r.COP) and not getattr(r, "bad", False):
                o = (r.COP, r)
            else:
                o = (-math.inf, r)
        except SOLVE_ERRORS:
            pass
        state["memo"][k] = o
        return o
    return f
