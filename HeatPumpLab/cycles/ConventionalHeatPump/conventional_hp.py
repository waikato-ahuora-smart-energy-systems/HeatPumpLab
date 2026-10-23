"""Conventional single-stage vapour-compression heat pump: evaporator, internal heat exchanger
(IHX), compressor, condenser, subcooler and expansion valve.

Flow path (state numbers used in ``states``)::

    1e evap_out    evaporator outlet, saturated vapour (dew point) at P_L
    1 suction      IHX vapour outlet = compressor suction, dT_superheat above the dew point at P_L
    2 discharge    compressor outlet at P_H, after the compressor heat loss (see below)
    2d cond_dew    saturated vapour at P_H (end of desuperheating, inside the condenser)
    3 cond_out     condenser outlet at P_H: the bubble point, or subcooled (preheat_subcool_share)
    4 subcool_out  subcooler outlet, dT_subcool below the bubble point at P_H
    4i ihx_out     IHX liquid outlet = valve inlet; the liquid gives the suction gas its superheat
    5 evap_in      expansion valve outlet (isenthalpic) at P_L

The IHX always sits after the subcooler (3 -> 4 -> 4i -> valve) and gives all of the suction
superheat: it is in use whenever dT_superheat > 0 (with 0 the evaporator outlet goes straight to the
compressor). It is checked for a minimum temperature approach ``dTmin_ihx`` (counterflow).

Inputs are the heat-source temperature, the steam temperature, the hot-water inlet and outlet
temperatures, the suction superheat and a compressor curve (see ``compressor_curves`` and
``compressor_data``). The heat source is taken as isothermal: the evaporator outlet (dew point)
sits ``dTmin_source`` below it, which fixes the evaporating pressure. The condenser (with its desuperheater) raises steam: the refrigerant bubble
point at the condenser outlet sits ``dTmin_steam`` above the steam temperature, which fixes the
condensing pressure. At or above the critical temperature there is no bubble point: the gas-cooler exit sits at
T_steam + dTmin_steam on the critical isochore (rho = rho_crit), which sets the high pressure, and the compressor map
reads the critical-isochore state at the discharge pressure (``critical_isochore``). The subcooler heats water from ``T_water_in`` to ``T_water_out`` counter-
currently: the refrigerant is subcooled as far as the minimum approach ``dTmin_water`` allows (for a
liquid whose heat capacity rises with temperature the pinch is at the cold end, so the liquid leaves
at T_water_in + dTmin_water), and the water flow follows from the duty (constant water heat
capacity). ``T_evap`` / ``T_cond`` / ``dT_subcool`` set the evaporating or condensing temperature or
the subcooling directly instead, bypassing the source / steam / hot-water specification. For a zeotropic
mixture the saturation temperatures are taken at the dew point by default (the EN 12900
convention of manufacturer data); ``evap_basis`` / ``cond_basis`` can be "bubble" or "mean".
Compressor maps are always evaluated at the dew-point temperatures. Pressure drops are neglected:
each side of the cycle is at one pressure.

Compressor heat loss (``loss_mode``):
    "UA"                  Q_loss = UA (T_wall - T_ambient), with T_wall the mean of the suction and
                          discharge temperatures (``loss_dT_basis="mean"``, default) or the
                          discharge temperature ("discharge").
                          The discharge enthalpy solves h2 = h1 + w_el - Q_loss / m, with a fixed
                          UA (``loss_UA_kWK``, default 0.015 kW/K = 15 W/K). Needs a sized
                          compressor (mass flow).
    "discharge_superheat" the discharge is fixed at dT_discharge_superheat above the dew point and
                          the loss closes the balance.
    "curve"               the compressor curve's own discharge enthalpy (heat_loss_fraction);
                          used when loss_UA_kWK is None.
The loss is the part of the electrical input that does not end up in the gas.

Exchanger heat losses: Q = UA (T_mean - T_ambient) for the condenser (steam generator), subcooler and
evaporator, with T_mean the mean of the refrigerant inlet and outlet temperatures of that exchanger
(``loss_UA_cond_kWK``, ``loss_UA_subcool_kWK``, ``loss_UA_evap_kWK``, default 0). The refrigerant states
are unchanged: the condenser and subcooler losses come off the heat delivered to the steam and the hot
water, and the heat source supplies the evaporator duty plus the evaporator loss (a negative loss, below
ambient, is a gain). They need a mass flow (a sized compressor, or Q_heat_kW).

The model is direct (no iteration on the cycle apart from the 1-D heat-loss balance): the pressures
follow from the saturation temperatures, the compressor curve and heat loss give the discharge
state, and the subcooler outlet fixes the
valve. Results are per kg of refrigerant, plus absolute flows and duties when the compressor curve
has a displacement (or when ``Q_heat_kW`` scales the cycle).

Units follow HeatPumpLab: MPa, degC, kJ/kg, kg/m3, kg/s, kW.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field, replace
from typing import Iterable

from .compressor_curves import CompressorInlet
from .compressor_data import recip_r600a
from .compressor_window import limit_flags
from .water_steam import steam_from_feedwater

GAS_PHASES = frozenset({"saturated_vapour", "superheated_vapour", "supercritical"})
LIQUID_PHASES = frozenset({"liquid", "saturated_liquid"})
DENSE_PHASES = LIQUID_PHASES | {"supercritical"}     # liquid-like states on a supercritical high side


def critical_isochore(pr, P: float) -> dict:
    """State on the critical isochore (rho = rho_crit) at pressure P >= P_crit: the continuation of the high-side
    saturation line above the critical point (supercritical gas cooler)."""
    from scipy.optimize import brentq
    rc = getattr(pr, "rho_crit", float("nan"))
    if not rc > 0:
        raise ValueError("critical density not available for this refrigerant")

    def g(T):
        return pr.ref_TD(T, rc)["P"] - P
    lo, hi = pr.T_crit - 0.5, pr.T_crit + 20.0
    for _ in range(12):
        if g(hi) >= 0:
            break
        hi += 40.0
    T = brentq(g, lo, hi, xtol=1e-9)
    return dict(pr.ref_TD(T, rc)) | dict(P=P, phase="supercritical", q=-1.0)


# ----------------------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------------------


@dataclass
class ConventionalHPConfig:
    """Operating point. Temperatures in degC, temperature differences in K, pressures in MPa.

    Defaults: n-butane, 55 C heat source (evaporating dew point 50 C), 120 C steam (condenser outlet at 125 C), hot water 80 to
    95 C (refrigerant leaves the subcooler at 85 C), 15 K suction superheat, Recip-M R600a
    map (56.25 m3/h), compressor heat loss
    UA ((T_suction + T_discharge)/2 - T_ambient) with a fixed UA of 15 W/K.
    """

    refrigerant: str = "Butane"
    T_source: float = 55.0           # heat source (isothermal)
    dTmin_source: float = 5.0        # evaporator outlet (dew point) below the source, K
    T_steam: float = 120.0           # steam raised by the desuperheater and condenser (saturated)
    dTmin_steam: float = 5.0         # condenser bubble point above the steam temperature, K
    # Share (0-1) of the feedwater preheat (bled hot water heated to boiling in the steam generator) that
    # comes from subcooling the refrigerant liquid at the condenser exit; the rest is from condensing.
    # The subcooling heat goes to the steam instead of the hot water, so more steam.
    preheat_subcool_share: float = 0.0
    T_water_in: float = 80.0         # hot water into the subcooler
    T_water_out: float = 95.0        # hot water out of the subcooler
    dTmin_water: float = 5.0         # minimum subcooler approach to the water, K
    # Direct alternatives: set any to bypass the source / steam / hot-water specification.
    T_evap: float | None = None      # evaporating temperature (on evap_basis)
    T_cond: float | None = None      # condensing temperature (on cond_basis)
    dT_subcool: float | None = None  # subcooler outlet below the bubble point at P_H
    dT_superheat: float = 15.0       # compressor suction above the dew point at P_L
    dTmin_ihx: float = 5.0           # minimum IHX approach, K
    # Compressor heat loss, Q = UA (T_wall - T_ambient); see the module docstring.
    loss_UA_kWK: float | None = 0.015    # kW/K (15 W/K); None -> the curve's own discharge enthalpy
    T_ambient: float = 15.0
    loss_dT_basis: str = "mean"          # "mean" (of suction and discharge) | "discharge"
    # Exchanger heat losses, Q = UA (T_mean - T_ambient); see the module docstring.
    loss_UA_cond_kWK: float = 0.0
    loss_UA_subcool_kWK: float = 0.0
    loss_UA_evap_kWK: float = 0.0
    # Alternative: fix the discharge above the condensing dew point (takes precedence over UA).
    dT_discharge_superheat: float | None = None
    compressor: object = field(default_factory=lambda: recip_r600a("Recip-M"))
    evap_basis: str = "dew"          # "dew" | "bubble" | "mean"
    cond_basis: str = "dew"
    subcooler_heat_useful: bool = True   # count the subcooler duty as delivered heat
    Q_heat_kW: float | None = None   # scale a size-less cycle to this heat output
    dP_crit_min: float = 0.0         # MPa, minimum clearance of P_H below the critical pressure
    # Refrigerant pressure-loss coefficients: dP [kPa] = K (V [m3/h])^1.8, V at the mean specific volume of each
    # exchanger's inlet and outlet. Exits are pinned (condenser bubble point, evaporator dew point), so the
    # compressor discharges at P_H + dP_cond and the valve feeds the evaporator at P_L + dP_evap.
    dP_K_cond: float = 0.0
    dP_K_subcool: float = 0.0
    dP_K_evap: float = 0.0
    property_backend: str | None = None  # "refprop" | "coolprop"; None -> $HEATPUMPLAB_PROPERTY_BACKEND or CoolProp
    compressor_Hz: float | None = None   # compressor supply frequency; None keeps the curve's own (50 Hz)


# ----------------------------------------------------------------------------------------
# Property access (injectable, so the cycle logic can be tested without REFPROP)
# ----------------------------------------------------------------------------------------


class FluidProperties:
    """Thin wrapper over a HeatPumpLab property package for one refrigerant.

    backend: "refprop" (PropertyPackage) or "coolprop" (CoolPropPackage); None follows
    $HEATPUMPLAB_PROPERTY_BACKEND and defaults to REFPROP.
    """

    def __init__(self, refrigerant: str, backend: str | None = None):
        from ...utils.coolprop_properties import make_property_package, property_backend
        self.backend = property_backend(backend)
        self.ref = make_property_package(refrigerant, backend=self.backend)
        self.P_crit = float(self.ref.P_crit)
        self.T_crit = float(self.ref.T_crit)
        # critical density, kg/m3 on both backends (PropertyPackage converts REFPROP's mol/dm3)
        self.rho_crit = float(getattr(self.ref, "D_crit", float("nan")))

    def _r(self, pair, a, b):
        return self.ref.get_std_properties(pair, a, b, self.ref.fld_composition)

    def ref_PT(self, P, T):
        return self._r("PT", P, T)

    def ref_PH(self, P, h):
        return self._r("PH", P, h)

    def ref_PS(self, P, s):
        return self._r("PS", P, s)

    def ref_PQ(self, P, q):
        return self._r("PQ", P, q)

    def ref_TQ(self, T, q):
        return self._r("TQ", T, q)

    def ref_TD(self, T, D):
        return self._r("TD", T, D)


class RefpropProperties(FluidProperties):
    """REFPROP-backed properties (HeatPumpLab PropertyPackage)."""

    def __init__(self, refrigerant: str):
        super().__init__(refrigerant, backend="refprop")


class CoolPropProperties(FluidProperties):
    """CoolProp-backed properties (HeatPumpLab CoolPropPackage)."""

    def __init__(self, refrigerant: str):
        super().__init__(refrigerant, backend="coolprop")


# ----------------------------------------------------------------------------------------
# Result
# ----------------------------------------------------------------------------------------


@dataclass
class ConventionalHPResult:
    success: bool
    message: str
    COP: float | None = None          # heating COP on electrical input
    COP_cooling: float | None = None
    VHC: float | None = None          # MJ/m3, delivered heat per suction volume
    P_H: float | None = None          # MPa
    P_L: float | None = None
    report: dict = field(default_factory=dict)
    states: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def to_dict(self):
        return asdict(self)


_STATE_KEYS = ("T", "P", "h", "s", "D", "q", "phase")


def _state(d: dict) -> dict:
    return {k: d.get(k) for k in _STATE_KEYS}


# ----------------------------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------------------------


class ConventionalHeatPump:
    """Rating model for a single-stage heat pump; see the module docstring."""

    def __init__(self, config: ConventionalHPConfig, props=None):
        self.c = config
        self.props = props if props is not None else FluidProperties(config.refrigerant, config.property_backend)

    # ---- compressor heat loss --------------------------------------------------------

    @property
    def loss_mode(self) -> str:
        c = self.c
        if c.dT_discharge_superheat is not None:
            return "discharge_superheat"
        if c.loss_UA_kWK is not None:
            return "UA"
        return "curve"

    def _wall_temperature(self, T_suction: float, T_discharge: float) -> float:
        if self.c.loss_dT_basis == "discharge":
            return T_discharge
        if self.c.loss_dT_basis == "mean":
            return 0.5 * (T_suction + T_discharge)
        raise ValueError("loss_dT_basis must be 'discharge' or 'mean'")

    def loss_UA(self) -> float:
        """UA (kW/K) of the compressor heat loss (``loss_UA_kWK``)."""
        if self.c.loss_UA_kWK is None:
            raise ValueError("no UA: set loss_UA_kWK")
        return self.c.loss_UA_kWK

    WATER_CP = 4.19     # kJ/kg/K, liquid water at 60-110 C (within 1 %)

    def _water_subcooler(self, P_H: float, sc_in: dict, n: int = 20):
        """Subcooler outlet for hot water heated from T_water_in to T_water_out counter-currently.

        The refrigerant is cooled as far as the minimum approach dTmin_water allows: the outlet goes
        to T_water_in + dTmin_water unless the approach along the exchanger (water temperature linear
        in duty) is tighter somewhere else, in which case the outlet is raised until it is not.
        Returns (outlet state, dict(min_approach=K)) or (None, failure message).
        """
        c, pr = self.c, self.props
        T_in, h_in = sc_in["T"], sc_in["h"]
        dT = c.dTmin_water
        if T_in - c.T_water_out < dT - 1e-9:
            return None, (f"hot water outlet {c.T_water_out:g} C is less than {dT:g} K below the refrigerant "
                          f"entering the subcooler ({T_in:.1f} C)")

        def approach(T_out):
            """Smallest refrigerant-water temperature difference along the exchanger: n samples,
            then a golden-section refinement around the smallest interior one."""
            h_out = pr.ref_PT(P_H, T_out)["h"]

            def gap(x):                          # x: duty fraction from the cold end
                return (pr.ref_PH(P_H, h_out + (h_in - h_out) * x)["T"]
                        - (c.T_water_in + (c.T_water_out - c.T_water_in) * x))
            ends = min(T_out - c.T_water_in, T_in - c.T_water_out)
            vals = [(gap(i / n), i) for i in range(1, n)]
            g, i = min(vals)
            if g >= ends:
                return ends
            a, b, r = (i - 1) / n, (i + 1) / n, (math.sqrt(5) - 1) / 2
            for _ in range(40):
                x1, x2 = b - r * (b - a), a + r * (b - a)
                if gap(x1) < gap(x2):
                    b = x2
                else:
                    a = x1
            return min(g, gap(0.5 * (a + b)), ends)

        T_out = c.T_water_in + dT
        if T_out >= T_in:                        # no room to subcool
            return dict(sc_in), dict(min_approach=T_in - c.T_water_in)
        if approach(T_out) < dT - 1e-6:
            from scipy.optimize import brentq
            T_out = brentq(lambda T: approach(T) - dT, T_out, T_in - 1e-6, xtol=1e-6)
        out = pr.ref_PT(P_H, T_out)
        if out.get("phase") not in DENSE_PHASES:
            return None, f"subcooler outlet not liquid ({out.get('phase')})"
        return out, dict(min_approach=approach(T_out))

    def _ua_discharge(self, P_H: float, suction: dict, w_el: float, m_dot: float) -> dict:
        """Discharge state with Q_loss = UA (T_wall - T_ambient): h2 = h1 + w_el - Q_loss / m."""
        from scipy.optimize import brentq
        UA, T_amb, h1, T1 = self.loss_UA(), self.c.T_ambient, suction["h"], suction["T"]

        def g(h2):
            T2 = self.props.ref_PH(P_H, h2)["T"]
            return (h2 - h1) - w_el + UA * (self._wall_temperature(T1, T2) - T_amb) / m_dot

        h_adiabatic = h1 + w_el
        if UA == 0.0:
            return self.props.ref_PH(P_H, h_adiabatic)
        g0 = g(h_adiabatic)
        if abs(g0) < 1e-12:
            return self.props.ref_PH(P_H, h_adiabatic)
        step = abs(g0) + 1.0
        a = b = h_adiabatic
        for _ in range(60):                       # g increases with h2: bracket the root
            if g0 > 0:
                a = b - step
                if g(a) <= 0:
                    break
                b = a
            else:
                b = a + step
                if g(b) >= 0:
                    break
                a = b
            step *= 2
        else:
            raise ValueError("compressor heat-loss balance not bracketed")
        return self.props.ref_PH(P_H, brentq(g, a, b, xtol=1e-9))

    # ---- saturation ------------------------------------------------------------------

    def saturation_pressure(self, T: float, basis: str = "dew") -> float:
        """Pressure (MPa) at which the dew, bubble or mean saturation temperature equals T."""
        pr = self.props
        if basis == "dew":
            return pr.ref_TQ(T, 1.0)["P"]
        if basis == "bubble":
            return pr.ref_TQ(T, 0.0)["P"]
        if basis != "mean":
            raise ValueError(f"saturation basis must be 'dew', 'bubble' or 'mean', not '{basis}'")
        P_dew, P_bub = pr.ref_TQ(T, 1.0)["P"], pr.ref_TQ(T, 0.0)["P"]
        if abs(P_bub / P_dew - 1.0) < 1e-9:      # pure fluid or azeotrope
            return P_dew

        def residual(lnP):
            P = math.exp(lnP)
            return 0.5 * (pr.ref_PQ(P, 1.0)["T"] + pr.ref_PQ(P, 0.0)["T"]) - T

        from scipy.optimize import brentq
        return math.exp(brentq(residual, math.log(P_dew), math.log(P_bub), xtol=1e-10))

    # ---- one operating point ---------------------------------------------------------

    def solve(self) -> ConventionalHPResult:
        """Rate the cycle at the configured operating point."""
        try:
            from ...utils.property_cache import active_property_call_cache, property_call_cache
        except ImportError:                       # pragma: no cover - standalone use
            return self._solve()
        if active_property_call_cache() is None:
            with property_call_cache():
                return self._solve()
        return self._solve()

    def _solve(self) -> ConventionalHPResult:
        # Sized by Q_heat_kW, the mass flow is only known after the cycle; the preheat share's condenser loss
        # per kg needs it, so re-solve with the last mass flow until it settles.
        self._m_dot_hint = None
        res = self._solve_dp()
        c = self.c
        if not (c.loss_UA_cond_kWK and c.preheat_subcool_share > 0):
            return res
        for _ in range(20):
            if not res.success or res.report.get("sized_by") != "Q_heat_kW":
                break
            m = res.report["m_dot_kgs"]
            if self._m_dot_hint is not None and abs(m - self._m_dot_hint) <= 1e-10 * m:
                break
            self._m_dot_hint = m
            res = self._solve_dp()
        return res

    def _solve_dp(self) -> ConventionalHPResult:
        c = self.c
        K = (c.dP_K_cond or 0.0, c.dP_K_subcool or 0.0, c.dP_K_evap or 0.0)
        dP = (0.0, 0.0, 0.0)
        res = self._solve_once(dP)
        if not any(K) or not res.success or res.report.get("Vmean_cond_m3h") is None:
            return res
        for _ in range(15):                     # the drops depend on the flow: re-solve until they settle
            r = res.report
            nxt = tuple(k * v ** 1.8 / 1000.0 for k, v in
                        zip(K, (r["Vmean_cond_m3h"], r["Vmean_subcool_m3h"], r["Vmean_evap_m3h"])))   # MPa
            err = sum(abs(a - b) for a, b in zip(nxt, dP))
            dP = nxt
            res = self._solve_once(dP)
            if not res.success or err < 1e-7:
                break
        return res

    def _solve_once(self, dP) -> ConventionalHPResult:
        c, pr = self.c, self.props
        dP_cond, dP_sub, dP_evap = dP
        warnings: list[str] = []
        try:
            if (c.dT_subcool is not None and c.dT_subcool < 0) or c.dT_superheat < 0:
                return ConventionalHPResult(False, "subcooling and superheat must be >= 0")
            if c.dT_subcool is None and not c.T_water_out > c.T_water_in:
                return ConventionalHPResult(False, "hot water outlet must be above its inlet")
            T_crit = getattr(pr, "T_crit", None)
            T_hot_set = c.T_cond if c.T_cond is not None else c.T_steam + c.dTmin_steam
            # at or above the critical temperature the high side is supercritical: the gas-cooler exit (where the
            # bubble point would be) at T_hot_set on the critical isochore sets the pressure
            sc = T_crit is not None and math.isfinite(T_crit) and T_hot_set >= T_crit
            if sc and not getattr(pr, "rho_crit", float("nan")) > 0:
                return ConventionalHPResult(False, f"condensing at {T_hot_set:.1f} C is at or above the critical "
                                                   f"temperature of this refrigerant ({T_crit:.1f} C), and its "
                                                   f"critical density is not available")
            liquid_ok = DENSE_PHASES if sc else LIQUID_PHASES
            if c.T_evap is not None:
                P_L = self.saturation_pressure(c.T_evap, c.evap_basis)
            else:                                   # evaporator outlet (dew) dTmin_source below the source
                P_L = self.saturation_pressure(c.T_source - c.dTmin_source, "dew")
            if sc:
                gc_out = pr.ref_TD(T_hot_set, pr.rho_crit)
                P_H = gc_out["P"]
            elif c.T_cond is not None:
                P_H = self.saturation_pressure(c.T_cond, c.cond_basis)
            else:                                   # condenser outlet (bubble) dTmin_steam above the steam
                P_H = self.saturation_pressure(c.T_steam + c.dTmin_steam, "bubble")
            if not P_H > P_L:
                return ConventionalHPResult(False, "condensing pressure not above evaporating pressure",
                                            P_H=P_H, P_L=P_L)
            P_crit = getattr(pr, "P_crit", None)
            if not sc and P_crit is not None and P_H > P_crit - c.dP_crit_min:
                return ConventionalHPResult(
                    False, f"condensing pressure {P_H * 10:.2f} bar within {c.dP_crit_min * 10:g} bar "
                           f"of the critical pressure {P_crit * 10:.2f} bar", P_H=P_H, P_L=P_L)

            # pressure drops (exits pinned): compressor discharge and evaporator inlet pressures
            P_dis, P_evin, P_valve = P_H + dP_cond, P_L + dP_evap, P_H - dP_sub
            if not sc and dP_cond and P_crit is not None and P_dis > P_crit - c.dP_crit_min:
                return ConventionalHPResult(
                    False, f"discharge pressure {P_dis * 10:.2f} bar (with the condenser pressure drop) is above "
                           f"the critical pressure {P_crit * 10:.2f} bar", P_H=P_H, P_L=P_L)
            if not P_valve > P_evin:
                return ConventionalHPResult(False, "pressure drops leave no pressure difference across the valve",
                                            P_H=P_H, P_L=P_L)
            evap_dew, evap_bub = pr.ref_PQ(P_L, 1.0), pr.ref_PQ(P_L, 0.0)
            # cond_dew: the dew point at the condenser inlet (compressor discharge pressure); supercritical: the
            # critical-isochore states at the discharge pressure and at the gas-cooler exit (T_hot_set)
            if sc:
                cond_dew = critical_isochore(pr, P_dis)
                cond_bub = dict(gc_out) | dict(P=P_H, phase="supercritical", q=-1.0)
            else:
                cond_dew, cond_bub = pr.ref_PQ(P_dis, 1.0), pr.ref_PQ(P_H, 0.0)
            cond_out = cond_bub                     # 3 condenser exit (subcooled below, after the compressor)
            share = min(max(c.preheat_subcool_share, 0.0), 1.0)
            # steam from bled hot water: needs the hot-water subcooler (dT_subcool None), as the steam bookkeeping does
            steam_mode = (c.T_cond is None and c.dT_subcool is None and c.T_water_out is not None
                          and c.T_water_out < c.T_steam)

            # 1 suction (IHX vapour outlet), 1e evaporator outlet
            suction = pr.ref_PT(P_L, evap_dew["T"] + c.dT_superheat) if c.dT_superheat > 0 else evap_dew
            if suction.get("phase") not in GAS_PHASES:
                return ConventionalHPResult(False, f"suction not gas ({suction.get('phase')})", P_H=P_H, P_L=P_L)
            evap_out = evap_dew                     # saturated: the IHX gives all the suction superheat
            q_ihx = suction["h"] - evap_out["h"]
            ihx = c.dT_superheat > 0

            # 2 compressor
            isentropic = pr.ref_PS(P_dis, suction["s"])
            h_is = isentropic["h"]

            def rated_suction(T):
                st = dict(pr.ref_PT(P_L, T))
                st["h_isentropic"] = pr.ref_PS(P_dis, st["s"])["h"]
                return st

            inlet = CompressorInlet(suction=suction, P_discharge=P_dis, h_isentropic=h_is,
                                    T_evap_dew=evap_dew["T"], T_cond_dew=cond_dew["T"],
                                    rated_suction=rated_suction,
                                    rho_evap_dew=evap_dew.get("D"), rho_cond_dew=cond_dew.get("D"))
            compressor = c.compressor
            if c.compressor_Hz is not None:
                if not hasattr(compressor, "frequency_Hz"):
                    return ConventionalHPResult(False, "this compressor curve has no speed setting (compressor_Hz)",
                                                P_H=P_H, P_L=P_L)
                compressor = replace(compressor, frequency_Hz=c.compressor_Hz)
            comp = compressor.evaluate(inlet)
            warnings += comp.warnings
            mode = self.loss_mode
            UA = None
            if mode == "UA":
                if comp.m_dot is None:
                    return ConventionalHPResult(False, "UA heat loss needs a sized compressor (displacement)",
                                                P_H=P_H, P_L=P_L, warnings=warnings)
                UA = self.loss_UA()
                discharge = self._ua_discharge(P_dis, suction, comp.w_electric, comp.m_dot)
                q_loss = 0.0 if UA == 0.0 else comp.w_electric - (discharge["h"] - suction["h"])
            elif mode == "discharge_superheat":
                discharge = pr.ref_PT(P_dis, cond_dew["T"] + c.dT_discharge_superheat)
                q_loss = comp.w_electric - (discharge["h"] - suction["h"])
                if q_loss < 0:
                    return ConventionalHPResult(
                        False, f"discharge superheat {c.dT_discharge_superheat:g} K needs "
                               f"{-q_loss:.1f} kJ/kg more than the compressor's electrical input",
                        P_H=P_H, P_L=P_L, warnings=warnings)
            else:
                discharge = pr.ref_PH(P_dis, comp.h_discharge)
                q_loss = comp.w_electric - (comp.h_discharge - suction["h"])
            if discharge.get("phase") not in GAS_PHASES:
                return ConventionalHPResult(False, f"wet compression (discharge {discharge.get('phase')})",
                                            P_H=P_H, P_L=P_L, warnings=warnings)

            # 3 condenser (steam generator) exit. With preheat_subcool_share > 0 the liquid is subcooled so that it
            # supplies that share of the feedwater preheat: per kg of refrigerant, q_x = share fp (q_c0 + q_x - l_c),
            # with fp = preheat / (steam - feed) enthalpy rise of the water and l_c the condenser loss per kg
            if share > 0 and steam_mode:
                fp = steam_from_feedwater(1.0, c.T_steam, c.T_water_out)
                fp = fp["preheat_kJkg"] / fp["dh_kJkg"]
                q_c0, q_x = discharge["h"] - cond_bub["h"], 0.0
                # mass flow from the compressor, or (sized by Q_heat_kW) from the previous pass of _solve
                m_c = comp.m_dot or getattr(self, "_m_dot_hint", None)
                for _ in range(4):                 # the condenser loss depends a little on the exit temperature
                    l_c = (c.loss_UA_cond_kWK * ((discharge["T"] + cond_out["T"]) / 2 - c.T_ambient) / m_c
                           if m_c else 0.0)
                    q_x = share * fp * (q_c0 - l_c) / (1 - share * fp)
                    cond_out = pr.ref_PH(P_H, cond_bub["h"] - q_x)
                if cond_out.get("phase") not in liquid_ok:
                    return ConventionalHPResult(False, f"condenser exit not liquid ({cond_out.get('phase')})",
                                                P_H=P_H, P_L=P_L, warnings=warnings)
            elif share > 0:
                warnings.append("preheat from subcooling needs the steam and hot-water specification; ignored")
            # 4 subcooler outlet
            h3 = cond_out["h"]
            water = None
            if c.dT_subcool is not None:
                sc_out = pr.ref_PT(P_H, cond_out["T"] - c.dT_subcool) if c.dT_subcool > 0 else cond_out
                if c.dT_subcool > 0 and sc_out.get("phase") not in liquid_ok:
                    return ConventionalHPResult(False, f"subcooler outlet not liquid ({sc_out.get('phase')})",
                                                P_H=P_H, P_L=P_L, warnings=warnings)
            else:
                sc_out, water = self._water_subcooler(P_H, cond_out)
                if sc_out is None:
                    return ConventionalHPResult(False, water, P_H=P_H, P_L=P_L, warnings=warnings)

            # 4i IHX liquid side (after the subcooler; its outlet is the valve inlet)
            h1, h2, h4 = suction["h"], discharge["h"], sc_out["h"]
            ihx_min_approach = None
            ihx_in, h_ihx_out = sc_out, h4 - q_ihx
            h_valve = h_ihx_out
            ihx_out = pr.ref_PH(P_H, h_ihx_out) if q_ihx > 0 else ihx_in
            if ihx:
                n = 6
                ihx_min_approach = min(
                    pr.ref_PH(P_H, h_ihx_out + q_ihx * i / n)["T"] - pr.ref_PH(P_L, evap_out["h"] + q_ihx * i / n)["T"]
                    for i in range(n + 1))
                if ihx_min_approach < c.dTmin_ihx - 1e-6:
                    return ConventionalHPResult(
                        False, f"IHX approach {ihx_min_approach:.1f} K below {c.dTmin_ihx:g} K: the liquid "
                               f"({ihx_in['T']:.1f} C) is too cold to superheat the suction gas to {suction['T']:.1f} C",
                        P_H=P_H, P_L=P_L, warnings=warnings)

            # 5 valve outlet
            evap_in = pr.ref_PH(P_evin, h_valve)

            q_desuperheat = h2 - cond_dew["h"]
            q_cond = h2 - h3
            q_sc = h3 - h4
            q_evap = evap_out["h"] - h_valve
            q_heat = q_cond + (q_sc if c.subcooler_heat_useful else 0.0)
            w_el, w_sh = comp.w_electric, comp.w_shaft
            balance = q_evap + (h2 - h1) - (q_cond + q_sc)

            cop = q_heat / w_el
            cop_c = q_evap / w_el
            v1 = 1.0 / suction["D"]
            vhc = q_heat / v1 / 1000.0

            # exchanger heat losses to ambient (kW), from the mean refrigerant temperature in each exchanger
            L_cond = c.loss_UA_cond_kWK * ((discharge["T"] + cond_out["T"]) / 2 - c.T_ambient)
            L_sc = c.loss_UA_subcool_kWK * ((cond_out["T"] + sc_out["T"]) / 2 - c.T_ambient)
            L_evap = c.loss_UA_evap_kWK * ((evap_in["T"] + evap_out["T"]) / 2 - c.T_ambient)
            L_heat = L_cond + (L_sc if c.subcooler_heat_useful else 0.0)
            m_dot, scaled = comp.m_dot, False
            if m_dot is None and c.Q_heat_kW is not None:
                m_dot, scaled = (c.Q_heat_kW + L_heat) / q_heat, True
            if m_dot is None and (c.loss_UA_cond_kWK or c.loss_UA_subcool_kWK or c.loss_UA_evap_kWK):
                warnings.append("exchanger heat losses need a mass flow (a sized compressor or Q_heat_kW); ignored")
            if m_dot is not None:          # delivered heat, net of the exchanger losses
                cop = (m_dot * q_heat - L_heat) / (m_dot * w_el)
                vhc = (m_dot * q_heat - L_heat) / (m_dot * v1) / 1000.0

            T_hot = c.T_cond if c.T_cond is not None else cond_bub["T"]   # Carnot between the condensing ...
            T_low = c.T_evap if c.T_evap is not None else evap_dew["T"]  # ... and evaporating temperatures
            report = dict(
                T_steam=None if c.T_cond is not None else c.T_steam,
                T_water_in=None if water is None else c.T_water_in, T_water_out=None if water is None else c.T_water_out,
                subcooler_min_approach_K=None if water is None else water["min_approach"],
                COP=cop, COP_cooling=cop_c, SEC=1 / cop, VHC=vhc, rho_suction_kgm3=suction["D"],
                COP_carnot=(T_hot + 273.15) / (T_hot - T_low) if T_hot > T_low else None,
                P_evap_bar=P_L * 10, P_cond_bar=P_H * 10, PR=P_dis / P_L,
                supercritical=sc, rho_crit_kgm3=pr.rho_crit if sc else None,
                P_discharge_bar=P_dis * 10, P_evap_in_bar=P_evin * 10, P_valve_in_bar=P_valve * 10,
                dP_cond_kPa=dP_cond * 1000, dP_subcool_kPa=dP_sub * 1000, dP_evap_kPa=dP_evap * 1000,
                Vmean_cond_m3h=None if m_dot is None else m_dot * (1 / discharge["D"] + 1 / cond_out["D"]) / 2 * 3600,
                Vmean_subcool_m3h=None if m_dot is None else m_dot * (1 / cond_out["D"] + 1 / sc_out["D"]) / 2 * 3600,
                Vmean_evap_m3h=None if m_dot is None else m_dot * (1 / evap_in["D"] + 1 / evap_out["D"]) / 2 * 3600,
                T_source=None if c.T_evap is not None else c.T_source,
                T_evap=T_low, T_cond=T_hot, evap_basis=c.evap_basis if c.T_evap is not None else "dew",
                cond_basis=c.cond_basis if c.T_cond is not None else "bubble",
                T_evap_dew=evap_dew["T"], T_evap_bubble=evap_bub["T"], evap_glide_K=evap_dew["T"] - evap_bub["T"],
                T_cond_dew=cond_dew["T"], T_cond_bubble=cond_bub["T"], cond_glide_K=cond_dew["T"] - cond_bub["T"],
                T_cond_out=cond_out["T"], dT_subcool_cond=cond_bub["T"] - cond_out["T"],
                preheat_subcool_share=share if steam_mode else 0.0,
                T_evap_in=evap_in["T"], x_evap_in=evap_in.get("q"), T_evap_out=evap_out["T"],
                T_suction=suction["T"], T_discharge=discharge["T"], T_subcool_out=sc_out["T"],
                T_valve_in=pr.ref_PH(P_H, h_valve)["T"] if h_valve != h4 else sc_out["T"],
                dT_superheat=c.dT_superheat, dT_subcool=cond_out["T"] - sc_out["T"],
                ihx=ihx, q_ihx_kJkg=q_ihx,
                ihx_min_approach_K=ihx_min_approach,
                eta_is_shaft=comp.eta_is, eta_is_overall=(h_is - h1) / w_el, eta_vol=comp.eta_vol,
                map_to=comp.to_map, map_tc=comp.tc_map,
                compressor_Hz=getattr(compressor, "frequency_Hz", None), compressor_rpm=getattr(compressor, "rpm", None),
                property_backend=getattr(pr, "backend", None),
                q_heat_kJkg=q_heat, q_cond_kJkg=q_cond, q_desuperheat_kJkg=q_desuperheat,
                q_subcool_kJkg=q_sc, q_evap_kJkg=q_evap, w_electric_kJkg=w_el, w_shaft_kJkg=w_sh,
                w_to_gas_kJkg=h2 - h1, q_compressor_loss_kJkg=q_loss,
                compressor_loss_fraction=q_loss / w_el, dT_discharge_superheat=discharge["T"] - cond_dew["T"],
                loss_mode=mode, UA_loss_kWK=UA, T_ambient=c.T_ambient,
                loss_dT_K=self._wall_temperature(suction["T"], discharge["T"]) - c.T_ambient,
                energy_balance_residual_kJkg=balance, v_suction_m3kg=v1,
                in_envelope=not any("envelope" in w for w in warnings),
            )
            # beyond the window: where the maps were read, the density ratio and how far past a ratio limit it is
            win = getattr(comp, "window", None)
            report.update(
                window_edge=None if win is None or win["inside"] else win["edge_name"],
                window_edge_kind=None if win is None or win["inside"] else win["edge_kind"],
                to_read=None if win is None else win["read"][0], tc_read=None if win is None else win["read"][1],
                density_ratio=None if win is None else win["ratio"],
                density_ratio_beyond=0.0 if win is None or win["inside"] else win["ratio_beyond"])
            swept = getattr(compressor, "swept_volume_m3h", None)
            report["w_swept_kJm3"] = (comp.m_dot * w_el / (swept / 3600.0)
                                      if comp.m_dot is not None and swept else None)
            # mechanical limits (seals, materials, motor) on the real cycle: flagged, not enforced
            disp, rpm0 = getattr(compressor, "displacement_m3h", None), getattr(compressor, "rated_rpm", None)
            report["limit_flags"] = limit_flags(getattr(compressor, "limits", None), P_dis * 10, discharge["T"],
                                                report["w_swept_kJm3"], P_L * 10,
                                                disp / (rpm0 * 60.0) if disp and rpm0 else None, w_sh / w_el)
            if m_dot is not None:
                report.update(
                    sized_by=("Q_heat_kW" if scaled else "compressor"),
                    m_dot_kgs=m_dot, Q_heat_kW=m_dot * q_heat - L_heat, Q_cond_kW=m_dot * q_cond - L_cond,
                    Q_subcool_kW=m_dot * q_sc - L_sc, Q_evap_kW=m_dot * q_evap, Q_ihx_kW=m_dot * q_ihx,
                    Q_cond_refrigerant_kW=m_dot * q_cond, Q_subcool_refrigerant_kW=m_dot * q_sc,
                    Q_source_kW=m_dot * q_evap + L_evap,
                    Q_loss_cond_kW=L_cond, Q_loss_subcool_kW=L_sc, Q_loss_evap_kW=L_evap,
                    UA_loss_cond_kWK=c.loss_UA_cond_kWK, UA_loss_subcool_kWK=c.loss_UA_subcool_kWK,
                    UA_loss_evap_kWK=c.loss_UA_evap_kWK, COP=cop, SEC=1 / cop, VHC=vhc,
                    W_electric_kW=m_dot * w_el, W_shaft_kW=m_dot * w_sh, Q_compressor_loss_kW=m_dot * q_loss,
                    V_suction_m3h=m_dot * v1 * 3600.0, V_discharge_m3h=m_dot / discharge["D"] * 3600.0,
                    m_water_kgs=(None if water is None else
                                 (m_dot * q_sc - L_sc) / (self.WATER_CP * (c.T_water_out - c.T_water_in))),
                )
            # steam from hot water bled at the subcooler outlet: the steam generator (desuperheating and
            # condensing) heats it from T_water_out to saturated steam at T_steam (water_steam.py)
            report.update(m_steam_kgs=None, steam_feed_T=None, steam_bleed_fraction=None, m_hot_water_net_kgs=None)
            if m_dot is not None and water is not None and c.T_cond is None and c.T_water_out < c.T_steam:
                sf = steam_from_feedwater(report["Q_cond_kW"], c.T_steam, c.T_water_out)
                report.update(m_steam_kgs=sf["m_kgs"], steam_feed_T=c.T_water_out, steam_dh_kJkg=sf["dh_kJkg"],
                              steam_preheat_fraction=sf["preheat_kJkg"] / sf["dh_kJkg"],
                              P_steam_bar=sf["P_steam_MPa"] * 10, steam_bleed_fraction=sf["m_kgs"] / report["m_water_kgs"],
                              m_hot_water_net_kgs=report["m_water_kgs"] - sf["m_kgs"])
                # steam generator approach, counter-current: the feedwater is heated to boiling over the first
                # preheat fraction of the duty from the refrigerant's cold (exit) end, then boils at T_steam. With a
                # subcooled exit the liquid must not fall within dTmin_steam of the water anywhere.
                report["steam_gen_min_approach_K"] = cond_bub["T"] - c.T_steam
                if h3 < cond_bub["h"] - 1e-9:
                    fp = report["steam_preheat_fraction"]
                    mn = math.inf
                    for f in sorted([i / 40 for i in range(41)] + [fp]):
                        Tr = pr.ref_PH(P_H + dP_cond * f, h3 + (h2 - h3) * f)["T"]
                        Tw = c.T_water_out + (c.T_steam - c.T_water_out) * f / fp if f < fp else c.T_steam
                        mn = min(mn, Tr - Tw)
                    report["steam_gen_min_approach_K"] = mn
                    if mn < c.dTmin_steam - 1e-6:
                        return ConventionalHPResult(
                            False, f"preheat from subcooling leaves the steam generator approach at {mn:.1f} K, below "
                                   f"{c.dTmin_steam:g} K", P_H=P_H, P_L=P_L, warnings=warnings)
            if report.get("m_water_kgs") is not None and report["m_water_kgs"] <= 0.0:
                return ConventionalHPResult(
                    False, f"the subcooler delivers no hot water ({report['Q_subcool_kW']:.1f} kW net of its heat loss)",
                    P_H=P_H, P_L=P_L, warnings=warnings)
            if report["m_steam_kgs"] is not None and report["m_steam_kgs"] > report["m_water_kgs"]:
                return ConventionalHPResult(
                    False, f"steam {report['m_steam_kgs']:.3f} kg/s needs more feedwater than the subcooler heats "
                           f"({report['m_water_kgs']:.3f} kg/s)", P_H=P_H, P_L=P_L, warnings=warnings)
            states = dict(suction=_state(suction), discharge=_state(discharge), cond_dew=_state(cond_dew),
                          cond_bub=_state(cond_bub),
                          cond_out=_state(cond_out), subcool_out=_state(sc_out), evap_in=_state(evap_in),
                          evap_dew=_state(evap_dew), evap_out=_state(evap_out), ihx_out=_state(ihx_out))
            states["isentropic_discharge"] = _state(isentropic)
            return ConventionalHPResult(True, "success", cop, cop_c, vhc, P_H, P_L, report, states, warnings)
        except (ValueError, ZeroDivisionError) as exc:     # RefpropError subclasses ValueError
            return ConventionalHPResult(False, f"{type(exc).__name__}: {exc}", warnings=warnings)

    # ---- maps ------------------------------------------------------------------------

    def rate(self, **overrides) -> ConventionalHPResult:
        """Solve at this model's config with fields overridden (shares the property package)."""
        return ConventionalHeatPump(replace(self.c, **overrides), props=self.props).solve()

    def rating_map(self, T_lows: Iterable[float], T_highs: Iterable[float], key: str = "T_steam",
                   low_key: str = "T_source", **overrides) -> list[dict]:
        """Flat records over a grid of source (``low_key`` "T_source", default) or evaporating
        ("T_evap") temperatures and steam (``key`` "T_steam", default) or condensing ("T_cond")
        temperatures.

        Each record holds the inputs, success/message, the scalar report entries and any warnings;
        ``pandas.DataFrame(records)`` gives a table.
        """
        if key not in ("T_steam", "T_cond") or low_key not in ("T_source", "T_evap"):
            raise ValueError("key must be 'T_steam' or 'T_cond'; low_key 'T_source' or 'T_evap'")
        records = []
        T_highs = list(T_highs)
        for T_l in T_lows:
            for T_h in T_highs:
                res = self.rate(**{low_key: float(T_l), key: float(T_h)}, **overrides)
                rec = dict(success=res.success, message=res.message, warnings="; ".join(res.warnings))
                rec.update({k: v for k, v in res.report.items() if not isinstance(v, (dict, list))})
                rec[low_key], rec[key] = float(T_l), float(T_h)
                records.append(rec)
        return records
