"""Supercritical CO2 reverse Brayton heat pump with a split compression and a preheater, after architectures 1 and 2
of Dalal, Klöppel, Rose and Walmsley (Euro sCO2 2027, "Evaluating sCO2 heat pump architectures for large-scale
industrial air drying above 150 °C").

Low side (suction pressure P_L, usually above the critical pressure): turbine (expander) -> LTHX (heated by the
source) -> IHX cold side -> split. A fraction x (the secondary stream) goes to compressor 1; its hot discharge
preheats the main stream (1 - x) in the preheater before compressor 2.

- Architecture 1: compressor 2 heats the air in the HTHX; the secondary stream, after the preheater, rejoins the main
  stream after the HTHX; together they give heat back in the IHX hot side, then expand in the turbine.
- Architecture 2 (with Held's split stream): the secondary stream, after the preheater, heats the cooler part of the
  air in HTHX 2 and the main stream the hotter part in HTHX 1; only the main stream passes the IHX hot side; the two
  rejoin before the turbine. The air temperature between the HTHXs is the lowest at which the main stream still
  keeps its approach (both exchangers at their pinch).

The recuperators (IHX, preheater) are set by an effectiveness: the duty is that fraction of the largest the
exchanger could pass (found along the profiles, so the pinch can sit inside near the pseudo-critical point). The
LTHX heats the CO2 as far as the source allows with the minimum approach dT_src (optionally capped at q_lthx_max per
kg). The HTHXs cool the CO2 to the coldest outlet that keeps the air approach dT_GC_min (they heat the dryer air
directly). The cycle closes by a damped fixed-point iteration on the suction enthalpy and the LTHX outlet.

The discharge pressure P_H is an input, or with ``P_H_min`` the lowest at which compressor 2's discharge reaches the
air target plus the approach. ``optimise`` maximises the COP over the suction pressure, the discharge pressure and
the split. Pressure drops: dP_hp per exchanger on the high side, dP_lp on the low side. Shaft work: compressors /
eta_m less turbine x eta_m; electrical: shaft / motor_eff. Duties per kg of total CO2; the swept volume (VHC) is the
total suction volume of both compressors.

Port of the "Advanced dryer heat pumps" web app (Ahuora Centre, University of Waikato): DRYCASC.braytonAt,
braytonMinP, brayton and braytonOptimiseSteps. Units follow HeatPumpLab: MPa, C, kJ/kg, kg/m3.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field, replace
from typing import Callable

from ..common.dryer import (
    MOTOR, SOLVE_ERRORS, AirModel, Fluid, air_profile, comp_path, cooler_outlet, iso_path, lin, max_duty,
    min_gap, point, source_T, stage,
)

# ----------------------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------------------


@dataclass
class BraytonHPConfig:
    """Design inputs (the web app's defaults). Pressures in MPa, temperatures in C, approaches in K."""

    arch: int = 1                    # 1 or 2 (see the module docstring)
    fluid: str = "CarbonDioxide"
    P_L: float = 8.0                 # suction pressure
    P_H: float = 16.5                # discharge pressure (ignored with P_H_min)
    P_H_min: bool = False            # use the lowest discharge pressure that reaches the air target
    x: float = 0.37                  # secondary (compressor 1) mass fraction, 0.01-0.95
    eps_ihx: float = 0.85
    eps_pre: float = 0.85
    eta_c: float = 0.85              # compressors' isentropic efficiency
    eta_t: float = 0.85              # turbine
    eta_m: float = 0.99              # mechanical (compressors and turbine)
    dP_hp: float = 0.05              # per exchanger on the high side (50 kPa)
    dP_lp: float = 0.01              # per exchanger on the low side (10 kPa)
    Pmax: float = 30.0               # upper bound of the discharge-pressure search and the optimiser
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0
    dT_GC_min: float = 20.0          # minimum approach to the air
    T_source_in: float = 40.0
    T_source_out: float = 30.0
    dT_src: float = 5.0              # minimum LTHX approach to the source
    source_curve: Callable[[float], float] | None = None   # source T at a duty fraction from its cold end (else water)
    source_name: str = "source water"
    q_lthx_max: float | None = None  # optional cap on the LTHX duty, kJ per kg of CO2
    motor_eff: float = MOTOR
    air_RH: float = 60.0             # dryer air inlet relative humidity, %
    P_atm: float = 101.325           # kPa
    backend: str | None = None       # property backend: "coolprop", "refprop" or None (the package default)


@dataclass
class BraytonHPResult:
    success: bool
    message: str
    COP: float | None = None
    VHC: float | None = None          # MJ/m3, air duty per total suction volume of both compressors
    q_air: float | None = None        # kJ per kg of CO2
    w_shaft: float | None = None
    q_source: float | None = None
    P_H: float | None = None          # MPa
    report: dict = field(default_factory=dict)
    states: dict = field(default_factory=dict)
    paths: list = field(default_factory=list)       # [{name, pts: [[T, P bar, h, s]]}]
    profiles: list = field(default_factory=list)    # exchanger profiles [{name, hot: [[q, T]], cold: [[q, T]]}]
    warnings: list = field(default_factory=list)

    def to_dict(self):
        return asdict(self)


@dataclass
class BraytonOptimum:
    success: bool
    message: str = ""
    COP: float | None = None
    P_L: float | None = None          # MPa
    P_H: float | None = None
    x: float | None = None
    start: float | None = None        # COP of the design as set (-inf when it fails)
    start_ok: bool = False
    n_eval: int = 0
    at_bound: dict = field(default_factory=dict)
    budget_hit: bool = False
    result: BraytonHPResult | None = None


def _fail(message: str, **kw) -> BraytonHPResult:
    return BraytonHPResult(False, message, **kw)


def _state(st: dict) -> dict:
    return {k: st.get(k) for k in ("T", "P", "h", "s", "D", "q", "phase")}


# ----------------------------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------------------------


class BraytonHeatPump:
    """sCO2 reverse Brayton heat pump (architectures 1 and 2); see the module docstring."""

    def __init__(self, config: BraytonHPConfig | None = None, fluid: Fluid | None = None, air: AirModel | None = None):
        self.c = config if config is not None else BraytonHPConfig()
        self.pr = fluid if fluid is not None else Fluid(self.c.fluid, self.c.backend)
        self.air = air if air is not None else AirModel.humid(self.c.T_air_inlet, self.c.air_RH, self.c.P_atm * 1000)

    def _Tph(self, P, h):
        return self.pr.T_ph_fast(P, h)

    def _src_T(self, f):
        c = self.c
        return source_T(c.T_source_in, c.T_source_out, f, c.source_curve)

    # ---- the cycle at a discharge pressure (braytonAt) --------------------------------------

    def cycle_at(self, P_H: float) -> dict:
        c, pr, air = self.c, self.pr, self.air
        x = min(max(float(c.x), 0.01), 0.95)
        arch = 2 if int(c.arch) == 2 else 1
        Tain, Ttar, a = c.T_air_inlet, c.T_air_tar, c.dT_GC_min
        dPh, dPl, PL = max(c.dP_hp or 0.0, 0.0), max(c.dP_lp or 0.0, 0.0), c.P_L
        Pto, Pl1, Pc2, Ph1, Pti = PL + 2 * dPl, PL + dPl, PL - dPl, P_H - dPh, P_H - 2 * dPh
        Tph = self._Tph
        dTs = c.dT_src
        cap = max(float(c.q_lthx_max), 0.0) if c.q_lthx_max is not None and math.isfinite(c.q_lthx_max) else math.inf

        def lthx_approach(ht):
            # the CO2 is heated as far as the source allows while keeping the source approach everywhere along the
            # exchanger (counterflow; the source profile is a fixed shape in the duty fraction)
            h_max = pr.PT(Pl1, c.T_source_in - dTs)["h"]
            if not h_max > ht:
                return ht

            def ok(h1):
                for f in lin(0.0, 1.0, 16):
                    if Tph(Pl1, ht + f * (h1 - ht)) > self._src_T(f) - dTs + 1e-6:
                        return False
                return True

            if ok(h_max):
                return h_max
            if not ok(ht + 1e-6):
                return ht
            lo, hi, k = ht, h_max, 0
            while k < 40 and hi - lo > 0.01:
                m = (lo + hi) / 2
                if ok(m):
                    lo = m
                else:
                    hi = m
                k += 1
            return lo

        def lthx_out(ht):
            return min(lthx_approach(ht), ht + cap)

        # start: the suction at the source inlet temperature less 20 K
        hci = pr.PT(PL, c.T_source_in - 20)["h"]
        S, d = None, math.inf
        for it in range(200):
            ci = pr.PH(PL, hci)
            s1 = stage(pr, PL, ci, P_H, c.eta_c)
            # preheater: secondary (x, hot, after compressor 1) against the main stream (1 - x, cold) before
            # compressor 2
            mh = x / (1 - x)
            qp = c.eps_pre * max_duty(pr, P_H, s1["out"]["h"], mh, pr, Pc2, hci, 0.0, fast=True)
            hc2 = hci + qp
            hpho = s1["out"]["h"] - qp / mh
            s2 = stage(pr, Pc2, pr.PH(Pc2, hc2), P_H, c.eta_c)
            h2, T2 = s2["out"]["h"], s2["out"]["T"]
            # HTHX(s): the coldest outlets that keep the air approach. Below the target plus the approach there is no
            # design; the air is then taken only to the discharge less the approach, so the discharge temperature
            # stays continuous in the pressure and the pressure search brackets the target
            ok_hot = T2 >= Ttar + a - 1e-9
            Tt = max(min(Ttar, T2 - a), Tain + 1)

            def cool(P, h_in, Ta_in, Ta_out, dT):
                try:
                    return cooler_outlet(pr, air, P, h_in, Ta_in, Ta_out, dT, fast=True)
                except SOLVE_ERRORS:
                    if ok_hot:
                        raise
                    return h_in

            h5, Tmid = None, None
            if arch == 1:
                h3 = cool(P_H, h2, Tain, Tt, min(a, T2 - Tt))
                h_hot_in, m_hot = (1 - x) * h3 + x * hpho, 1.0
            else:
                Tpho = Tph(Ph1, hpho)

                def sec(Tm):
                    hh5 = cool(Ph1, hpho, Tain, Tm, a)
                    return hh5, x * (hpho - hh5) / (air.h(Tm) - air.h(Tain))

                def main(Tm, m_air):                    # (1 - x) h3 from the balance
                    return (1 - x) * h2 - m_air * (air.h(Tt) - air.h(Tm))

                if Tpho - a > Tain + 0.5 and Tt > Tain + 1:
                    # the lowest air temperature between the two HTHXs at which the main stream still keeps its
                    # approach: both exchangers then sit at their pinch (the most heat to the air)
                    aa = min(a, T2 - Tt)

                    def G(Tm):
                        _, m_air = sec(Tm)
                        return main(Tm, m_air) / (1 - x) - cool(P_H, h2, Tm, Tt, aa)

                    lo, hi = Tain + 0.5, min(Tpho - a, Tt - 0.5)
                    if G(hi) < 0:
                        Tmid = hi
                    elif G(lo) >= 0:
                        Tmid = lo
                    else:
                        k = 0
                        while k < 40 and hi - lo > 0.02:
                            m = (lo + hi) / 2
                            if G(m) >= 0:
                                hi = m
                            else:
                                lo = m
                            k += 1
                        Tmid = hi
                    h5, m_air = sec(Tmid)
                    h3 = main(Tmid, m_air) / (1 - x)
                else:
                    Tmid, h5 = Tain, hpho
                    h3 = cool(P_H, h2, Tain, Tt, min(a, T2 - Tt))
                h_hot_in, m_hot = h3, 1 - x
            # IHX: the hot side (all the flow in architecture 1, the main stream in 2) heats the low side after the LTHX
            P_hot = Ph1
            h1 = S["h1"] if S else pr.PT(Pl1, c.T_source_in - 10)["h"]
            qi = c.eps_ihx * max_duty(pr, P_hot, h_hot_in, m_hot, pr, Pl1, h1, 0.0, fast=True)
            hci2 = h1 + qi
            hho = h_hot_in - qi / m_hot
            hti = hho if arch == 1 else (1 - x) * hho + x * h5
            ti = pr.PH(Pti, hti)
            ht = hti - c.eta_t * (hti - pr.PS(Pto, ti["s"])["h"])
            h1n = lthx_out(ht)
            S = dict(ci=ci, s1=s1, s2=s2, hc2=hc2, hpho=hpho, h2=h2, T2=T2, ok_hot=ok_hot, h3=h3, h5=h5, Tmid=Tmid,
                     h_hot_in=h_hot_in, m_hot=m_hot, h1=h1n, qi=qi, hho=hho, hti=hti, ht=ht, qp=qp, x=x, arch=arch,
                     P_hot=P_hot)
            d = abs(hci2 - hci) + abs(h1n - h1)
            if d < 1e-3 and it > 2:
                break
            hci = hci2 if it < 3 else hci + (0.5 if it < 40 else 0.2) * (hci2 - hci)   # damped: the recuperators feed back
        S["hci"] = hci
        S["P"] = dict(PL=PL, Pto=Pto, Pl1=Pl1, Pc2=Pc2, P_H=P_H, Ph1=Ph1, Pti=Pti)
        S["conv"] = d < 1e-2
        return S

    def min_pressure_cycle(self) -> dict:
        """The cycle at the lowest discharge pressure at which compressor 2 reaches the air target plus the approach
        (bisection; the discharge temperature is continuous in the pressure, see ``cycle_at``)."""
        c = self.c
        target = c.T_air_tar + c.dT_GC_min

        def g(P):
            S = self.cycle_at(P)
            return S, S["T2"] - target

        lo = c.P_L * 1.15
        hi = max(c.Pmax or 30.0, lo * 1.01)
        Sl, vl = g(lo)
        Sh, vh = g(hi)
        if vh < 0:
            return Sh
        if vl >= 0:
            return Sl
        S, k = Sh, 0
        while k < 30 and hi - lo > 1e-3:
            m = (lo + hi) / 2
            Sm, vm = g(m)
            if vm >= 0:
                hi, S = m, Sm
            else:
                lo = m
            k += 1
        return S

    # ---- solve (brayton) -------------------------------------------------------------------

    def solve(self, quick: bool = False) -> BraytonHPResult:
        """The design at the configured inputs. ``quick`` skips the lowest-pressure search that a failure below the
        target would otherwise report (as the optimiser does)."""
        try:
            return self._solve(quick)
        except SOLVE_ERRORS as exc:
            return _fail(f"{exc}")

    def _solve(self, quick: bool) -> BraytonHPResult:
        c, pr, air = self.c, self.pr, self.air
        warnings: list[str] = []
        a, dTs = c.dT_GC_min, c.dT_src
        Tph = self._Tph
        if not c.P_L > 0:
            return _fail("set the suction pressure")
        if c.P_L < pr.P_crit:
            warnings.append(f"suction {c.P_L * 10:.0f} bar is below the critical pressure ({pr.P_crit * 10:.1f} bar): "
                            "the low side runs subcritical and the turbine may exhaust two-phase")

        def src_check(S):
            Tto, Tlim = Tph(S["P"]["Pl1"], S["ht"]), self._src_T(0.0) - dTs
            if Tto > Tlim + 0.05:
                return _fail(f"LTHX: the CO₂ leaves the turbine at {Tto:.1f} °C, above the {c.source_name} outlet "
                             f"({self._src_T(0.0):.1f} °C) less the {dTs:g} K source approach, so the source cannot "
                             "heat it: raise the source outlet temperature, lower the source approach, or lower the "
                             "suction pressure")
            return None

        if c.P_H_min:
            S = self.min_pressure_cycle()
            sc = src_check(S)
            if sc:
                return sc
            if not S["ok_hot"]:
                return _fail(f"compressor 2 reaches only {S['T2']:.0f} °C at {S['P']['P_H'] * 10:.0f} bar, short of the "
                             f"{c.T_air_tar:g} °C target plus the {a:g} K approach: raise the max pressure or the "
                             "preheater/IHX effectiveness")
        else:
            if not c.P_H > c.P_L * 1.05:
                return _fail("the discharge pressure must be above the suction pressure")
            S = self.cycle_at(c.P_H)
            sc = src_check(S)
            if sc:
                return sc
            if not S["ok_hot"]:
                m = {"ok_hot": False} if quick else self.min_pressure_cycle()
                at_least = f" to at least {m['P']['P_H'] * 10:.1f} bar" if m["ok_hot"] else ""
                return _fail(f"compressor 2 discharges at {S['T2']:.1f} °C, below the {c.T_air_tar:g} °C target plus "
                             f"the {a:g} K approach: raise the discharge pressure{at_least}", P_H=c.P_H)
        x, arch, P = S["x"], S["arch"], S["P"]
        Tain, Ttar = c.T_air_inlet, c.T_air_tar
        # duties and work per kg of total CO2
        qA1 = (1 - x) * (S["h2"] - S["h3"])
        qA2 = x * (S["hpho"] - S["h5"]) if arch == 2 else 0.0
        q_air = qA1 + qA2
        wc = x * S["s1"]["w"] + (1 - x) * S["s2"]["w"]
        wt = S["hti"] - S["ht"]
        w_shaft = wc / c.eta_m - wt * c.eta_m
        w_el = w_shaft / c.motor_eff
        cop = q_air / w_el
        qE = S["h1"] - S["ht"]
        if not w_shaft > 0:
            return _fail("the turbine gives back more than the compressors take")
        # the loop has settled when the cycle's energy balance closes (source + compressors = air + turbine)
        bal = qE + wc - wt - q_air
        if abs(bal) > 0.005 * q_air + 0.2:
            return _fail(f"the cycle does not settle (energy balance off by {bal:.1f} kJ/kg): the IHX and preheater feed "
                         "heat back to the suction without a steady state at this pressure and split; raise the "
                         "discharge pressure or lower the split")

        # LTHX against the source: no temperature cross
        def cold_T(f):
            return Tph(P["Pl1"], S["ht"] + f * qE)

        gapE = min(self._src_T(f) - cold_T(f) for f in lin(0.0, 1.0, 16))
        if gapE < dTs - 0.1:
            return _fail(f"LTHX approach {gapE:.1f} K, below the {dTs:g} K source approach")
        # gaps on the air
        pA1 = air_profile(pr, air, P["P_H"], S["h2"], S["h3"], S["Tmid"] if arch == 2 else Tain, Ttar, 16, fast=True)
        gap1 = min_gap(pA1)
        pA2 = air_profile(pr, air, P["Ph1"], S["hpho"], S["h5"], Tain, S["Tmid"], 16, fast=True) \
            if arch == 2 and qA2 > 1e-9 else None
        gap2 = min_gap(pA2) if pA2 else None
        if gap1 < a - 0.1 or (gap2 is not None and gap2 < a - 0.1):
            warnings.append(f"air approach {min(gap1, gap2 if gap2 is not None else math.inf):.1f} K, below the {a:g} K set")
        c2_in = pr.PH(P["Pc2"], S["hc2"])
        v_suc = x / S["ci"]["D"] + (1 - x) / c2_in["D"]
        st = dict(
            suction=S["ci"], c1_out=S["s1"]["out"], c2_in=c2_in, c2_out=S["s2"]["out"], hthx_out=pr.PH(P["P_H"], S["h3"]),
            pre_hot_out=pr.PH(P["Ph1"], S["hpho"]), ihx_hot_in=pr.PH(S["P_hot"], S["h_hot_in"]),
            ihx_hot_out=pr.PH(S["P_hot"], S["hho"]), turb_in=pr.PH(P["Pti"], S["hti"]), turb_out=pr.PH(P["Pto"], S["ht"]),
            lthx_out=pr.PH(P["Pl1"], S["h1"]),
        )
        if arch == 2:
            st["hthx2_out"] = pr.PH(P["Ph1"], S["h5"])
        paths = [
            {"name": "Turbine (expander)", "pts": comp_path(pr, P["Pti"], S["hti"], P["Pto"], S["ht"])},
            {"name": "LTHX (source)", "pts": iso_path(pr, P["Pto"], S["ht"], S["h1"])},
            {"name": "Internal heat exchanger (IHX)", "pts": iso_path(pr, P["Pl1"], S["h1"], S["hci"])},
            {"name": "Compressor 1 (secondary)", "pts": comp_path(pr, P["PL"], S["hci"], P["P_H"], S["s1"]["out"]["h"])},
            {"name": "Preheater recuperator", "pts": iso_path(pr, P["Pc2"], S["hci"], S["hc2"])},
            {"name": "Compressor 2 (main)", "pts": comp_path(pr, P["Pc2"], S["hc2"], P["P_H"], S["h2"])},
            {"name": "HTHX 1 gas cooler" if arch == 2 else "HTHX gas cooler", "pts": iso_path(pr, P["P_H"], S["h2"], S["h3"])},
            {"name": "Preheater recuperator", "pts": iso_path(pr, P["P_H"], S["s1"]["out"]["h"], S["hpho"])},
        ]
        if arch == 2:
            paths.append({"name": "HTHX 2 gas cooler", "pts": iso_path(pr, P["Ph1"], S["hpho"], S["h5"])})
        else:
            paths += [{"name": "Mixing", "pts": [point(st["pre_hot_out"]), point(st["ihx_hot_in"])]},
                      {"name": "Mixing", "pts": [point(st["hthx_out"]), point(st["ihx_hot_in"])]}]
        paths.append({"name": "Internal heat exchanger (IHX)", "pts": iso_path(pr, S["P_hot"], S["h_hot_in"], S["hho"])})
        if arch == 2:
            paths += [{"name": "Mixing", "pts": [point(st["ihx_hot_out"]), point(st["turb_in"])]},
                      {"name": "Mixing", "pts": [point(st["hthx2_out"]), point(st["turb_in"])]}]
        qI, qP = S["qi"], (1 - x) * S["qp"]
        profiles = [
            {"name": "HTHX 1: main stream to air" if arch == 2 else "HTHX: CO2 to air",
             "hot": [[p[0] * (1 - x), p[1]] for p in pA1], "cold": [[p[0] * (1 - x), p[2]] for p in pA1]},
        ]
        if pA2:
            profiles.append({"name": "HTHX 2: secondary stream to air", "hot": [[p[0] * x, p[1]] for p in pA2],
                             "cold": [[p[0] * x, p[2]] for p in pA2]})
        fs = lin(0.0, 1.0, 16)
        profiles += [
            {"name": "LTHX (source)", "hot": [[f * qE, self._src_T(f)] for f in fs], "cold": [[f * qE, cold_T(f)] for f in fs]},
            {"name": "IHX", "hot": [[f * qI, Tph(S["P_hot"], S["hho"] + f * qI / S["m_hot"])] for f in lin(0, 1, 14)],
             "cold": [[f * qI, Tph(P["Pl1"], S["h1"] + f * qI)] for f in lin(0, 1, 14)]},
            {"name": "Preheater", "hot": [[f * qP, Tph(P["Ph1"], S["hpho"] + f * qP / x)] for f in lin(0, 1, 14)],
             "cold": [[f * qP, Tph(P["Pc2"], S["hci"] + f * S["qp"])] for f in lin(0, 1, 14)]},
        ]
        report = dict(
            COP=cop, COP_shaft=q_air / w_shaft, arch=arch, x=x, P_L_bar=P["PL"] * 10, P_H_bar=P["P_H"] * 10,
            PR=P["P_H"] / P["PL"], T_suction=S["ci"]["T"], T_c1_out=S["s1"]["out"]["T"], T_c2_in=c2_in["T"],
            T_c2_out=S["T2"], T_hthx_out=st["hthx_out"]["T"], T_hthx2_out=st["hthx2_out"]["T"] if arch == 2 else None,
            T_air_mid=S["Tmid"] if arch == 2 else None, T_turb_in=st["turb_in"]["T"], T_turb_out=st["turb_out"]["T"],
            T_lthx_out=st["lthx_out"]["T"], balance_kJkg=bal, w_c1=x * S["s1"]["w"], w_c2=(1 - x) * S["s2"]["w"],
            w_turb=wt, w_shaft=w_shaft, w_el_kJkg=w_el, q_air_1=qA1, q_air_2=qA2, q_ihx=qI, q_pre=qP, q_lthx=qE,
            gap_hthx=gap1, gap_hthx2=gap2, gap_lthx=gapE, eps_ihx=c.eps_ihx, eps_pre=c.eps_pre, dT_src=dTs,
            q_lthx_cap=float(c.q_lthx_max) if c.q_lthx_max is not None and math.isfinite(c.q_lthx_max) else None,
            eta_c=c.eta_c, eta_t=c.eta_t, eta_m=c.eta_m, rho_suction=1 / v_suc, Pmax_bar=P["P_H"] * 10,
            T_dis_max=max(S["s1"]["out"]["T"], S["T2"]), converged=S["conv"],
        )
        return BraytonHPResult(True, "success", cop, q_air / v_suc / 1000, q_air, w_shaft, qE, P["P_H"], report,
                               {k: _state(v) for k, v in st.items()}, paths, profiles, warnings)

    # ---- optimiser (braytonOptimiseSteps) --------------------------------------------------

    def _with(self, **kw) -> "BraytonHeatPump":
        return BraytonHeatPump(replace(self.c, **kw), self.pr, self.air)

    def _PH_x_search(self, PL, v0, step0, cnt, budget, tol=(0.01, 0.004)):
        """Compass search over (P_H, x) at a fixed suction pressure; infeasible points lose."""
        Pmax = self.c.Pmax or 30.0
        lo, hi = [PL * 1.1, 0.02], [Pmax, 0.9]

        def f(v):
            r = self._with(P_L=PL, P_H_min=False, P_H=v[0], x=v[1]).solve(quick=True)
            cnt["n"] += 1
            return r.COP if r.success else -math.inf

        if not lo[0] < hi[0]:
            return {"COP": -math.inf}
        v = [min(max(v0[0], lo[0]), hi[0]), min(max(v0[1], lo[1]), hi[1])]
        best = f(v)
        start = best
        # an infeasible start (pressure too low): walk the pressure up until it solves
        while not math.isfinite(best) and v[0] < hi[0] and cnt["n"] < budget:
            v = [min(v[0] * 1.08, hi[0]), v[1]]
            best = f(v)
        if not math.isfinite(best):
            return {"COP": -math.inf, "start": start}
        step = list(step0)
        while (step[0] > tol[0] or step[1] > tol[1]) and cnt["n"] < budget:
            moved = False
            for d in (0, 1):
                if moved:
                    break
                for sg in (1, -1):
                    w = list(v)
                    w[d] = min(max(w[d] + sg * step[d], lo[d]), hi[d])
                    if w[d] == v[d]:
                        continue
                    fw = f(w)
                    if fw > best + 1e-6:
                        best, v, moved = fw, w, True
                        break
            if not moved:
                step = [t / 2 for t in step]
        return {"COP": best, "P_H": v[0], "x": v[1], "start": start,
                "at_bound": {"P_H": v[0] >= hi[0] - 1e-6 or v[0] <= lo[0] + 1e-6,
                             "x": v[1] <= lo[1] + 1e-6 or v[1] >= hi[1] - 1e-6}}

    def optimise(self, opt_P_L: bool = True, budget: int | None = None) -> BraytonOptimum:
        """Maximise the COP over the suction pressure P_L, the discharge pressure P_H and the split x.

        Inner: a compass search over P_H (from just above the suction to Pmax) and x (0.02-0.9) at a fixed P_L.
        Outer (``opt_P_L``): a compass search on P_L (steps 10, 5, 2.5, 1.25 bar), each trial re-optimising P_H and x
        from the best design so far (P_H scaled with P_L) to a coarser step, then a finish at the chosen P_L. A solve
        budget (600, or 300 for architecture 2) caps the run. P_L runs from 6 bar to Pmax / 1.15. ``result`` is the
        design at the optimum."""
        c = self.c
        Pmax = c.Pmax or 30.0
        cnt = {"n": 0}
        budget = budget or (300 if int(c.arch) == 2 else 600)
        PL_lo, PL_hi = 0.6, Pmax / 1.15
        PL0 = min(max(c.P_L, PL_lo), PL_hi)
        PL = PL0
        v0 = [(c.P_H or 17.0) if c.P_H_min else c.P_H, c.x]
        r0 = self._with(P_L=PL, P_H_min=False, P_H=v0[0], x=v0[1]).solve(quick=True)
        cnt["n"] += 1
        start = r0.COP if r0.success else -math.inf
        cur = self._PH_x_search(PL, v0, [1.0, 0.08], cnt, budget)
        # no design at the set suction pressure (the max pressure too low for it): step the suction pressure down
        k = 0
        while opt_P_L and not math.isfinite(cur["COP"]) and PL > PL_lo + 1e-9 and k < 12 and cnt["n"] < budget:
            PL = max(PL - 1.0, PL_lo)
            cur = self._PH_x_search(PL, [min(v0[0], Pmax), v0[1]], [1.0, 0.08], cnt, budget)
            k += 1
        if not math.isfinite(cur["COP"]):
            span = f" at a suction pressure from {PL * 10:.0f} to {PL0 * 10:.0f} bar" if opt_P_L else ""
            return BraytonOptimum(False, f"no design up to the max pressure{span}", n_eval=cnt["n"], start=start)
        cur["P_L"] = PL
        if opt_P_L:
            s, direction = 1.0, 1
            while s >= 0.125 - 1e-9 and cnt["n"] < budget:
                moved = False
                for sg in (direction, -direction):
                    PLn = min(max(cur["P_L"] + sg * s, PL_lo), PL_hi)
                    if abs(PLn - cur["P_L"]) < 1e-9:
                        continue
                    tr = self._PH_x_search(PLn, [cur["P_H"] * PLn / cur["P_L"], cur["x"]], [0.25, 0.02], cnt, budget,
                                           (0.04, 0.008))
                    if tr["COP"] > cur["COP"] + 1e-5:
                        cur = {**tr, "P_L": PLn}
                        direction, moved = sg, True
                        break
                if not moved:
                    s /= 2
            # the P_L trials stop P_H and x at a coarser step: finish them at the chosen P_L
            fin = self._PH_x_search(cur["P_L"], [cur["P_H"], cur["x"]], [0.04, 0.008], cnt, budget + 80)
            if fin["COP"] > cur["COP"]:
                cur = {**fin, "P_L": cur["P_L"]}
        result = self._with(P_L=cur["P_L"], P_H_min=False, P_H=cur["P_H"], x=cur["x"]).solve()
        at_bound = {"P_L": opt_P_L and (cur["P_L"] <= PL_lo + 1e-6 or cur["P_L"] >= PL_hi - 1e-6),
                    "P_H": cur["at_bound"]["P_H"], "x": cur["at_bound"]["x"]}
        return BraytonOptimum(True, "success", cur["COP"], cur["P_L"], cur["P_H"], cur["x"], start, math.isfinite(start),
                              cnt["n"], at_bound, cnt["n"] >= budget, result)
