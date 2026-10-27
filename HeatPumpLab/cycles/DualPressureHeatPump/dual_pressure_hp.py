"""Dual-pressure condensation heat pump for air drying: one evaporator, two condensing pressures in series on the air,
after the dual-pressure condensation high-temperature heat pumps of Dai et al. (Energy Conversion and Management 218
(2020) 112997), which approximate a gliding sink better than a single condensing pressure.

Flow path, per kg of refrigerant through the evaporator::

    evaporator -> low-stage compressor (P_e -> P_L) -> split:
        1 - a  -> low-pressure condenser (P_L)
        a      -> high-stage compressor (P_L -> P_H) -> high-pressure condenser (P_H) -> valve to P_L -> joins the
                  low-pressure condenser inlet
    low-pressure condenser outlet (all the flow) -> valve -> evaporator

The air passes the low-pressure condenser first, from the inlet to T_mid, then the high-pressure condenser, from T_mid
to the target. The high-pressure liquid is throttled back to P_L and mixed into the low-pressure condenser, where its
flash vapour condenses and its liquid is subcooled with the rest (one of the arrangements Dai et al. compare; this
model does not reproduce their water-sink cases). Each condenser's outlet is the coldest that keeps the air approach.

Unknowns: the split a follows from the air balance (the high-pressure condenser heats the air from T_mid to the target,
the low-pressure one from the inlet to T_mid); T_mid is the highest the low-pressure condenser allows (its pinch at
the low condensing temperature); P_H is the lowest that heats the air to the target. P_L is set by its saturation
temperature ``T_cond_low`` (C), and ``optimise()`` searches it for the best COP. Suction superheat is at least
``superheat`` (K), raised if needed for a dry discharge 5 K above the dew point.

Units: MPa, C, kJ/kg, per kg of evaporator flow.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

from ...utils.property_cache import active_property_call_cache, property_call_cache
from ..common.dryer import (
    DEFAULT_MOTOR_EFF,
    AirModel,
    Fluid,
    LoopResult,
    SOLVE_ERRORS,
    air_profile,
    air_secondary,
    comp_path,
    compass_search,
    cooler_outlet,
    dry_air,
    dry_superheat,
    humid_air,
    iso_path,
    js_round,
    min_gap,
    source_T,
    stage,
)

DualPressureResult = LoopResult


@dataclass
class DualPressureConfig:
    """Design inputs. Temperatures in C, approaches in K."""

    fluid: str = "Cyclopentane"
    T_cond_low: float = 110.0         # saturation (dew) temperature at the low condensing pressure P_L
    superheat: float = 5.0            # minimum evaporator superheat
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0
    dT_GC_min: float = 20.0           # air approach in both condensers
    T_source_in: float = 60.0
    T_source_out: float = 50.0
    dT_src: float = 5.0               # evaporating temperature below the source water outlet
    air_ratio: float = 0.0            # secondary / primary air (heated with the primary from the inlet to T_mt)
    T_mt: float = 102.0
    P_max_bar: float = 60.0           # highest high-side pressure searched
    comp_model: int = 5               # 0 ideal, 1 constant eta_const, 2 polynomial, 3 Brendel, 5 Pareto front
    eta_const: float = 0.73
    motor_eff: float = DEFAULT_MOTOR_EFF
    air_RH: float = 60.0
    P_air_kPa: float = 101.325
    property_backend: str | None = None   # "refprop" | "coolprop"; None -> $HEATPUMPLAB_PROPERTY_BACKEND


def _fail(message: str, **extra) -> LoopResult:
    return LoopResult(False, message, **extra)


class DualPressureHeatPump:
    """Dual-pressure condensation heat pump; see the module docstring."""

    def __init__(self, config: DualPressureConfig | None = None, fluid: Fluid | None = None, air: AirModel | None = None):
        self.c = config if config is not None else DualPressureConfig()
        c = self.c
        self.fluid = fluid if fluid is not None else Fluid(c.fluid, c.property_backend)
        if air is None:
            air = humid_air(c.T_air_inlet, c.air_RH, c.P_air_kPa) if c.air_RH > 0 else \
                dry_air(c.P_air_kPa, c.property_backend)
            air = air_secondary(air, c.T_air_inlet, c.T_mt, c.air_ratio, c.T_air_tar)
        self.air = air

    # ---- public ------------------------------------------------------------------------

    def solve(self, **overrides) -> LoopResult:
        """The design at the configuration (fields can be overridden by keyword, except the fluid and air)."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._solve(c))

    def optimise(self, max_eval: int = 120, **overrides) -> LoopResult:
        """Best-COP low condensing temperature (compass search from the configured value and a coarse grid), rounded
        to 0.1 K. The result is the full solution there, with the search in report["optimiser"]."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._optimise(c, max_eval))

    @staticmethod
    def _cached(fn):
        if active_property_call_cache() is None:
            with property_call_cache():
                return fn()
        return fn()

    # ---- model -------------------------------------------------------------------------

    def _at(self, c: DualPressureConfig, P_H: float, Pe: float, Te: float, P_L: float, low: tuple) -> dict:
        """All states at a high pressure P_H (raises ValueError when the air cannot be heated to the target)."""
        fl, air = self.fluid, self.air
        sh, suc, lo = low
        h1 = lo["out"]["h"]                                  # low-stage discharge at P_L
        hi = stage(fl, P_L, lo["out"], P_H, None, c.comp_model, c.eta_const)
        h2 = hi["out"]["h"]                                  # high-stage discharge at P_H
        Tain, Ttar = c.T_air_inlet, c.T_air_tar
        ha_in, ha_tar = air.h(Tain), air.h(Ttar)
        T_dew_L = fl.PQ(P_L, 1.0)["T"]

        def sections(T_mid):
            """Both condensers for the air leaving the low-pressure one at T_mid; the split from the air balance."""
            h3H = cooler_outlet(fl, air, P_H, h2, T_mid, Ttar, c.dT_GC_min)
            qH = h2 - h3H                                    # per kg through the high stage
            a = 0.5
            for _ in range(60):
                h_inL = (1 - a) * h1 + a * h3H               # low-stage discharge + throttled high-pressure liquid
                h3L = cooler_outlet(fl, air, P_L, h_inL, Tain, T_mid, c.dT_GC_min)
                qL = h_inL - h3L                             # per kg of evaporator flow
                # air balance: a qH heats it T_mid -> target, qL heats it inlet -> T_mid (the same air flow)
                a_new = qL * (ha_tar - air.h(T_mid)) / (qH * (air.h(T_mid) - ha_in))
                if not 0 < a_new < 1:
                    raise ValueError("no split of the flow balances the two condensers")
                if abs(a_new - a) < 1e-9:
                    a = a_new
                    break
                a = a + 0.8 * (a_new - a)
            h_inL = (1 - a) * h1 + a * h3H
            h3L = cooler_outlet(fl, air, P_L, h_inL, Tain, T_mid, c.dT_GC_min)
            return dict(a=a, h3H=h3H, qH=qH, h_inL=h_inL, h3L=h3L, qL=h_inL - h3L)

        # the highest T_mid the low-pressure condenser allows: its pinch is at its condensing temperature, so scan
        # down from there and bisect on feasibility
        top = T_dew_L - c.dT_GC_min
        if top <= Tain + 1:
            raise ValueError(f"the low condensing temperature ({T_dew_L:.0f} C) is too close to the air inlet")
        best, lo_T, hi_T = None, Tain + 1.0, None
        for k in range(21):
            T = top - (top - Tain - 1.0) * k / 20
            try:
                S = sections(T)
                best, lo_T = (T, S), T
                break
            except SOLVE_ERRORS:
                hi_T = T
        if best is None:
            raise ValueError("the low-pressure condenser cannot take any of the air heating")
        if hi_T is not None:
            for _ in range(25):
                if hi_T - lo_T < 0.02:
                    break
                m = (lo_T + hi_T) / 2
                try:
                    best, lo_T = (m, sections(m)), m
                except SOLVE_ERRORS:
                    hi_T = m
        T_mid, S = best
        a = S["a"]
        he = S["h3L"]
        q_air = a * S["qH"] + S["qL"]
        w = lo["w"] + a * hi["w"]
        q_src = suc["h"] - he
        gapH = min_gap(air_profile(fl, air, P_H, h2, S["h3H"], T_mid, Ttar, 18))
        gapL = min_gap(air_profile(fl, air, P_L, S["h_inL"], S["h3L"], Tain, T_mid, 18))
        return dict(P_H=P_H, P_L=P_L, Pe=Pe, Te=Te, sh=sh, suc=suc, lo=lo, hi=hi, T_mid=T_mid, h1=h1, h2=h2,
                    he=he, q_air=q_air, w=w, q_src=q_src, gap=min(gapH, gapL), gapH=gapH, gapL=gapL, **S)

    def _solve(self, c: DualPressureConfig) -> LoopResult:
        fl, warnings = self.fluid, []
        try:
            Te = (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_src
            Pe = fl.TQ(Te, 1.0)["P"]
            if math.isfinite(fl.T_crit) and c.T_cond_low >= fl.T_crit - 1:
                return _fail(f"the low condensing temperature must be below {fl.fluid}'s critical temperature")
            if c.T_cond_low <= Te + 5:
                return _fail("the low condensing temperature must be above the evaporating temperature")
            P_L = fl.TQ(c.T_cond_low, 1.0)["P"]
            sh = dry_superheat(fl, Pe, Te, P_L, c.superheat, c.comp_model, c.eta_const)
            suc = fl.PT(Pe, Te + sh)
            low = (sh, suc, stage(fl, Pe, suc, P_L, None, c.comp_model, c.eta_const))
            P_max = c.P_max_bar / 10.0
            P_lo = P_L * 1.05
            if P_lo >= P_max:
                return _fail(f"{fl.fluid} needs more than {c.P_max_bar:g} bar")

            def at(P_H):
                try:
                    return self._at(c, P_H, Pe, Te, P_L, low)
                except SOLVE_ERRORS as exc:
                    return exc

            prev, found, last = P_lo, None, None
            N = 40
            for i in range(N + 1):
                P = P_lo * (P_max / P_lo) ** (i / N)
                r = at(P)
                if isinstance(r, dict) and r["gap"] >= c.dT_GC_min - 0.05:
                    found = (P, r)
                    break
                last = r
                prev = P
            if found is None:
                why = f" ({last})" if isinstance(last, Exception) else ""
                return _fail(f"no high pressure up to {c.P_max_bar:g} bar heats the air to {c.T_air_tar:g} C with a "
                             f"{c.dT_GC_min:g} K approach{why}")
            lo, (hi, S) = prev, found
            for _ in range(30):
                if hi / lo - 1 < 1e-5:
                    break
                m = math.sqrt(lo * hi)
                r = at(m)
                if isinstance(r, dict) and r["gap"] >= c.dT_GC_min - 0.05:
                    hi, S = m, r
                else:
                    lo = m
            return self._result(c, S, warnings)
        except SOLVE_ERRORS as exc:
            return _fail(str(exc), warnings=warnings)

    def _result(self, c: DualPressureConfig, S: dict, warnings: list) -> LoopResult:
        fl, air = self.fluid, self.air
        wEl = S["w"] / c.motor_eff
        cop = S["q_air"] / wEl
        P_H, P_L, Pe = S["P_H"], S["P_L"], S["Pe"]
        supercritical = math.isfinite(fl.P_crit) and P_H >= fl.P_crit
        prs = [S["lo"]["pr"], S["hi"]["pr"]]
        if max(prs) > 20:
            warnings.append(f"stage pressure ratio {max(prs):.0f}, beyond the compressor efficiency correlations "
                            "(fitted up to 20): treat the efficiency, and the COP, with caution")
        dis = [S["lo"]["out"]["T"], S["hi"]["out"]["T"]]
        if max(dis) > 250:
            warnings.append(f"discharge at {max(dis):.0f} C: check the refrigerant's and the oil's stability limits")
        a = S["a"]
        vs = 1.0 / S["suc"]["D"] + a / S["lo"]["out"]["D"]
        VHC = S["q_air"] / vs / 1000
        st = lambda P, h: (lambda x: dict(T=x["T"], P=x["P"], h=x["h"], s=x["s"]))(fl.PH(P, h))
        states = dict(evap_in=st(Pe, S["he"]), suction=st(Pe, S["suc"]["h"]), low_discharge=st(P_L, S["h1"]),
                      high_discharge=st(P_H, S["h2"]), high_condenser_out=st(P_H, S["h3H"]),
                      low_condenser_in=st(P_L, S["h_inL"]), low_condenser_out=st(P_L, S["h3L"]))
        paths = {"Cycle": [
            dict(name="Expansion valve", pts=comp_path(fl, P_L, S["h3L"], Pe, S["he"])),
            dict(name="Evaporator", pts=iso_path(fl, Pe, S["he"], S["suc"]["h"])),
            dict(name="Low-stage compressor", pts=comp_path(fl, Pe, S["suc"]["h"], P_L, S["h1"])),
            dict(name="High-stage compressor", pts=comp_path(fl, P_L, S["h1"], P_H, S["h2"])),
            dict(name="High-pressure condenser" + (" / gas cooler" if supercritical else ""),
                 pts=iso_path(fl, P_H, S["h2"], S["h3H"])),
            dict(name="High-pressure valve", pts=comp_path(fl, P_H, S["h3H"], P_L, S["h3H"])),
            dict(name="Low-pressure condenser", pts=iso_path(fl, P_L, S["h_inL"], S["h3L"]))]}
        pH = air_profile(fl, air, P_H, S["h2"], S["h3H"], S["T_mid"], c.T_air_tar, 18)
        pL = air_profile(fl, air, P_L, S["h_inL"], S["h3L"], c.T_air_inlet, S["T_mid"], 18)
        exchangers = [
            dict(name="Low-pressure condenser to air", hot=[[p[0], p[1]] for p in pL], cold=[[p[0], p[2]] for p in pL],
                 other="dryer air", gap=S["gapL"]),
            dict(name="High-pressure condenser to air", hot=[[a * p[0], p[1]] for p in pH],
                 cold=[[a * p[0], p[2]] for p in pH], other="dryer air", gap=S["gapH"]),
            dict(name="Evaporator",
                 hot=[[0.0, source_T(c.T_source_in, c.T_source_out, 0.0)],
                      [S["q_src"], source_T(c.T_source_in, c.T_source_out, 1.0)]],
                 cold=[[p[2] - S["he"], p[0]] for p in iso_path(fl, Pe, S["he"], S["suc"]["h"])],
                 other="source water")]
        report = dict(COP=cop, VHC=VHC, fluid=fl.fluid, P_high_bar=P_H * 10, P_low_cond_bar=P_L * 10,
                      P_evap_bar=Pe * 10, T_cond_low=c.T_cond_low, T_evap=S["Te"], T_mid=S["T_mid"],
                      split_high=a, superheat=S["sh"], supercritical=supercritical, w_el_kJkg=wEl,
                      q_air=S["q_air"], q_air_low=S["qL"], q_air_high=a * S["qH"], q_evap=S["q_src"],
                      gap_air=S["gap"], gap_low=S["gapL"], gap_high=S["gapH"], T_dis_low=dis[0], T_dis_high=dis[1],
                      air_per_ref=S["q_air"] / (air.h(c.T_air_tar) - air.h(c.T_air_inlet)), air_W=air.W,
                      balance=S["q_air"] - S["w"] - S["q_src"])
        return LoopResult(True, "success", COP=cop, VHC=VHC, Pmax_bar=P_H * 10, T_dis_max=max(dis), q_air=S["q_air"],
                          w_shaft=S["w"], q_source=S["q_src"], report=report, states=states, paths=paths,
                          exchangers=exchangers, warnings=warnings)

    # ---- optimiser ---------------------------------------------------------------------

    def _optimise(self, c: DualPressureConfig, max_eval: int) -> LoopResult:
        fl = self.fluid
        state = dict(nEval=0)
        memo: dict = {}
        Te = (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_src
        top = min(c.T_air_tar, (fl.T_crit - 2) if math.isfinite(fl.T_crit) else c.T_air_tar)
        lo, hi = [max(Te + 10.0, c.T_air_inlet + c.dT_GC_min + 5.0)], [top - 5.0]
        if lo[0] >= hi[0]:
            return _fail("no room for a low condensing temperature between the source and the air target")

        def f(x):
            key = f"{x[0]:.4f}"
            if key in memo:
                return memo[key]
            state["nEval"] += 1
            r = self._solve(replace(c, T_cond_low=float(x[0])))
            o = dict(y=r.COP if r.success else -math.inf, r=r)
            memo[key] = o
            return o

        x0 = [min(max(c.T_cond_low, lo[0]), hi[0])]
        b0 = f(x0)
        xg, bg = None, dict(y=-math.inf)
        for frac in (0.2, 0.35, 0.5, 0.65, 0.8):
            z = [lo[0] + frac * (hi[0] - lo[0])]
            o = f(z)
            if o["y"] > bg["y"]:
                bg, xg = o, z
        x, best, converged = compass_search(f, [(x0, b0), (xg, bg)], lo, hi, [8.0], [0.05], [False], state, max_eval,
                                            score=lambda o: o["y"], empty={"y": -math.inf})
        if x is None:
            return _fail("no low condensing temperature gives a feasible design",
                         report=dict(optimiser=dict(nEval=state["nEval"])))
        xr = [js_round(x[0], 1)]
        o = f(xr)
        if math.isfinite(o["y"]):
            x, best = xr, o
        final = best["r"]
        final.report["optimiser"] = dict(COP=best["y"], start=b0["y"], T_cond_low=float(x[0]), nEval=state["nEval"],
                                         converged=converged, lo=lo, hi=hi)
        return final
