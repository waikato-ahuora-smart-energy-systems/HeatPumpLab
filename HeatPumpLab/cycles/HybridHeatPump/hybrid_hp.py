"""Hybrid absorption-compression (Osenbrueck) heat pump with ammonia-water for air drying.

The working fluid is an ammonia-water solution. Per kg of solution circulating through the absorber::

    desorber (P_L): the solution, throttled from the absorber, is heated by the source to T_des_out, where a vapour
        fraction q_des has boiled off
    separator: vapour (ammonia-rich, q_des kg) -> compressor -> P_H
               lean liquid (1 - q_des kg) -> solution pump -> P_H, less any part injected at the compressor suction
               (wet compression, m_inj kg)
    absorber (P_H): vapour and liquid mix and the vapour is absorbed while the solution heats the dryer air
    valve: absorber outlet -> P_L -> desorber

Because the solution keeps one overall composition (``w_rich``, the ammonia mass fraction) through the absorber and the
desorber, both behave like a condenser and an evaporator for a zeotropic fluid with a large temperature glide, set by
the composition and the degassing ``q_des`` (Jensen, Markussen, Reinholdt and Elmegaard, Int. J. Refrigeration 58
(2015) 79-89). Compression of the separated vapour keeps the pressure ratio and the discharge temperature down, but
the high pressure limits the supply temperature: with a 170 C discharge limit Jensen et al. found about 111, 129 and
147 C maximum supply temperatures for 28, 50 and 140 bar components.

Design: the desorber outlet is T_source_in - dT_src (counterflow, against the hottest source water), which with
q_des sets the low pressure; the desorber must keep dT_src from the source water everywhere (a failed design
otherwise). The absorber outlet is the coldest that keeps the air approach dT_GC_min, and the high pressure is the
lowest that heats the air to the target (scanned up to P_max_bar, then bisected). Dry compression of the ammonia-rich
vapour runs very hot at these pressure ratios, so with ``wet_compression`` (the default) the smallest amount of lean
liquid that holds the discharge at T_dis_max (Jensen et al.'s 170 C, for standard ammonia compressors and oil) is
injected at the suction and the two-phase mixture is compressed (the compressor model's efficiency is applied to it
unchanged, an optimistic simplification for wet compression). Without it, a discharge above T_dis_max is a warning,
or a failure with ``limit_discharge``.
``optimise()`` searches w_rich and q_des for the best COP.

Properties: the IAPWS 2001 ammonia-water formulation through teqp (``HeatPumpLab.utils.ammonia_water``); dryer air as
the other cycles. Units: MPa, C, kJ/kg, per kg of solution through the absorber. Needs the ``ammonia-water`` extra.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, replace

from ...utils.ammonia_water import AmmoniaWater, AmmoniaWaterStream
from ..common.dryer import (
    DEFAULT_MOTOR_EFF,
    AirModel,
    LoopResult,
    SOLVE_ERRORS,
    air_profile,
    air_secondary,
    compass_search,
    cooler_outlet,
    dry_air,
    humid_air,
    js_round,
    lin,
    min_gap,
    source_T,
    stage,
)

HybridResult = LoopResult


@dataclass
class HybridConfig:
    """Design inputs. Temperatures in C, approaches in K."""

    w_rich: float = 0.5               # ammonia mass fraction of the circulating solution
    q_des: float = 0.4                # vapour mass fraction boiled off in the desorber (the degassing)
    T_air_inlet: float = 15.0
    T_air_tar: float = 120.0
    dT_GC_min: float = 10.0           # air approach in the absorber
    T_source_in: float = 60.0
    T_source_out: float = 40.0
    dT_src: float = 5.0               # desorber approach to the source water
    air_ratio: float = 0.0            # secondary / primary air (heated with the primary from the inlet to T_mt)
    T_mt: float = 80.0
    P_max_bar: float = 50.0           # component pressure class (Jensen et al.: 28, 50 or 140 bar)
    T_dis_max: float = 170.0          # discharge temperature limit (ammonia compressors and oil)
    wet_compression: bool = True      # inject lean liquid at the suction to hold the discharge at T_dis_max
    limit_discharge: bool = False     # dry compression only: a discharge above T_dis_max fails (else a warning)
    comp_model: int = 5               # 0 ideal, 1 constant eta_const, 2 polynomial, 3 Brendel, 5 Pareto front
    eta_const: float = 0.73
    eta_pump: float = 0.6
    motor_eff: float = DEFAULT_MOTOR_EFF
    air_RH: float = 60.0
    P_air_kPa: float = 101.325
    property_backend: str | None = None   # only for dry air (air_RH = 0)


def _fail(message: str, **extra) -> LoopResult:
    return LoopResult(False, message, **extra)


class HybridHeatPump:
    """Hybrid absorption-compression heat pump; see the module docstring."""

    def __init__(self, config: HybridConfig | None = None, air: AirModel | None = None, model: AmmoniaWater | None = None):
        self.c = config if config is not None else HybridConfig()
        c = self.c
        self.model = model if model is not None else AmmoniaWater()
        if air is None:
            air = humid_air(c.T_air_inlet, c.air_RH, c.P_air_kPa) if c.air_RH > 0 else \
                dry_air(c.P_air_kPa, c.property_backend)
            air = air_secondary(air, c.T_air_inlet, c.T_mt, c.air_ratio, c.T_air_tar)
        self.air = air
        self._streams: dict = {}

    def stream(self, w: float) -> AmmoniaWaterStream:
        """One stream object (with its caches) per composition."""
        key = round(w, 12)
        if key not in self._streams:
            self._streams[key] = AmmoniaWaterStream(w, self.model)
        return self._streams[key]

    # ---- public ------------------------------------------------------------------------

    def solve(self, **overrides) -> LoopResult:
        c = replace(self.c, **overrides) if overrides else self.c
        return self._solve(c)

    def optimise(self, max_eval: int = 80, **overrides) -> LoopResult:
        """Best-COP solution composition and degassing (compass search from the configured design)."""
        c = replace(self.c, **overrides) if overrides else self.c
        state = dict(nEval=0)
        memo: dict = {}
        lo, hi = [0.2, 0.05], [0.9, 0.9]

        def f(x):
            key = f"{x[0]:.4f},{x[1]:.4f}"
            if key in memo:
                return memo[key]
            state["nEval"] += 1
            r = self._solve(replace(c, w_rich=float(x[0]), q_des=float(x[1])))
            ok = r.success and not (c.limit_discharge and r.T_dis_max > c.T_dis_max)
            o = dict(y=r.COP if ok else -math.inf, r=r)
            memo[key] = o
            return o

        x0 = [min(max(c.w_rich, lo[0]), hi[0]), min(max(c.q_des, lo[1]), hi[1])]
        b0 = f(x0)
        starts = [(x0, b0)]
        if not math.isfinite(b0["y"]):
            for w in (0.35, 0.5, 0.65):
                for q in (0.2, 0.4, 0.6):
                    o = f([w, q])
                    if math.isfinite(o["y"]):
                        starts.append(([w, q], o))
        x, best, converged = compass_search(f, starts, lo, hi, [0.08, 0.08], [0.002, 0.002], [False, False], state,
                                            max_eval, score=lambda o: o["y"], empty={"y": -math.inf})
        if x is None:
            return _fail("no solution composition and degassing gives a feasible design",
                         report=dict(optimiser=dict(nEval=state["nEval"])))
        xr = [js_round(x[0], 3), js_round(x[1], 3)]
        o = f(xr)
        if math.isfinite(o["y"]):
            x, best = xr, o
        final = best["r"]
        final.report["optimiser"] = dict(COP=best["y"], start=b0["y"], w_rich=float(x[0]), q_des=float(x[1]),
                                         nEval=state["nEval"], converged=converged)
        return final

    # ---- model -------------------------------------------------------------------------

    def _solve(self, c: HybridConfig) -> LoopResult:
        warnings: list = []
        if not 0 < c.w_rich < 1 or not 0 < c.q_des < 1:
            return _fail("w_rich and q_des must lie between 0 and 1")
        air = self.air
        try:
            R = self.stream(c.w_rich)
            T_do = c.T_source_in - c.dT_src
            # low pressure: the solution at T_do with q_des boiled off
            des_out = R.TQ(T_do, c.q_des)
            P_L = des_out["P"]
            V, L = self.stream(des_out["w_vap"]), self.stream(des_out["w_liq"])
            vap = V.PQ(P_L, 1.0)                 # separator vapour (dew at P_L)
            liq = L.PQ(P_L, 0.0)                 # lean liquid (bubble at P_L)
            if abs(vap["T"] - T_do) > 0.05 or abs(liq["T"] - T_do) > 0.05:
                raise ValueError("the separator phases are not at the desorber outlet temperature")
            P_max = c.P_max_bar / 10.0
            if P_L * 1.2 >= P_max:
                return _fail(f"the low pressure ({P_L * 10:.1f} bar) leaves no room below {c.P_max_bar:g} bar")

            w_v, w_l = des_out["w_vap"], des_out["w_liq"]

            def compression(P_H, m):
                """Compress the separator vapour (q_des kg) with m kg of lean liquid injected at the suction."""
                if m <= 0:
                    return stage(V, P_L, vap, P_H, None, c.comp_model, c.eta_const), V, 0.0
                n = c.q_des + m
                Cm = self.stream((c.q_des * w_v + m * w_l) / n)
                suc = Cm.PH(P_L, (c.q_des * vap["h"] + m * liq["h"]) / n)
                return stage(Cm, P_L, suc, P_H, None, c.comp_model, c.eta_const), Cm, m

            def at(P_H):
                try:
                    comp, _, m_inj = compression(P_H, 0.0)
                    if c.wet_compression and comp["out"]["T"] > c.T_dis_max:
                        m_hi = 1.0 - c.q_des - 1e-6
                        if compression(P_H, m_hi)[0]["out"]["T"] > c.T_dis_max:
                            raise ValueError("even injecting all the lean liquid cannot hold the discharge at "
                                             f"{c.T_dis_max:g} C")
                        # regula falsi (Illinois) on the discharge temperature, which falls smoothly with m
                        a_m, b_m = 0.0, m_hi
                        ga = comp["out"]["T"] - c.T_dis_max
                        gb = compression(P_H, b_m)[0]["out"]["T"] - c.T_dis_max
                        side = 0
                        best = b_m
                        for _ in range(40):
                            m = b_m - gb * (b_m - a_m) / (gb - ga)
                            gm = compression(P_H, m)[0]["out"]["T"] - c.T_dis_max
                            if gm <= 0:
                                best = m
                            if abs(gm) < 0.02:
                                break
                            if (gm > 0) == (ga > 0):
                                a_m, ga = m, gm
                                if side == -1:
                                    gb /= 2
                                side = -1
                            else:
                                b_m, gb = m, gm
                                if side == 1:
                                    ga /= 2
                                side = 1
                        comp, _, m_inj = compression(P_H, m if abs(gm) < 0.02 else best)
                    m_c = c.q_des + m_inj                                       # through the compressor
                    m_p = 1.0 - m_c                                             # through the pump
                    w_pump = (P_H - P_L) * 1e3 / liq["D"] / c.eta_pump        # kJ/kg (MPa * 1e3 / kg m-3)
                    h_pump = liq["h"] + w_pump
                    h_in = m_c * comp["out"]["h"] + m_p * h_pump                # absorber inlet, per kg solution
                    h_out = cooler_outlet(R, air, P_H, h_in, c.T_air_inlet, c.T_air_tar, c.dT_GC_min, fast=True)
                    gap = min_gap(air_profile(R, air, P_H, h_in, h_out, c.T_air_inlet, c.T_air_tar, 18, fast=True))
                    return dict(P_H=P_H, comp=comp, m_inj=m_inj, m_c=m_c, m_p=m_p, w_pump=w_pump, h_pump=h_pump,
                                h_in=h_in, h_out=h_out, gap=gap)
                except SOLVE_ERRORS as exc:
                    return exc

            lo_P = P_L * 1.2
            prev, found, last = lo_P, None, None
            N = 24
            for i in range(N + 1):
                P = lo_P * (P_max / lo_P) ** (i / N)
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
            for _ in range(20):
                if hi / lo - 1 < 2e-4:
                    break
                m = math.sqrt(lo * hi)
                r = at(m)
                if isinstance(r, dict) and r["gap"] >= c.dT_GC_min - 0.05:
                    hi, S = m, r
                else:
                    lo = m
            P_H = S["P_H"]
            # desorber: throttled absorber outlet heated to the outlet state; it must keep dT_src from the source
            h_di, h_do = S["h_out"], des_out["h"]
            q_src = h_do - h_di
            pts = [[f * q_src, R.T_ph_fast(P_L, h_di + f * q_src),
                    source_T(c.T_source_in, c.T_source_out, f)] for f in lin(0.0, 1.0, 16)]
            gap_src = min(p[2] - p[1] for p in pts)
            if gap_src < c.dT_src - 0.05:
                return _fail(f"the desorber comes within {gap_src:.1f} K of the source water (needs {c.dT_src:g} K): "
                             "raise the degassing or the ammonia fraction, or lower the source outlet")
            w = S["m_c"] * S["comp"]["w"] + S["m_p"] * S["w_pump"]
            q_air = S["h_in"] - S["h_out"]
            wEl = w / c.motor_eff
            cop = q_air / wEl
            T_dis = S["comp"]["out"]["T"]
            if T_dis > c.T_dis_max + 0.05:
                msg = f"compressor discharge at {T_dis:.0f} C, above the {c.T_dis_max:g} C limit"
                if c.limit_discharge:
                    return _fail(msg)
                warnings.append(msg)
            return self._result(c, S, des_out, vap, liq, P_L, P_H, q_air, w, wEl, cop, q_src, gap_src, pts, warnings)
        except SOLVE_ERRORS as exc:
            return _fail(str(exc), warnings=warnings)

    def _result(self, c, S, des_out, vap, liq, P_L, P_H, q_air, w, wEl, cop, q_src, gap_src, src_pts, warnings):
        R, air = self.stream(c.w_rich), self.air
        st = lambda s: dict(T=s["T"], P=s["P"], h=s["h"], s=s["s"])
        states = dict(desorber_out=st(des_out), vapour=st(vap), lean_liquid=st(liq), discharge=st(S["comp"]["out"]),
                      absorber_in=st(R.PH(P_H, S["h_in"])), absorber_out=st(R.PH(P_H, S["h_out"])),
                      desorber_in=st(R.PH(P_L, S["h_out"])))
        Tb, Td, _, _ = R.bubble_dew(P_H)
        prof = air_profile(R, air, P_H, S["h_in"], S["h_out"], c.T_air_inlet, c.T_air_tar, 18, fast=True)
        exchangers = [
            dict(name="Absorber to air", hot=[[p[0], p[1]] for p in prof], cold=[[p[0], p[2]] for p in prof],
                 other="dryer air", gap=S["gap"]),
            dict(name="Desorber", hot=[[p[0], p[2]] for p in src_pts], cold=[[p[0], p[1]] for p in src_pts],
                 other="source water", gap=gap_src)]
        vs = S["m_c"] / R.PH(P_L, S["comp"]["out"]["h"] - S["comp"]["w"])["D"] if S["m_inj"] > 0 else c.q_des / vap["D"]
        report = dict(COP=cop, VHC=q_air / vs / 1000, w_rich=c.w_rich, q_des=c.q_des, w_vapour=des_out["w_vap"],
                      w_lean=des_out["w_liq"], P_low_bar=P_L * 10, P_high_bar=P_H * 10, PR=P_H / P_L,
                      T_des_out=des_out["T"], T_des_in=states["desorber_in"]["T"], T_abs_in=states["absorber_in"]["T"],
                      T_abs_out=states["absorber_out"]["T"], T_bubble_high=Tb - 273.15, T_dew_high=Td - 273.15,
                      T_dis=S["comp"]["out"]["T"], eta_comp=S["comp"]["eta"], w_comp=S["m_c"] * S["comp"]["w"],
                      w_pump=S["m_p"] * S["w_pump"], m_injected=S["m_inj"], wet_compression=S["m_inj"] > 0,
                      w_el_kJkg=wEl, q_air=q_air, q_source=q_src,
                      gap_air=S["gap"], gap_source=gap_src,
                      air_per_solution=q_air / (air.h(c.T_air_tar) - air.h(c.T_air_inlet)), air_W=air.W,
                      balance=q_air - w - q_src, properties="IAPWS G4-01 (Tillner-Roth and Friend) via teqp")
        return LoopResult(True, "success", COP=cop, VHC=report["VHC"], Pmax_bar=P_H * 10, T_dis_max=S["comp"]["out"]["T"],
                          q_air=q_air, w_shaft=w, q_source=q_src, report=report, states=states,
                          exchangers=exchangers, warnings=warnings)
