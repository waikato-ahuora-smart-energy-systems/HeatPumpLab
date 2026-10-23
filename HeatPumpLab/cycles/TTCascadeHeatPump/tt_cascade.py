"""Transcritical-transcritical CO2 cascade heat pump for air drying, after Kong (2026, PhD thesis, University of
Waikato, https://hdl.handle.net/10289/18056), cycle TT-4 with the additional gas cooler and secondary air (Fig. 58);
see also Kong et al., Applied Thermal Engineering 238 (2024) 122005, and, for expansion work recovery, Kong et al.,
Energies 19 (2026) 3767. A port of ``DRYCASC.kong`` and
``kongOptimise`` from the "Advanced dryer heat pumps" web app.

Flow path::

    bottom: evaporator -> IHX (cold) -> compressor -> cascade gas cooler (CGC, hot side: its only heat sink)
            -> IHX (hot) -> expansion valve (or expander) -> evaporator
    top:    evaporator -> IHX (cold) -> CGC (cold side: the bottom heat superheats the top suction) -> compressor
            -> HT gas cooler (primary air) -> MT gas cooler (primary + secondary air) -> IHX (hot) -> valve

Air: primary and secondary air pass the MT gas cooler together (inlet -> T_mt); the secondary leaves at T_mt and
the primary goes on through the HT gas cooler to the dryer target. The two gas coolers are one counterflow exchanger
against that composite air curve; the top gas-cooler outlet is the coldest that keeps the approaches (dT_GC_min in
the MT gas cooler, dT_ht in the HT one), raised by gc_off when that pays through the top IHX. The air flow follows
from the duty. air_ratio = secondary / primary air mass flow (0: no secondary air); r = bottom / top refrigerant mass
flow. The bottom cycle may run subcritical.

Closure. For a top gas-cooler outlet temperature, the top IHX duty fixes the top CGC cold inlet; the bottom cycle is
then solved for its CGC outlet temperature (its own IHX raises its discharge and so the CGC outlet again: the lowest
self-consistent outlet is taken), the CGC passes the largest duty that keeps dT_cgc, and the top compression and
gas-cooler outlet follow. The outlet is iterated to 0.01 K. Each IHX passes the largest duty both its ends allow
(approach dT_ihx at each end). Both evaporators sit on the source water in series (top first); their approaches
dT_ev to the water leaving each are checked (warnings).

Kong's reference (Table 9, section 6.2.4): eta_s 0.70, 150 bar, 20 K CO2-air and 5 K CO2-CO2 approaches,
evaporating 10 / 20 C, a third of the air to the fluidised bed (air_ratio 0.5): COP 2.55. The defaults here are the
web app's (170 C dryer air, no secondary air, 60 % RH inlet air, source water 40 -> 30 C, Pareto-front compressor).

Duties are per kg of top-cycle refrigerant. Units: MPa, C, kJ/kg.
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
    air_profile,
    comp_path,
    compass_search,
    dry_air,
    expander_electric,
    expander_warning,
    expand_to,
    humid_air,
    iso_path,
    js_num,
    js_round,
    lin,
    max_duty,
    min_gap,
    plain_state,
    source_T,
    stage,
)


@dataclass
class TTCascadeConfig:
    """Design inputs (the web app's defaults). Temperatures in C, pressures in MPa, approaches in K."""

    bottom_fluid: str = "CO2"
    top_fluid: str = "CO2"
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0         # dryer (primary) air target
    dT_GC_min: float = 20.0          # gas-cooler approach to the air (MT gas cooler; HT too unless dT_ht is set)
    dT_ht: float | None = None       # HT gas-cooler approach (None: dT_GC_min)
    T_source_in: float = 40.0        # source water, through the top then the bottom evaporator
    T_source_out: float = 30.0
    T_ev_b: float = 10.0             # bottom evaporating temperature (dew)
    T_ev_t: float = 20.0             # top evaporating temperature (dew)
    P_b: float = 15.0                # bottom discharge pressure, MPa
    P_t: float = 15.0                # top discharge pressure, MPa
    r: float = 0.43                  # bottom / top refrigerant mass flow
    air_ratio: float = 0.0           # secondary / primary air mass flow
    T_mt: float | None = None        # secondary air target (None: midway between the air inlet and target)
    dT_cgc: float = 5.0              # cascade gas-cooler approach
    dT_ihx_b: float = 5.0
    dT_ihx_t: float = 5.0
    dT_ev: float = 5.0               # evaporator approach to the water leaving it (checked)
    gc_off: float = 0.0              # top gas-cooler outlet raised above its coldest feasible value, K
    eta: float | None = None         # fixed isentropic efficiency of both compressors (None: comp_model)
    comp_model: int = 5              # 0 ideal, 1 constant eta_const, 2 polynomial, 3 Brendel, 5 Pareto front
    eta_const: float = 0.73
    motor_eff: float = DEFAULT_MOTOR_EFF
    expander_eta: float | None = None    # expanders in place of both valves (isentropic efficiency), None: valves
    air_RH: float = 60.0             # inlet air relative humidity, % (0: dry air through the property package)
    P_air_kPa: float = 101.325
    P_max: float = 15.0              # optimiser bound on both discharge pressures, MPa
    property_backend: str | None = None  # "refprop" | "coolprop"; None -> $HEATPUMPLAB_PROPERTY_BACKEND


TTCascadeResult = LoopResult


def _fail(message: str, **extra) -> LoopResult:
    return LoopResult(False, message, **extra)


class TTCascadeHeatPump:
    """TT CO2 cascade (Kong TT-4); see the module docstring."""

    def __init__(self, config: TTCascadeConfig | None = None, bottom: Fluid | None = None,
                 top: Fluid | None = None, air: AirModel | None = None):
        self.c = config if config is not None else TTCascadeConfig()
        c = self.c
        self.bottom = bottom if bottom is not None else Fluid(c.bottom_fluid, c.property_backend)
        if top is not None:
            self.top = top
        elif c.top_fluid == c.bottom_fluid:
            self.top = self.bottom
        else:
            self.top = Fluid(c.top_fluid, c.property_backend)
        if air is None:
            air = humid_air(c.T_air_inlet, c.air_RH, c.P_air_kPa) if c.air_RH > 0 else \
                dry_air(c.P_air_kPa, c.property_backend)
        self.air = air

    # ---- public ------------------------------------------------------------------------

    def solve(self, **overrides) -> LoopResult:
        """The design at the configuration (fields can be overridden by keyword)."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._kong(c))

    def optimise(self, max_eval: int = 500, **overrides) -> LoopResult:
        """Best-COP bottom and top discharge pressures, mass-flow ratio and gas-cooler outlet raise (the app's
        kongOptimise): a coarse grid and a compass search from the configured design and from the best grid point,
        then rounded as the app shows them (0.01 MPa, 0.001 in r, 0.1 K), r stepped up until the rounded design solves.
        The result is the full solution there, with the search in report["optimiser"]."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._optimise(c, max_eval))

    @staticmethod
    def _cached(fn):
        if active_property_call_cache() is None:
            with property_call_cache():
                return fn()
        return fn()

    # ---- model -------------------------------------------------------------------------

    def _kong(self, c: TTCascadeConfig, lite: bool = False) -> LoopResult:
        fb, ft, air = self.bottom.fast_view(), self.top.fast_view(), self.air
        warnings: list = []
        eta = c.eta or None

        def stg(fl, P_in, st, P_out):
            return stage(fl, P_in, st, P_out, eta, c.comp_model, c.eta_const)

        try:
            r, Tain, phi = c.r, c.T_air_inlet, max(0.0, c.air_ratio or 0.0)
            T_mt = c.T_mt if c.T_mt is not None else (Tain + c.T_air_tar) / 2
            if not r > 0:
                return _fail("the bottom/top mass flow must be positive")
            if not (Tain < T_mt < c.T_air_tar):
                return _fail(f"the secondary air target must lie between the air inlet ({js_num(Tain)} C) and the dryer "
                             f"inlet target ({js_num(c.T_air_tar)} C)")
            # composite air curve per kg primary air: [0, H1] both streams inlet -> T_mt, [H1, H] primary only
            ha, hm, ho = air.h(Tain), air.h(T_mt), air.h(c.T_air_tar)
            H1, H2 = (1 + phi) * (hm - ha), ho - hm
            H = H1 + H2
            f1 = H1 / H

            def T_comp(x):
                return air.T(ha + x / (1 + phi)) if x <= H1 else air.T(hm + (x - H1))

            eb, et = fb.TQ(c.T_ev_b, 1.0), ft.TQ(c.T_ev_t, 1.0)
            Peb, Pet = eb["P"], et["P"]
            if c.P_t <= Pet or c.P_b <= Peb:
                return _fail("discharge pressure below the evaporating pressure")

            def ihx(pr, Pc, cin, Ph, Thot, dT):
                """cold side from the evaporator outlet, hot side from Thot (same flow): the largest duty both ends allow"""
                if Thot - dT <= cin["T"]:
                    return 0.0
                return max(0.0, min(pr.PT(Ph, Thot)["h"] - pr.PT(Ph, cin["T"] + dT)["h"],
                                    pr.PT(Pc, Thot - dT)["h"] - cin["h"]))

            def bottom_at(ctIn, Tb):
                qib = ihx(fb, Peb, eb, c.P_b, Tb, c.dT_ihx_b)
                cbIn = fb.PH(Peb, eb["h"] + qib)
                cb = stg(fb, Peb, cbIn, c.P_b)
                cg = {"q": max_duty(fb, c.P_b, cb["out"]["h"], r, ft, Pet, ctIn["h"], c.dT_cgc)}
                hbo = cb["out"]["h"] - cg["q"] / r
                return dict(qib=qib, cbIn=cbIn, cb=cb, cg=cg, hbo=hbo, f=fb.Tph(c.P_b, hbo) - Tb)

            def solve_bottom(ctIn):
                """the lowest self-consistent bottom CGC outlet (scan up from no IHX, then bisect); None: the bottom
                cycle brings more heat than the top suction can take"""
                lo = eb["T"] + c.dT_ihx_b
                A = bottom_at(ctIn, lo)
                if A["f"] <= 0:
                    return A
                hi = None
                T = lo + 5
                while T <= 400:
                    B = bottom_at(ctIn, T)
                    if B["f"] <= 0:
                        hi = T
                        break
                    lo, A = T, B
                    T += 5
                if hi is None:
                    return None
                k = 0
                while k < 30 and hi - lo > 0.005:
                    m = (lo + hi) / 2
                    Mm = bottom_at(ctIn, m)
                    if Mm["f"] > 0:
                        lo, A = m, Mm
                    else:
                        hi = m
                    k += 1
                return bottom_at(ctIn, hi)

            dTht = c.dT_ht if c.dT_ht is not None else c.dT_GC_min

            def composite_outlet(P, hIn, dT, dT2, n=32):
                """coldest gas-cooler outlet against the composite air curve that keeps the approaches"""
                def gap(h):
                    m = math.inf
                    for i in range(n + 1):
                        f = i / n
                        m = min(m, ft.Tph(P, h + f * (hIn - h)) - T_comp(f * H) - (dT if f <= f1 else dT2))
                    if f1 < 1:
                        m = min(m, ft.Tph(P, h + f1 * (hIn - h)) - T_comp(H1) - min(dT, dT2))
                    return m
                try:
                    hLow = ft.PT(P, Tain + dT)["h"]
                except Exception:
                    hLow = hIn - 1
                if gap(hLow) >= -1e-6:
                    return dict(h=hLow, ok=True)
                lo, hi = hLow, hIn - 1e-6 * abs(hIn)
                if gap(hi) < 0:
                    return dict(h=hi, ok=False)
                N = 40
                for i in range(1, N + 1):
                    h = hLow + (hi - hLow) * i / N
                    if gap(h) >= 0:
                        lo, hi = hLow + (hi - hLow) * (i - 1) / N, h
                        break
                k = 0
                while k < 40 and hi - lo > 0.01:
                    m = (lo + hi) / 2
                    if gap(m) >= 0:
                        hi = m
                    else:
                        lo = m
                    k += 1
                return dict(h=hi, ok=True)

            def gc_out(ctp):
                g = composite_outlet(c.P_t, ctp["out"]["h"], c.dT_GC_min, dTht)
                if not g["ok"] or not c.gc_off > 0:
                    return g
                T = ft.Tph(c.P_t, g["h"]) + c.gc_off
                if T >= ctp["out"]["T"] - 1:
                    return dict(h=g["h"], ok=False)
                return dict(h=ft.PT(c.P_t, T)["h"], ok=True)

            # outer: the top gas-cooler outlet sets the top IHX duty
            Tg, S, conv = Tain + c.dT_GC_min, None, False
            for it in range(40):
                qit = ihx(ft, Pet, et, c.P_t, Tg, c.dT_ihx_t)
                ctIn = ft.PH(Pet, et["h"] + qit)
                B = solve_bottom(ctIn)
                if B is None:
                    return _fail("the bottom cycle brings more heat than the cascade gas cooler can pass to the top "
                                 "suction (lower the bottom/top mass flow or the bottom pressure)")
                ts = ft.PH(Pet, ctIn["h"] + B["cg"]["q"])
                ctp = stg(ft, Pet, ts, c.P_t)
                go = gc_out(ctp)
                if not go["ok"]:
                    return _fail(f"the top discharge ({ctp['out']['T']:.0f} C) cannot heat the primary air to "
                                 f"{js_num(c.T_air_tar)} C with a {js_num(dTht)} K approach (raise the top pressure or the bottom/top "
                                 "mass flow)")
                nTg = ft.Tph(c.P_t, go["h"])
                S = dict(ctIn=ctIn, cbIn=B["cbIn"], cb=B["cb"], cg=B["cg"], ts=ts, hbo=B["hbo"], ctp=ctp, go=go)
                done = abs(nTg - Tg) < 0.01
                Tg += (1 if it < 3 else 0.5) * (nTg - Tg)
                if done:
                    conv = True
                    break
            if not conv:
                warnings.append("the top gas-cooler outlet did not settle to 0.01 K")
            ctIn, cbIn, cb, cg, ts, hbo, ctp, go = (S[k] for k in ("ctIn", "cbIn", "cb", "cg", "ts", "hbo", "ctp", "go"))
            if not fb.is_gas(cbIn) or not ft.is_gas(ts):
                return _fail("suction not gas")
            # IHX with the converged outlets
            hOut = go["h"]
            Tout = ft.Tph(c.P_t, hOut)
            qit = ihx(ft, Pet, et, c.P_t, Tout, c.dT_ihx_t)
            qib = ihx(fb, Peb, eb, c.P_b, fb.Tph(c.P_b, hbo), c.dT_ihx_b)
            hvt, hvb = hOut - qit, hbo - qib
            # valves, or expanders: the evaporator inlets
            het = expand_to(self.top, c.P_t, hvt, Pet, c.expander_eta)
            heb = expand_to(self.bottom, c.P_b, hvb, Peb, c.expander_eta)
            wx = (hvt - het) + r * (hvb - heb)
            qAir = ctp["out"]["h"] - hOut
            mP = qAir / H
            hMid = hOut + mP * H1
            qMT, qHT = mP * H1, mP * H2
            qPri = mP * (air.h(c.T_air_tar) - air.h(Tain))
            qSec = phi * mP * (air.h(T_mt) - air.h(Tain))
            w = ctp["w"] + r * cb["w"]
            wEl = expander_electric(w, wx, c.motor_eff)
            cop = qAir / wEl
            qEt, qEb = et["h"] - het, r * (eb["h"] - heb)
            xq = [expander_warning(self.top, Pet, het, warnings, "top expander"),
                  expander_warning(self.bottom, Peb, heb, warnings, "bottom expander")] if c.expander_eta else None
            if cg["q"] <= 1e-6:
                warnings.append("the cascade gas cooler carries no heat: the bottom discharge is too cold for the top "
                                "suction")
            if max(ctp["out"]["T"], cb["out"]["T"]) > 320:
                warnings.append(f"discharge at {max(ctp['out']['T'], cb['out']['T']):.0f} C: the bottom cycle brings "
                                "more heat than the top suction can usefully take (lower the bottom/top mass flow)")
            # source water in series: the top evaporator first, then the bottom one
            Tw_out = c.T_source_out if c.T_source_out is not None else c.T_source_in
            fMid = qEb / (qEt + qEb)
            Tw_mid = source_T(c.T_source_in, c.T_source_out, fMid)
            gapEt, gapEb = Tw_mid - c.T_ev_t, Tw_out - c.T_ev_b
            if gapEt < c.dT_ev - 1e-6:
                warnings.append(f"top evaporator: the source water leaves it at {Tw_mid:.1f} C, {gapEt:.1f} K above "
                                f"the {js_num(c.T_ev_t)} C evaporating temperature (approach {js_num(c.dT_ev)} K); lower the top "
                                "evaporating temperature or raise the source outlet")
            if gapEb < c.dT_ev - 1e-6:
                warnings.append(f"bottom evaporator: the source water leaves at {Tw_out:.1f} C, {gapEb:.1f} K above "
                                f"the {js_num(c.T_ev_b)} C evaporating temperature (approach {js_num(c.dT_ev)} K); lower the bottom "
                                "evaporating temperature or raise the source outlet")
            pHT = air_profile(ft, air, c.P_t, ctp["out"]["h"], hMid, T_mt, c.T_air_tar, 16)
            pMT = air_profile(ft, air, c.P_t, hMid, hOut, Tain, T_mt, 16)
            gap_HT, gap_MT = min_gap(pHT), min_gap(pMT)
            if lite:
                return LoopResult(True, "success", COP=cop, warnings=warnings,
                                  report=dict(r=r, air_ratio=phi, gap_HT=gap_HT, gap_MT=gap_MT, T0=[Tg]))
            vt, vb = 1 / ts["D"], 1 / cbIn["D"]
            top, bot = self.top, self.bottom
            states = dict(
                top=dict(evap_out=plain_state(top, Pet, et["h"]), ihx_cold_out=plain_state(top, Pet, ctIn["h"]),
                         suction=plain_state(top, Pet, ts["h"]), discharge=plain_state(top, c.P_t, ctp["out"]["h"]),
                         htgc_out=plain_state(top, c.P_t, hMid), gc_out=plain_state(top, c.P_t, hOut),
                         valve_in=plain_state(top, c.P_t, hvt), valve_out=plain_state(top, Pet, het)),
                bot=dict(evap_out=plain_state(bot, Peb, eb["h"]), suction=plain_state(bot, Peb, cbIn["h"]),
                         discharge=plain_state(bot, c.P_b, cb["out"]["h"]), cgc_out=plain_state(bot, c.P_b, hbo),
                         valve_in=plain_state(bot, c.P_b, hvb), valve_out=plain_state(bot, Peb, heb)),
                r=r, air_primary=mP, air_secondary=phi * mP,
            )
            exp_name = "Expander" if c.expander_eta else "Expansion valve"
            ihx_top = [dict(name="Internal heat exchanger", pts=iso_path(top, Pet, et["h"], ctIn["h"]))] if qit > 0 else []
            ihx_top_hot = [dict(name="Internal heat exchanger", pts=iso_path(top, c.P_t, hOut, hvt))] if qit > 0 else []
            ihx_bot = [dict(name="Internal heat exchanger", pts=iso_path(bot, Peb, eb["h"], cbIn["h"]))] if qib > 0 else []
            ihx_bot_hot = [dict(name="Internal heat exchanger", pts=iso_path(bot, c.P_b, hbo, hvb))] if qib > 0 else []
            paths = {
                "Top cycle": [
                    dict(name=exp_name, pts=comp_path(top, c.P_t, hvt, Pet, het)),
                    dict(name="Evaporator", pts=iso_path(top, Pet, het, et["h"])),
                    *ihx_top,
                    dict(name="Cascade gas cooler (cold)", pts=iso_path(top, Pet, ctIn["h"], ts["h"])),
                    dict(name="Compressor", pts=comp_path(top, Pet, ts["h"], c.P_t, ctp["out"]["h"])),
                    dict(name="HT gas cooler", pts=iso_path(top, c.P_t, ctp["out"]["h"], hMid)),
                    dict(name="MT gas cooler", pts=iso_path(top, c.P_t, hMid, hOut)),
                    *ihx_top_hot,
                ],
                "Bottom cycle": [
                    dict(name=exp_name, pts=comp_path(bot, c.P_b, hvb, Peb, heb)),
                    dict(name="Evaporator", pts=iso_path(bot, Peb, heb, eb["h"])),
                    *ihx_bot,
                    dict(name="Compressor", pts=comp_path(bot, Peb, cbIn["h"], c.P_b, cb["out"]["h"])),
                    dict(name="Cascade gas cooler (hot)", pts=iso_path(bot, c.P_b, cb["out"]["h"], hbo)),
                    *ihx_bot_hot,
                ],
            }
            # cascade gas cooler: both sides per kg top flow from the cold end
            Qc = ts["h"] - ctIn["h"]
            cgc_hot = [[f * Qc, bot.PH(c.P_b, hbo + f * (cb["out"]["h"] - hbo))["T"]] for f in lin(0, 1, 14)]
            cgc_cold = [[f * Qc, top.PH(Pet, ctIn["h"] + f * Qc)["T"]] for f in lin(0, 1, 14)]
            water = lambda Q, q0, fa, fb_: [[q0, source_T(c.T_source_in, c.T_source_out, fa)],
                                            [q0 + Q, source_T(c.T_source_in, c.T_source_out, fb_)]]
            exchangers = [
                dict(name="HT gas cooler", hot=[[qMT + p[0], p[1]] for p in pHT], cold=[[qMT + p[0], p[2]] for p in pHT],
                     other="primary air", gap=gap_HT),
                dict(name="MT gas cooler", hot=[[p[0], p[1]] for p in pMT], cold=[[p[0], p[2]] for p in pMT],
                     other="primary + secondary air", gap=gap_MT),
                dict(name="Cascade gas cooler", hot=cgc_hot, cold=cgc_cold, other="top-cycle suction"),
                dict(name="Top evaporator", hot=water(qEt, qEb, fMid, 1.0),
                     cold=[[qEb + (p[2] - het), p[0]] for p in iso_path(top, Pet, het, et["h"])], other="source water"),
                dict(name="Bottom evaporator", hot=water(qEb, 0.0, 0.0, fMid),
                     cold=[[r * (p[2] - heb), p[0]] for p in iso_path(bot, Peb, heb, eb["h"])], other="source water"),
            ]
            VHC = qAir / (vt + r * vb) / 1000
            report = dict(
                COP=cop, VHC=VHC, P_top_bar=c.P_t * 10, P_bot_bar=c.P_b * 10, P_evap_top_bar=Pet * 10,
                P_evap_bot_bar=Peb * 10, PR_top=c.P_t / Pet, PR_bot=c.P_b / Peb, eta_top=ctp["eta"], eta_bot=cb["eta"],
                T_dis_top=ctp["out"]["T"], T_dis_bot=cb["out"]["T"], T_suc_top=ts["T"], q_cgc=cg["q"], q_HT=qHT,
                q_MT=qMT, q_primary=qPri, q_secondary=qSec, share_sec=qSec / qAir, q_ihx_top=qit, q_ihx_bot=r * qib,
                gap_HT=gap_HT, gap_MT=gap_MT, T_gc_mid=ft.Tph(c.P_t, hMid), T_gc_out=Tout,
                T_bot_cgc_out=fb.Tph(c.P_b, hbo), T_water_mid=Tw_mid, gap_evap_top=gapEt, gap_evap_bot=gapEb,
                q_evap_top=qEt, q_evap_bot=qEb, gc_off=c.gc_off or 0.0, T_bot_valve=fb.Tph(c.P_b, hvb),
                primary_per_top=mP, bottom_supercritical=(c.P_b > bot.P_crit) if math.isfinite(bot.P_crit) else None,
                w_el_kJkg=wEl, r=r, air_ratio=phi, T_mt=T_mt, expander=bool(c.expander_eta), w_exp_kJkg=wx,
                exp_out_quality=xq, air_W=air.W,
                # basis for sizing: the bottom compressor (as the app): suction density, flow per kg top, PR
                ref=dict(rho=cbIn["D"], flow=r, PR=c.P_b / Peb),
            )
            return LoopResult(True, "success", COP=cop, VHC=VHC, Pmax_bar=max(c.P_t, c.P_b) * 10,
                              T_dis_max=max(ctp["out"]["T"], cb["out"]["T"]), q_air=qAir, w_shaft=w,
                              q_source=qEt + qEb, report=report, states=states, paths=paths, exchangers=exchangers,
                              warnings=warnings)
        except Exception as exc:  # property failures and numerical edge cases, as the app reports them
            return _fail(f"property calculation failed: {exc}", warnings=warnings)

    # ---- optimiser ---------------------------------------------------------------------

    def _optimise(self, c: TTCascadeConfig, max_eval: int) -> LoopResult:
        Pmax = c.P_max
        Peb = self.bottom.TQ(c.T_ev_b, 1.0)["P"]
        Pet = self.top.TQ(c.T_ev_t, 1.0)["P"]
        lo = [Peb * 1.05, Pet * 1.05, math.log(0.05), 0.0]
        hi = [Pmax, Pmax, math.log(5.0), 80.0]

        def clip(x):
            return [min(max(v, lo[i]), hi[i]) for i, v in enumerate(x)]

        state = dict(nEval=0)
        memo: dict = {}

        def at(x):
            return replace(c, P_b=x[0], P_t=x[1], r=math.exp(x[2]), gc_off=x[3])

        def f(x):
            k = ",".join(f"{v:.6f}" for v in x)
            if k in memo:
                return memo[k]
            state["nEval"] += 1
            y, res = -math.inf, None
            try:
                res = self._kong(at(x), lite=True)
                if res.success and math.isfinite(res.COP):
                    y = res.COP
            except Exception:
                pass
            o = dict(y=y, r=res)
            memo[k] = o
            return o

        x0 = clip([c.P_b, c.P_t, math.log(c.r), c.gc_off or 0.0])
        b0 = f(x0)
        start = b0["y"]
        # coarse grid: the COP surface has several local optima (bottom cycle sub- or transcritical, the HT approach
        # edge), so the compass search starts from the configured design and from the best grid point
        xg, bg = None, dict(y=-math.inf)
        for a in (0.3, 0.5, 0.7, 0.9):
            for b in (0.6, 0.8, 1.0):
                for r0 in (0.2, 0.35, 0.5, 0.8):
                    z = clip([lo[0] + a * (hi[0] - lo[0]), lo[1] + b * (hi[1] - lo[1]), math.log(r0), 0.0])
                    o = f(z)
                    if o["y"] > bg["y"]:
                        bg, xg = o, z
        x, best, converged = compass_search(f, [(x0, b0), (xg, bg)], lo, hi, [1.0, 1.0, 0.25, 8.0],
                                            [0.002, 0.002, 0.002, 0.05], [False] * 4, state, max_eval,
                                            score=lambda o: o["y"], empty={"y": -math.inf})
        if x is None:
            return _fail("no feasible design found within the bounds", report=dict(optimiser=dict(nEval=state["nEval"])))
        # the optimum often sits on the edge of feasibility (the HT gas-cooler approach): round as the app shows the
        # values and step r up until the rounded design solves
        P_b, P_t = js_round(x[0], 2), js_round(x[1], 2)
        rr, off = js_round(math.exp(x[2]), 3), js_round(x[3], 1)
        chk = None
        for _ in range(60):
            chk = self._kong(replace(c, P_b=P_b, P_t=P_t, r=rr, gc_off=off))
            if chk.success:
                break
            rr = js_round(rr + 0.001, 3)
        if chk is not None and chk.success:
            x = [P_b, P_t, math.log(rr), off]
            final = chk
        else:
            final = self._kong(at(x))

        def near(v, b):
            return abs(v - b) < 1e-3

        final.report["optimiser"] = dict(
            COP=final.COP, start=start, P_b=x[0], P_t=x[1], r=math.exp(x[2]), gc_off=x[3], nEval=state["nEval"],
            converged=converged, best_search_COP=best["y"],
            atBound=dict(P_b=near(x[0], hi[0]) or near(x[0], lo[0]), P_t=near(x[1], hi[1]) or near(x[1], lo[1]),
                         r=near(x[2], lo[2]) or near(x[2], hi[2])))
        return final

