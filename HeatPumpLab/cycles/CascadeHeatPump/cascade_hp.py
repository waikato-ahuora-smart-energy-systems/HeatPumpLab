"""n-level cascade heat pump for air drying, after Hoang et al.: 2 to 4 vapour-compression levels, level 1 evaporating
on the source water and level k+1 in level k's cascade condenser. A port of ``DRYCASC.cascade`` and
``cascadeOptimise`` from the "Advanced dryer heat pumps" web app.

Flow path of level k::

    evaporator -> [IHX cold side] -> compressor -> condenser (top level: desuperheats, condenses and subcools into
    the air; lower levels: a fraction f_k of the desuperheating and condensing heat to the air, the rest to the
    cascade condenser) -> [subcooler to the air] -> [IHX hot side] -> expansion valve (or expander)

Air passes the levels in series, 1 to n, each level's subcooler then its condenser. A lower level with cond_air on
condenses into the air until the air leaves it at its condensing temperature less the approach (that fixes f_k); with
sub_air on, its liquid is subcooled into the air to the air inlet + approach (but no lower than its IHX needs, and
capped so the air leaving its subcooler stays an approach below its condensing temperature). The lower levels'
condensing temperatures T_cond are inputs; the top level's is the lowest (scan in 5 K steps, then bisect to 0.02 K)
whose air approach holds. Level k+1 evaporates dT_casc below level k's condensing temperature; level 1 evaporates
dT_src below the source water outlet.

Suction superheat: at least ``superheat`` (K), and enough for a dry discharge 5 K above the condensing dew point;
any superheat above ``superheat`` is taken from the level's liquid (an internal heat exchanger). With ihx on for a
level, the suction is heated up to the liquid entering the IHX less dT_ihx.

Duties are per kg of level-1 refrigerant; each level's flow is in report["levels"]. Units: MPa, C, kJ/kg.

Defaults are the web app's: isobutane / cyclopentane, level 1 condensing at 100 C, 170 C dryer air from 15 C
(60 % RH), source water 40 -> 30 C, 20 K air and 5 K cascade approaches, Pareto-front compressors.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace

from ...utils.property_cache import active_property_call_cache, property_call_cache
from ..common.dryer_loops import (
    DEFAULT_MOTOR_EFF,
    AirModel,
    LoopFluid,
    LoopResult,
    air_profile,
    air_secondary,
    comp_path,
    compass_search,
    dry_air,
    exp_el,
    exp_warn,
    expand_to,
    humid_air,
    iso_path,
    js_num,
    js_round,
    lin,
    min_gap,
    source_T,
    stage,
)

MAX_LEVELS = 4


@dataclass
class CascadeConfig:
    """Design inputs (the web app's defaults). Temperatures in C, approaches in K. Lists are per level, bottom first;
    entries past the number of levels are ignored."""

    fluids: list = field(default_factory=lambda: ["Isobutane", "Cyclopentane", "Cyclopentane", "Cyclopentane"])
    n: int = 2                        # levels, 2-4
    T_cond: list = field(default_factory=lambda: [100.0, 150.0, 180.0])   # condensing (bubble) T of levels 1..n-1
    cond_air: list = field(default_factory=lambda: [True, True, True])    # lower level condenses partly into the air
    sub_air: list = field(default_factory=lambda: [True, True, True, True])   # level subcools its liquid into the air
    ihx: list = field(default_factory=lambda: [False, False, False, False])  # IHX on the level
    dT_ihx: float = 5.0
    dT_casc: float = 5.0              # cascade condenser: level k+1 evaporates this far below level k's condensing T
    superheat: float = 5.0            # minimum suction superheat
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0
    dT_GC_min: float = 20.0           # air approach
    T_source_in: float = 40.0
    T_source_out: float = 30.0
    dT_src: float = 5.0               # level-1 evaporating temperature below the source water outlet
    air_ratio: float = 0.0            # secondary / primary air (heated with the primary from the inlet to T_mt)
    T_mt: float = 102.0               # secondary air target
    comp_model: int = 5               # 0 ideal, 1 constant eta_const, 2 polynomial, 3 Brendel, 4/5 Pareto front v1/v2
    eta_const: float = 0.73
    motor_eff: float = DEFAULT_MOTOR_EFF
    expander_eta: float | None = None     # expanders in place of the valves (isentropic efficiency), None: valves
    air_RH: float = 60.0              # inlet air relative humidity, % (0: dry air through the property package)
    P_air_kPa: float = 101.325
    property_backend: str | None = None   # "refprop" | "coolprop"; None -> $HEATPUMPLAB_PROPERTY_BACKEND


CascadeResult = LoopResult


def _fail(message: str, **extra) -> LoopResult:
    return LoopResult(False, message, **extra)


class CascadeHeatPump:
    """n-level cascade (Hoang et al.); see the module docstring."""

    def __init__(self, config: CascadeConfig | None = None, fluids: list | None = None, air: AirModel | None = None):
        self.c = config if config is not None else CascadeConfig()
        c = self.c
        if fluids is None:
            made: dict = {}
            fluids = []
            for name in c.fluids[:max(c.n, 0)]:
                if name not in made:
                    made[name] = LoopFluid(name, c.property_backend)
                fluids.append(made[name])
        self.fluids = fluids
        if air is None:
            air = humid_air(c.T_air_inlet, c.air_RH, c.P_air_kPa) if c.air_RH > 0 else \
                dry_air(c.P_air_kPa, c.property_backend)
            air = air_secondary(air, c.T_air_inlet, c.T_mt, c.air_ratio, c.T_air_tar)
        self.air = air

    # ---- public ------------------------------------------------------------------------

    def solve(self, **overrides) -> LoopResult:
        """The design at the configuration (fields can be overridden by keyword, except the fluids and air)."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._cascade(c))

    def optimise(self, max_eval: int = 400, **overrides) -> LoopResult:
        """Best-COP condensing temperatures of the lower levels (the app's cascadeOptimise): a grid of evenly spread
        lifts and a compass search from it and from the configured design; a design counts only if it keeps every air
        approach and every IHX closes. Rounded to 0.1 K (nudged up until the rounded design still counts). The
        result is the full solution there, with the search in report["optimiser"]."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._optimise(c, max_eval))

    @staticmethod
    def _cached(fn):
        if active_property_call_cache() is None:
            with property_call_cache():
                return fn()
        return fn()

    # ---- model -------------------------------------------------------------------------

    def _cascade(self, c: CascadeConfig) -> LoopResult:
        n, fls, air, warnings = c.n, self.fluids, self.air, []
        if not 2 <= n <= MAX_LEVELS:
            return _fail(f"a cascade has 2 to {MAX_LEVELS} levels")
        if len(fls) < n:
            return _fail("a refrigerant is needed for every level")

        def stg(fl, P_in, st, P_out):
            return stage(fl, P_in, st, P_out, None, c.comp_model, c.eta_const)

        def flag(lst, k):
            return bool(lst[k]) if k < len(lst) else False

        try:
            Tain = c.T_air_inlet
            ev1 = (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_src

            def level(pr, Te, Tc, sh, Tihx=None):
                ev = pr.TQ(Te, 1.0)
                Pe = ev["P"]
                bub, dew = pr.TQ(Tc, 0.0), pr.TQ(Tc, 1.0)
                Pc = max(bub["P"], dew["P"])
                if Pc <= Pe:
                    raise ValueError(f"condensing {Tc:.1f} C not above evaporating {Te:.1f} C")

                def dry(s):    # discharge above the dew point less 5 K
                    return stg(pr, Pe, pr.PT(Pe, Te + s), Pc)["out"]["T"] - dew["T"] - 5

                shN = max(sh, 0.5)
                if dry(shN) < 0:
                    hi = shN
                    while dry(hi) < 0 and hi < Tc - Te - 5:
                        hi += 10
                    if dry(hi) < 0:
                        raise ValueError(f"level with {Te:.0f} -> {Tc:.0f} C cannot compress dry (raise the superheat "
                                         "or change the fluid)")
                    lo, k = shN, 0
                    while k < 30 and hi - lo > 0.05:
                        m = (lo + hi) / 2
                        if dry(m) >= 0:
                            hi = m
                        else:
                            lo = m
                        k += 1
                    shN = hi
                if Tihx is not None and Tihx - Te > shN:
                    shN = Tihx - Te
                suc = pr.PT(Pe, Te + shN)
                qIhx = max(0.0, suc["h"] - pr.PT(Pe, Te + min(sh, shN))["h"])
                cp = stg(pr, Pe, suc, Pc)
                lb = pr.PQ(Pc, 0.0)
                return dict(ev=ev, Pe=Pe, Pc=Pc, suc=suc, cp=cp, hb=lb["h"], Tb=lb["T"], Td=dew["T"], sh=shN, qIhx=qIhx)

            Tcs = [c.T_cond[k] if k < len(c.T_cond) and c.T_cond[k] is not None else 100.0 for k in range(n - 1)]
            # state carried across solves (as the app): condensing-to-air fractions, subcooling caps, IHX liquid T
            fK = [0.2] * n
            subMax = [math.inf] * n
            Tliq = [None] * n

            def fA(k):
                return fK[k] if flag(c.cond_air, k) else 0.0

            def solve_at(Tcn, airIn0):
                """for a top condensing temperature: fixed point on the air temperatures entering each level"""
                airIn = list(airIn0)
                S = None
                for _ in range(60):
                    L = []
                    Te = ev1
                    for k in range(n):
                        Tc = Tcs[k] if k < n - 1 else Tcn
                        Tihx = Tliq[k] - c.dT_ihx if flag(c.ihx, k) and Tliq[k] is not None else None
                        L.append(level(fls[k], Te, Tc, c.superheat, Tihx))
                        Te = Tc - c.dT_casc
                    for k in range(n):
                        l, pr, top = L[k], fls[k], k == n - 1
                        # subcooled into the air down to the air inlet + approach, but not below what the IHX needs
                        Tfloor = l["ev"]["T"] + l["sh"] + 5
                        if flag(c.sub_air, k):
                            T_sub = min(max(airIn[k] + c.dT_GC_min, Tfloor), l["Tb"] - 0.01)
                            l["hsub"] = min(l["hb"], pr.PT(l["Pc"], T_sub)["h"])
                        else:
                            l["hsub"] = l["hb"]
                        if not top and math.isfinite(subMax[k]):
                            l["hsub"] = max(l["hsub"], l["hb"] - max(subMax[k], 0.0))
                        hd = l["cp"]["out"]["h"]
                        l["qa"] = hd - l["hsub"] if top else fA(k) * (hd - l["hb"]) + (l["hb"] - l["hsub"])
                        l["qc"] = 0.0 if top else (1 - fA(k)) * (hd - l["hb"])
                        l["hv"] = l["hsub"] - l["qIhx"]
                        l["he"] = expand_to(pr, l["Pc"], l["hv"], l["Pe"], c.expander_eta)
                        l["ihxShort"] = pr.PH(l["Pc"], l["hv"])["T"] < l["ev"]["T"] + 2
                        if flag(c.ihx, k):
                            t = pr.PH(l["Pc"], l["hsub"])["T"]
                            l["dTliq"] = math.inf if Tliq[k] is None else abs(t - Tliq[k])
                            Tliq[k] = t
                        else:
                            l["dTliq"] = 0.0
                    m = [1.0]
                    for k in range(n - 1):
                        m.append(m[k] * L[k]["qc"] / (L[k + 1]["suc"]["h"] - L[k + 1]["qIhx"] - L[k + 1]["he"]))
                    Q = [m[k] * l["qa"] for k, l in enumerate(L)]
                    Qt = sum(Q)
                    ha0 = air.h(Tain)
                    mA = Qt / (air.h(c.T_air_tar) - ha0)
                    nIn, ha = [], ha0
                    for k in range(n):
                        nIn.append(air.T(ha))
                        ha += Q[k] / mA
                    dfMax = 0.0
                    for k in range(n - 1):
                        if not flag(c.cond_air, k):
                            continue
                        l = L[k]
                        want = mA * (air.h(l["Tb"] - c.dT_GC_min) - air.h(nIn[k])) / m[k] - (l["hb"] - l["hsub"])
                        f = min(max(want / (l["cp"]["out"]["h"] - l["hb"]), 0.0), 0.95)
                        dfMax = max(dfMax, abs(f - fK[k]))
                        fK[k] = fK[k] + 0.7 * (f - fK[k])
                    dsMax = 0.0
                    for k in range(n - 1):
                        if not flag(c.sub_air, k):
                            continue
                        room = max(0.0, mA * (air.h(L[k]["Tb"] - c.dT_GC_min) - air.h(nIn[k])) / m[k])
                        cur = L[k]["hb"] - L[k]["hsub"]
                        if cur > room + 1e-9 or math.isfinite(subMax[k]):
                            nv = subMax[k] + 0.7 * (room - subMax[k]) if math.isfinite(subMax[k]) else room
                            dsMax = max(dsMax, abs(nv - (subMax[k] if math.isfinite(subMax[k]) else cur)))
                            subMax[k] = nv
                    conv = (max(abs(v - airIn[k]) for k, v in enumerate(nIn)) < 0.01 and dfMax < 1e-4
                            and dsMax < 1e-3 and max(l["dTliq"] for l in L) < 0.02)
                    airIn = nIn
                    S = dict(L=L, m=m, Q=Q, mA=mA, airIn=airIn, f=list(fK))
                    if conv:
                        break
                t = S["L"][n - 1]
                S["gap"] = min_gap(air_profile(fls[n - 1], air, t["Pc"], t["cp"]["out"]["h"], t["hsub"],
                                               S["airIn"][n - 1], c.T_air_tar, 18))
                return S

            prTop = fls[n - 1]
            TcritTop = prTop.T_crit - 1 if math.isfinite(prTop.T_crit) else c.T_air_tar + 60
            lo = (Tcs[n - 2] if n > 1 else ev1) + c.dT_casc + 2
            hi = min(TcritTop, c.T_air_tar + 40)
            airIn = [Tain] * n
            # scan up for the first top condensing temperature that keeps the air approach, then bisect
            prev, found, Sf, lastErr = lo, None, None, ""
            T = lo
            while T <= hi + 1e-9:
                St = None
                try:
                    St = solve_at(T, airIn)
                except Exception as exc:
                    lastErr = str(exc)
                if St is not None:
                    airIn = St["airIn"]
                    if St["gap"] >= c.dT_GC_min:
                        found, Sf = T, St
                        break
                prev = T
                T += 5
            if found is None:
                raise ValueError(f"no top condensing temperature up to {hi:.0f} C heats the air to {js_num(c.T_air_tar)} C "
                                 f"with a {js_num(c.dT_GC_min)} K approach" + (f" ({lastErr})" if lastErr else ""))
            lo, hi, S = prev, found, Sf
            k = 0
            while k < 25 and hi - lo > 0.02:
                mid = (lo + hi) / 2
                try:
                    Sm = solve_at(mid, airIn)
                except Exception:
                    Sm = None
                if Sm is not None and Sm["gap"] >= c.dT_GC_min:
                    hi, S, airIn = mid, Sm, Sm["airIn"]
                else:
                    lo = mid
                k += 1
            Tcn = hi
            L, m, Q, mA = S["L"], S["m"], S["Q"], S["mA"]
            airIn = S["airIn"]
            for k, l in enumerate(L):
                if l["ihxShort"]:
                    warnings.append(f"level {k + 1}: its IHX needs more liquid subcooling than is left")
            w = sum(m[k] * l["cp"]["w"] for k, l in enumerate(L))
            wx = sum(m[k] * (l["hv"] - l["he"]) for k, l in enumerate(L))
            wEl = exp_el(w, wx, c.motor_eff)
            qAir = sum(Q)
            cop = qAir / wEl
            xq = [exp_warn(fls[k], l["Pe"], l["he"], warnings, f"level {k + 1} expander") for k, l in enumerate(L)] \
                if c.expander_eta else None
            paths, exchangers = {}, []
            hamb = air.h(Tain)
            exp_name = "Expander" if c.expander_eta else "Expansion valve"
            for k, l in enumerate(L):
                pr, top = fls[k], k == n - 1
                hv, hval, hEvOut, hd = l["hsub"], l["hv"], l["suc"]["h"] - l["qIhx"], l["cp"]["out"]["h"]
                pp = [dict(name=exp_name, pts=comp_path(pr, l["Pc"], hval, l["Pe"], l["he"])),
                      dict(name="Evaporator", pts=iso_path(pr, l["Pe"], l["he"], hEvOut))]
                if l["qIhx"] > 1e-6:
                    pp.append(dict(name="Internal heat exchanger", pts=iso_path(pr, l["Pe"], hEvOut, l["suc"]["h"])))
                pp.append(dict(name="Compressor", pts=comp_path(pr, l["Pe"], l["suc"]["h"], l["Pc"], hd)))
                if top:
                    pp.append(dict(name="Air condenser", pts=iso_path(pr, l["Pc"], hd, hv)))
                else:
                    pp.append(dict(name="Cascade condenser", pts=iso_path(pr, l["Pc"], hd, l["hb"])))
                    if l["hb"] - hv > 1e-6:
                        pp.append(dict(name="Subcooler", pts=iso_path(pr, l["Pc"], l["hb"], hv)))
                if l["qIhx"] > 1e-6:
                    pp.append(dict(name="Internal heat exchanger", pts=iso_path(pr, l["Pc"], hv, hval)))
                paths[f"Level {k + 1}"] = pp
                # air side of this level: the refrigerant duty to air from the coldest point (subcooling, then the
                # fraction f of desuperheating + condensing)
                Ta0, Ta1 = air.T(hamb), air.T(hamb + Q[k] / mA)
                hamb += Q[k] / mA
                l["gap_air"] = None
                if Q[k] > 1e-9:
                    pts = []
                    qs = l["hb"] - hv
                    for f in lin(0.0, 1.0, 16):
                        q = f * l["qa"]
                        if q <= qs or top:
                            T = pr.PH(l["Pc"], hv + q)["T"]
                        else:
                            fr = (q - qs) / fA(k) if fA(k) > 1e-9 else 0.0
                            T = pr.PH(l["Pc"], min(l["hb"] + fr, hd))["T"]
                        pts.append([m[k] * q, T, air.T(air.h(Ta0) + f * Q[k] / mA)])
                    gap = min_gap(pts)
                    l["gap_air"] = gap
                    if not top and gap < c.dT_GC_min - 0.05:
                        warnings.append(f"level {k + 1} air approach {gap:.1f} K, below the {js_num(c.dT_GC_min)} K minimum: "
                                        "lower its condensing-to-air fraction or raise its condensing temperature")
                    if gap < 0:
                        raise ValueError(f"level {k + 1} cannot heat the air from {Ta0:.0f} to {Ta1:.0f} C "
                                         "(temperature cross)")
                    exchangers.append(dict(name=f"Level {k + 1} to air", hot=[[p[0], p[1]] for p in pts],
                                           cold=[[p[0], p[2]] for p in pts], other="dryer air", gap=gap))
                if k < n - 1:
                    nx, prn, qc = L[k + 1], fls[k + 1], m[k] * l["qc"]
                    hot = [[f * qc, pr.PH(l["Pc"], l["hb"] + f * (hd - l["hb"]))["T"]] for f in lin(0, 1, 14)]
                    hc0, hc1 = nx["he"], nx["suc"]["h"] - nx["qIhx"]
                    cold = [[f * qc, prn.PH(nx["Pe"], hc0 + f * (hc1 - hc0))["T"]] for f in lin(0, 1, 14)]
                    exchangers.append(dict(name=f"Cascade {k + 1}->{k + 2}", hot=hot, cold=cold,
                                           other=f"level {k + 2} evaporator"))
                l["gap_ihx"] = None
                if l["qIhx"] > 1e-6:
                    Qi = m[k] * l["qIhx"]
                    hot = [[f * Qi, pr.PH(l["Pc"], hval + f * l["qIhx"])["T"]] for f in lin(0, 1, 12)]
                    cold = [[f * Qi, pr.PH(l["Pe"], hEvOut + f * l["qIhx"])["T"]] for f in lin(0, 1, 12)]
                    gap = min(p[1] - cold[i][1] for i, p in enumerate(hot))
                    l["gap_ihx"] = gap
                    if flag(c.ihx, k) and gap < c.dT_ihx - 0.05:
                        warnings.append(f"level {k + 1} IHX approach {gap:.1f} K, below its {js_num(c.dT_ihx)} K")
                    exchangers.append(dict(name=f"Level {k + 1} IHX", hot=hot, cold=cold,
                                           other=f"level {k + 1} suction", gap=gap))
            l1 = L[0]
            qE = l1["suc"]["h"] - l1["qIhx"] - l1["he"]
            exchangers.append(dict(
                name="Level 1 evaporator",
                hot=[[0.0, source_T(c.T_source_in, c.T_source_out, 0.0)], [qE, source_T(c.T_source_in, c.T_source_out, 1.0)]],
                cold=[[p[2] - l1["he"], p[0]] for p in iso_path(fls[0], l1["Pe"], l1["he"], l1["suc"]["h"] - l1["qIhx"])],
                other="source water"))
            vs = sum(m[k] / l["suc"]["D"] for k, l in enumerate(L))
            VHC = qAir / vs / 1000
            states = {}
            for k, l in enumerate(L):
                pr = fls[k]
                st = lambda P, h: (lambda x: dict(T=x["T"], P=x["P"], h=x["h"], s=x["s"]))(pr.PH(P, h))
                states[f"level_{k + 1}"] = dict(
                    evap_in=st(l["Pe"], l["he"]), evap_out=st(l["Pe"], l["suc"]["h"] - l["qIhx"]),
                    suction=st(l["Pe"], l["suc"]["h"]), discharge=st(l["Pc"], l["cp"]["out"]["h"]),
                    bubble=st(l["Pc"], l["hb"]), subcooler_out=st(l["Pc"], l["hsub"]), valve_in=st(l["Pc"], l["hv"]))
            levels = [dict(Te=l["ev"]["T"], Tc=l["Tb"], Pe_bar=l["Pe"] * 10, Pc_bar=l["Pc"] * 10, PR=l["Pc"] / l["Pe"],
                           eta=l["cp"]["eta"], Tdis=l["cp"]["out"]["T"], flow=m[k], q_air=Q[k], air_in=airIn[k],
                           superheat=l["sh"], q_ihx=l["qIhx"], f_air=fA(k) if k < n - 1 else 1.0,
                           gap_air=l["gap_air"], gap_ihx=l["gap_ihx"], fluid=fls[k].fluid)
                      for k, l in enumerate(L)]
            report = dict(COP=cop, VHC=VHC, levels=levels, w_el_kJkg=wEl, T_cond_top=Tcn, air_per_ref=mA,
                          expander=bool(c.expander_eta), w_exp_kJkg=wx, exp_out_quality=xq, air_W=air.W,
                          q_evap=qE, ref=dict(rho=l1["suc"]["D"], flow=1.0, PR=l1["Pc"] / l1["Pe"]))
            return LoopResult(True, "success", COP=cop, VHC=VHC, Pmax_bar=max(l["Pc"] for l in L) * 10,
                              T_dis_max=max(l["cp"]["out"]["T"] for l in L), q_air=qAir, w_shaft=w, q_source=qE,
                              report=report, states=states, paths=paths, exchangers=exchangers, warnings=warnings)
        except Exception as exc:  # property failures and infeasible designs, as the app reports them
            return _fail(str(exc), warnings=warnings)

    # ---- optimiser ---------------------------------------------------------------------

    def _optimise(self, c: CascadeConfig, max_eval: int) -> LoopResult:
        n, fls = c.n, self.fluids
        dim = n - 1
        if not 2 <= n <= MAX_LEVELS:
            return _fail(f"a cascade has 2 to {MAX_LEVELS} levels")
        if len(fls) < n:
            return _fail("a refrigerant is needed for every level")
        state = dict(nEval=0)
        memo: dict = {}
        ev1 = (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_src
        top = c.T_air_tar - c.dT_GC_min + 30

        def tc(k):
            return fls[k].T_crit - 2 if math.isfinite(fls[k].T_crit) else top

        lo = [ev1 + 10 + k * (c.dT_casc + 10) for k in range(dim)]
        hi = [min(tc(k), top - (dim - 1 - k) * (c.dT_casc + 10)) for k in range(dim)]
        if any(v >= hi[k] for k, v in enumerate(lo)):
            return _fail("no room between the source and the air target for the lower levels")

        def bad(r):
            return r.success and any("approach" in w or "IHX" in w for w in r.warnings)

        def f(x):
            key = ",".join(f"{v:.4f}" for v in x)
            if key in memo:
                return memo[key]
            state["nEval"] += 1
            o = dict(y=-math.inf, r=None)
            try:
                r = self._cascade(replace(c, T_cond=[float(v) for v in x]))
                o = dict(y=r.COP, r=r) if r.success and math.isfinite(r.COP) and not bad(r) else dict(y=-math.inf, r=r)
            except Exception:
                pass
            memo[key] = o
            return o

        def clip(x):
            return [min(max(v, lo[i]), hi[i]) for i, v in enumerate(x)]

        x0 = clip([c.T_cond[k] if k < len(c.T_cond) and c.T_cond[k] is not None else lo[k] for k in range(dim)])
        b0 = f(x0)
        start = b0["r"].COP if b0["r"] is not None and b0["r"].success else float("nan")
        # grid: level k condensing a fraction of the way from the source to the air target, lifts spread evenly
        xg, bg = None, dict(y=-math.inf)
        span = top - ev1
        fr = (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8) if dim == 1 else (0.8, 0.9, 1, 1.1, 1.25)
        for a in fr:
            z = clip([ev1 + a * span if dim == 1 else ev1 + span * ((k + 1) / n) ** a for k in range(dim)])
            o = f(z)
            if o["y"] > bg["y"]:
                bg, xg = o, z
        x, best, converged = compass_search(f, [(x0, b0), (xg, bg)], lo, hi, [8.0] * dim, [0.05] * dim, state,
                                            max_eval)
        if x is None:
            return _fail("no design within the bounds keeps every approach",
                         report=dict(optimiser=dict(nEval=state["nEval"])))
        # round as the app shows them (0.1 K), nudging up until the rounded design still keeps its approaches
        xr = [js_round(v, 1) for v in x]
        for _ in range(20):
            o = f(xr)
            if math.isfinite(o["y"]):
                x, best = xr, o
                break
            xr = [js_round(v + 0.1, 1) for v in xr]
        final = best["r"]

        def near(v, b):
            return abs(v - b) < 0.05

        final.report["optimiser"] = dict(
            COP=best["y"], start=start, startOk=math.isfinite(b0["y"]), T_cond=[float(v) for v in x],
            T_cond_top=final.report["T_cond_top"], nEval=state["nEval"], converged=converged,
            atBound=[near(v, lo[i]) or near(v, hi[i]) for i, v in enumerate(x)], lo=lo, hi=hi)
        return final
