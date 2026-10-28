"""Multi-pressure condensation heat pump for air drying: one evaporator and n condensing pressures in series on the
air. This generalises the dual-pressure condensation heat pumps of Dai et al. (Energy Conversion and Management 218
(2020) 112997), which approximate a gliding sink better than a single condensing pressure, to n pressures.

Flow path, per kg of refrigerant through the evaporator (pressures P_1 < P_2 < ... < P_n)::

    evaporator -> compressor 1 (P_e -> P_1) -> at each P_k (k < n) the flow splits:
        m_k - m_(k+1)  -> condenser k (P_k)
        m_(k+1)        -> compressor k+1 (P_k -> P_(k+1))
    compressor n discharges into condenser n (P_n); the liquid of every condenser k > 1 is throttled into condenser
    k - 1, where its flash vapour condenses and its liquid is subcooled with the rest; condenser 1's outlet (all the
    flow) -> valve -> evaporator

with m_1 = 1 the flow through compressor 1 and m_k the flow through compressor k (and condenser k). The air passes
the condensers in series, 1 to n, condenser k heating it from T_(k-1) to T_k (T_0 the inlet, T_n the target). Each
condenser's outlet is the coldest that keeps the air approach.

Unknowns: the flows m_k follow from the air balances (every condenser heats the same air over its own range); the
air temperatures between the condensers are the highest the condensers allow (each lower condenser's pinch is at its
condensing temperature: they are lowered together, by one fraction of the way to the air inlet, until every condenser
is feasible); P_n is the lowest that heats the air to the target. The lower pressures P_1 ... P_(n-1) are set by
their saturation (dew) temperatures ``T_cond`` (C, rising), and ``optimise()`` searches them for the best COP. Suction
superheat is at least ``superheat`` (K), raised if needed for a dry discharge 5 K above the dew point. With n = 2 this
is the dual-pressure cycle (``DualPressureHeatPump``).

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

MultiPressureResult = LoopResult
MAX_PRESSURES = 5


@dataclass
class MultiPressureConfig:
    """Design inputs. Temperatures in C, approaches in K."""

    fluid: str = "Cyclopentane"
    n_pressures: int = 3              # condensing pressures (and compressors), 2 to MAX_PRESSURES
    T_cond: list | None = None        # dew temperatures at P_1 ... P_(n-1), rising; None: equal steps (default_T_cond)
    superheat: float = 5.0            # minimum evaporator superheat
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0
    dT_GC_min: float = 20.0           # air approach in every condenser
    T_source_in: float = 60.0
    T_source_out: float = 50.0
    dT_src: float = 5.0               # evaporating temperature below the source water outlet
    air_ratio: float = 0.0            # secondary / primary air (heated with the primary from the inlet to T_mt)
    T_mt: float = 102.0
    P_max_bar: float = 60.0           # highest top pressure searched
    comp_model: int = 5               # 0 ideal, 1 constant eta_const, 2 polynomial, 3 Brendel, 5 Pareto front
    eta_const: float = 0.73
    motor_eff: float = DEFAULT_MOTOR_EFF
    air_RH: float = 60.0
    P_air_kPa: float = 101.325
    property_backend: str | None = None   # "refprop" | "coolprop"; None -> $HEATPUMPLAB_PROPERTY_BACKEND


def _fail(message: str, **extra) -> LoopResult:
    return LoopResult(False, message, **extra)


def _split(g) -> float:
    """The flow ratio a in (0, 1) where g(a) = 0: g falls with a (more liquid from above, less vapour stopping at this
    pressure), so a bracket on (0, 1) and regula falsi (Illinois). Fails when even all the flow from above (a -> 1)
    gives the condenser more heat than its air takes."""
    lo, hi = 1e-6, 1.0 - 1e-6
    g_lo, g_hi = g(lo), g(hi)
    if not (g_lo > 0 > g_hi):
        raise ValueError("no split of the flow balances the condensers")
    side, a = 0, lo
    for _ in range(80):
        a = hi - g_hi * (hi - lo) / (g_hi - g_lo)
        ga = g(a)
        if abs(ga) < 1e-11 or hi - lo < 1e-13:
            break
        if ga > 0:
            lo, g_lo = a, ga
            if side == 1:
                g_hi /= 2
            side = 1
        else:
            hi, g_hi = a, ga
            if side == -1:
                g_lo /= 2
            side = -1
    return a


def evap_T(c) -> float:
    return (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_src


def default_T_cond(c, T_crit: float = math.inf) -> list:
    """Lower condensing temperatures in equal steps from the evaporating temperature to the air target (or 2 K below
    the critical temperature): T_e + k (top - T_e) / n for k = 1 ... n - 1."""
    Te, n = evap_T(c), int(c.n_pressures)
    top = min(c.T_air_tar, T_crit - 2.0)
    return [Te + (top - Te) * k / n for k in range(1, n)]


class MultiPressureHeatPump:
    """Multi-pressure condensation heat pump; see the module docstring."""

    def __init__(self, config: MultiPressureConfig | None = None, fluid: Fluid | None = None, air: AirModel | None = None):
        self.c = config if config is not None else MultiPressureConfig()
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
        """Best-COP lower condensing temperatures (compass search from the configured values and a coarse grid),
        rounded to 0.1 K. The result is the full solution there, with the search in report["optimiser"]."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._optimise(c, max_eval))

    @staticmethod
    def _cached(fn):
        if active_property_call_cache() is None:
            with property_call_cache():
                return fn()
        return fn()

    # ---- model -------------------------------------------------------------------------

    def _at(self, c: MultiPressureConfig, P_H: float, Pe: float, Te: float, Ps: list, low: dict) -> dict:
        """All states at a top pressure P_H (raises ValueError when the air cannot be heated to the target)."""
        fl, air = self.fluid, self.air
        n = len(Ps) + 1
        P = list(Ps) + [P_H]
        comps = low["comps"] + [stage(fl, Ps[-1], low["comps"][-1]["out"], P_H, None, c.comp_model, c.eta_const)]
        h_dis = [s["out"]["h"] for s in comps]               # discharge of compressor k at P_k
        Tain, Ttar, dT = c.T_air_inlet, c.T_air_tar, c.dT_GC_min

        def sections(Ta: list) -> dict:
            """Every condenser for the air leaving condenser k at Ta[k] (k < n - 1); flows from the air balances,
            from the top down (scaled to m_n = 1, then to m_1 = 1)."""
            T = [Tain] + list(Ta) + [Ttar]
            h_in, h_out, r = [0.0] * n, [0.0] * n, [1.0] * n
            h_in[-1] = h_dis[-1]
            h_out[-1] = cooler_outlet(fl, air, P[-1], h_in[-1], T[-2], Ttar, dT)
            m = [0.0] * n
            m[-1] = 1.0
            G = (h_in[-1] - h_out[-1]) / (air.h(Ttar) - air.h(T[-2]))     # air per kg through compressor n
            for k in range(n - 2, -1, -1):
                dh = air.h(T[k + 1]) - air.h(T[k])

                def g(a):
                    """Air balance at a = m_(k+1) / m_k: (m_(k+1) / a) (h_in - h_out) heats G kg of air over T[k] ->
                    T[k+1]; as the ratio it implies less a."""
                    hi_ = (1 - a) * h_dis[k] + a * h_out[k + 1]   # vapour stopping at P_k + throttled liquid from above
                    return m[k + 1] * (hi_ - cooler_outlet(fl, air, P[k], hi_, T[k], T[k + 1], dT)) / (G * dh) - a

                a = _split(g)
                h_in[k] = (1 - a) * h_dis[k] + a * h_out[k + 1]
                h_out[k] = cooler_outlet(fl, air, P[k], h_in[k], T[k], T[k + 1], dT)
                r[k] = a
                m[k] = m[k + 1] / a
            m = [v / m[0] for v in m]
            q = [m[k] * (h_in[k] - h_out[k]) for k in range(n)]   # per kg of evaporator flow
            return dict(T=T, m=m, r=r[:-1], h_in=h_in, h_out=h_out, q=q)

        # the air between the condensers: each lower condenser's pinch sits at its condensing temperature, so start
        # from those less the approach and lower all of them together (a fraction lam of the way to the air inlet)
        # until every condenser is feasible: a scan, then bisection on lam
        T_dew = [fl.PQ(Pk, 1.0)["T"] for Pk in Ps]
        top = [t - dT for t in T_dew]
        if top[0] <= Tain + 1:
            raise ValueError(f"the lowest condensing temperature ({T_dew[0]:.0f} C) is too close to the air inlet")
        if any(t >= Ttar for t in top):
            raise ValueError("a lower condensing temperature is too close to the air target")
        floor = [Tain + 1.0 + k for k in range(n - 1)]

        def temps(lam):
            return [t - lam * (t - f) for t, f in zip(top, floor)]

        def try_sections(lam):
            Ta = temps(lam)
            if any(b <= a_ for a_, b in zip(Ta, Ta[1:])):
                raise ValueError("the condensing temperatures must rise")
            return sections(Ta)

        best, lo_l, hi_l = None, 1.0, None
        for k in range(21):
            lam = k / 20
            try:
                best, lo_l = (lam, try_sections(lam)), lam
                break
            except SOLVE_ERRORS:
                hi_l = lam
        if best is None:
            raise ValueError("the lower condensers cannot take any of the air heating")
        span = max(t - f for t, f in zip(top, floor))
        if hi_l is not None:
            for _ in range(25):
                if (lo_l - hi_l) * span < 0.02:
                    break
                mid = (lo_l + hi_l) / 2
                try:
                    best, lo_l = (mid, try_sections(mid)), mid
                except SOLVE_ERRORS:
                    hi_l = mid
        S = best[1]
        q_air = sum(S["q"])
        w = sum(S["m"][k] * comps[k]["w"] for k in range(n))
        q_src = low["suc"]["h"] - S["h_out"][0]
        gaps = [min_gap(air_profile(fl, air, P[k], S["h_in"][k], S["h_out"][k], S["T"][k], S["T"][k + 1], 18))
                for k in range(n)]
        return dict(P_H=P_H, P=P, Pe=Pe, Te=Te, sh=low["sh"], suc=low["suc"], comps=comps, h_dis=h_dis,
                    he=S["h_out"][0], q_air=q_air, w=w, q_src=q_src, gap=min(gaps), gaps=gaps, **S)

    def _check(self, c: MultiPressureConfig):
        """The lower condensing temperatures, or a failure message."""
        fl = self.fluid
        n = int(c.n_pressures)
        if not 2 <= n <= MAX_PRESSURES:
            return f"2 to {MAX_PRESSURES} condensing pressures"
        Tc = list(c.T_cond) if c.T_cond is not None else default_T_cond(c, fl.T_crit if math.isfinite(fl.T_crit) else math.inf)
        if len(Tc) != n - 1:
            return f"{n} condensing pressures need {n - 1} lower condensing temperatures"
        if any(b <= a for a, b in zip(Tc, Tc[1:])):
            return "the lower condensing temperatures must rise"
        if math.isfinite(fl.T_crit) and Tc[-1] >= fl.T_crit - 1:
            return f"the condensing temperatures must be below {fl.fluid}'s critical temperature"
        if Tc[0] <= evap_T(c) + 5:
            return "the condensing temperatures must be above the evaporating temperature"
        return [float(t) for t in Tc]

    def _solve(self, c: MultiPressureConfig) -> LoopResult:
        fl, warnings = self.fluid, []
        try:
            Tc = self._check(c)
            if isinstance(Tc, str):
                return _fail(Tc)
            Te = evap_T(c)
            Pe = fl.TQ(Te, 1.0)["P"]
            Ps = [fl.TQ(t, 1.0)["P"] for t in Tc]
            sh = dry_superheat(fl, Pe, Te, Ps[0], c.superheat, c.comp_model, c.eta_const)
            suc = fl.PT(Pe, Te + sh)
            comps = [stage(fl, Pe, suc, Ps[0], None, c.comp_model, c.eta_const)]
            for k in range(1, len(Ps)):
                comps.append(stage(fl, Ps[k - 1], comps[-1]["out"], Ps[k], None, c.comp_model, c.eta_const))
            low = dict(sh=sh, suc=suc, comps=comps)
            P_max = c.P_max_bar / 10.0
            P_lo = Ps[-1] * 1.05
            if P_lo >= P_max:
                return _fail(f"{fl.fluid} needs more than {c.P_max_bar:g} bar")

            def at(P_H):
                try:
                    return self._at(c, P_H, Pe, Te, Ps, low)
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
            return self._result(c, Tc, S, warnings)
        except SOLVE_ERRORS as exc:
            return _fail(str(exc), warnings=warnings)

    def _result(self, c: MultiPressureConfig, Tc: list, S: dict, warnings: list) -> LoopResult:
        fl, air = self.fluid, self.air
        n = len(S["P"])
        wEl = S["w"] / c.motor_eff
        cop = S["q_air"] / wEl
        P_H, Pe, P = S["P_H"], S["Pe"], S["P"]
        supercritical = math.isfinite(fl.P_crit) and P_H >= fl.P_crit
        prs = [s["pr"] for s in S["comps"]]
        if max(prs) > 20:
            warnings.append(f"stage pressure ratio {max(prs):.0f}, beyond the compressor efficiency correlations "
                            "(fitted up to 20): treat the efficiency, and the COP, with caution")
        dis = [s["out"]["T"] for s in S["comps"]]
        if max(dis) > 250:
            warnings.append(f"discharge at {max(dis):.0f} C: check the refrigerant's and the oil's stability limits")
        m = S["m"]
        # suction volume: compressor 1 takes the evaporator vapour, compressor k the discharge of k - 1
        vs = 1.0 / S["suc"]["D"] + sum(m[k] / S["comps"][k - 1]["out"]["D"] for k in range(1, n))
        VHC = S["q_air"] / vs / 1000
        st = lambda P_, h: (lambda x: dict(T=x["T"], P=x["P"], h=x["h"], s=x["s"]))(fl.PH(P_, h))
        states = dict(evap_in=st(Pe, S["he"]), suction=st(Pe, S["suc"]["h"]))
        for k in range(n):
            states[f"discharge_{k + 1}"] = st(P[k], S["h_dis"][k])
            states[f"condenser_{k + 1}_in"] = st(P[k], S["h_in"][k])
            states[f"condenser_{k + 1}_out"] = st(P[k], S["h_out"][k])
        top_name = f"Condenser {n}" + (" / gas cooler" if supercritical else "")
        cname = [f"Condenser {k + 1}" for k in range(n - 1)] + [top_name]
        paths = [dict(name="Expansion valve", pts=comp_path(fl, P[0], S["h_out"][0], Pe, S["he"])),
                 dict(name="Evaporator", pts=iso_path(fl, Pe, S["he"], S["suc"]["h"])),
                 dict(name="Compressor 1", pts=comp_path(fl, Pe, S["suc"]["h"], P[0], S["h_dis"][0]))]
        for k in range(1, n):
            paths.append(dict(name=f"Compressor {k + 1}", pts=comp_path(fl, P[k - 1], S["h_dis"][k - 1], P[k], S["h_dis"][k])))
        for k in range(n - 1, -1, -1):
            paths.append(dict(name=cname[k], pts=iso_path(fl, P[k], S["h_in"][k], S["h_out"][k])))
            if k > 0:
                paths.append(dict(name=f"Valve {k + 1} to {k}", pts=comp_path(fl, P[k], S["h_out"][k], P[k - 1], S["h_out"][k])))
        exchangers = []
        for k in range(n):
            prof = air_profile(fl, air, P[k], S["h_in"][k], S["h_out"][k], S["T"][k], S["T"][k + 1], 18)
            exchangers.append(dict(name=f"{cname[k]} to air", hot=[[m[k] * p[0], p[1]] for p in prof],
                                   cold=[[m[k] * p[0], p[2]] for p in prof], other="dryer air", gap=S["gaps"][k]))
        exchangers.append(dict(
            name="Evaporator",
            hot=[[0.0, source_T(c.T_source_in, c.T_source_out, 0.0)],
                 [S["q_src"], source_T(c.T_source_in, c.T_source_out, 1.0)]],
            cold=[[p[2] - S["he"], p[0]] for p in iso_path(fl, Pe, S["he"], S["suc"]["h"])],
            other="source water"))
        report = dict(COP=cop, VHC=VHC, fluid=fl.fluid, n_pressures=n, T_cond=list(Tc), P_cond_bar=[p * 10 for p in P],
                      P_high_bar=P_H * 10, P_evap_bar=Pe * 10, T_evap=S["Te"], T_air_between=S["T"][1:-1],
                      flows=list(m), q_air_levels=list(S["q"]), gaps=list(S["gaps"]), T_dis=dis, PR=prs,
                      superheat=S["sh"], supercritical=supercritical, w_el_kJkg=wEl, q_air=S["q_air"],
                      q_evap=S["q_src"], gap_air=S["gap"], air_per_ref=S["q_air"] / (air.h(c.T_air_tar) - air.h(c.T_air_inlet)),
                      air_W=air.W, balance=S["q_air"] - S["w"] - S["q_src"])
        return LoopResult(True, "success", COP=cop, VHC=VHC, Pmax_bar=P_H * 10, T_dis_max=max(dis), q_air=S["q_air"],
                          w_shaft=S["w"], q_source=S["q_src"], report=report, states=states, paths={"Cycle": paths},
                          exchangers=exchangers, warnings=warnings)

    # ---- optimiser ---------------------------------------------------------------------

    def _optimise(self, c: MultiPressureConfig, max_eval: int) -> LoopResult:
        fl = self.fluid
        n = int(c.n_pressures)
        if not 2 <= n <= MAX_PRESSURES:
            return _fail(f"2 to {MAX_PRESSURES} condensing pressures")
        dim = n - 1
        state = dict(nEval=0)
        memo: dict = {}
        Te = evap_T(c)
        top = min(c.T_air_tar, (fl.T_crit - 2) if math.isfinite(fl.T_crit) else c.T_air_tar)
        lo, hi = [max(Te + 10.0, c.T_air_inlet + c.dT_GC_min + 5.0)] * dim, [top - 5.0] * dim
        if lo[0] >= hi[0]:
            return _fail("no room for the lower condensing temperatures between the source and the air target")

        def f(x):
            key = ",".join(f"{v:.4f}" for v in x)
            if key in memo:
                return memo[key]
            state["nEval"] += 1
            if any(b <= a for a, b in zip(x, x[1:])):
                o = dict(y=-math.inf, r=None)
            else:
                r = self._solve(replace(c, T_cond=[float(v) for v in x]))
                o = dict(y=r.COP if r.success else -math.inf, r=r)
            memo[key] = o
            return o

        T0 = list(c.T_cond) if c.T_cond is not None else default_T_cond(c, fl.T_crit if math.isfinite(fl.T_crit) else math.inf)
        x0 = [min(max(v, lo[0]), hi[0]) for v in T0] if len(T0) == dim else None
        b0 = f(x0) if x0 is not None else dict(y=-math.inf, r=None)
        # grid: level k at a fraction of the way up, a ** ((dim - k) / dim), so the levels rise
        xg, bg = None, dict(y=-math.inf)
        for frac in (0.2, 0.35, 0.5, 0.65, 0.8):
            z = [lo[0] + frac ** ((dim - k) / dim) * (hi[0] - lo[0]) for k in range(dim)]
            o = f(z)
            if o["y"] > bg["y"]:
                bg, xg = o, z
        x, best, converged = compass_search(f, [(x0, b0), (xg, bg)], lo, hi, [8.0] * dim, [0.05] * dim, [False] * dim,
                                            state, max_eval, score=lambda o: o["y"], empty={"y": -math.inf})
        if x is None:
            return _fail("no set of lower condensing temperatures gives a feasible design",
                         report=dict(optimiser=dict(nEval=state["nEval"])))
        xr = [js_round(v, 1) for v in x]
        o = f(xr)
        if math.isfinite(o["y"]):
            x, best = xr, o
        final = best["r"]
        final.report["optimiser"] = dict(COP=best["y"], start=b0["y"], T_cond=[float(v) for v in x],
                                         nEval=state["nEval"], converged=converged, lo=lo, hi=hi)
        return final
