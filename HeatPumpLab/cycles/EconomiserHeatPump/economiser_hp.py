"""Two-stage vapour-compression heat pump for air drying with an intermediate pressure: a flash tank (open
intercooler), a closed economiser (vapour injection), or parallel compression, with a single-stage cycle for
comparison. Standard high-temperature heat pump enhancements (see Arpagaus et al., Energy 152 (2018) 985-1010 and
Adamson et al., Renewable and Sustainable Energy Reviews 167 (2022) 112798).

One refrigerant evaporates on the source water and condenses (or, above its critical pressure, cools) into the dryer
air in one condenser and subcooler, from the air inlet to the target. Per kg of the high-pressure flow:

``flash``
    The condenser liquid expands to the intermediate pressure P_i into a flash tank. The flash vapour (fraction x)
    mixes with the low-stage discharge at the high-stage suction; the saturated liquid expands to the evaporator.
``economiser``
    A fraction y of the condenser liquid expands to P_i and evaporates (to ``superheat_inj`` above saturation) against
    the rest, which is subcooled to dT_econ above the saturation temperature at P_i; the injected vapour mixes with the
    low-stage discharge. The subcooled main liquid expands to the evaporator.
``parallel``
    As ``flash``, but the flash vapour (heated ``superheat_inj`` above saturation by the high-pressure liquid if
    needed for a dry compression) is compressed by its own compressor from P_i straight to the high pressure; the
    low-stage compressor takes the evaporator vapour from the evaporating pressure to the high pressure in one step.
``single``
    One compressor and one expansion valve (no intermediate pressure).

The high pressure is the lowest that heats the air to the target with the air approach dT_GC_min everywhere (the
condenser outlet is the coldest that keeps it): scanned up in pressure, then bisected. Subcritical and transcritical
operation follow from the fluid and the target. P_i is set by its saturation temperature ``T_int`` (C), or, when
None, by the geometric mean of the evaporating and high pressures (the usual starting point); ``optimise()``
searches it for the best COP. The suction superheat is at least ``superheat`` (K), raised if needed for a dry
discharge 5 K above the dew point (overhanging fluids).

Units: MPa, C, kJ/kg, per kg of high-pressure flow (report["flows"] gives the others).
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

ARRANGEMENTS = ("flash", "economiser", "parallel", "single")
EconomiserResult = LoopResult


@dataclass
class EconomiserConfig:
    """Design inputs. Temperatures in C, approaches in K, pressures in bar where named so."""

    fluid: str = "Cyclopentane"
    arrangement: str = "flash"        # "flash" | "economiser" | "parallel" | "single"
    T_int: float | None = None        # saturation temperature at the intermediate pressure; None: geometric mean
    dT_econ: float = 5.0              # economiser: main liquid leaves this far above the saturation T at P_i
    superheat_inj: float = 5.0        # economiser / parallel: injected or parallel-compressor vapour above saturation
    superheat: float = 5.0            # minimum evaporator superheat
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0
    dT_GC_min: float = 20.0           # air approach in the condenser / gas cooler
    T_source_in: float = 60.0
    T_source_out: float = 50.0
    dT_src: float = 5.0               # evaporating temperature below the source water outlet
    air_ratio: float = 0.0            # secondary / primary air (heated with the primary from the inlet to T_mt)
    T_mt: float = 102.0
    P_max_bar: float = 60.0           # highest high-side pressure searched
    comp_model: int = 5               # 0 ideal, 1 constant eta_const, 2 polynomial, 3 Brendel, 5 Pareto front
    eta_const: float = 0.73
    motor_eff: float = DEFAULT_MOTOR_EFF
    air_RH: float = 60.0              # inlet air relative humidity, % (0: dry air through the property package)
    P_air_kPa: float = 101.325
    property_backend: str | None = None   # "refprop" | "coolprop"; None -> $HEATPUMPLAB_PROPERTY_BACKEND


def _fail(message: str, **extra) -> LoopResult:
    return LoopResult(False, message, **extra)


class EconomiserHeatPump:
    """Two-stage heat pump with a flash tank, economiser or parallel compression; see the module docstring."""

    def __init__(self, config: EconomiserConfig | None = None, fluid: Fluid | None = None, air: AirModel | None = None):
        self.c = config if config is not None else EconomiserConfig()
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
        """Best-COP intermediate saturation temperature (compass search from the geometric-mean pressure and a coarse
        grid), rounded to 0.1 K. The result is the full solution there, with the search in report["optimiser"]."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._cached(lambda: self._optimise(c, max_eval))

    @staticmethod
    def _cached(fn):
        if active_property_call_cache() is None:
            with property_call_cache():
                return fn()
        return fn()

    # ---- model -------------------------------------------------------------------------

    def _low_stage(self, c, Pe, Te, P_out):
        """Evaporator outlet with the superheat a dry discharge needs, and its compression to P_out."""
        fl = self.fluid
        sh = self._dry_superheat(c, Pe, Te, P_out, max(c.superheat, 0.5))
        return sh, fl.PT(Pe, Te + sh), stage(fl, Pe, fl.PT(Pe, Te + sh), P_out, None, c.comp_model, c.eta_const)

    def _dry_superheat(self, c, P_in, T_sat, P_out, sh_min):
        return dry_superheat(self.fluid, P_in, T_sat, P_out, sh_min, c.comp_model, c.eta_const)

    def _cycle(self, c: EconomiserConfig, P_H: float, Pe: float, Te: float, P_i: float | None) -> dict:
        """All states at a high pressure P_H (raises ValueError when the air cannot be heated to the target)."""
        fl, air = self.fluid, self.air
        arr = c.arrangement
        two_stage = arr != "single"
        if two_stage and not (Pe < P_i < P_H):
            raise ValueError("the intermediate pressure must lie between the evaporating and high pressures")
        lo_out = P_i if arr in ("flash", "economiser") else P_H
        sh, suc, low = self._low_stage(c, Pe, Te, lo_out)
        if two_stage:
            lq, vq = fl.PQ(P_i, 0.0), fl.PQ(P_i, 1.0)
            if arr == "parallel":
                # the parallel compressor's suction is heated (at least superheat_inj above saturation, more if a dry
                # discharge 5 K above the dew point needs it) by the high-pressure liquid before it expands
                sh_inj = self._dry_superheat(c, P_i, vq["T"], P_H, c.superheat_inj)
            else:
                sh_inj = max(c.superheat_inj, 0.0) if arr == "economiser" else 0.0
            h_inj = fl.PT(P_i, vq["T"] + sh_inj)["h"] if sh_inj > 0 else vq["h"]
            par = stage(fl, P_i, fl.PH(P_i, h_inj), P_H, None, c.comp_model, c.eta_const) if arr == "parallel" else None
        else:
            h_inj = par = None

        def evaluate(v):
            """The cycle at an injected (flash, economiser or parallel) fraction v of the high-pressure flow: the
            condenser inlet, its outlet, and the fraction that outlet implies."""
            hi = None
            if arr == "single":
                h_in = low["out"]["h"]
            elif arr in ("flash", "economiser"):
                st_mix = fl.PH(P_i, (1 - v) * low["out"]["h"] + v * h_inj)
                if not fl.is_gas(st_mix):
                    raise ValueError("the high-stage suction is wet (raise the superheat or the intermediate pressure)")
                hi = stage(fl, P_i, st_mix, P_H, None, c.comp_model, c.eta_const)
                h_in = hi["out"]["h"]
            else:
                h_in = (1 - v) * low["out"]["h"] + v * par["out"]["h"]
            h3 = cooler_outlet(fl, air, P_H, h_in, c.T_air_inlet, c.T_air_tar, c.dT_GC_min)
            if arr == "single":
                v_new = 0.0
            elif arr == "economiser":
                h_sub = min(h3, fl.PT(P_H, vq["T"] + c.dT_econ)["h"])
                v_new = (h3 - h_sub) / (h_inj - h_sub) if h3 > h_sub else 0.0
            else:
                # flash (dh = 0) or parallel (the liquid gives v * dh to the parallel suction before it flashes):
                # v = (h3 - h_l) / (h_g - h_l + dh)
                dh = h_inj - vq["h"]
                v_new = min(max((h3 - lq["h"]) / (vq["h"] - lq["h"] + dh), 0.0), 1.0)
            return dict(h_in=h_in, h3=h3, hi=hi, v_new=v_new)

        # the outlet depends on the condenser inlet, which depends on the injected fraction, which depends on the
        # outlet. Iterate directly; where that does not settle (more injection cools the condenser inlet and can
        # lift the pinch-limited outlet), take the smallest self-consistent fraction by a scan and bisection from 0.
        v, E = 0.0, evaluate(0.0)
        settled = not two_stage or E["v_new"] < 1e-12
        for _ in range(25):
            if settled:
                break
            v_next = E["v_new"]
            E_next = evaluate(v_next)
            if abs(E_next["v_new"] - v_next) < 1e-8:
                v, E, settled = v_next, E_next, True
            else:
                v, E = v_next, E_next
        if not settled:
            def g(u):
                return evaluate(u)["v_new"] - u
            lo_v, hi_v, found = 0.0, None, False
            for k in range(1, 41):
                u = k / 40.0
                if g(u) <= 0:
                    hi_v, found = u, True
                    break
                lo_v = u
            if not found:
                raise ValueError("no self-consistent injected fraction")
            for _ in range(40):
                if hi_v - lo_v < 1e-8:
                    break
                m = (lo_v + hi_v) / 2
                if g(m) > 0:
                    lo_v = m
                else:
                    hi_v = m
            v = hi_v
            E = evaluate(v)
            v = E["v_new"] if abs(E["v_new"] - v) < 1e-5 else v
            E = evaluate(v)
        h_in, h3, hi = E["h_in"], E["h3"], E["hi"]
        x = v if arr in ("flash", "parallel") else 0.0
        y = v if arr == "economiser" else 0.0
        # the evaporator and its feed, per kg of high-pressure flow
        if arr == "single":
            m_low, h_valve = 1.0, h3
        elif arr in ("flash", "parallel"):
            # saturated liquid from the flash tank (with no flash vapour, the liquid itself)
            m_low = 1.0 - x
            h_valve = lq["h"] if x > 0 else h3
        else:
            # 1 - y through the evaporator and the low stage, y injected; the main liquid after the economiser
            m_low = 1.0 - y
            h_valve = h3 - y / (1 - y) * (h_inj - h3) if y > 0 else h3
        he = h_valve
        q_air = h_in - h3
        w = m_low * low["w"] + (hi["w"] if hi is not None else 0.0) + (x * par["w"] if par is not None else 0.0)
        q_src = m_low * (suc["h"] - he)
        gap = min_gap(air_profile(fl, air, P_H, h_in, h3, c.T_air_inlet, c.T_air_tar, 18))
        return dict(P_H=P_H, P_i=P_i, Pe=Pe, Te=Te, sh=sh, suc=suc, low=low, hi=hi, par=par, x=x, y=y, h_inj=h_inj,
                    h_in=h_in, h3=h3, h_valve=h_valve, he=he, m_low=m_low, q_air=q_air, w=w, q_src=q_src, gap=gap)

    def _solve(self, c: EconomiserConfig) -> LoopResult:
        fl, warnings = self.fluid, []
        if c.arrangement not in ARRANGEMENTS:
            return _fail(f"arrangement must be one of {', '.join(ARRANGEMENTS)}")
        try:
            Te = (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_src
            Pe = fl.TQ(Te, 1.0)["P"]
            P_max = c.P_max_bar / 10.0
            Tc_crit = fl.T_crit if math.isfinite(fl.T_crit) else math.inf
            # floor of the pressure scan: condensing well below the air target (desuperheating can heat the top of
            # the air range), and below the critical point
            T0 = min(max(Te + 10.0, c.T_air_tar - 60.0), Tc_crit - 1.0)
            P_lo = fl.TQ(T0, 1.0)["P"] if T0 > Te + 5 else Pe * 2
            if P_lo >= P_max:
                return _fail(f"{fl.fluid} needs more than {c.P_max_bar:g} bar to reach {c.T_air_tar:g} C")

            def P_int(P_H):
                if c.arrangement == "single":
                    return None
                if c.T_int is not None:
                    return fl.TQ(c.T_int, 1.0)["P"]
                return math.sqrt(Pe * P_H)

            def at(P_H):
                try:
                    return self._cycle(c, P_H, Pe, Te, P_int(P_H))
                except SOLVE_ERRORS as exc:
                    return exc

            # scan up in pressure (40 log steps to P_max) for the first feasible design, then bisect in log P
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

    def _result(self, c: EconomiserConfig, S: dict, warnings: list) -> LoopResult:
        fl, air = self.fluid, self.air
        wEl = S["w"] / c.motor_eff
        cop = S["q_air"] / wEl
        P_H, P_i, Pe = S["P_H"], S["P_i"], S["Pe"]
        supercritical = math.isfinite(fl.P_crit) and P_H >= fl.P_crit
        # suction volumes: the largest flow sets the compressor size
        flows = dict(high=1.0, low=S["m_low"], injected=S["x"] if c.arrangement in ("flash", "parallel") else S["y"])
        vs = S["m_low"] / S["suc"]["D"]
        if S["hi"] is not None:
            vs += 1.0 / fl.PH(P_i, S["hi"]["out"]["h"] - S["hi"]["w"])["D"]
        if S["par"] is not None:
            vs += S["x"] / fl.PH(P_i, S["h_inj"])["D"]
        VHC = S["q_air"] / vs / 1000
        st = lambda P, h: (lambda x: dict(T=x["T"], P=x["P"], h=x["h"], s=x["s"]))(fl.PH(P, h))
        states = dict(evap_in=st(Pe, S["he"]), suction=st(Pe, S["suc"]["h"]), low_discharge=st(S["low"]["out"]["P"],
                      S["low"]["out"]["h"]), condenser_in=st(P_H, S["h_in"]), condenser_out=st(P_H, S["h3"]))
        if S["hi"] is not None:
            states["high_suction"] = st(P_i, S["hi"]["out"]["h"] - S["hi"]["w"])
            states["high_discharge"] = st(P_H, S["hi"]["out"]["h"])
        if S["par"] is not None:
            states["parallel_suction"] = st(P_i, S["h_inj"])
            states["parallel_discharge"] = st(P_H, S["par"]["out"]["h"])
        if P_i is not None:
            states["injected"] = st(P_i, S["h_inj"])
        paths = {"Cycle": [
            dict(name="Expansion valve", pts=comp_path(fl, P_H if c.arrangement in ("single",) else P_i, S["h_valve"],
                                                         Pe, S["he"])),
            dict(name="Evaporator", pts=iso_path(fl, Pe, S["he"], S["suc"]["h"])),
            dict(name="Low-stage compressor", pts=comp_path(fl, Pe, S["suc"]["h"], S["low"]["out"]["P"],
                                                           S["low"]["out"]["h"])),
            dict(name="Condenser" + (" / gas cooler" if supercritical else ""),
                 pts=iso_path(fl, P_H, S["h_in"], S["h3"]))]}
        if S["hi"] is not None:
            paths["Cycle"].append(dict(name="High-stage compressor", pts=comp_path(
                fl, P_i, S["hi"]["out"]["h"] - S["hi"]["w"], P_H, S["hi"]["out"]["h"])))
        if S["par"] is not None:
            paths["Cycle"].append(dict(name="Parallel compressor", pts=comp_path(fl, P_i, S["h_inj"], P_H,
                                                                                  S["par"]["out"]["h"])))
        Q = S["q_air"]
        prof = air_profile(fl, air, P_H, S["h_in"], S["h3"], c.T_air_inlet, c.T_air_tar, 18)
        exchangers = [dict(name="Condenser to air", hot=[[p[0], p[1]] for p in prof], cold=[[p[0], p[2]] for p in prof],
                           other="dryer air", gap=S["gap"])]
        qE = S["q_src"]
        exchangers.append(dict(
            name="Evaporator",
            hot=[[0.0, source_T(c.T_source_in, c.T_source_out, 0.0)], [qE, source_T(c.T_source_in, c.T_source_out, 1.0)]],
            cold=[[S["m_low"] * (p[2] - S["he"]), p[0]] for p in iso_path(fl, Pe, S["he"], S["suc"]["h"])],
            other="source water"))
        dis = [S["low"]["out"]["T"]] + ([S["hi"]["out"]["T"]] if S["hi"] else []) + ([S["par"]["out"]["T"]] if S["par"] else [])
        prs = [S["low"]["pr"]] + ([S["hi"]["pr"]] if S["hi"] else []) + ([S["par"]["pr"]] if S["par"] else [])
        if max(prs) > 20:
            warnings.append(f"stage pressure ratio {max(prs):.0f}, beyond the compressor efficiency correlations "
                            "(fitted up to 20): treat the efficiency, and the COP, with caution")
        if max(dis) > 250:
            warnings.append(f"discharge at {max(dis):.0f} C: check the refrigerant's and the oil's stability limits")
        report = dict(COP=cop, VHC=VHC, arrangement=c.arrangement, fluid=fl.fluid, P_high_bar=P_H * 10,
                      P_int_bar=P_i * 10 if P_i is not None else None, P_evap_bar=Pe * 10,
                      T_int=fl.PQ(P_i, 1.0)["T"] if P_i is not None else None, T_evap=S["Te"], superheat=S["sh"],
                      supercritical=supercritical, flows=flows, x_flash=S["x"], y_econ=S["y"], w_el_kJkg=wEl,
                      q_air=Q, q_evap=qE, gap_air=S["gap"], T_dis_low=S["low"]["out"]["T"],
                      T_dis_high=S["hi"]["out"]["T"] if S["hi"] else None,
                      T_dis_parallel=S["par"]["out"]["T"] if S["par"] else None,
                      air_per_ref=Q / (air.h(c.T_air_tar) - air.h(c.T_air_inlet)), air_W=air.W,
                      balance=Q - S["w"] - qE)
        return LoopResult(True, "success", COP=cop, VHC=VHC, Pmax_bar=P_H * 10, T_dis_max=max(dis), q_air=Q,
                          w_shaft=S["w"], q_source=qE, report=report, states=states, paths=paths,
                          exchangers=exchangers, warnings=warnings)

    # ---- optimiser ---------------------------------------------------------------------

    def _optimise(self, c: EconomiserConfig, max_eval: int) -> LoopResult:
        if c.arrangement == "single":
            return self._solve(c)
        fl = self.fluid
        state = dict(nEval=0)
        memo: dict = {}
        Te = (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_src
        top = min(c.T_air_tar, (fl.T_crit - 2) if math.isfinite(fl.T_crit) else c.T_air_tar)
        lo, hi = [Te + 5.0], [top - 5.0]
        if lo[0] >= hi[0]:
            return _fail("no room for an intermediate pressure between the evaporator and the condenser")

        def f(x):
            key = f"{x[0]:.4f}"
            if key in memo:
                return memo[key]
            state["nEval"] += 1
            r = self._solve(replace(c, T_int=float(x[0])))
            o = dict(y=r.COP if r.success else -math.inf, r=r)
            memo[key] = o
            return o

        base = self._solve(replace(c, T_int=None))
        x0 = [min(max(base.report["T_int"], lo[0]), hi[0])] if base.success else None
        b0 = f(x0) if x0 is not None else dict(y=-math.inf, r=base)
        xg, bg = None, dict(y=-math.inf)
        for a in (0.2, 0.35, 0.5, 0.65, 0.8):
            z = [lo[0] + a * (hi[0] - lo[0])]
            o = f(z)
            if o["y"] > bg["y"]:
                bg, xg = o, z
        x, best, converged = compass_search(f, [(x0, b0), (xg, bg)], lo, hi, [8.0], [0.05], [False], state, max_eval,
                                            score=lambda o: o["y"], empty={"y": -math.inf})
        if x is None:
            return _fail("no intermediate pressure gives a feasible design", report=dict(optimiser=dict(nEval=state["nEval"])))
        xr = [js_round(x[0], 1)]
        o = f(xr)
        if math.isfinite(o["y"]):
            x, best = xr, o
        final = best["r"]
        final.report["optimiser"] = dict(COP=best["y"], start=base.COP if base.success else float("nan"),
                                         T_int=float(x[0]), nEval=state["nEval"], converged=converged, lo=lo, hi=hi)
        return final
