"""Heat pump + MVR steam ("HP + MVR"): a vapour-compression bottom cycle heats the dryer air first and raises
low-pressure steam with the rest of its condenser heat; a chain of mechanical-vapour-recompression (MVR) fans then
compresses the steam, which desuperheats and condenses into the air up to the target. After
Walmsley et al., "Steam cogeneration heat pumps for the food, dairy and meat industries", 5th High-Temperature
Heat Pump Symposium 2026, Book of Presentations, pp. 270-275.

Bottom cycle (ammonia by default)::

    evaporator (source water, dew T_source_out - dT_src, superheat to at most T_source_in - dT_src)
      -> [IHX cold side] -> compressor -> condenser at T_cond: steam generator (raises steam at T_cond - dT_sg)
      -> air condenser (air T_air_inlet -> T_cond - dT_GC_min) -> subcooler into the air (sub_air)
      -> [IHX hot side] -> expansion valve or expander -> evaporator

Steam: raised at T_cond - dT_sg (sub-atmospheric below 100 C); n_mvr fans at a pressure ratio PR_mvr each. After
each fan but the last, water is injected to sh_inject above saturation (inject), then a share ``split[i]`` of the
exhaust goes to the process (condenses into the air at that pressure) and the rest to the next fan; the last fan's
exhaust all goes to the process. The condensate flashes back to the steam generator. With ``split_fill`` each
fan's process share instead heats the air to that fan's saturation temperature less the approach (every cap
lowered by a common margin when the next stream would miss its approach). The fans have their own efficiency
model (``fan_model``, ``fan_eta`` for the constant model) and motor (``fan_motor``).

Air: heated from T_air_inlet by the bottom cycle (subcooler and condenser) to T_cond - dT_GC_min, then by the
process steam streams in turn to T_air_tar. Closure: per kg of bottom refrigerant, the condenser duty is the air
section plus the steam raised; the steam (per kg raised) heats the rest of the air. When the subcooler pinches
inside (a liquid whose cp climbs steeply, like propane near its critical point) the liquid outlet is raised until
the approach holds, or the bottom cycle heats the air less (``t1_drop``), whichever gives the higher COP.

Duties are per kg of bottom refrigerant. COP = air duty / electrical input (compressor / motor_eff, fans /
fan_motor, less any expander work x motor_eff).

Port of the "Advanced dryer heat pumps" web app (Ahuora Centre, University of Waikato): DRYCASC.booster,
boosterFill, boosterAt and boosterOptimiseSteps. Units follow HeatPumpLab: MPa, C, kJ/kg, kg/m3.
"""
from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field, replace
from typing import Sequence

from ..common.dryer import (
    MOTOR, SOLVE_ERRORS, AirModel, Fluid, comp_path, compass_search, eta_model, expand_to, expander_electric,
    expander_warning, is_gas, iso_path, js_round, lin, memo_objective, min_gap, point, stage,
)

# ----------------------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------------------


@dataclass
class MVRHPConfig:
    """Design inputs (the web app's defaults). Temperatures in C, approaches in K."""

    refrigerant: str = "Ammonia"     # bottom cycle
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0
    dT_GC_min: float = 20.0          # minimum approach to the air (bottom cycle and steam)
    T_source_in: float = 40.0        # source water in
    T_source_out: float = 30.0       # source water out
    dT_src: float = 5.0              # evaporator dew below the source outlet; superheat limit below its inlet
    T_cond: float = 90.0             # bottom condensing temperature
    dT_sg: float = 5.0               # steam raised this far below T_cond
    n_mvr: int = 9                   # MVR fans
    PR_mvr: float = 1.5              # pressure ratio per fan
    inject: bool = True              # water injection after each fan but the last
    sh_inject: float = 5.0           # K above saturation after injection
    split: Sequence[float] = ()      # share of fan i's exhaust to the process (all but the last fan; missing = 0)
    split_fill: bool = False         # set the shares to fill the air to each fan's saturation less the approach
    superheat: float = 5.0           # evaporator outlet superheat
    sub_air: bool = True             # subcool the liquid into the air (to T_air_inlet + dT_GC_min)
    ihx: bool = False                # bottom IHX: liquid heats the suction to its temperature less dT_ihx
    dT_ihx: float = 5.0
    comp_model: int = 5              # bottom compressor: 0 ideal, 1 constant, 2 polynomial, 3 Brendel, 5 Pareto front
    comp_eta: float = 0.73           # isentropic efficiency for the constant model
    motor_eff: float = MOTOR
    fan_model: int | None = 5        # MVR fans' efficiency model (None: the compressor's)
    fan_eta: float = 0.8             # constant model
    fan_motor: float | None = MOTOR  # None: motor_eff
    expander_eta: float | None = None    # bottom-cycle expander instead of the valve (isentropic efficiency)
    air_RH: float = 60.0             # dryer air inlet relative humidity, %
    P_atm: float = 101.325           # kPa
    backend: str | None = None       # property backend: "coolprop", "refprop" or None (the package default)


@dataclass
class MVRHPResult:
    success: bool
    message: str
    COP: float | None = None
    VHC: float | None = None          # MJ/m3, air duty per bottom compressor suction volume
    q_air: float | None = None        # kJ per kg of bottom refrigerant
    w_shaft: float | None = None
    q_source: float | None = None
    report: dict = field(default_factory=dict)
    states: dict = field(default_factory=dict)
    paths: list = field(default_factory=list)       # [{name, fluid, flow, paths: [{name, pts: [[T, P bar, h, s]]}]}]
    profiles: list = field(default_factory=list)    # exchanger profiles [{name, hot: [[q, T]], cold: [[q, T]]}]
    warnings: list = field(default_factory=list)

    def to_dict(self):
        return asdict(self)


@dataclass
class MVROptimum:
    success: bool
    message: str = ""
    COP: float | None = None
    T_cond: float | None = None
    n_mvr: int | None = None
    splits: list | None = None
    start: float | None = None        # COP of the design as set (NaN when it fails)
    start_ok: bool = False
    n_eval: int = 0
    converged: bool = False
    at_bound: dict = field(default_factory=dict)
    bounds: tuple = ()
    result: MVRHPResult | None = None


def _fail(message: str, **kw) -> MVRHPResult:
    return MVRHPResult(False, message, **kw)


def _state(st: dict) -> dict:
    return {k: st.get(k) for k in ("T", "P", "h", "s", "D", "q", "phase")}


# ----------------------------------------------------------------------------------------
# Model
# ----------------------------------------------------------------------------------------


class MVRHeatPump:
    """Heat pump + MVR steam; see the module docstring."""

    def __init__(self, config: MVRHPConfig | None = None, fluid: Fluid | None = None, water: Fluid | None = None,
                 air: AirModel | None = None):
        self.c = config if config is not None else MVRHPConfig()
        self.fb = fluid if fluid is not None else Fluid(self.c.refrigerant, self.c.backend)
        self.water = water if water is not None else Fluid("Water", self.c.backend)
        self.air = air if air is not None else AirModel.humid(self.c.T_air_inlet, self.c.air_RH, self.c.P_atm * 1000)

    def _eta_comp(self, pr: float, P_in: float) -> float:
        return eta_model(self.c.comp_model, self.c.comp_eta, pr, P_in)

    def _eta_fan(self, P_in: float) -> float:
        c = self.c
        if c.fan_model is None:
            return self._eta_comp(c.PR_mvr, P_in)
        return eta_model(c.fan_model, c.fan_eta, c.PR_mvr, P_in)

    # ---- solve (booster) ---------------------------------------------------------------

    def solve(self) -> MVRHPResult:
        """The design at the configured inputs (the page's ``booster``)."""
        r = self._fill()
        if not r.success or not (r.report["sub_raise"] > 0.05):
            return r
        # the subcooler pinched inside and raising the liquid outlet fixed it; the other way is to keep the liquid at
        # the air inlet plus the approach and let the bottom cycle heat the air less (the steam takes the rest).
        # Keep whichever gives the higher COP.
        dT = self.c.dT_GC_min

        def at(d):
            return self._fill(t1_drop=d, no_sub_raise=True)

        def ok(q):
            return q.success and q.report["gap_bottom"] >= dT - 0.02

        hi = 2.0
        rh = at(hi)
        while not ok(rh) and hi < 64:
            hi *= 2
            rh = at(hi)
        if not ok(rh):
            return r
        lo, k = 0.0, 0
        while k < 12 and hi - lo > 0.05:
            m = (lo + hi) / 2
            rm = at(m)
            if ok(rm):
                hi, rh = m, rm
            else:
                lo = m
            k += 1
        return rh if rh.COP > r.COP else r

    def _fill(self, t1_drop: float = 0.0, no_sub_raise: bool = False) -> MVRHPResult:
        """split_fill: lower every fan's air cap by the smallest common margin that keeps the steam approach."""
        c = self.c
        if not c.split_fill:
            return self.evaluate(t1_drop=t1_drop, no_sub_raise=no_sub_raise)
        dT = c.dT_GC_min

        def at(m):
            return self.evaluate(t1_drop=t1_drop, no_sub_raise=no_sub_raise, fill_margin=m)

        def ok(r):
            return r.success and r.report["gap_steam"] >= dT - 0.05

        r0 = at(0.0)
        if not r0.success or ok(r0):
            return r0
        hi = 2.0
        rh = at(hi)
        while not ok(rh) and hi < 64:
            hi *= 2
            rh = at(hi)
        if not ok(rh):
            return r0
        lo, k = 0.0, 0
        while k < 18 and hi - lo > 0.01:
            m = (lo + hi) / 2
            rm = at(m)
            if ok(rm):
                hi, rh = m, rm
            else:
                lo = m
            k += 1
        rh.report["fill_margin"] = hi
        return rh

    # ---- one design (boosterAt) ----------------------------------------------------------

    def evaluate(self, t1_drop: float = 0.0, no_sub_raise: bool = False, fill_margin: float | None = None) -> MVRHPResult:
        try:
            return self._evaluate(t1_drop, no_sub_raise, fill_margin)
        except SOLVE_ERRORS as exc:
            return _fail(f"{exc}")

    def _evaluate(self, t1_drop, no_sub_raise, fill_margin) -> MVRHPResult:
        c, fb, water, air = self.c, self.fb, self.water, self.air
        warnings: list[str] = []
        motor = c.motor_eff
        mot_f = c.fan_motor if c.fan_motor is not None else motor
        dT = c.dT_GC_min
        Tain = c.T_air_inlet
        Te = c.T_source_out - c.dT_src
        ev = fb.TQ(Te, 1.0)
        Pe = ev["P"]
        cb = fb.TQ(c.T_cond, 0.0)
        Pc = max(cb["P"], fb.TQ(c.T_cond, 1.0)["P"])
        ev_out = fb.PT(Pe, min(Te + c.superheat, c.T_source_in - c.dT_src)) if c.superheat > 0 else ev

        def comp_at(hsub):
            # suction: the evaporator outlet, or with the IHX on, heated by the liquid (after any subcooling to the
            # air) to its temperature less the IHX approach
            suc = ev_out
            if c.ihx:
                Ts_ = fb.PH(Pc, hsub)["T"] - c.dT_ihx
                if Ts_ > ev_out["T"] + 0.01:
                    suc = fb.PT(Pe, Ts_)
            cp = stage(fb, Pe, suc, Pc, self._eta_comp(Pc / Pe, Pe))
            return suc, cp, suc["h"] - ev_out["h"]

        hb = fb.PQ(Pc, 0.0)["h"]
        # steam: the fan states do not depend on the splits, only the masses do, so the chain is built once per kg
        # of steam raised
        Ts = c.T_cond - c.dT_sg
        sv, sl = water.TQ(Ts, 1.0), water.TQ(Ts, 0.0)
        nS = max(1, js_round(c.n_mvr))
        T1 = c.T_cond - dT - (t1_drop or 0.0)        # air after the bottom cycle
        if T1 <= Tain:
            return _fail("condensing temperature too low to heat the air")
        ha0, ha1, ha2 = air.h(Tain), air.h(T1), air.h(c.T_air_tar)
        P, st = sv["P"], sv
        fans, spaths = [], []
        for i in range(nS):
            P2 = P * c.PR_mvr
            e = self._eta_fan(P)
            s = stage(water, P, st, P2, e)
            last = i == nS - 1
            spaths.append({"name": f"MVR {i + 1}", "pts": comp_path(water, P, st["h"], P2, s["out"]["h"], 5)})
            out, inj = s["out"], 0.0
            if c.inject and not last:
                tgt = water.PT(P2, water.PQ(P2, 1.0)["T"] + max(c.sh_inject, 0.1))
                if out["h"] > tgt["h"]:
                    inj = (out["h"] - tgt["h"]) / (tgt["h"] - sl["h"])
                    spaths.append({"name": "Water injection", "pts": [point(out), point(tgt)]})
                    out = tgt
            lf, sat = water.PQ(P2, 0.0), water.PQ(P2, 1.0)
            fans.append({"P_in": P, "P_out": P2, "T_out": s["out"]["T"], "w": s["w"], "eta": e, "inj": inj, "out": out,
                         "lf": lf, "Tsat": sat["T"], "q": out["h"] - lf["h"]})
            P, st = P2, out

        def chain(x):
            """Masses per kg raised for process fractions x (the last forced to 1)."""
            M, J, wS, E = 1.0, 0.0, 0.0, []
            for i, f in enumerate(fans):
                wS += M * f["w"]
                J += M * f["inj"]
                Mo = M * (1 + f["inj"])
                xi = 1.0 if i == nS - 1 else min(max(x[i], 0.0), 1.0)
                f["M"] = M
                E.append(Mo * xi)
                M = Mo * (1 - xi)
            qS = sum(e * fans[i]["q"] for i, e in enumerate(E))
            g = sv["h"] + J * sl["h"] - sum(e * fans[i]["lf"]["h"] for i, e in enumerate(E))   # generator, per kg raised
            return E, wS, qS, g

        def split_at(i):
            try:
                v = float(c.split[i])
            except (IndexError, TypeError, ValueError):
                return 0.0
            return v if math.isfinite(v) else 0.0

        x = [1.0 if i == nS - 1 else min(max(split_at(i), 0.0), 1.0) for i in range(nS)]
        if c.split_fill and nS > 1:
            def fill(a):
                M, ha, xs = 1.0, ha1, []
                for i, f in enumerate(fans):
                    Mo = M * (1 + f["inj"])
                    if i == nS - 1:
                        xs.append(1.0)
                        ha += Mo * f["q"] / a
                        continue
                    cap = air.h(min(f["Tsat"] - dT - (fill_margin or 0.0), c.T_air_tar))
                    Ei = min(Mo, a * (cap - ha) / f["q"]) if cap > ha else 0.0
                    xs.append(Ei / Mo if Mo > 0 else 0.0)
                    ha += Ei * f["q"] / a
                    M = Mo - Ei
                return xs, ha - ha2

            lo, hi = math.log(1e-3), math.log(1e5)
            if fill(math.exp(lo))[1] > 0 and fill(math.exp(hi))[1] < 0:
                for _ in range(80):
                    m = (lo + hi) / 2
                    if fill(math.exp(m))[1] > 0:
                        lo = m
                    else:
                        hi = m
                x = fill(math.exp(lo))[0]
        E, wS, qS, g = chain(x)
        last_f = fans[-1]
        Pf, stF = last_f["P_out"], last_f["out"]
        if last_f["Tsat"] < c.T_air_tar and stF["T"] < c.T_air_tar + dT:
            warnings.append("final steam is cooler than the air target plus the approach: add MVR stages")

        def bottom(sub_off):
            # subcooling to the air, the air section up to T1, the rest raises steam. The subcooled liquid's T(h)
            # bends (its cp rises with temperature), so the pinch can sit a little inside the subcooler: the
            # caller raises the liquid outlet until it holds
            hsub = fb.PT(Pc, min(Tain + dT + sub_off, cb["T"] - 0.01))["h"] if c.sub_air else hb
            suc, cp, q_ihx = comp_at(hsub)
            q_cond = cp["out"]["h"] - hsub
            # per kg bottom refrigerant: mA (ha1 - ha0) + raised g = q_cond ; mA (ha2 - ha1) = raised qS
            k2 = (ha2 - ha1) / qS
            mA = q_cond / ((ha1 - ha0) + k2 * g)
            raised = mA * k2
            qA1 = mA * (ha1 - ha0)
            p1 = [[q, fb.PH(Pc, hsub + q)["T"], air.T(ha0 + q / mA)] for q in (f * qA1 for f in lin(0.0, 1.0, 24))]
            return dict(hsub=hsub, mA=mA, raised=raised, qA1=qA1, p1=p1, g1=min_gap(p1), suc=suc, cp=cp, q_ihx=q_ihx)

        sub_off = 0.0
        B = bottom(0.0)
        k = 0
        while k < 25 and c.sub_air and not no_sub_raise and B["g1"] < dT - 0.01 and B["mA"] > 0:
            sub_off += (dT - B["g1"]) * 1.5 + 0.005
            B = bottom(sub_off)
            k += 1
        hsub, mA, raised, qA1, p1, g1 = B["hsub"], B["mA"], B["raised"], B["qA1"], B["p1"], B["g1"]
        suc, cp, q_ihx = B["suc"], B["cp"], B["q_ihx"]
        hv = hsub - q_ihx
        he = expand_to(fb, Pc, hv, Pe, c.expander_eta)
        wx = hv - he
        if not is_gas(cp["out"]):
            return _fail("wet compression in the bottom cycle")
        g_ihx = None
        if q_ihx > 1e-6:
            g_ihx = fb.PH(Pc, hv)["T"] - ev_out["T"]
            if g_ihx < c.dT_ihx - 0.05:
                warnings.append(f"bottom IHX cold-end approach {g_ihx:.1f} K, below its {c.dT_ihx} K")
        if not raised > 0 or not mA > 0:
            return _fail("no steam split balances the condenser")
        qA2 = mA * (ha2 - ha1)
        q_air = qA1 + qA2
        w = cp["w"] + raised * wS
        w_el = expander_electric(cp["w"], wx, motor) + raised * wS / mot_f
        cop = q_air / w_el
        xq = expander_warning(fb, Pe, he, warnings) if c.expander_eta else None
        # approach checks: the bottom cycle (above), then each process steam stream in turn (counterflow within each)
        p2, q_acc, ha_acc = [], qA1, ha1
        for i, f in enumerate(fans):
            Q = raised * E[i] * f["q"]
            if not Q > 1e-9:
                continue
            for fr in lin(0.0, 1.0, 12):
                p2.append([q_acc + fr * Q, water.PH(f["P_out"], f["lf"]["h"] + fr * f["q"])["T"],
                           air.T(ha_acc + fr * Q / mA)])
            q_acc += Q
            ha_acc += Q / mA
        g2 = min_gap(p2)
        if g2 < 0:
            return _fail(f"the process steam cannot heat the air to {c.T_air_tar:g} C (temperature cross): add MVR "
                         "stages, send less exhaust to the process early, or raise the bottom condensing temperature")
        if g2 < dT - 0.05:
            warnings.append(f"steam-to-air approach {g2:.1f} K, below the {dT:g} K minimum: add MVR stages or send less "
                            "exhaust to the process early")
        if g1 < dT - 0.05:
            warnings.append(f"bottom-cycle air approach {g1:.1f} K")
        if sv["P"] < 0.1013:
            warnings.append(f"steam raised at {sv['P'] * 10:.2f} bar (vacuum)")
        q_sg = raised * g
        q_e = ev_out["h"] - he
        m_s = raised * sum(E)
        stages = [dict(P_in=f["P_in"], P_out=f["P_out"], T_out=f["T_out"], T_sat=f["Tsat"], inj=f["inj"], mass=f["M"],
                       eta=f["eta"], split=x[i], to_process=raised * E[i], q_air=raised * E[i] * f["q"])
                  for i, f in enumerate(fans)]

        cond_paths = []
        for i, f in enumerate(fans):
            if not E[i] > 1e-9:
                continue
            cond_paths.append({"name": "Steam condenser" if i == nS - 1 else f"Process steam {i + 1}",
                               "pts": iso_path(water, f["P_out"], f["out"]["h"], f["lf"]["h"])})
            cond_paths.append({"name": "Flash", "pts": [point(f["lf"]), point(water.PH(sv["P"], f["lf"]["h"]))]})
        bottom_paths = [
            {"name": "Expander" if c.expander_eta else "Expansion valve", "pts": comp_path(fb, Pc, hv, Pe, he)},
            {"name": "Evaporator", "pts": iso_path(fb, Pe, he, ev_out["h"])},
        ]
        if q_ihx > 1e-6:
            bottom_paths += [{"name": "Internal heat exchanger", "pts": iso_path(fb, Pe, ev_out["h"], suc["h"])},
                             {"name": "Internal heat exchanger", "pts": iso_path(fb, Pc, hsub, hv)}]
        bottom_paths += [
            {"name": "Compressor", "pts": comp_path(fb, Pe, suc["h"], Pc, cp["out"]["h"])},
            {"name": "Steam generator", "pts": iso_path(fb, Pc, cp["out"]["h"], hsub + qA1)},
            {"name": "Air condenser", "pts": iso_path(fb, Pc, hsub + qA1, hsub)},
        ]
        paths = [
            {"name": "Bottom cycle", "fluid": c.refrigerant, "flow": 1.0, "paths": bottom_paths},
            {"name": "Steam", "fluid": "Water", "flow": m_s,
             "paths": [{"name": "Steam generator", "pts": [point(sl), point(sv)]}, *spaths, *cond_paths]},
        ]
        profiles = [
            {"name": "Bottom cycle to air", "hot": [[p[0], p[1]] for p in p1], "cold": [[p[0], p[2]] for p in p1]},
            {"name": "Steam to air", "hot": [[p[0], p[1]] for p in p2], "cold": [[p[0], p[2]] for p in p2]},
            {"name": "Steam generator",
             "hot": [[f * q_sg, fb.PH(Pc, hsub + qA1 + f * (cp["out"]["h"] - hsub - qA1))["T"]] for f in lin(0, 1, 14)],
             "cold": [[0.0, Ts], [q_sg, Ts]]},
        ]
        states = dict(
            evap_in=_state(fb.PH(Pe, he)), evap_out=_state(ev_out), suction=_state(suc), discharge=_state(cp["out"]),
            cond_bubble=_state(cb), subcool_out=_state(fb.PH(Pc, hsub)), valve_in=_state(fb.PH(Pc, hv)),
            steam_raised=_state(sv), steam_final=_state(stF),
        )
        report = dict(
            COP=cop, P_evap_bar=Pe * 10, P_cond_bar=Pc * 10, PR_bottom=Pc / Pe, eta_bottom=cp["eta"],
            T_dis_bottom=cp["out"]["T"], T_steam=Ts, P_steam_bar=sv["P"] * 10, P_final_bar=Pf * 10,
            T_sat_final=last_f["Tsat"], T_steam_final=st["T"], eta_mvr=self._eta_fan(sv["P"]), PR_mvr=c.PR_mvr,
            fan_motor=mot_f, steam_per_kg=raised, steam_final_per_kg=m_s, w_bottom=cp["w"], w_mvr=raised * wS,
            q_air_bottom=qA1, q_air_steam=qA2, q_sg=q_sg, T_air_mid=T1, gap_bottom=g1, gap_steam=g2, stages=stages,
            w_el_kJkg=w_el, splits=x[:-1], w_el_fans=raised * wS / mot_f, sub_raise=sub_off, t1_drop=t1_drop or 0.0,
            q_ihx=q_ihx, T_suction=suc["T"], gap_ihx=g_ihx, expander=bool(c.expander_eta), w_exp_kJkg=wx,
            exp_out_quality=xq, Pmax_bar=max(Pc, Pf) * 10, T_dis_max=max([cp["out"]["T"]] + [s["T_out"] for s in stages]),
            rho_suction=suc["D"],
        )
        return MVRHPResult(True, "success", cop, q_air / (1 / suc["D"]) / 1000, q_air, w, q_e, report, states, paths,
                           profiles, warnings)

    # ---- optimiser (boosterOptimiseSteps) --------------------------------------------------

    def optimise(self, max_eval: int = 1200, n_max: int = 20) -> MVROptimum:
        """Maximise the COP over the bottom condensing temperature and the number of MVR stages. A design counts
        only if both the bottom-cycle and the steam air approaches hold. Scans every stage count on a 5 K condensing
        grid, then refines the condensing temperature for the best three counts and the current design's; the
        optimum is rounded to 0.1 K as the page shows it."""
        c = self.c
        state = {"nEval": 0, "memo": {}}
        Te = c.T_source_out - c.dT_src
        lo = [max(c.T_air_inlet + c.dT_GC_min + 2, Te + 10), 1]
        T_crit = self.fb.T_crit
        hi = [min(T_crit - 2 if math.isfinite(T_crit) else c.T_air_tar, c.T_air_tar + c.dT_GC_min), n_max]
        if lo[0] >= hi[0]:
            return MVROptimum(False, "no room for the bottom condensing temperature")

        def bad(r):
            return r.report["gap_steam"] < c.dT_GC_min - 0.05 or r.report["gap_bottom"] < c.dT_GC_min - 0.05

        def solve(x):
            m = MVRHeatPump(replace(c, T_cond=x[0], n_mvr=int(x[1])), self.fb, self.water, self.air)
            r = m.solve()
            if r.success:
                r.bad = bad(r)
            return r

        f = memo_objective(state, solve)

        def clip(x):
            return [min(max(x[0], lo[0]), hi[0]), js_round(min(max(x[1], lo[1]), hi[1]))]

        x0 = clip([c.T_cond, c.n_mvr])
        b0 = f(x0)
        start = b0[1].COP if b0[1] is not None and b0[1].success else float("nan")
        start_ok = math.isfinite(b0[0])
        best_by_n: dict = {}
        for m in range(lo[1], hi[1] + 1):
            T = lo[0]
            while T <= hi[0] + 1e-9:
                z = clip([T, m])
                o = f(z)
                cur = best_by_n.get(m)
                if cur is None or o[0] > cur[1][0]:
                    best_by_n[m] = (z, o)
                T += 5
        picks = sorted((v for v in best_by_n.values() if math.isfinite(v[1][0])), key=lambda v: -v[1][0])[:3]
        if math.isfinite(b0[0]) and not any(v[0][1] == x0[1] for v in picks):
            picks.append((x0, b0))
        X, Bst, conv = None, (-math.inf, None), True
        for z, o in picks:
            m = z[1]

            def f1(zz, m=m):
                return f([zz[0], m])

            starts = [([z[0]], o)]
            if z[1] == x0[1] and z is not x0 and math.isfinite(b0[0]):
                starts.append(([x0[0]], b0))
            xs, bs, cv = compass_search(f1, starts, [lo[0]], [hi[0]], [4], [0.05], [False], state, max_eval)
            if bs[0] > Bst[0]:
                Bst, X, conv = bs, [xs[0], m], cv
        if X is None:
            return MVROptimum(False, "no design within the bounds keeps both air approaches", n_eval=state["nEval"],
                              start=start, start_ok=start_ok)
        x, best = X, Bst
        xr = [js_round(x[0] * 10) / 10, x[1]]
        for _ in range(20):
            o = f(xr)
            if math.isfinite(o[0]):
                x, best = xr, o
                break
            xr = [js_round((xr[0] + 0.1) * 10) / 10, xr[1]]

        def near(v, b):
            return abs(v - b) < 0.05

        return MVROptimum(True, "success", best[0], x[0], x[1], best[1].report["splits"], start, start_ok, state["nEval"],
                          conv, {"T_cond": near(x[0], lo[0]) or near(x[0], hi[0]), "n_mvr": x[1] == hi[1]},
                          (tuple(lo), tuple(hi)), best[1])
