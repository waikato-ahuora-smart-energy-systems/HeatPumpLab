"""Split-stream heat pump: the suction-gas recuperator (RHX1) runs in parallel with the
low-temperature gas cooler instead of in series with it (Design A).

Flow path (per kg of refrigerant)::

    evaporator (Tref_evap_inlet -> dew) -> superheater (to T_source_in - dTmin_rhx)
      -> RHX1 cold side (lift dT_rhx1_cold) -> compressor stage 1 -> stage 2 (no intercooling)
      -> high-temperature gas cooler (air T_mid -> T_air_tar)
      -> split at h_split:  fraction 1 - a -> low-temperature gas cooler, cooled as far as practical,
                               to T_air_inlet + dT_GC_min (air T_air_inlet -> T_mid)
                            fraction a     -> RHX1 hot side, cooled to T_superheat + dT_rhx1_cold_end
      -> branches mix (different temperatures allowed) -> expansion valve -> evaporator

dT_rhx1_cold_end sets the RHX1 size: the parallel arrangement gives nearly balanced RHX1
streams, so the LMTD is close to this cold-end difference. Design A's RHX1 runs at an LMTD of
about 19-29 K (5 K at its hot end, 46-88 K at its cold end), so 25 K is the like-for-like default.

Closure. For a gas-cooler pressure P_H and RHX1 branch fraction a, both branch outlet states are
known, so the mixed valve inlet fixes the evaporator pressure P_L (refrigerant enters the
evaporator at Tref_evap_inlet). The suction state, RHX1 lift and compressor power split fix the
compression. The RHX1 duty then fixes the split enthalpy, h_split = h_rhx_out + q_rhx1 / a, and
the energy balances fix the air flow and the air temperature between the gas coolers. P_H and a
are optimised for COP subject to the minimum approaches (gas coolers dT_GC_min, RHX1 dTmin_rhx
anywhere along it) and an optional minimum RHX1 LMTD.

The model is independent of the RecuperativeHP class; it shares the property package and the
compressor efficiency correlations. Units follow HeatPumpLab: MPa, C, kJ/kg, kg/m3.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, asdict
from typing import Callable

import numpy as np

from ...components.compressor import pareto_front_efficiency
from ...utils.properties_package import COMPRESSOR_GAS_PHASES, RefpropError
from ...utils.coolprop_properties import make_property_package
from ...utils.property_cache import active_property_call_cache, property_call_cache

# ----------------------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------------------


@dataclass
class SplitHPConfig:
    """Design inputs. Temperatures in C, pressures in MPa, approaches in K."""

    refrigerant: str
    T_air_inlet: float = 15.0
    T_air_tar: float = 200.0
    T_source_in: float = 40.0
    Tref_evap_inlet: float = 10.0
    dT_GC_min: float = 20.0          # minimum approach in both gas coolers (and their cold ends)
    dTmin_rhx: float = 5.0           # minimum approach in RHX1 and the source superheater
    dT_rhx1_cold: float = 80.0       # RHX1 cold-side temperature lift (design variable)
    comp2_power_frac: float = 0.5    # stage-2 share of compressor shaft power (design variable)
    n_comp: int = 2
    comp_model_index: int = 5        # 0 ideal, 1 constant 0.729, 4 Pareto front v1, 5 Pareto front v2
    motor_eff: float = 0.93
    dT_rhx1_cold_end: float = 25.0   # RHX1 hot-branch outlet minus suction-gas inlet (sets RHX1 size)
    rhx1_lmtd_min: float | None = None   # optional minimum RHX1 LMTD (K)
    n_alpha_scan: int = 12           # RHX1 branch-fraction scan points per gas-cooler pressure
    P_air: float = 0.101             # MPa
    P_H_bounds: tuple[float, float] = (0.5, 6.0)   # MPa, search range for the gas-cooler pressure
    Pmax_cap: float | None = None    # MPa, optional hard cap on P_H
    n_profile: int = 41              # points per exchanger profile check
    n_PH_scan: int = 28              # coarse gas-cooler pressure scan points


# ----------------------------------------------------------------------------------------
# Property access (injectable, so the cycle logic can be tested without REFPROP)
# ----------------------------------------------------------------------------------------


class RefpropProperties:
    """Thin wrapper over HeatPumpLab's property package (CoolProp by default) for one refrigerant and air."""

    def __init__(self, refrigerant: str, sink: str = "Air"):
        self.ref = make_property_package(refrigerant)
        self.air = make_property_package(sink, "snk")

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

    def air_h(self, P, T):
        return self.air.get_std_properties("PT", P, T, self.air.fld_composition)["h"]

    def air_T(self, P, h):
        return self.air.get_std_properties("PH", P, h, self.air.fld_composition)["T"]


# ----------------------------------------------------------------------------------------
# Compressor helpers (same correlations as components/compressor.py)
# ----------------------------------------------------------------------------------------


def isentropic_efficiency(pressure_ratio: float, model_index: int, P_in: float | None = None) -> float:
    pr = max(float(pressure_ratio), 1.0)
    if model_index == 0:
        eff = 1.0
    elif model_index == 1:
        eff = 0.7 / 0.96
    elif model_index == 2:
        eff = (0.65 + 0.015 * pr - 0.0015 * pr ** 2) / 0.96
    elif model_index == 4:
        eff = pareto_front_efficiency(pr, "isentropic")
    elif model_index == 5:
        eff = pareto_front_efficiency(pr, "isentropic", version=2)
    else:
        raise ValueError(f"compressor model {model_index} is not supported by SplitHeatPump")
    return min(max(eff, 0.4), 1.0)


def volumetric_efficiency(pressure_ratio: float) -> float:
    """Brendel et al. correlation, as in components/compressor.py."""
    return 1 - 0.08244 * pow(max(pressure_ratio, 1.0) - 1, 0.72773)


# ----------------------------------------------------------------------------------------
# Pure exchanger-profile helpers (no property calls)
# ----------------------------------------------------------------------------------------


def counterflow_min_approach(hot_T: Callable[[np.ndarray], np.ndarray], cold_T: Callable[[np.ndarray], np.ndarray],
                             n: int = 41) -> float:
    """Minimum (T_hot - T_cold) along a counterflow exchanger.

    ``hot_T(f)`` and ``cold_T(f)`` give temperatures at duty fraction f in [0, 1], where f = 0 is the
    cold end (hot-side outlet, cold-side inlet) and f = 1 the hot end.
    """
    f = np.linspace(0.0, 1.0, n)
    return float(np.min(hot_T(f) - cold_T(f)))


class _Table:
    """Monotone T(h) table along an isobar, built once and interpolated."""

    def __init__(self, h: np.ndarray, T: np.ndarray):
        order = np.argsort(h)
        self.h, self.T = np.asarray(h)[order], np.asarray(T)[order]

    def __call__(self, h):
        return np.interp(h, self.h, self.T)


# ----------------------------------------------------------------------------------------
# Cycle
# ----------------------------------------------------------------------------------------


@dataclass
class SplitHPResult:
    success: bool
    message: str
    COP: float | None = None
    VHC: float | None = None              # MJ/m3, heat to air per total suction volume of the stages
    P_H: float | None = None              # MPa
    P_I: float | None = None
    P_L: float | None = None
    report: dict = field(default_factory=dict)
    states: dict = field(default_factory=dict)
    # coarse gas-cooler pressure scan behind the optimum: [P_H MPa, COP, VHC, T_comp_discharge C, Pmax bar]
    # (COP None where infeasible); lets a pressure limit be applied after the fact
    scan: list = field(default_factory=list)

    def to_dict(self):
        return asdict(self)


class SplitHeatPump:
    """Split-stream heat pump model; see module docstring."""

    def __init__(self, config: SplitHPConfig, props=None):
        self.c = config
        self.props = props if props is not None else RefpropProperties(config.refrigerant)

    # ---- building blocks ------------------------------------------------------------

    @property
    def T_superheat(self) -> float:
        return self.c.T_source_in - self.c.dTmin_rhx

    @property
    def T_LT_out(self) -> float:
        return self.c.T_air_inlet + self.c.dT_GC_min

    @property
    def T_RX_out(self) -> float:
        return self.T_superheat + self.c.dT_rhx1_cold_end

    def _evaporator_pressure(self, h_valve: float, P_H: float) -> float:
        """P_L such that the valve outlet (isenthalpic, h_valve) sits at Tref_evap_inlet."""
        target = self.c.Tref_evap_inlet
        lo, hi = 0.02, 0.98 * P_H

        def T_at(P):
            return self.props.ref_PH(P, h_valve)["T"]

        T_lo, T_hi = T_at(lo), T_at(hi)
        if not T_lo < target < T_hi:
            raise ValueError(f"evaporator inlet {target} C not bracketed ({T_lo:.1f} .. {T_hi:.1f} C)")
        for _ in range(50):
            mid = math.sqrt(lo * hi)
            if T_at(mid) < target:
                lo = mid
            else:
                hi = mid
            if hi / lo - 1 < 2e-6:
                break
        return math.sqrt(lo * hi)

    def _stage(self, P_in: float, s_in: float, h_in: float, P_out: float):
        eta = isentropic_efficiency(P_out / P_in, self.c.comp_model_index, P_in)
        h_is = self.props.ref_PS(P_out, s_in)["h"]
        w = (h_is - h_in) / eta
        out = self.props.ref_PH(P_out, h_in + w)
        return w, eta, out

    def _compress(self, suction: dict, P_L: float, P_H: float):
        """Two stages (or one) with the stage-2 shaft-power share comp2_power_frac."""
        c = self.c
        if c.n_comp == 1 or c.comp2_power_frac <= 0.0:
            w, eta, out = self._stage(P_L, suction["s"], suction["h"], P_H)
            return dict(P_I=None, w1=w, w2=0.0, eta1=eta, eta2=None, out1=out, in2=None, out2=out)
        target = c.comp2_power_frac
        lo, hi = math.log(P_L * 1.0001), math.log(P_H / 1.0001)

        def evaluate(lnP):
            P_I = math.exp(lnP)
            w1, e1, o1 = self._stage(P_L, suction["s"], suction["h"], P_I)
            w2, e2, o2 = self._stage(P_I, o1["s"], o1["h"], P_H)
            return w2 / (w1 + w2), (P_I, w1, w2, e1, e2, o1, o2)

        frac_lo, _ = evaluate(lo)
        frac_hi, _ = evaluate(hi)
        if not (frac_hi <= target <= frac_lo):
            raise ValueError("compressor power split not reachable")
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            frac, pack = evaluate(mid)
            if frac > target:
                lo = mid
            else:
                hi = mid
            if abs(frac - target) < 2e-4:
                break
        P_I, w1, w2, e1, e2, o1, o2 = pack
        return dict(P_I=P_I, w1=w1, w2=w2, eta1=e1, eta2=e2, out1=o1, in2=o1, out2=o2)

    # ---- one design point (gas-cooler pressure, RHX1 branch fraction) ----------------

    def evaluate(self, P_H: float, alpha: float) -> SplitHPResult:
        """alpha = mass fraction of refrigerant sent through the RHX1 hot side."""
        c, pr = self.c, self.props
        try:
            if c.Pmax_cap is not None and P_H > c.Pmax_cap + 1e-12:
                return SplitHPResult(False, "above pressure cap", P_H=P_H)
            if c.dT_rhx1_cold <= 0.0:
                alpha = 0.0
            elif not 0.0 < alpha < 1.0:
                return SplitHPResult(False, "RHX1 branch fraction outside (0, 1)", P_H=P_H)
            lt = pr.ref_PT(P_H, self.T_LT_out)
            rx = pr.ref_PT(P_H, self.T_RX_out) if alpha > 0 else lt
            h_mix = (1 - alpha) * lt["h"] + alpha * rx["h"]
            P_L = self._evaporator_pressure(h_mix, P_H)
            dew = pr.ref_PQ(P_L, 1.0)
            if dew["T"] > self.T_superheat + 1e-9:
                return SplitHPResult(False, f"evaporator dew {dew['T']:.1f} C above superheater outlet "
                                            f"{self.T_superheat:.1f} C", P_H=P_H, P_L=P_L)
            sup = pr.ref_PT(P_L, self.T_superheat)
            T_suc = self.T_superheat + c.dT_rhx1_cold
            suction = pr.ref_PT(P_L, T_suc) if c.dT_rhx1_cold > 0 else sup
            if suction.get("phase") not in COMPRESSOR_GAS_PHASES:
                return SplitHPResult(False, "suction not gas", P_H=P_H, P_L=P_L)
            comp = self._compress(suction, P_L, P_H)
            for key in ("out1", "out2"):
                if comp[key].get("phase") not in COMPRESSOR_GAS_PHASES:
                    return SplitHPResult(False, f"wet compression ({key})", P_H=P_H, P_L=P_L)
            h_dis, T_dis = comp["out2"]["h"], comp["out2"]["T"]
            q_rhx = suction["h"] - sup["h"]
            w = comp["w1"] + comp["w2"]
            if alpha > 0:
                h_split = rx["h"] + q_rhx / alpha
                if h_split > h_dis:
                    return SplitHPResult(False, "RHX1 duty exceeds the branch (split above discharge)", P_H=P_H, P_L=P_L)
                if h_split < lt["h"]:
                    return SplitHPResult(False, "split below the low-temperature gas-cooler outlet", P_H=P_H, P_L=P_L)
            else:
                h_split = None
            q_air = h_dis - h_mix - q_rhx
            ha_in, ha_tar = pr.air_h(c.P_air, c.T_air_inlet), pr.air_h(c.P_air, c.T_air_tar)
            r = q_air / (ha_tar - ha_in)
            q_LT = (1 - alpha) * ((h_split if h_split is not None else h_dis) - lt["h"])
            q_HT = q_air - q_LT
            ha_mid = ha_in + q_LT / r

            # Exchanger profiles
            n = c.n_profile
            h_lo = min(lt["h"], rx["h"])
            h_tab = np.linspace(h_lo, h_dis, 90)
            hp_tab = _Table(h_tab, np.array([pr.ref_PH(P_H, h)["T"] for h in h_tab]))
            T_air_grid = np.linspace(c.T_air_inlet, c.T_air_tar, 60)
            air_tab = _Table(np.array([pr.air_h(c.P_air, T) for T in T_air_grid]), T_air_grid)
            h_top = h_split if h_split is not None else h_dis
            approaches = {}
            approaches["LT_gas_cooler"] = counterflow_min_approach(
                lambda f: hp_tab(lt["h"] + f * (h_top - lt["h"])),
                lambda f: air_tab(ha_in + f * (ha_mid - ha_in)), n)
            if h_split is not None:
                approaches["HT_gas_cooler"] = counterflow_min_approach(
                    lambda f: hp_tab(h_split + f * (h_dis - h_split)),
                    lambda f: air_tab(ha_mid + f * (ha_tar - ha_mid)), n)
                h_lp = np.linspace(sup["h"], suction["h"], 30)
                lp_tab = _Table(h_lp, np.array([pr.ref_PH(P_L, h)["T"] for h in h_lp]))
                approaches["RHX1"] = counterflow_min_approach(
                    lambda f: hp_tab(rx["h"] + f * (h_split - rx["h"])),
                    lambda f: lp_tab(sup["h"] + f * (suction["h"] - sup["h"])), n)
                d_hot = float(hp_tab(h_split)) - T_suc
                d_cold = self.T_RX_out - self.T_superheat
                lmtd = d_hot if abs(d_hot - d_cold) < 1e-6 else (d_hot - d_cold) / math.log(d_hot / d_cold) \
                    if d_hot > 0 and d_cold > 0 else float("nan")
            else:
                lmtd = None
            margins = [approaches["LT_gas_cooler"] - c.dT_GC_min]
            if "HT_gas_cooler" in approaches:
                margins += [approaches["HT_gas_cooler"] - c.dT_GC_min, approaches["RHX1"] - c.dTmin_rhx]
                if c.rhx1_lmtd_min is not None:
                    margins.append((lmtd if lmtd == lmtd else -1e9) - c.rhx1_lmtd_min)
            margin = min(margins)
            if margin < -1e-6:
                return SplitHPResult(False, f"approach violated by {-margin:.2f} K", P_H=P_H, P_L=P_L,
                                     report=dict(approaches_K=approaches))

            v1 = 1.0 / suction["D"]
            v2 = 1.0 / comp["in2"]["D"] if comp["in2"] is not None else 0.0
            cop = q_air / (w / c.motor_eff)
            vhc = q_air / (v1 + v2) / 1000.0
            P_I = comp["P_I"]
            pr1 = (P_I if P_I else P_H) / P_L
            pr2 = P_H / P_I if P_I else None
            T_split = float(hp_tab(h_split)) if h_split is not None else None
            report = dict(
                COP=cop, VHC=vhc, SEC=1 / cop,
                Pmax_bar=P_H * 10, P_evap_bar=P_L * 10, P_int_bar=(P_I * 10 if P_I else None),
                PR_1=pr1, PR_2=pr2, Overall_PR=P_H / P_L,
                eta_is_1=comp["eta1"], eta_is_2=comp["eta2"],
                n_vol_1=volumetric_efficiency(pr1), n_vol_2=(volumetric_efficiency(pr2) if pr2 else None),
                T_evap_inlet=c.Tref_evap_inlet, T_dew_evap=dew["T"], T_superheat_out=self.T_superheat,
                T_suction=T_suc, Comp_1_T_out=comp["out1"]["T"], T_comp_discharge=T_dis,
                T_split=T_split, T_LT_gc_out=self.T_LT_out, T_RHX1_hot_out=(self.T_RX_out if alpha > 0 else None),
                T_valve_in=pr.ref_PH(P_H, h_mix)["T"], T_air_mid=float(air_tab(ha_mid)),
                alpha_RHX1=alpha, alpha_LT_gas_cooler=1 - alpha,
                q_air_kJkg=q_air, q_rhx1_kJkg=q_rhx, w_shaft_kJkg=w, q_evap_kJkg=sup["h"] - h_mix,
                q_HT_gc_kJkg=q_HT, q_LT_gc_kJkg=q_LT, air_per_ref_kgkg=r,
                rhx1_LMTD_K=lmtd, min_approach_margin_K=margin, approaches_K=approaches,
                v_suction_1=v1, v_suction_2=v2,
            )
            states = dict(
                LT_gc_out=dict(P=P_H, T=self.T_LT_out, h=lt["h"]),
                RHX1_hot_out=dict(P=P_H, T=self.T_RX_out, h=rx["h"]) if alpha > 0 else None,
                valve_in=dict(P=P_H, h=h_mix), evap_in=dict(P=P_L, T=c.Tref_evap_inlet, h=h_mix),
                evap_dew=dict(P=P_L, T=dew["T"], h=dew["h"]),
                superheat_out=dict(P=P_L, T=self.T_superheat, h=sup["h"]),
                suction=dict(P=P_L, T=T_suc, h=suction["h"], s=suction["s"]),
                stage1_out=dict(P=comp["out1"]["P"], T=comp["out1"]["T"], h=comp["out1"]["h"]),
                discharge=dict(P=P_H, T=T_dis, h=h_dis),
                split=dict(P=P_H, T=T_split, h=h_split) if h_split is not None else None,
            )
            return SplitHPResult(True, "success", cop, vhc, P_H, P_I, P_L, report, states)
        except (RefpropError, ValueError, ZeroDivisionError) as exc:
            return SplitHPResult(False, f"{type(exc).__name__}: {exc}", P_H=P_H)

    # ---- optimisation over (P_H, alpha) ---------------------------------------------

    @staticmethod
    def _golden(f, a, b, iters, start=None):
        g = (math.sqrt(5) - 1) / 2
        x1, x2 = b - g * (b - a), a + g * (b - a)
        r1, r2 = f(x1), f(x2)
        score = lambda r: r.COP if r.success else -1.0
        best = start
        for _ in range(iters):
            if score(r1) > score(r2):
                b, x2, r2 = x2, x1, r1
                x1 = b - g * (b - a)
                r1 = f(x1)
            else:
                a, x1, r1 = x1, x2, r2
                x2 = a + g * (b - a)
                r2 = f(x2)
        for cand in (r1, r2):
            if cand.success and (best is None or not best.success or cand.COP > best.COP):
                best = cand
        return best

    def _best_alpha(self, P_H: float) -> SplitHPResult:
        if self.c.dT_rhx1_cold <= 0.0:
            return self.evaluate(P_H, 0.0)
        grid = np.linspace(0.04, 0.96, self.c.n_alpha_scan)
        res = [self.evaluate(P_H, float(a)) for a in grid]
        ok = [(r.COP, i) for i, r in enumerate(res) if r.success]
        if not ok:
            # report the most informative failure
            return res[len(res) // 2]
        _, i = max(ok)
        lo, hi = float(grid[max(i - 1, 0)]), float(grid[min(i + 1, len(grid) - 1)])
        best = self._golden(lambda a: self.evaluate(P_H, a), lo, hi, 10, res[i])
        return best

    def solve(self) -> SplitHPResult:
        """Best-COP feasible gas-cooler pressure and RHX1 branch fraction."""
        if active_property_call_cache() is None:
            with property_call_cache():
                return self._solve()
        return self._solve()

    def _solve(self) -> SplitHPResult:
        c = self.c
        lo, hi = c.P_H_bounds
        if c.Pmax_cap is not None:
            hi = min(hi, c.Pmax_cap)
        grid = np.geomspace(lo, hi, c.n_PH_scan)
        results = [self._best_alpha(float(P)) for P in grid]
        scan = [[float(P), r.COP if r.success else None, r.VHC if r.success else None,
                 r.report.get('T_comp_discharge') if r.success else None,
                 r.report.get('Pmax_bar') if r.success else None] for P, r in zip(grid, results)]
        feasible = [(res.COP, i) for i, res in enumerate(results) if res.success]
        if not feasible:
            reasons = {}
            for res in results:
                key = res.message.split(':')[0][:60]
                reasons[key] = reasons.get(key, 0) + 1
            return SplitHPResult(False, f"no feasible design in {lo}-{hi} MPa: {reasons}", scan=scan)
        _, i = max(feasible)
        a = float(grid[max(i - 1, 0)])
        b = float(grid[min(i + 1, len(grid) - 1)])
        best = self._golden(self._best_alpha, a, b, 12, results[i])
        best.scan = scan
        return best
