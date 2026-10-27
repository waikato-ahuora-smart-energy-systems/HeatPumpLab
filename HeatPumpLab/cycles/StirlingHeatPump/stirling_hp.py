"""Reverse-Stirling heat pumps heating dryer air: a semi-empirical block, not a physical Stirling model.

Industrial reverse-Stirling heat pumps (helium, about 400-720 kW: Olvondo HighLift, Enerin HoegTemp) deliver heat at
150-200 C, but their published performance is thin: Hoeg et al. (2026, IEA Heat Pump Conference) frame it as a
fraction of the Carnot COP against the temperature ratio, from about 15,000 h of field data at 134-199 C sinks, without
tables; Haggqvist, Tveit and Zevenhoven (Case Studies in Thermal Engineering 47 (2023) 103066) give a measured
500 kW unit for model validation; and the supplier-reported RFC Maasdam dairy case gives COP about 2 for 155 C steam
from waste heat at about 75 C. A physical (adiabatic or second-order) Stirling model would need machine geometry, so
this block uses the published framing instead:

    COP_unit = eta_carnot * T_hot / (T_hot - T_cold)            (temperatures in K)

with eta_carnot a constant (default 0.37, from the RFC Maasdam point: COP 2 at about 428 K / 348 K) or a function of
the temperature ratio T_hot / T_cold supplied by the user (for example a fit to Hoeg et al. 2026 once its data are
available). Treat results as indicative.

A Stirling unit takes and gives heat at nearly constant temperature, so with a gliding air sink each unit's hot side
sits ``dT_hot`` above the air temperature it delivers (its pinch, through an intermediate water or steam loop), and its
cold side ``dT_cold`` below the source outlet. ``n_units`` units in series on the air each heat one band of it; their
boundaries are equal temperature steps, or set by ``boundaries``, and ``optimise()`` searches them for the best COP.

Duties are per kg of dry air: q_air is the air's enthalpy rise, W the electrical input (COP is heat to air over
electrical input, as the unit COPs are reported), q_source = q_air - W.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Callable

from ..common.dryer import AirModel, LoopResult, air_secondary, compass_search, dry_air, humid_air, js_round

StirlingResult = LoopResult
DEMONSTRATED_T_HOT = 200.0          # C, about the highest sink temperature in the published field data


@dataclass
class StirlingConfig:
    """Design inputs. Temperatures in C, approaches in K."""

    n_units: int = 1                   # units in series on the air, 1-4
    boundaries: list | None = None     # air temperatures between units (n_units - 1 of them); None: equal steps
    eta_carnot: float = 0.37           # fraction of the Carnot COP (used when eta_of_ratio is None)
    eta_of_ratio: Callable[[float], float] | None = field(default=None, repr=False)   # eta(T_hot / T_cold), K / K
    dT_hot: float = 20.0               # unit hot side above the air it delivers (through the intermediate loop)
    dT_cold: float = 5.0               # unit cold side below the source water outlet
    T_air_inlet: float = 15.0
    T_air_tar: float = 170.0
    T_source_in: float = 60.0
    T_source_out: float = 50.0
    air_ratio: float = 0.0             # secondary / primary air (heated with the primary from the inlet to T_mt)
    T_mt: float = 102.0
    air_RH: float = 60.0
    P_air_kPa: float = 101.325
    property_backend: str | None = None   # only for dry air (air_RH = 0)


def _fail(message: str, **extra) -> LoopResult:
    return LoopResult(False, message, **extra)


class StirlingHeatPump:
    """Reverse-Stirling heat pump units heating dryer air; see the module docstring."""

    def __init__(self, config: StirlingConfig | None = None, air: AirModel | None = None):
        self.c = config if config is not None else StirlingConfig()
        c = self.c
        if air is None:
            air = humid_air(c.T_air_inlet, c.air_RH, c.P_air_kPa) if c.air_RH > 0 else \
                dry_air(c.P_air_kPa, c.property_backend)
            air = air_secondary(air, c.T_air_inlet, c.T_mt, c.air_ratio, c.T_air_tar)
        self.air = air

    def eta(self, c: StirlingConfig, T_hot_K: float, T_cold_K: float) -> float:
        return c.eta_of_ratio(T_hot_K / T_cold_K) if c.eta_of_ratio is not None else c.eta_carnot

    def solve(self, **overrides) -> LoopResult:
        """The units at the configuration (fields can be overridden by keyword, except the air)."""
        c = replace(self.c, **overrides) if overrides else self.c
        return self._solve(c)

    def optimise(self, max_eval: int = 200, **overrides) -> LoopResult:
        """Best-COP air temperatures between the units (compass search from equal steps), rounded to 0.1 K."""
        c = replace(self.c, **overrides) if overrides else self.c
        n = c.n_units
        if n <= 1:
            return self._solve(c)
        state = dict(nEval=0)
        lo = [c.T_air_inlet + 1.0] * (n - 1)
        hi = [c.T_air_tar - 1.0] * (n - 1)

        def f(x):
            state["nEval"] += 1
            if any(b <= a for a, b in zip(x, x[1:])):
                return dict(y=-math.inf, r=None)
            r = self._solve(replace(c, boundaries=[float(v) for v in x]))
            return dict(y=r.COP if r.success else -math.inf, r=r)

        x0 = self._boundaries(c)
        b0 = f(x0)
        x, best, converged = compass_search(f, [(x0, b0)], lo, hi, [10.0] * (n - 1), [0.05] * (n - 1),
                                            [False] * (n - 1), state, max_eval, score=lambda o: o["y"],
                                            empty={"y": -math.inf})
        if x is None:
            return _fail("no set of unit boundaries gives a feasible design")
        xr = [js_round(v, 1) for v in x]
        o = f(xr)
        if math.isfinite(o["y"]):
            x, best = xr, o
        final = best["r"]
        final.report["optimiser"] = dict(COP=best["y"], start=b0["y"], boundaries=[float(v) for v in x],
                                         nEval=state["nEval"], converged=converged)
        return final

    # ---- model -------------------------------------------------------------------------

    @staticmethod
    def _boundaries(c: StirlingConfig) -> list:
        n = c.n_units
        if c.boundaries is not None:
            return [float(v) for v in c.boundaries]
        return [c.T_air_inlet + (c.T_air_tar - c.T_air_inlet) * k / n for k in range(1, n)]

    def _solve(self, c: StirlingConfig) -> LoopResult:
        air, warnings = self.air, []
        n = c.n_units
        if not 1 <= n <= 4:
            return _fail("1 to 4 Stirling units in series")
        bounds = self._boundaries(c)
        if len(bounds) != n - 1:
            return _fail(f"{n} units need {n - 1} boundary temperatures")
        T = [c.T_air_inlet] + bounds + [c.T_air_tar]
        if any(b <= a for a, b in zip(T, T[1:])):
            return _fail("the unit boundaries must rise from the air inlet to the target")
        T_cold = (c.T_source_out if c.T_source_out is not None else c.T_source_in) - c.dT_cold
        units, q_air, W = [], 0.0, 0.0
        for k in range(n):
            T_hot = T[k + 1] + c.dT_hot
            Th, Tc = T_hot + 273.15, T_cold + 273.15
            if Th <= Tc:
                return _fail(f"unit {k + 1}: its hot side ({T_hot:.0f} C) is not above its cold side ({T_cold:.0f} C)")
            eta = self.eta(c, Th, Tc)
            if not 0 < eta <= 1:
                return _fail(f"unit {k + 1}: Carnot fraction {eta} outside (0, 1]")
            carnot = Th / (Th - Tc)
            cop = eta * carnot
            q = air.h(T[k + 1]) - air.h(T[k])
            w = q / cop
            q_air += q
            W += w
            units.append(dict(T_air_in=T[k], T_air_out=T[k + 1], T_hot=T_hot, T_cold=T_cold, carnot_COP=carnot,
                              eta_carnot=eta, COP=cop, q_air=q, w_el=w, q_source=q - w, lift=T_hot - T_cold))
            if T_hot > DEMONSTRATED_T_HOT:
                warnings.append(f"unit {k + 1} hot side at {T_hot:.0f} C, above the about {DEMONSTRATED_T_HOT:.0f} C "
                                "of the published field data")
        cop = q_air / W
        report = dict(COP=cop, units=units, n_units=n, boundaries=bounds, T_cold=T_cold, w_el_kJkg=W,
                      q_air=q_air, q_source=q_air - W, air_W=air.W,
                      eta_model="eta_of_ratio" if c.eta_of_ratio is not None else f"constant {c.eta_carnot:g}",
                      balance=0.0)
        exchangers = [dict(name=f"Unit {k + 1} to air", hot=[[0.0, u["T_hot"]], [u["q_air"], u["T_hot"]]],
                           cold=[[0.0, u["T_air_in"]], [u["q_air"], u["T_air_out"]]], other="dryer air",
                           gap=u["T_hot"] - u["T_air_out"]) for k, u in enumerate(units)]
        return LoopResult(True, "success", COP=cop, q_air=q_air, w_shaft=W, q_source=q_air - W,
                          T_dis_max=max(u["T_hot"] for u in units), report=report, exchangers=exchangers,
                          warnings=warnings)
