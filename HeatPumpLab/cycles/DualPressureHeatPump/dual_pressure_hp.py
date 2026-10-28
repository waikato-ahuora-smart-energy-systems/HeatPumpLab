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

This is MultiPressureHeatPump with n_pressures = 2; the inputs and the report keep the two-pressure names
(T_cond_low, T_mid, split_high, ...).

Units: MPa, C, kJ/kg, per kg of evaporator flow.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields

from ..common.dryer import DEFAULT_MOTOR_EFF, AirModel, Fluid, LoopResult
from ..MultiPressureHeatPump.multi_pressure_hp import MultiPressureConfig, MultiPressureHeatPump

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


def _multi(c: DualPressureConfig) -> MultiPressureConfig:
    d = asdict(c)
    T = d.pop("T_cond_low")
    return MultiPressureConfig(n_pressures=2, T_cond=[T], **d)


_DUAL_FIELDS = {f.name for f in fields(DualPressureConfig)}


def _overrides(kw: dict) -> dict:
    kw = dict(kw)
    bad = set(kw) - _DUAL_FIELDS - {"max_eval"}
    if bad:
        raise TypeError(f"unknown inputs: {', '.join(sorted(bad))}")
    if "T_cond_low" in kw:
        kw["T_cond"] = [kw.pop("T_cond_low")]
    return kw


def _legacy(r: LoopResult) -> LoopResult:
    """Add the two-pressure report names."""
    if not r.success:
        return r
    p = r.report
    q, g, m = p["q_air_levels"], p["gaps"], p["flows"]
    p.update(P_low_cond_bar=p["P_cond_bar"][0], T_cond_low=p["T_cond"][0], T_mid=p["T_air_between"][0],
             split_high=m[1], q_air_low=q[0], q_air_high=q[1], gap_low=g[0], gap_high=g[1],
             T_dis_low=p["T_dis"][0], T_dis_high=p["T_dis"][1])
    if "optimiser" in p and "T_cond" in p["optimiser"]:
        p["optimiser"]["T_cond_low"] = p["optimiser"]["T_cond"][0]
    return r


class DualPressureHeatPump(MultiPressureHeatPump):
    """Dual-pressure condensation heat pump; see the module docstring."""

    def __init__(self, config: DualPressureConfig | None = None, fluid: Fluid | None = None, air: AirModel | None = None):
        self.dual = config if config is not None else DualPressureConfig()
        super().__init__(_multi(self.dual), fluid, air)

    def solve(self, **overrides) -> LoopResult:
        """The design at the configuration (fields can be overridden by keyword, except the fluid and air)."""
        return _legacy(super().solve(**_overrides(overrides)))

    def optimise(self, max_eval: int = 120, **overrides) -> LoopResult:
        """Best-COP low condensing temperature (compass search from the configured value and a coarse grid), rounded
        to 0.1 K. The result is the full solution there, with the search in report["optimiser"]."""
        return _legacy(super().optimise(max_eval, **_overrides(overrides)))
