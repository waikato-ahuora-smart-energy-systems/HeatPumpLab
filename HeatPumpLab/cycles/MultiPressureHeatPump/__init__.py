"""Multi-pressure condensation heat pump for air drying: one evaporator and n condensing pressures in series on the
air (the dual-pressure cycle of Dai et al. 2020, generalised)."""

from .multi_pressure_hp import MAX_PRESSURES, MultiPressureConfig, MultiPressureHeatPump, MultiPressureResult

__all__ = ["MAX_PRESSURES", "MultiPressureConfig", "MultiPressureHeatPump", "MultiPressureResult"]
