"""Dual-pressure condensation heat pump for air drying (after Dai et al. 2020): one evaporator, two condensing
pressures in series on the air."""

from .dual_pressure_hp import DualPressureConfig, DualPressureHeatPump, DualPressureResult

__all__ = ["DualPressureConfig", "DualPressureHeatPump", "DualPressureResult"]
