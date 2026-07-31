"""Expose common component base classes and physical-domain errors."""

from .base_heat_exchanger import BaseHeatExchanger, TemperatureCrossError

__all__ = ["BaseHeatExchanger", "TemperatureCrossError"]
