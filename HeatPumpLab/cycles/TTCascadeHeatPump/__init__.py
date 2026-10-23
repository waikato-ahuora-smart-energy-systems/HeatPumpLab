"""Transcritical-transcritical CO2 cascade heat pump for air drying (Kong 2026, PhD thesis, University of Waikato,
cycle TT-4 with the additional gas cooler and secondary air); ported from the Advanced dryer heat pumps web app."""

from .tt_cascade import TTCascadeConfig, TTCascadeHeatPump, TTCascadeResult

__all__ = ["TTCascadeConfig", "TTCascadeHeatPump", "TTCascadeResult"]
