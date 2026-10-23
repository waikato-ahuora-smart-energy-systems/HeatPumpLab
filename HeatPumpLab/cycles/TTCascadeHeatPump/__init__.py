"""Transcritical-transcritical CO2 cascade heat pump for air drying (Kong 2026, PhD thesis, University of Waikato,
cycle TT-4 with the additional gas cooler and secondary air; see also Kong et al., Applied Thermal Engineering 238
(2024) 122005, and Kong et al., Energies 19 (2026) 3767, for the expanders); ported from the Advanced dryer heat pumps
web app."""

from .tt_cascade import TTCascadeConfig, TTCascadeHeatPump, TTCascadeResult

__all__ = ["TTCascadeConfig", "TTCascadeHeatPump", "TTCascadeResult"]
