"""Hybrid absorption-compression (Osenbrueck) heat pump with ammonia-water for air drying; needs the ammonia-water
extra (teqp)."""

from .hybrid_hp import HybridConfig, HybridHeatPump, HybridResult

__all__ = ["HybridConfig", "HybridHeatPump", "HybridResult"]
