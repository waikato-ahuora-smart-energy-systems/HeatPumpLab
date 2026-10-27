"""Two-stage vapour-compression heat pump for air drying with a flash tank, a closed economiser (vapour injection) or
parallel compression, and a single-stage cycle for comparison."""

from .economiser_hp import ARRANGEMENTS, EconomiserConfig, EconomiserHeatPump, EconomiserResult

__all__ = ["ARRANGEMENTS", "EconomiserConfig", "EconomiserHeatPump", "EconomiserResult"]
