"""n-level (2-4) cascade heat pump for air drying (after Hoang et al., Applied Thermal Engineering 279 (2025) 127901), with air condensers and subcoolers per level
and optional IHXs; ported from the Advanced dryer heat pumps web app."""

from .cascade_hp import CascadeConfig, CascadeHeatPump, CascadeResult

__all__ = ["CascadeConfig", "CascadeHeatPump", "CascadeResult"]
