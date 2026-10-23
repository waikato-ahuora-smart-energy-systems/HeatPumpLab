"""Which pressure-targeting and retry path a working fluid uses."""

from __future__ import annotations

__all__ = ["LEGACY_PURE_SOLVER_FLAG", "uses_band_pressure_targeting"]

# Case/model flag that restores the pre-v21 pure-fluid path for historical replays.
LEGACY_PURE_SOLVER_FLAG = "legacy_pure_solver_path"


def uses_band_pressure_targeting(hp) -> bool:
    """Return whether this model uses the band-target pressure search and retries.

    Mixtures always do. Since v21, pure fluids do too, unless the legacy flag is
    set: the old pure path accepted feasible endpoints above the minimum
    discharge pressure and depended on the initial pressure guess.
    """
    if ";" in getattr(hp, "refrigerant", ""):
        return True
    return not bool(getattr(hp, LEGACY_PURE_SOLVER_FLAG, False))
