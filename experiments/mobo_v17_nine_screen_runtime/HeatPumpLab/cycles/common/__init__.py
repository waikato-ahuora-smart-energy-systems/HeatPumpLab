"""Expose shared cycle base classes."""

from .heatPump import BaseHeatPump
from .scalar_search import (
    ScalarSearchResult,
    find_lowest_log_bracketed_target,
    find_sampled_log_target,
)
from .solver import (
    CycleFailure,
    calculate_cop,
    calculate_sec,
    cycle_convergence_failure,
    cycle_has_converged,
)

__all__ = [
    "BaseHeatPump",
    "CycleFailure",
    "ScalarSearchResult",
    "calculate_cop",
    "calculate_sec",
    "cycle_convergence_failure",
    "cycle_has_converged",
    "find_lowest_log_bracketed_target",
    "find_sampled_log_target",
]
