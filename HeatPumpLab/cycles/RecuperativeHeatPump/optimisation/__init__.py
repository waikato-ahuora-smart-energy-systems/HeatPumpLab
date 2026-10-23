"""Outer-design optimisation tools for recuperative heat-pump cases."""

from importlib import import_module
from typing import TYPE_CHECKING

from ._types import (
    OPTIMISER_VERSION,
    DesignEvaluation,
    HistoricalIncumbent,
    IncumbentReplay,
    OuterDesignOptimisationConfig,
    OuterDesignOptimisationResult,
    PhaseEvaluation,
)

if TYPE_CHECKING:
    from ._batch import (
        optimise_recuperative_heat_pump_case,
        optimise_recuperative_heat_pump_cases,
    )
    from .outer_design import optimise_recuperative_heat_pump

# The optimiser entry points pull in torch/BoTorch. Load them on first use so that
# isolated exact-evaluation workers, which only need ``_objective``, start quickly.
_LAZY_EXPORTS = {
    "optimise_recuperative_heat_pump_case": "._batch",
    "optimise_recuperative_heat_pump_cases": "._batch",
    "optimise_recuperative_heat_pump": ".outer_design",
}


def __getattr__(name: str):
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(module_name, __name__), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(_LAZY_EXPORTS))


__all__ = [
    "OPTIMISER_VERSION",
    "DesignEvaluation",
    "HistoricalIncumbent",
    "IncumbentReplay",
    "OuterDesignOptimisationConfig",
    "OuterDesignOptimisationResult",
    "PhaseEvaluation",
    "optimise_recuperative_heat_pump_case",
    "optimise_recuperative_heat_pump_cases",
    "optimise_recuperative_heat_pump",
]
