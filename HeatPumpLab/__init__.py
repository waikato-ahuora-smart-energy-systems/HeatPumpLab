"""Recuperative Heat Pump public API helpers."""

import warnings

# linear_operator 0.6.1 (via GPyTorch/BoTorch) decorates two CG helpers with @torch.jit.script at
# import time, and recent torch releases warn that torch.jit.script is deprecated. Harmless; silence only that
# message, before anything imports BoTorch. Child processes inherit nothing, but they import this
# package first as well.
warnings.filterwarnings(
    "ignore", message=r"`torch\.jit\.script` is deprecated", category=FutureWarning
)

from .utils import get_parameter_sweep  # noqa: E402

__all__ = [
    "get_parameter_sweep",
]
