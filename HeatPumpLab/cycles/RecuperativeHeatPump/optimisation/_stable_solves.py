"""Triangular solves whose rounding does not depend on where a buffer sits in memory.

Apple Accelerate's ``dtrsm`` handles an unpaired right-hand-side column with a kernel whose
rounding depends on the address of that column (seen on an M3 Max, for 64+ rows). Identical inputs
then give results that differ in the last bits. The GP's ill-conditioned kernel matrix amplifies
this to about 1e-10 in the acquisition, enough to change BO proposals between two same-seed
processes, whose heaps differ. Paired columns are not affected, so an odd column count is padded
with a zero column. Results are unchanged beyond rounding, and other BLAS libraries are unaffected.
"""

import torch
from torch.overrides import TorchFunctionMode

__all__ = ["AddressIndependentSolves"]


def _solve_padded(rhs, call):
    """Return ``call(rhs with an extra zero column)`` without that column."""
    padded = torch.cat([rhs, torch.zeros_like(rhs[..., :1])], -1)
    return call(padded)[..., :-1]


class AddressIndependentSolves(TorchFunctionMode):
    """Within this context, solve odd-column right-hand sides with one extra zero column."""

    def __torch_function__(self, func, types, args=(), kwargs=None):
        kwargs = kwargs or {}
        if func is torch.linalg.solve_triangular:
            A, B = args
            if (
                kwargs.get("left", True)
                and kwargs.get("out") is None
                and type(B) is torch.Tensor
                and B.shape[-1] % 2
            ):
                return _solve_padded(B, lambda rhs: func(A, rhs, **kwargs))
        elif func in (torch.cholesky_solve, torch.Tensor.cholesky_solve):
            B, *rest = args
            if kwargs.get("out") is None and type(B) is torch.Tensor and B.shape[-1] % 2:
                return _solve_padded(B, lambda rhs: func(rhs, *rest, **kwargs))
        return func(*args, **kwargs)
