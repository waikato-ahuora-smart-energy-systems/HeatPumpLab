"""Make the venv's ``ninja`` executable visible to PyTorch's C++ extension loader.

BoTorch JIT-compiles a fused qLogEHVI kernel (~3x faster) with torch.utils.cpp_extension, which
runs ``ninja --version`` from PATH. Calling ``.venv/bin/python`` without activating the venv leaves
``.venv/bin`` off PATH, so the pip-installed ninja is missed and BoTorch silently falls back to the
pure-Python path. Results are identical either way; only speed differs.
"""
from __future__ import annotations

import os
import shutil
import sys

__all__ = ["ensure_ninja_on_path"]


def ensure_ninja_on_path() -> bool:
    """Prepend the interpreter's bin directory to PATH when it holds ninja; True if ninja is found."""
    if shutil.which("ninja"):
        return True
    bin_dir = os.path.dirname(os.path.abspath(sys.executable))
    if os.path.isfile(os.path.join(bin_dir, "ninja")) or os.path.isfile(os.path.join(bin_dir, "ninja.exe")):
        os.environ["PATH"] = bin_dir+os.pathsep+os.environ.get("PATH", "")
        return True
    return False
