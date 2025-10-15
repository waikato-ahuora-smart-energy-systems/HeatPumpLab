"""Ensure every module, class, and function in RecuperativeHeatPump has a docstring."""

from __future__ import annotations

import ast
from pathlib import Path


def _collect_missing_docstrings(root: Path) -> list[str]:
    missing: list[str] = []
    for path in root.rglob("*.py"):
        if "__pycache__" in path.parts:
            continue
        tree = ast.parse(path.read_text())
        module_doc = ast.get_docstring(tree)
        if module_doc is None:
            missing.append(f"{path}:module")
        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                if ast.get_docstring(node) is None:
                    missing.append(f"{path}:{node.name}:{node.lineno}")
            elif isinstance(node, ast.AsyncFunctionDef):
                if ast.get_docstring(node) is None:
                    missing.append(f"{path}:{node.name}:{node.lineno}")
            elif isinstance(node, ast.ClassDef):
                if ast.get_docstring(node) is None:
                    missing.append(f"{path}:{node.name}:{node.lineno}")
                for sub in node.body:
                    if isinstance(sub, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        if ast.get_docstring(sub) is None:
                            missing.append(
                                f"{path}:{node.name}.{sub.name}:{sub.lineno}"
                            )
    return missing


def test_recuperative_heat_pump_docstrings_present() -> None:
    """Fail when a module/class/function in RecuperativeHeatPump lacks a docstring."""
    project_root = Path(__file__).resolve().parents[1]
    target_root = project_root / "RecuperativeHeatPump"
    missing = _collect_missing_docstrings(target_root)
    assert not missing, f"Missing docstrings:\n" + "\n".join(missing)
