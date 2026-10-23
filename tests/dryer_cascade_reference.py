"""Shared checks of the multi-loop dryer cycles against the web app's JS models.

tests/data/dryer_cascades_js_reference.json holds, per case, the Python config overrides and the results of the
app's DRYCASC models run in Node over CoolProp 8 (WebAssembly); see its "source" field. The Python ports run on the
CoolProp backend, so they should agree to round-off (the same equations, property library and iterations).
"""
import json
import math
from pathlib import Path

import pytest

DATA = Path(__file__).parent / "data" / "dryer_cascades_js_reference.json"
REL, ABS = 1e-6, 1e-5


def reference_cases(model: str) -> list:
    return [c for c in json.loads(DATA.read_text())["cases"] if c["model"] == model]


def coolprop_or_skip():
    pytest.importorskip("CoolProp")


def close(a, b, what):
    if isinstance(b, bool) or b is None or isinstance(b, str):
        assert a == b, what
    else:
        assert a == pytest.approx(b, rel=REL, abs=ABS), what


def compare_scalars(got: dict, want: dict, where: str, skip=()):
    for k, v in want.items():
        if k in skip or k not in got or isinstance(v, (list, dict)):
            continue
        if isinstance(v, float) and not math.isfinite(v):
            continue
        close(got[k], v, f"{where}.{k}")


def compare_result(res, exp: dict):
    assert res.success == exp["success"], res.message
    if not exp["success"]:
        assert res.message == exp["message"]
        return
    for k in ("COP", "VHC", "Pmax_bar", "T_dis_max", "q_air", "w_shaft", "q_source"):
        close(getattr(res, k), exp[k], k)
    assert res.warnings == exp["warnings"]
    compare_scalars(res.report, exp["report"], "report")
