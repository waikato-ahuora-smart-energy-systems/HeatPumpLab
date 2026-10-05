import json
import math
from pathlib import Path

import pytest

FIT_JSON = (
    Path(__file__).resolve().parents[1]
    / "paper/compressor_research/efficiency-pr-front-20261004"
    / "v2-manufacturer-software/front_coefficients.json"
)

needs_fit = pytest.mark.skipif(
    not FIT_JSON.exists(), reason="v2 fit data under paper/ (not in git)"
)


@pytest.fixture(scope="module")
def comp_mod():
    try:
        from HeatPumpLab.components import compressor
    except Exception as exc:
        pytest.skip(f"HeatPumpLab components unavailable: {exc}")
    return compressor


@needs_fit
def test_v2_coefficients_match_fit(comp_mod):
    fit = json.loads(FIT_JSON.read_text())
    for basis, key in (("isentropic", "eta_shaft"), ("overall", "eta_elec")):
        c = fit[key]["design_fit"]
        assert comp_mod.PARETO_FRONT_V2_COEFFS[basis] == pytest.approx(
            tuple(c[k] for k in ("eta_1", "a", "b", "d", "L0", "w"))
        )


@needs_fit
def test_v2_tracks_consensus_front(comp_mod):
    fit = json.loads(FIT_JSON.read_text())
    for basis, key in (("isentropic", "eta_shaft"), ("overall", "eta_elec")):
        lo, hi = fit[key]["design_fit"]["pr_range"]
        for row in fit[key]["table"]:
            if row["design"] is None or not lo <= row["PR"] <= hi:
                continue
            eta = comp_mod.pareto_front_efficiency(row["PR"], basis, version=2)
            # 0.075: the smooth formula sits under the front; the largest gap is the single-source
            # peak near PR 4-5 (Mayekawa J NH3, 0.89)
            assert eta == pytest.approx(row["design"], abs=0.075)
            assert eta <= row["design"] + 1e-9  # lower envelope: never above the design front


@pytest.mark.parametrize("basis", ["isentropic", "overall"])
def test_v2_non_increasing_and_bounded_to_pr_20(comp_mod, basis):
    prs = [1.0 + 0.1 * i for i in range(191)]
    etas = [comp_mod.pareto_front_efficiency(pr, basis, version=2) for pr in prs]
    assert all(e2 <= e1 + 1e-12 for e1, e2 in zip(etas, etas[1:]))
    assert all(0.35 < e < 0.88 for e in etas)


def test_v2_overall_is_motor_efficiency_times_isentropic(comp_mod):
    # both fronts pool shaft and electrical data with a 0.93 motor efficiency
    f = comp_mod.pareto_front_efficiency
    for pr in (1.5, 3.0, 6.0, 10.0, 15.0):
        assert f(pr, "overall", version=2) == pytest.approx(0.93 * f(pr, "isentropic", version=2), abs=2e-3)


def test_v2_values(comp_mod):
    f = comp_mod.pareto_front_efficiency
    assert f(1.0, version=2) == pytest.approx(0.8695, abs=1e-4)  # Piller MVR level at PR 1
    assert f(2.0, version=2) == pytest.approx(0.8547, abs=1e-3)  # no early drop
    assert f(1.5, version=2) == pytest.approx(0.8628, abs=1e-3)
    assert f(3.0, version=2) == pytest.approx(0.8403, abs=1e-3)
    assert f(10.0, version=2) == pytest.approx(0.7718, abs=1e-3)
    assert f(15.0, version=2) == pytest.approx(0.6929, abs=1e-3)
    assert f(20.0, version=2) == pytest.approx(0.5783, abs=1e-3)
    assert f(3.0) == pytest.approx(0.846, abs=1e-3)  # v1 default unchanged


def test_compressor_model_index_5_uses_v2(hpl):
    compressor = hpl.Compressor(fluid="Butane", model_index=5, motor_eff=0.96)
    compressor.flow_in = hpl.Stream(
        fluid="Butane", temperature=10.0, quality=1.0, mass_flow=1.0
    )
    compressor.flow_out = hpl.Stream(
        fluid="Butane", pressure=compressor.flow_in.pressure * 3.0, mass_flow=1.0
    )
    compressor.solve()

    from HeatPumpLab.components.compressor import pareto_front_efficiency

    assert compressor.isentropic_eff == pytest.approx(
        pareto_front_efficiency(3.0, version=2)
    )
    assert math.isfinite(compressor.power) and compressor.power > 0.0
