import json
import math
from pathlib import Path

import pytest

FIT_JSON = (
    Path(__file__).resolve().parents[1]
    / "paper/compressor_research/efficiency-pr-front-20261004/front_coefficients.json"
)


needs_fit = pytest.mark.skipif(not FIT_JSON.exists(), reason="fit data under paper/ (not in git)")


@pytest.fixture(scope="module")
def comp_mod():
    try:
        from HeatPumpLab.components import compressor
    except Exception as exc:
        pytest.skip(f"HeatPumpLab components unavailable: {exc}")
    return compressor


@needs_fit
def test_coefficients_match_published_fit(comp_mod):
    fit = json.loads(FIT_JSON.read_text())
    for basis, key in (("isentropic", "eta_is"), ("overall", "eta_overall")):
        eta_1, a, b = comp_mod.PARETO_FRONT_COEFFS[basis]
        assert (eta_1, a, b) == pytest.approx(
            (fit[key]["eta_1"], fit[key]["a"], fit[key]["b"])
        )


@pytest.mark.parametrize("basis", ["isentropic", "overall"])
def test_front_is_non_increasing_and_bounded(comp_mod, basis):
    prs = [1.0 + 0.1 * i for i in range(91)]
    etas = [comp_mod.pareto_front_efficiency(pr, basis) for pr in prs]
    assert all(e2 <= e1 + 1e-12 for e1, e2 in zip(etas, etas[1:]))
    assert all(0.6 < e < 0.9 for e in etas)


@needs_fit
def test_front_passes_close_to_its_defining_points(comp_mod):
    fit = json.loads(FIT_JSON.read_text())
    for basis, key in (("isentropic", "eta_is"), ("overall", "eta_overall")):
        for p in fit[key]["front_points"]:
            assert comp_mod.pareto_front_efficiency(
                p["pressure_ratio"], basis
            ) == pytest.approx(p["eta"], abs=0.01)


def test_overall_front_sits_below_isentropic_front(comp_mod):
    for pr in (1.0, 2.0, 4.0, 6.0, 10.0):
        iso = comp_mod.pareto_front_efficiency(pr, "isentropic")
        ovr = comp_mod.pareto_front_efficiency(pr, "overall")
        assert 0.85 < ovr / iso < 0.97


def test_pressure_ratio_below_one_is_clamped(comp_mod):
    assert comp_mod.pareto_front_efficiency(0.8) == comp_mod.pareto_front_efficiency(
        1.0
    )


def test_compressor_model_index_4_uses_front(hpl):
    compressor = hpl.Compressor(fluid="Butane", model_index=4, motor_eff=0.96)
    compressor.flow_in = hpl.Stream(
        fluid="Butane", temperature=10.0, quality=1.0, mass_flow=1.0
    )
    compressor.flow_out = hpl.Stream(
        fluid="Butane", pressure=compressor.flow_in.pressure * 3.0, mass_flow=1.0
    )
    compressor.solve()

    from HeatPumpLab.components.compressor import pareto_front_efficiency

    assert compressor.isentropic_eff == pytest.approx(pareto_front_efficiency(3.0))
    assert compressor.overall_eff == pytest.approx(compressor.isentropic_eff * 0.96)
    assert math.isfinite(compressor.power) and compressor.power > 0.0


def test_production_default_uses_front_with_96_percent_motor():
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import DEFAULT_BASE_CASE
    assert DEFAULT_BASE_CASE["comp_model_index"] == 4
    assert DEFAULT_BASE_CASE["motor_eff"] == pytest.approx(0.96)


def test_rhp_stages_use_front_at_their_own_pressure_ratio(tmp_path):
    from HeatPumpLab.components.compressor import pareto_front_efficiency
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
    objective = _OuterDesignObjective(max_evals=1, refrigerant="Propane;Pentane|0.5;0.5|mass",
                                      results_dir=tmp_path)
    try:
        result = objective.evaluate([0.0, 60.0, 0.5, 0.4, 0.6])
    finally:
        objective.close()
    report = result.report
    assert report["comp_model_index"] == 4 and report["motor_eff"] == pytest.approx(0.96)
    for stage in (1, 2):
        assert report[f"n_is_{stage}"] == pytest.approx(pareto_front_efficiency(report[f"PR_{stage}"]))


def test_projection_to_pr_20_stays_monotonic_and_above_the_floor(comp_mod):
    assert comp_mod.PARETO_FRONT_PR_RANGE == (1.0, 20.0)
    prs = [10.0 + 0.25 * i for i in range(41)]
    etas = [comp_mod.pareto_front_efficiency(pr) for pr in prs]
    assert all(e2 < e1 for e1, e2 in zip(etas, etas[1:]))
    assert etas[0] == pytest.approx(0.692, abs=1e-3) and etas[-1] == pytest.approx(0.513, abs=1e-3)
    assert min(etas) > 0.4                       # the compressor's 0.4 floor is not reached by PR 20


def test_search_limit_is_stage_pr_20():
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import DEFAULT_MAX_STAGE_PRESSURE_RATIO
    assert DEFAULT_MAX_STAGE_PRESSURE_RATIO == 20.0
