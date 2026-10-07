"""Split-stream heat pump (RHX1 in parallel with the low-temperature gas cooler)."""
import numpy as np
import pytest

from HeatPumpLab.components.compressor import pareto_front_efficiency


def split_module():
    from HeatPumpLab.cycles.SplitHeatPump import split_hp
    return split_hp


# ---------------------------------------------------------------- pure helpers (no REFPROP)

def test_counterflow_min_approach_finds_internal_pinch():
    sh = split_module()
    hot = lambda f: 100 + 50 * f
    cold = lambda f: 20 + 120 * f ** 2
    expected = min(100 + 50 * f - 20 - 120 * f ** 2 for f in np.linspace(0, 1, 201))
    assert sh.counterflow_min_approach(hot, cold, 201) == pytest.approx(expected)


def test_compressor_efficiency_matches_shared_correlation():
    sh = split_module()
    for pr in (1.5, 3.0, 9.5):
        assert sh.isentropic_efficiency(pr, 5) == pytest.approx(pareto_front_efficiency(pr, "isentropic", version=2))
    assert sh.isentropic_efficiency(3.0, 1) == pytest.approx(0.7 / 0.96)
    assert sh.volumetric_efficiency(1.0) == pytest.approx(1.0)


def test_branch_outlet_temperatures_follow_the_config():
    sh = split_module()
    hp = sh.SplitHeatPump(sh.SplitHPConfig(refrigerant="X", dT_rhx1_cold_end=25.0), props=object())
    assert hp.T_superheat == pytest.approx(35.0)
    assert hp.T_LT_out == pytest.approx(35.0)      # 15 C air + 20 K: cooled as far as practical
    assert hp.T_RX_out == pytest.approx(60.0)      # suction-gas inlet + 25 K


# ---------------------------------------------------------------- REFPROP-backed cycle

FLUID = "Propane;Isopentane|0.225;0.775|mass"


@pytest.fixture
def solved(hpl):
    sh = split_module()

    def _solve(**kw):
        cfg = sh.SplitHPConfig(refrigerant=FLUID, **({"T_air_tar": 150.0, "dT_rhx1_cold": 60.0} | kw))
        return sh.SplitHeatPump(cfg).solve()
    return _solve


def test_split_cycle_solves_and_closes_energy_balance(solved):
    res = solved()
    assert res.success, res.message
    r, st = res.report, res.states
    assert 2.0 < res.COP < 5.0
    assert r["q_HT_gc_kJkg"] + r["q_LT_gc_kJkg"] == pytest.approx(r["q_air_kJkg"], rel=1e-6)
    a = r["alpha_RHX1"]
    assert 0 < a < 1
    assert a * (st["split"]["h"] - st["RHX1_hot_out"]["h"]) == pytest.approx(r["q_rhx1_kJkg"], rel=1e-6)
    assert st["valve_in"]["h"] == pytest.approx((1 - a) * st["LT_gc_out"]["h"] + a * st["RHX1_hot_out"]["h"])
    assert st["discharge"]["h"] - st["valve_in"]["h"] == pytest.approx(r["q_air_kJkg"] + r["q_rhx1_kJkg"], rel=1e-6)
    assert res.COP == pytest.approx(r["q_air_kJkg"] / (r["w_shaft_kJkg"] / 0.93), rel=1e-9)
    assert r["T_LT_gc_out"] == pytest.approx(35.0) and r["T_RHX1_hot_out"] == pytest.approx(60.0)
    for name, need in (("HT_gas_cooler", 20.0), ("LT_gas_cooler", 20.0), ("RHX1", 5.0)):
        assert r["approaches_K"][name] >= need - 0.05
    assert r["rhx1_LMTD_K"] > 5.0


def test_zero_lift_has_no_rhx1_branch(solved):
    # Without RHX1 the 22.5/77.5 mixture evaporating from 10 C ends its glide above the 35 C
    # superheater outlet, so evaporate from 5 C. Any remaining failure must be a physical one.
    res = solved(dT_rhx1_cold=0.0, Tref_evap_inlet=5.0)
    if not res.success:
        assert "wet compression" in res.message or "evaporator dew" in res.message, res.message
        return
    assert res.report["alpha_RHX1"] == pytest.approx(0.0)
    assert res.report["q_rhx1_kJkg"] == pytest.approx(0.0, abs=1e-9)


def test_lmtd_floor_is_respected(solved):
    res = solved(rhx1_lmtd_min=25.0)
    if res.success:
        assert res.report["rhx1_LMTD_K"] >= 25.0 - 1e-6
