"""Hybrid absorption-compression (ammonia-water) heat pump (HeatPumpLab/cycles/HybridHeatPump)."""
import pytest

pytest.importorskip("teqp")

from HeatPumpLab.cycles.HybridHeatPump import HybridConfig, HybridHeatPump  # noqa: E402


@pytest.fixture(scope="module")
def hp():
    return HybridHeatPump(HybridConfig())


@pytest.fixture(scope="module")
def base(hp):
    return hp.solve()


def test_balances_and_approaches(base):
    r = base
    assert r.success, r.message
    rp = r.report
    # per kg of solution: heat to air = compressor + pump work + heat from the source
    assert r.q_air == pytest.approx(r.w_shaft + r.q_source, abs=1e-6)
    assert r.w_shaft == pytest.approx(rp["w_comp"] + rp["w_pump"], rel=1e-12)
    assert rp["gap_air"] >= 10.0 - 0.06 and rp["gap_source"] >= 5.0 - 0.06
    assert rp["T_des_out"] == pytest.approx(55.0, abs=0.01)
    # the separator splits the solution into ammonia-rich vapour and lean liquid
    assert rp["w_lean"] < 0.5 < rp["w_vapour"]
    assert 0.4 * rp["w_vapour"] + 0.6 * rp["w_lean"] == pytest.approx(0.5, rel=1e-6)
    assert rp["P_low_bar"] < rp["P_high_bar"] <= 50.0


def test_wet_compression_holds_the_discharge_limit(hp, base):
    rp = base.report
    assert rp["wet_compression"] and 0 < rp["m_injected"] < 0.6
    assert rp["T_dis"] == pytest.approx(170.0, abs=0.05)
    dry = hp.solve(wet_compression=False)
    assert dry.success and dry.report["T_dis"] > 250 and any("170 C limit" in w for w in dry.warnings)
    assert base.COP > dry.COP                      # the injected liquid cools the compression
    assert not hp.solve(wet_compression=False, limit_discharge=True).success


def test_higher_supply_needs_more_pressure_and_lowers_cop(hp, base):
    hot = hp.solve(T_air_tar=150.0)
    assert hot.success, hot.message
    assert hot.report["P_high_bar"] > base.report["P_high_bar"] and hot.COP < base.COP
    capped = hp.solve(T_air_tar=150.0, P_max_bar=28.0)
    assert not capped.success and "28 bar" in capped.message


def test_bad_inputs(hp):
    assert not hp.solve(q_des=1.2).success
    assert not hp.solve(w_rich=0.0).success
