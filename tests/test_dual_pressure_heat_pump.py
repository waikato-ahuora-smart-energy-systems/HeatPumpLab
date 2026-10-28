"""Dual-pressure condensation heat pump (HeatPumpLab/cycles/DualPressureHeatPump)."""
import pytest

from HeatPumpLab.cycles.DualPressureHeatPump import DualPressureConfig, DualPressureHeatPump
from HeatPumpLab.cycles.EconomiserHeatPump import EconomiserConfig, EconomiserHeatPump

pytest.importorskip("CoolProp")


@pytest.fixture(scope="module")
def hp():
    return DualPressureHeatPump(DualPressureConfig(property_backend="coolprop"))


@pytest.mark.parametrize("T_cond_low", [100.0, 120.0, 140.0])
def test_balances_and_approaches(hp, T_cond_low):
    r = hp.solve(T_cond_low=T_cond_low)
    assert r.success, r.message
    rp = r.report
    assert r.q_air == pytest.approx(r.w_shaft + r.q_source, abs=1e-6)           # per kg of evaporator flow
    assert rp["q_air"] == pytest.approx(rp["q_air_low"] + rp["q_air_high"], rel=1e-12)
    assert 0.0 < rp["split_high"] < 1.0
    assert min(rp["gap_low"], rp["gap_high"]) >= 20.0 - 0.06
    # the air leaves the low-pressure condenser at its pinch: the low condensing temperature less the approach
    assert rp["T_mid"] == pytest.approx(T_cond_low - 20.0, abs=0.05)
    # the two condensers share one air flow: each heats its own range with the duty it carries
    air = hp.air
    mA = rp["air_per_ref"]
    assert rp["q_air_low"] == pytest.approx(mA * (air.h(rp["T_mid"]) - air.h(15.0)), rel=1e-6)
    assert rp["q_air_high"] == pytest.approx(mA * (air.h(170.0) - air.h(rp["T_mid"])), rel=1e-6)


def test_beats_a_single_condensing_pressure():
    """At the same conditions and compressors, splitting the condensing duty over two pressures beats the two-stage
    flash cycle that condenses everything at the high pressure (Dai et al. report 3-9 % with water sinks)."""
    dual = DualPressureHeatPump(DualPressureConfig(property_backend="coolprop")).optimise(max_eval=40)
    flash = EconomiserHeatPump(EconomiserConfig(property_backend="coolprop")).optimise(max_eval=40)
    assert dual.success and flash.success
    assert dual.COP > flash.COP
    assert dual.COP / flash.COP - 1 < 0.3                   # a gain, not a different machine


def test_optimiser_is_no_worse_than_the_start(hp):
    start = hp.solve()
    best = hp.optimise(max_eval=40)
    opt = best.report["optimiser"]
    assert best.COP >= start.COP - 1e-9 and opt["COP"] == pytest.approx(best.COP)
    assert opt["lo"][0] <= opt["T_cond_low"] <= opt["hi"][0]


def test_infeasible_designs_fail_with_a_reason(hp):
    near_source = hp.solve(T_cond_low=50.0)
    assert not near_source.success and "above the evaporating temperature" in near_source.message
    supercrit = hp.solve(T_cond_low=260.0)
    assert not supercrit.success and "critical temperature" in supercrit.message
    hot_inlet = hp.solve(T_cond_low=110.0, T_air_inlet=95.0)
    assert not hot_inlet.success and "too close to the air inlet" in hot_inlet.message


def test_low_condensing_temperature_well_below_the_target(hp):
    """A low condenser at 80 C still balances (most of the flow goes on to the high stage): the split is solved by
    a bracketed root, where the old damped iteration stopped at its first overshoot."""
    r = hp.solve(T_cond_low=80.0)
    assert r.success, r.message
    rp = r.report
    assert 0.75 < rp["split_high"] < 0.9 and abs(rp["balance"]) < 1e-6
    assert min(rp["gap_low"], rp["gap_high"]) >= 20.0 - 0.06
