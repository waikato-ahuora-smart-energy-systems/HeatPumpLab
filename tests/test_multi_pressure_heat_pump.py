"""Multi-pressure condensation heat pump (HeatPumpLab/cycles/MultiPressureHeatPump)."""
import pytest

from HeatPumpLab.cycles.DualPressureHeatPump import DualPressureConfig, DualPressureHeatPump
from HeatPumpLab.cycles.MultiPressureHeatPump import MultiPressureConfig, MultiPressureHeatPump

pytest.importorskip("CoolProp")


@pytest.fixture(scope="module")
def hp():
    return MultiPressureHeatPump(MultiPressureConfig(property_backend="coolprop"))


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_balances_flows_and_approaches(hp, n):
    r = hp.solve(n_pressures=n)
    assert r.success, r.message
    rp = r.report
    assert r.q_air == pytest.approx(r.w_shaft + r.q_source, abs=1e-6)            # per kg of evaporator flow
    assert rp["q_air"] == pytest.approx(sum(rp["q_air_levels"]), rel=1e-12)
    m = rp["flows"]
    assert m[0] == 1.0 and all(b < a for a, b in zip(m, m[1:])) and m[-1] > 0     # each compressor takes less
    assert min(rp["gaps"]) >= 20.0 - 0.06
    # one air flow through the condensers in series, each over its own range
    T = [15.0] + rp["T_air_between"] + [170.0]
    for k, q in enumerate(rp["q_air_levels"]):
        assert q == pytest.approx(rp["air_per_ref"] * (hp.air.h(T[k + 1]) - hp.air.h(T[k])), rel=1e-6)
    # each lower condenser's air leaves at most at its pinch (its condensing temperature less the approach); with
    # equal steps up to four pressures they all reach it
    for Tc, Ta in zip(rp["T_cond"], rp["T_air_between"]):
        assert Ta <= Tc - 20.0 + 0.05
        if n <= 4:
            assert Ta == pytest.approx(Tc - 20.0, abs=0.05)


def test_more_pressures_raise_the_cop(hp):
    cops = [hp.solve(n_pressures=n).COP for n in (2, 3, 4)]
    assert cops[0] < cops[1] < cops[2]


def test_two_pressures_are_the_dual_pressure_cycle(hp):
    dual = DualPressureHeatPump(DualPressureConfig(property_backend="coolprop", T_cond_low=120.0)).solve()
    multi = hp.solve(n_pressures=2, T_cond=[120.0])
    assert dual.success and multi.success
    assert dual.COP == pytest.approx(multi.COP, rel=1e-12)
    assert dual.report["T_mid"] == pytest.approx(multi.report["T_air_between"][0], rel=1e-12)
    assert dual.report["split_high"] == pytest.approx(multi.report["flows"][1], rel=1e-12)


def test_optimiser(hp):
    start = hp.solve(n_pressures=3)
    best = hp.optimise(n_pressures=3, max_eval=60)
    o = best.report["optimiser"]
    assert best.success and best.COP >= start.COP - 1e-9 and o["COP"] == pytest.approx(best.COP)
    assert o["T_cond"][0] < o["T_cond"][1] and best.report["T_cond"] == o["T_cond"]


@pytest.mark.parametrize("over,msg", [
    (dict(n_pressures=6), "2 to 5"),
    (dict(n_pressures=3, T_cond=[120.0]), "need 2"),
    (dict(n_pressures=3, T_cond=[130.0, 110.0]), "must rise"),
    (dict(n_pressures=3, T_cond=[48.0, 110.0]), "above the evaporating temperature"),
    (dict(n_pressures=3, T_cond=[110.0, 260.0]), "critical temperature"),
])
def test_bad_inputs(hp, over, msg):
    r = hp.solve(**over)
    assert not r.success and msg in r.message
