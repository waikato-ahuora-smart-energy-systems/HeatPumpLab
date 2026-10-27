"""Two-stage economiser / flash-tank / parallel-compression heat pump (HeatPumpLab/cycles/EconomiserHeatPump)."""
import pytest

from HeatPumpLab.cycles.EconomiserHeatPump import ARRANGEMENTS, EconomiserConfig, EconomiserHeatPump

pytest.importorskip("CoolProp")

_MODELS: dict = {}


def solve(arrangement, T_air_inlet=15.0, **kw):
    """Solve with cyclopentane on CoolProp; warm (recirculated) inlet air is dry so its h(T) table is cheap."""
    key = (T_air_inlet,)
    if key not in _MODELS:
        _MODELS[key] = EconomiserHeatPump(EconomiserConfig(property_backend="coolprop", T_air_inlet=T_air_inlet,
                                                           air_RH=60.0 if T_air_inlet < 60 else 0.0))
    return _MODELS[key].solve(arrangement=arrangement, **kw)


@pytest.mark.parametrize("arrangement", ARRANGEMENTS)
@pytest.mark.parametrize("T_air_inlet", [15.0, 110.0])
def test_energy_balance_and_air_approach(arrangement, T_air_inlet):
    r = solve(arrangement, T_air_inlet)
    assert r.success, r.message
    # heat to air = compressor work + heat from the source, per kg of high-pressure flow
    assert r.q_air == pytest.approx(r.w_shaft + r.q_source, abs=1e-6)
    assert r.report["gap_air"] >= 20.0 - 0.06
    assert r.COP == pytest.approx(r.q_air / (r.w_shaft / 0.93), rel=1e-12)
    flows = r.report["flows"]
    assert flows["high"] == 1.0 and 0.0 <= flows["injected"] < 1.0


def test_cold_inlet_air_leaves_no_flash_gas():
    """With 15 C air the subcooler cools the liquid far below the intermediate saturation temperature, so the flash
    tank and economiser carry nothing: the two-stage gain is then the compressors' (lower pressure ratio per stage)."""
    flash, econ, single = solve("flash"), solve("economiser"), solve("single")
    assert flash.report["x_flash"] == 0.0 and econ.report["y_econ"] == 0.0
    assert flash.COP == pytest.approx(econ.COP, rel=1e-12)
    assert flash.COP > single.COP
    assert any("pressure ratio" in w for w in single.warnings) or single.report["P_high_bar"] / single.report["P_evap_bar"] < 20


def test_warm_inlet_air_flash_and_economiser_recover_throttling_loss():
    """With recirculated air entering at 110 C the liquid leaves the condenser hot: the single-stage valve loses a lot,
    and the flash tank and economiser carry vapour and recover much of it."""
    flash, econ, single = solve("flash", 110.0), solve("economiser", 110.0), solve("single", 110.0)
    assert flash.report["x_flash"] > 0.1 and econ.report["y_econ"] > 0.1
    assert flash.COP > 1.4 * single.COP and econ.COP > 1.4 * single.COP
    assert flash.COP >= econ.COP - 1e-9                 # the economiser's approach costs a little


def test_parallel_compression_superheats_for_a_dry_discharge():
    r = solve("parallel", 110.0)
    assert r.success and r.report["x_flash"] > 0.1
    st = r.states
    # cyclopentane overhangs: the parallel suction is superheated and its discharge is dry
    assert st["parallel_suction"]["T"] > r.report["T_int"]
    assert r.report["T_dis_parallel"] > st["condenser_in"]["T"] - 50


def test_optimised_intermediate_pressure_is_no_worse():
    hp = EconomiserHeatPump(EconomiserConfig(property_backend="coolprop", T_air_inlet=110.0, air_RH=0.0))
    base = hp.solve()
    best = hp.optimise(max_eval=60)
    assert best.success
    opt = best.report["optimiser"]
    assert best.COP >= base.COP - 1e-9 and opt["COP"] == pytest.approx(best.COP)
    assert opt["lo"][0] <= opt["T_int"] <= opt["hi"][0]
    assert best.report["T_int"] == pytest.approx(opt["T_int"], abs=0.05)


def test_fixed_intermediate_temperature_and_bad_inputs():
    r = solve("flash", 110.0, T_int=90.0)
    assert r.success and r.report["T_int"] == pytest.approx(90.0, abs=1e-6)
    bad = solve("triple")
    assert not bad.success and "arrangement" in bad.message
    unreachable = solve("single", T_air_tar=240.0, P_max_bar=5.0)
    assert not unreachable.success
