import pytest


def test_base_case_overrides_compressor_model():
    mobo = pytest.importorskip("main_RHPfixed_mobo")
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import DEFAULT_BASE_CASE

    default = mobo.base_case()
    assert default == dict(DEFAULT_BASE_CASE, T_air_tar=200.0)
    case = mobo.base_case(5, 0.93)
    assert case["comp_model_index"] == 5 and case["motor_eff"] == 0.93
    assert {k: v for k, v in case.items() if k not in ("comp_model_index", "motor_eff")} == {
        k: v for k, v in default.items() if k not in ("comp_model_index", "motor_eff")
    }
