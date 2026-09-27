from types import SimpleNamespace

from HeatPumpLab.cycles.RecuperativeHeatPump.solver._fluid_policy import (
    uses_band_pressure_targeting,
)


def test_mixtures_always_use_band_targeting():
    hp = SimpleNamespace(refrigerant="Propane;Butane|0;1|mass", legacy_pure_solver_path=True)
    assert uses_band_pressure_targeting(hp)


def test_pure_fluids_use_band_targeting_by_default():
    assert uses_band_pressure_targeting(SimpleNamespace(refrigerant="Butane"))


def test_legacy_flag_restores_old_pure_path():
    hp = SimpleNamespace(refrigerant="Butane", legacy_pure_solver_path=True)
    assert not uses_band_pressure_targeting(hp)
