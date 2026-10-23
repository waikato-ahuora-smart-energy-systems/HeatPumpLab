import math
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
for _p in (REPO_ROOT, REPO_ROOT / "scripts"):      # scripts/: the rhp_* research search modules
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


BASE_RHP_ARGS = {
    "T_source_in": 40.0,
    "Tref_evap_inlet": 10.0,
    "T_air_inlet": 15.0,
    "T_air_tar": 200.0,
    "dT_GC_min_tar": 20.0,
    "dTmin_rhx": 5.0,
    "dT_rhx1_cold": 0.0,
    "dT_rhx2_cold": 0.0,
    "dT_rhx3_cold": 0.0,
    "comp2_power_frac": 0.5,
    "comp3_power_frac": 0.0,
    "comp_PR_uniform": False,
    "comp_model_index": 1,
    "n_comp": 2,
    "ihx_on": True,
    "motor_eff": 0.96,
    "sink": "Air",
    "source": "water",
    "debug": False,
}


@pytest.fixture(scope="session")
def hpl():
    try:
        from HeatPumpLab.components import (
            Compressor,
            CustomRelationship,
            Economiser,
            Evaporator,
            ExpansionValve,
            GasCooler,
            Mixer,
            PlateHeatExchanger,
            SetRelationship,
            Splitter,
            Stream,
            Turbine,
        )
        from HeatPumpLab.components.heat_exchanger import HeatExchanger
        from HeatPumpLab.cycles.RecuperativeHeatPump import RecuperativeHP
        from HeatPumpLab.utils.properties_package import PropertyPackage

        PropertyPackage("Butane")
    except Exception as exc:
        pytest.skip(f"REFPROP-backed HeatPumpLab tests are unavailable: {exc}")

    return SimpleNamespace(
        Compressor=Compressor,
        CustomRelationship=CustomRelationship,
        Economiser=Economiser,
        Evaporator=Evaporator,
        ExpansionValve=ExpansionValve,
        GasCooler=GasCooler,
        HeatExchanger=HeatExchanger,
        Mixer=Mixer,
        PlateHeatExchanger=PlateHeatExchanger,
        PropertyPackage=PropertyPackage,
        RecuperativeHP=RecuperativeHP,
        SetRelationship=SetRelationship,
        Splitter=Splitter,
        Stream=Stream,
        Turbine=Turbine,
    )


@pytest.fixture
def base_rhp_args():
    return BASE_RHP_ARGS.copy()


@pytest.fixture
def solve_recuperative_case(hpl, base_rhp_args):
    def _solve(refrigerant, **overrides):
        hp = hpl.RecuperativeHP(
            **(base_rhp_args | overrides | {"refrigerant": refrigerant})
        )
        return hp.solve()

    return _solve


@pytest.fixture
def assert_successful_cycle():
    def _assert(hp):
        assert hp.feasible, getattr(hp, "solve_error", "")
        assert hp.converged, getattr(hp, "solve_error", "")
        assert math.isfinite(hp.COP)
        assert math.isfinite(hp.SEC)
        assert hp.COP > 1.0
        assert sum(c.power for c in hp.compressors) > 0.0
        assert sum(gc.duty for gc in hp.gas_coolers) > sum(
            c.power for c in hp.compressors
        )
        assert hp.evaporator.sink_in.temperature <= hp.Tref_evap_inlet + 1e-4
        assert hp.gas_coolers[0].sink_out.temperature == pytest.approx(
            hp.T_air_tar, abs=1e-6
        )

    return _assert
