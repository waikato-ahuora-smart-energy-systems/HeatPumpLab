from types import SimpleNamespace
import pytest
import numpy as np
from main_RHPmobo import PAIRS, metadata_for, serialise
from rhp_mobo_search import ParetoBOSearch


def test_pairs_and_fraction_semantics():
    assert len(PAIRS) == 9 and len({frozenset(v) for v in PAIRS.values()}) == 9
    for key, pair in PAIRS.items():
        meta = metadata_for(key, 800, 20260912)
        assert meta["refrigerant"] == ";".join(pair) + "|0.5;0.5|mass"
        assert meta["variables"][-1] == "first_component_mass_fraction"
        assert meta["property_model_status"]["estimated_mixture_parameters"] == (
            "Cyclopentane" in pair
        )
    u = np.array([0.1, 0.2, 0.3, 0.4])
    assert np.allclose(ParetoBOSearch.encode(ParetoBOSearch.decode(u)), u)
    assert np.allclose(ParetoBOSearch.decode(u)[3:], [0.4, 0.6])


def record(report):
    return SimpleNamespace(
        x=[0, 20, 0.3, 0.5, 0.5],
        success=True,
        sec=0.5,
        cop=2.0,
        pmax=5000.0,
        elapsed_seconds=1.0,
        message="success",
        constraint_errors={},
        report=report,
    )


def test_temperature_and_pressure_units_and_missing_fields():
    report = dict(
        VHC=2000.0,
        Comp_1_T_out=260.0,
        Comp_2_T_out=240.0,
        T_comp_discharge=240.0,
        T_gc_out=35.0,
        Tref_evap_outlet=20.0,
        PR_1=12.0,
        PR_2=2.0,
        Overall_PR=24.0,
    )
    r = serialise(record(report))
    assert r["Pmax_MPa"] == 5 and r["VHC_MJ_per_m3"] == 2
    assert (
        r["max_reported_refrigerant_C"] == 260
        and r["max_compressor_discharge_C"] == 260
    )
    assert r["PR_1"] == 12 and r["PR_2"] == 2
    with pytest.raises(ValueError):
        serialise(record({}))


def test_global_search_resume_reproduces_next_candidate():
    s = ParetoBOSearch(123.33333333333333, 20260912)
    for i in range(32):
        x = s.ask()
        s.tell(x, 2 + i * 0.001, 2 - i * 0.002, True)
    s.ask()  # fit and save warm-start parameters
    restored = ParetoBOSearch(s.max_lift, s.seed)
    for x, y, good in zip(s.x, s.y, s.success):
        restored.tell(s.decode(x), *y, good)
    restored.kernel_state = s.kernel_state
    restored.mean_state = s.mean_state
    assert np.allclose(restored.ask(), s.ask(), atol=1e-9, rtol=0)
