from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump import rhp_helpers


def hp_args(**overrides):
    args = {
        "T_air_tar": 200.0,
        "T_air_inlet": 15.0,
        "dT_GC_min_tar": 20.0,
        "dT_gc_out_min_offset": 0.0,
        "T_source_in": 40.0,
        "dT_rhx1_cold": 0.0,
        "dT_rhx2_cold": 0.0,
        "dT_rhx3_cold": 0.0,
        "comp_PR_uniform": False,
        "comp2_power_frac": 0.5,
        "comp3_power_frac": 0.0,
        "n_comp": 2,
        "ihx_on": 1,
    }
    return SimpleNamespace(**(args | overrides))


@pytest.mark.parametrize(
    ("overrides", "expected"),
    [
        ({}, True),
        ({"dT_rhx1_cold": 200.0}, False),
        ({"dT_gc_out_min_offset": -0.1}, False),
        ({"dT_gc_out_min_offset": 50.1}, False),
        ({"comp2_power_frac": 0.8, "comp3_power_frac": 0.3}, False),
    ],
)
def test_rhp_input_feasibility_rejects_impossible_design_knobs(overrides, expected):
    assert rhp_helpers._check_input_param_feasibility(hp_args(**overrides)) is expected


def test_rhp_structure_normalises_uniform_three_compressor_cycle():
    hp = hp_args(
        n_comp=3,
        comp_PR_uniform=True,
        comp2_power_frac=0.8,
        comp3_power_frac=0.7,
        dT_rhx1_cold=-5.0,
        dT_rhx2_cold=10.0,
        dT_rhx3_cold=15.0,
        T_source_in=40.0,
        T_air_inlet=15.0,
        dT_GC_min_tar=20.0,
    )

    rhp_helpers._get_RecuperativeHP_structure(hp)

    assert hp.comp2_power_frac == pytest.approx(1 / 3)
    assert hp.comp3_power_frac == pytest.approx(1 / 3)
    assert hp.comp1_power_frac == pytest.approx(1 / 3)
    assert hp.comp2_on == 1
    assert hp.comp3_on == 1
    assert hp.rhx1_on == 0
    assert hp.rhx2_on == 1
    assert hp.rhx3_on == 1
    assert hp.T_gc_out == pytest.approx(35.0)
    assert hp.ihx_on == 0


@pytest.mark.parametrize(
    ("initial_n_comp", "comp2_frac", "comp3_frac", "expected_n_comp"),
    [
        (3, 0.5, 0.0, 2),
        (2, 0.0, 0.0, 1),
    ],
)
def test_rhp_structure_reduces_unused_compressor_stages(
    initial_n_comp, comp2_frac, comp3_frac, expected_n_comp
):
    hp = hp_args(
        n_comp=initial_n_comp,
        comp2_power_frac=comp2_frac,
        comp3_power_frac=comp3_frac,
        dT_rhx2_cold=10.0,
        dT_rhx3_cold=10.0,
    )

    rhp_helpers._get_RecuperativeHP_structure(hp)

    assert hp.n_comp == expected_n_comp
    assert hp.comp2_on == (1 if expected_n_comp > 1 else 0)
    assert hp.comp3_on == (1 if expected_n_comp > 2 else 0)
    # an inactive stage keeps no RHX lift
    assert hp.dT_rhx2_cold == (10.0 if expected_n_comp > 1 else 0.0)
    assert hp.dT_rhx3_cold == (10.0 if expected_n_comp > 2 else 0.0)


def test_inactive_stage_lift_does_not_make_a_single_stage_case_infeasible():
    # 2/3 of the 50 K air rise is 33.3 K: RHX1's 10 K fits, the unused RHX2's 30 K must not count
    hp = hp_args(T_air_tar=65.0, comp2_power_frac=0.0, dT_rhx1_cold=10.0, dT_rhx2_cold=30.0)

    rhp_helpers._get_RecuperativeHP_structure(hp)

    assert hp.n_comp == 1 and hp.dT_rhx2_cold == 0.0
    assert rhp_helpers._check_input_param_feasibility(hp)


class FakeRelationship:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


def component(name):
    return SimpleNamespace(
        name=name,
        sink_in=SimpleNamespace(name=f"{name}.sink_in"),
        sink_out=SimpleNamespace(name=f"{name}.sink_out"),
        source_in=SimpleNamespace(name=f"{name}.source_in"),
        source_out=SimpleNamespace(name=f"{name}.source_out"),
        flow_in=SimpleNamespace(name=f"{name}.flow_in"),
        flow_out=SimpleNamespace(name=f"{name}.flow_out"),
    )


def logic_hp(n_comp=3):
    hp = SimpleNamespace(
        rhx1_on=1,
        rhx2_on=1,
        rhx3_on=1,
        n_comp=n_comp,
        dT_rhx1_cold=5.0,
        dT_rhx2_cold=6.0,
        dT_rhx3_cold=7.0,
        dTmin_rhx=4.0,
        comp2_power_frac=0.3,
        comp3_power_frac=0.2,
        tol=1e-6,
        rhx1=component("rhx1"),
        rhx2=component("rhx2"),
        rhx3=component("rhx3"),
        comp1=component("comp1"),
        comp2=component("comp2"),
        comp_discharge=SimpleNamespace(name="discharge"),
    )
    return hp


class FakeComponent:
    def __init__(self, **kwargs):
        self.kwargs = kwargs
        self.name = kwargs.get("name", self.__class__.__name__)


def test_rhp_components_include_optional_third_stage_and_recuperators(monkeypatch):
    constructors = {
        "Compressor": [],
        "GasCooler": [],
        "PlateHeatExchanger": [],
        "ExpansionValve": [],
        "Evaporator": [],
    }

    def factory(kind):
        def _create(**kwargs):
            component_obj = FakeComponent(**kwargs)
            constructors[kind].append(component_obj)
            return component_obj

        return _create

    for kind in constructors:
        monkeypatch.setattr(rhp_helpers, kind, factory(kind))

    hp = SimpleNamespace(
        refrigerant="Butane",
        sink="Air",
        source="water",
        motor_eff=0.96,
        comp_model_index=1,
        dT_GC_min_tar=20.0,
        dTmin_rhx=5.0,
        comp2_on=1,
        comp3_on=1,
        rhx1_on=1,
        rhx2_on=1,
        rhx3_on=1,
        dT_rhx1_cold=5.0,
        dT_rhx2_cold=6.0,
        dT_rhx3_cold=7.0,
    )

    rhp_helpers._get_RecuperativeHP_components(hp)

    assert len(constructors["Compressor"]) == 3
    assert len(constructors["GasCooler"]) == 4
    assert len(constructors["PlateHeatExchanger"]) == 5
    assert len(constructors["ExpansionValve"]) == 1
    assert len(constructors["Evaporator"]) == 1
    assert hp.comp3.kwargs["fluid"] == "Butane"
    assert hp.gas_cooler4.kwargs["sink_fluid"] == "Air"
    assert hp.rhx3.kwargs["dT_c"] == pytest.approx(7.0)


def test_rhp_components_omit_disabled_optional_stages(monkeypatch):
    compressors = []
    gas_coolers = []

    def compressor_factory(**kwargs):
        component_obj = FakeComponent(**kwargs)
        compressors.append(component_obj)
        return component_obj

    def gas_cooler_factory(**kwargs):
        component_obj = FakeComponent(**kwargs)
        gas_coolers.append(component_obj)
        return component_obj

    monkeypatch.setattr(rhp_helpers, "Compressor", compressor_factory)
    monkeypatch.setattr(rhp_helpers, "GasCooler", gas_cooler_factory)
    monkeypatch.setattr(rhp_helpers, "PlateHeatExchanger", lambda **kwargs: FakeComponent(**kwargs))
    monkeypatch.setattr(rhp_helpers, "ExpansionValve", lambda **kwargs: FakeComponent(**kwargs))
    monkeypatch.setattr(rhp_helpers, "Evaporator", lambda **kwargs: FakeComponent(**kwargs))

    hp = SimpleNamespace(
        refrigerant="Butane",
        sink="Air",
        source="water",
        motor_eff=0.96,
        comp_model_index=1,
        dT_GC_min_tar=20.0,
        dTmin_rhx=5.0,
        comp2_on=0,
        comp3_on=0,
        rhx1_on=0,
        rhx2_on=0,
        rhx3_on=0,
    )

    rhp_helpers._get_RecuperativeHP_components(hp)

    assert len(compressors) == 1
    assert len(gas_coolers) == 2
    assert not hasattr(hp, "comp2")
    assert not hasattr(hp, "gas_cooler3")
    assert not hasattr(hp, "rhx1")


def test_rhp_logic_blocks_wire_recuperator_and_pressure_relationships(monkeypatch):
    monkeypatch.setattr(rhp_helpers, "SetRelationship", FakeRelationship)
    monkeypatch.setattr(rhp_helpers, "CustomRelationship", FakeRelationship)
    hp = logic_hp(n_comp=3)

    rhp_helpers._get_RecuperativeHP_logic_blocks(hp)

    assert hp.set_rhx1_dT_cold.kwargs["c"] == pytest.approx(5.0)
    assert hp.set_gc2_hot_out.kwargs["c"] == pytest.approx(4.0)
    assert hp.set_rhx2_dT_cold.kwargs["c"] == pytest.approx(6.0)
    assert hp.set_gc3_hot_out.kwargs["c"] == pytest.approx(4.0)
    assert hp.set_rhx3_dT_cold.kwargs["c"] == pytest.approx(7.0)
    assert hp.set_gc4_hot_out.kwargs["c"] == pytest.approx(4.0)
    assert hp.set_P_ints.kwargs["relation"] == 1
    assert hp.set_P_ints.kwargs["target_obj"] == [hp.comp1.flow_out, hp.comp2.flow_out]


def test_rhp_logic_blocks_target_discharge_for_two_compressor_cycle(monkeypatch):
    monkeypatch.setattr(rhp_helpers, "SetRelationship", FakeRelationship)
    monkeypatch.setattr(rhp_helpers, "CustomRelationship", FakeRelationship)
    hp = logic_hp(n_comp=2)

    rhp_helpers._get_RecuperativeHP_logic_blocks(hp)

    assert hp.set_P_ints.kwargs["target_obj"] == [
        hp.comp1.flow_out,
        hp.comp_discharge,
    ]


def test_rhp_solver_pathway_orders_optional_operations():
    hp = SimpleNamespace(
        evap_sup="evap_sup",
        n_comp=3,
        set_P_ints="set_P_ints",
        ihx_on=1,
        ihx="ihx",
        rhx1_on=1,
        set_rhx1_dT_cold="set_rhx1",
        set_gc2_hot_out="set_gc2",
        comp1="comp1",
        rhx2_on=1,
        set_rhx2_dT_cold="set_rhx2",
        set_gc3_hot_out="set_gc3",
        comp2_on=1,
        comp2="comp2",
        rhx3_on=1,
        set_rhx3_dT_cold="set_rhx3",
        set_gc4_hot_out="set_gc4",
        comp3_on=1,
        comp3="comp3",
        gas_cooler4="gc4",
        rhx3="rhx3",
        gas_cooler3="gc3",
        rhx2="rhx2",
        gas_cooler2="gc2",
        rhx1="rhx1",
        gas_cooler1="gc1",
        exp_valve="expansion",
        evaporator="evaporator",
    )

    rhp_helpers._get_RecuperativeHP_solver_pathway(hp)

    assert hp.operations == [
        "evap_sup",
        "set_P_ints",
        "ihx",
        "set_rhx1",
        "set_gc2",
        "comp1",
        "set_rhx2",
        "set_gc3",
        "comp2",
        "set_rhx3",
        "set_gc4",
        "comp3",
        "gc4",
        "rhx3",
        "gc3",
        "rhx2",
        "gc2",
        "rhx1",
        "gc1",
        "ihx",
        "expansion",
        "evaporator",
    ]


def test_rhp_solver_pathway_omits_disabled_optional_operations():
    hp = SimpleNamespace(
        evap_sup="evap_sup",
        n_comp=1,
        ihx_on=0,
        rhx1_on=0,
        comp1="comp1",
        rhx2_on=0,
        comp2_on=0,
        rhx3_on=0,
        comp3_on=0,
        gas_cooler1="gc1",
        exp_valve="expansion",
        evaporator="evaporator",
    )

    rhp_helpers._get_RecuperativeHP_solver_pathway(hp)

    assert hp.operations == ["evap_sup", "comp1", "gc1", "expansion", "evaporator"]
