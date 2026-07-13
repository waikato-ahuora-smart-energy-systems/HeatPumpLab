from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump import rhp_base
from HeatPumpLab.cycles.RecuperativeHeatPump.rhp_base import RecuperativeHP
from HeatPumpLab.cycles.common import heatPump as heat_pump_module
from HeatPumpLab.cycles.common.heatPump import BaseHeatPump


def state(**kwargs):
    defaults = {
        "pressure": 1.0,
        "temperature": 50.0,
        "enthalpy": 100.0,
        "entropy": 1.0,
        "vol_flow": 0.1,
    }
    return SimpleNamespace(**(defaults | kwargs))


def component(converged=True, **kwargs):
    defaults = {
        "converged": converged,
        "pressure_ratio": 2.0,
        "isentropic_eff": 0.75,
        "flow_out": state(pressure=2.0, temperature=80.0),
    }
    return SimpleNamespace(**(defaults | kwargs))


def report_hp(n_comp):
    hp = object.__new__(RecuperativeHP)
    hp.cycle_type = "Transcritical"
    hp.COP = 3.1
    hp.SEC = 0.48
    hp.VHC = 12.0
    hp.eDTL = 45.0
    hp.SLE = 10.7
    hp.dQdW = 8.5
    hp.LMTD_gcs = 11.0
    hp.converged = True
    hp.refrigerant = "Butane"
    hp.n_comp = n_comp
    hp.T_source_in = 40.0
    hp.dT_GC_min_tar = 20.0
    hp.T_air_inlet = 15.0
    hp.T_air_tar = 200.0
    hp.comp_PR_uniform = False
    hp.Tref_evap_inlet = 10.0
    hp.dTmin_rhx = 5.0
    hp.ihx_on = 1
    hp.comp_model_index = 1
    hp.motor_eff = 0.96
    hp.sink = "Air"
    hp.source = "water"
    hp.dT_rhx1_cold = 1.0
    hp.dT_rhx2_cold = 2.0
    hp.dT_rhx3_cold = 3.0
    hp.comp1_power_frac = 0.5
    hp.comp2_power_frac = 0.3
    hp.comp3_power_frac = 0.2
    hp.dT_gc_out_min_offset = 7.0
    hp.T_gc_out = 42.0
    hp.Tref_evap_outlet = 13.0
    hp.comp_discharge = state(pressure=5.0, temperature=120.0)
    hp.evaporator = SimpleNamespace(sink_in=state(pressure=0.4))
    hp.comp1 = component(
        pressure_ratio=2.1,
        isentropic_eff=0.71,
        flow_out=state(pressure=1.5, temperature=70.0),
    )
    if n_comp > 1:
        hp.comp2 = component(
            pressure_ratio=2.2,
            isentropic_eff=0.72,
            flow_out=state(pressure=3.0, temperature=90.0),
        )
    if n_comp > 2:
        hp.comp3 = component(
            pressure_ratio=2.3,
            isentropic_eff=0.73,
            flow_out=state(pressure=5.0, temperature=120.0),
        )
    return hp


@pytest.mark.parametrize("n_comp", [1, 3])
def test_rhp_report_exposes_domain_parameters_and_performance(n_comp):
    hp = report_hp(n_comp)

    params = hp.get_model_parameters()
    report = hp.get_report(dt=2.5)

    assert params == {
        "refrigerant": "Butane",
        "refrigerant_component_count": 1,
        "dT_GC_min_tar": 20.0,
        "T_air_tar": 200.0,
        "T_air_inlet": 15.0,
        "T_source_in": 40.0,
        "n_comp": n_comp,
        "dT_rhx1_cold": 1.0,
        "dT_rhx2_cold": 2.0,
        "dT_rhx3_cold": 3.0,
        "comp2_power_frac": 0.3,
        "comp3_power_frac": 0.2,
        "comp_PR_uniform": False,
    }
    assert report["COP"] == pytest.approx(3.1)
    assert report["Pmax"] == pytest.approx(5000.0)
    assert report["dT_rhx_cold_sum"] == pytest.approx(6.0)
    assert report["Overall_PR"] == pytest.approx(12.5)
    assert report["Solve time"] == pytest.approx(2.5)

    if n_comp == 1:
        assert report["Comp_2_P_out"] == pytest.approx(5000.0)
        assert report["PR_2"] == pytest.approx(1.0)
        assert report["PR_3"] == pytest.approx(1.0)
        assert report["n_is_2"] == pytest.approx(1.0)
        assert report["n_is_3"] == pytest.approx(1.0)
    else:
        assert report["Comp_2_P_out"] == pytest.approx(3000.0)
        assert report["PR_2"] == pytest.approx(2.2)
        assert report["PR_3"] == pytest.approx(2.3)
        assert report["n_is_2"] == pytest.approx(0.72)
        assert report["n_is_3"] == pytest.approx(0.73)


def test_rhp_graph_component_sets_include_only_converged_available_components():
    hp = object.__new__(RecuperativeHP)
    hp.ihx = component(converged=True)
    hp.exp_valve = component(converged=True)
    hp.evaporator = component(converged=False)
    hp.gas_cooler1 = component(converged=True)

    cycle_components, sink_components = hp.get_components_sets_to_graph()

    assert cycle_components == [
        [hp.ihx, "source"],
        [hp.exp_valve, "siso"],
        [hp.ihx, "sink"],
        [hp.gas_cooler1, "source"],
    ]
    assert sink_components == [[hp.gas_cooler1, "sink"]]


def test_rhp_solve_skips_solver_when_case_is_infeasible(monkeypatch):
    hp = object.__new__(RecuperativeHP)
    hp.feasible = False

    monkeypatch.setattr(
        rhp_base,
        "solve_heat_pump_cycle",
        lambda solved_hp: pytest.fail("infeasible case should not be solved"),
    )

    assert hp.solve() is hp


def test_base_heat_pump_prints_current_case_report(monkeypatch):
    printed = []

    class PrintableHeatPump(BaseHeatPump):
        def get_report(self):
            return {"COP": 3.2}

    monkeypatch.setattr(
        heat_pump_module,
        "print_result_to_terminal",
        lambda report: printed.append(report),
    )

    PrintableHeatPump().print_currrent_case()

    assert printed == [{"COP": 3.2}]


def test_base_heat_pump_delegates_cycle_plot_to_graph_module(monkeypatch):
    plotted = []

    class GraphableHeatPump(BaseHeatPump):
        def get_components_sets_to_graph(self):
            return [["cycle"]], [["sink"]]

    monkeypatch.setattr(
        heat_pump_module,
        "plot_cycle",
        lambda **kwargs: plotted.append(kwargs) or "cycle-figure",
    )

    fig = GraphableHeatPump().get_cycle_plot(
        x_axis="s",
        y_axis="P",
        refrigerant="Butane",
        include_sink=False,
    )

    assert fig == "cycle-figure"
    assert plotted == [
        {
            "cycle_obj": [["cycle"]],
            "sink_obj": [["sink"]],
            "x_axis": "s",
            "y_axis": "P",
            "refrigerant": "Butane",
            "include_sink": False,
        }
    ]


def test_base_heat_pump_delegates_driving_force_plot_to_graph_module(monkeypatch):
    plotted = []

    monkeypatch.setattr(
        heat_pump_module,
        "plot_driving_force",
        lambda **kwargs: plotted.append(kwargs) or "driving-force-figure",
    )

    fig = BaseHeatPump().get_driving_force_plot(
        hx_obj_ls=["gc1"],
        dtmin=12.0,
        x_axis="T",
        refrigerant="Butane",
    )

    assert fig == "driving-force-figure"
    assert plotted == [
        {
            "hx_obj": ["gc1"],
            "x_axis": "T",
            "dtmin": 12.0,
            "refrigerant": "Butane",
        }
    ]
