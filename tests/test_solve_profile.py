from types import SimpleNamespace

import numpy as np
import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.solver import cycle_solver
from HeatPumpLab.utils import analysis as heat_analysis
from HeatPumpLab.utils.properties_package import PropertyPackage
from HeatPumpLab.utils.solve_profile import profile_solve


def test_profile_context_records_named_sections():
    with profile_solve() as profile:
        with profile.section("custom"):
            pass

    report = profile.as_dict()
    assert report["counts"]["custom"] == 1
    assert report["timings"]["custom"] >= 0.0


def test_profile_context_records_property_package_calls():
    class FakeRefprop:
        def ABFLSHdll(self, *args):
            return SimpleNamespace(P=1000.0, h=20.0, q=-1.0, T=300.15, s=4.0, Cp=6.0, D=2.0)

        def REFPROPdll(self, *args):
            return SimpleNamespace(Output=[1.0, 2.0, 3.0])

    package = object.__new__(PropertyPackage)
    package.RP = FakeRefprop()
    package.fld_name = "Butane"
    package.fld_composition = [1.0]
    package.UNITS = 0
    package.iMass = 0
    package.iFlag = 0
    package.molar_mass = 2.0
    package.TC_to_K = 273.15
    package.DEFAULT_VAL = -9999990.0

    with profile_solve() as profile:
        package.get_std_properties("PT", 1.0, 20.0, [1.0])
        package.get_specified_properties("PRANDTL;VIS;TCX", "PH", 1.0, 2.0)

    report = profile.as_dict()
    assert report["counts"]["refprop_std"] == 1
    assert report["counts"]["refprop_specified"] == 1
    assert report["timings"]["refprop_std"] >= 0.0
    assert report["timings"]["refprop_specified"] >= 0.0


def test_profile_context_records_recuperative_metric_timing():
    hp = SimpleNamespace(
        converged=True,
        feasible=True,
        gas_coolers=[
            SimpleNamespace(
                duty=10.0,
                delta_entropy_src=2.0,
                get_LMTD_ave=lambda: 5.0,
            )
        ],
        evaporators=[SimpleNamespace(duty=4.0, delta_entropy_snk=1.0)],
        compressors=[SimpleNamespace(power=2.0)],
        comp1=SimpleNamespace(flow_in=SimpleNamespace(vol_flow=0.5)),
        comp_discharge=SimpleNamespace(pressure=5.0, P_crit=4.0),
    )

    with profile_solve() as profile:
        cycle_solver._get_metrics(hp)

    report = profile.as_dict()
    assert report["counts"]["metrics"] == 1
    assert report["timings"]["metrics"] >= 0.0
    assert hp.SEC == pytest.approx(0.25)


def test_profile_context_records_pressure_search_timing(monkeypatch):
    hp = SimpleNamespace(comp_discharge=SimpleNamespace(pressure_bnds=None))

    monkeypatch.setattr(
        cycle_solver,
        "_get_discharge_pressure_bounds",
        lambda solved_hp: ((1.0, 2.0), 1.5, True, (1.0, -1.0)),
    )
    monkeypatch.setattr(
        cycle_solver,
        "_find_pressure_for_gc_target",
        lambda bounds, x0, dt_bounds, solved_hp: (1.5, True),
    )
    monkeypatch.setattr(
        cycle_solver,
        "_solve_cycle_at_discharge_pressure",
        lambda pressure, solved_hp: solved_hp,
    )

    with profile_solve() as profile:
        cycle_solver._find_hp_with_min_discharge_pressure(hp)

    report = profile.as_dict()
    assert report["counts"]["pressure_search"] == 1
    assert report["timings"]["pressure_search"] >= 0.0
    assert hp.comp_discharge.pressure_bnds == pytest.approx((1.0, 2.0))


def test_profile_context_records_reporting_timing(monkeypatch):
    reports = []
    hp = SimpleNamespace(output={}, get_report=lambda: hp.output.update({"COP": 3.2}))

    monkeypatch.setattr(heat_analysis, "get_unique_f_suffix", lambda param: "-Butane")
    monkeypatch.setattr(
        heat_analysis,
        "print_results",
        lambda **kwargs: reports.append(kwargs),
    )

    with profile_solve() as profile:
        heat_analysis.get_results_report(
            hp=hp,
            param={"additional_data": {}},
            result={"nfev": 2, "message": "ok", "success": True},
            obj=0.4,
            dt=1.0,
        )

    report = profile.as_dict()
    assert report["counts"]["reporting"] == 1
    assert report["timings"]["reporting"] >= 0.0
    assert reports[0]["results"]["COP"] == pytest.approx(3.2)


def test_profile_context_records_objective_evaluations(monkeypatch):
    def fake_solve_single_case_with_timer(hp_cls, args):
        return SimpleNamespace(SEC=args["candidate"]), 0.0

    def fake_minimize(fun, x0, bounds, method, tol):
        fun(np.array([0.75]))
        fun(np.array([0.5]))
        return SimpleNamespace(success=True)

    monkeypatch.setattr(
        heat_analysis,
        "solve_single_case_with_timer",
        fake_solve_single_case_with_timer,
    )
    monkeypatch.setattr(heat_analysis, "loc_min", fake_minimize)

    with profile_solve() as profile:
        heat_analysis.get_local_optimum(
            HP_fn=object,
            param={"additional_data": {"vars_input": {"candidate": {"xi": 0.0}}}},
            vars_input={"candidate": {"xi": 0.0, "bnds": (0.0, 1.0)}},
            tol=1e-4,
            f_obj=heat_analysis.compute_SEC_obj,
            method="COBYQA",
        )

    report = profile.as_dict()
    assert report["counts"]["objective_eval"] == 2
    assert report["timings"]["objective_eval"] >= 0.0
