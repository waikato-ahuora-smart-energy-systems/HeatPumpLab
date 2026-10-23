from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.solver import cycle_solver
from HeatPumpLab.utils import properties_package


def test_pure_fluid_state_can_be_recreated_from_temperature_entropy(hpl):
    saturated_vapour = hpl.Stream(fluid="Butane", temperature=10.0, quality=1.0)

    recreated = hpl.Stream(
        fluid="Butane",
        temperature=saturated_vapour.temperature,
        entropy=saturated_vapour.entropy,
    )

    assert recreated.pressure == pytest.approx(saturated_vapour.pressure, abs=1e-9)
    assert recreated.enthalpy == pytest.approx(saturated_vapour.enthalpy, abs=1e-9)
    assert recreated.quality == pytest.approx(1.0, abs=1e-9)


def test_binary_fluid_exposes_temperature_glide_at_saturation(hpl):
    fluid = "Propane;Butane|0.5;0.5"
    bubble = hpl.Stream(fluid=fluid, temperature=10.0, quality=0.0)
    dew = hpl.Stream(fluid=fluid, temperature=10.0, quality=1.0)
    bubble_at_dew_pressure = hpl.Stream(fluid=fluid, pressure=dew.pressure, quality=0.0)

    assert bubble.pressure > dew.pressure
    assert dew.enthalpy > bubble.enthalpy
    assert dew.temperature - bubble_at_dew_pressure.temperature > 8.0


def test_ternary_refrigerant_list_solves_as_equal_mass_mixture(
    solve_recuperative_case,
    assert_successful_cycle,
):
    hp = solve_recuperative_case(["Propane", "Butane", "Pentane"])

    assert_successful_cycle(hp)
    assert hp.refrigerant == (
        "Propane;Butane;Pentane|0.333333333333;0.333333333333;0.333333333333|mass"
    )
    assert hp.comp1.flow_in.PP.n_components == 3


def test_mass_and_mole_fraction_inputs_describe_different_mixtures(hpl):
    mole_basis = hpl.Stream(
        fluid="Propane;Butane|0.5;0.5|mol", temperature=10.0, quality=1.0
    )
    mass_basis = hpl.Stream(
        fluid="Propane;Butane|0.5;0.5|mass", temperature=10.0, quality=1.0
    )
    equivalent_mole_basis = hpl.Stream(
        fluid="Propane;Butane|0.5681983211;0.4318016789|mol",
        temperature=10.0,
        quality=1.0,
    )

    assert mass_basis.pressure > mole_basis.pressure
    assert mass_basis.pressure == pytest.approx(
        equivalent_mole_basis.pressure, abs=2e-4
    )


def test_mixture_critical_properties_are_mixture_values(hpl):
    fluid = hpl.PropertyPackage("Propane;Butane|0.5;0.5")

    assert 120.0 < fluid.T_crit < 135.0
    assert 4.0 < fluid.P_crit < 4.5


@pytest.mark.parametrize(
    "fluid",
    [
        "Propane;Butane|0.5|mol",
        "Propane;Butane|0.5;0.5|volume",
        "Propane;Butane|0.5;-0.5|mol",
        "Propane;Butane|0;0|mol",
    ],
)
def test_invalid_mixture_descriptions_are_rejected(hpl, fluid):
    with pytest.raises(ValueError):
        hpl.PropertyPackage(fluid)


def test_mixture_description_rejects_component_fraction_count_mismatch(hpl):
    with pytest.raises(ValueError, match="composition length"):
        hpl.PropertyPackage("Propane;Butane|0.5;0.25;0.25|mol")


@pytest.mark.parametrize(
    ("input_pair", "prop1", "prop2", "expected"),
    [
        ("HP", 10.0, 2.0, (20.0, 2000.0)),
        ("SP", 3.0, 2.0, (6.0, 2000.0)),
        ("TP", 20.0, 4.0, (293.15, 4000.0)),
        ("QP", 0.5, 4.0, (0.5, 4000.0)),
        ("QT", 0.5, 20.0, (0.5, 293.15)),
        ("ST", 3.0, 20.0, (6.0, 293.15)),
        ("DT", 8.0, 20.0, (4.0, 293.15)),
        ("DP", 8.0, 4.0, (4.0, 4000.0)),
    ],
)
def test_property_package_converts_supported_si_input_pairs(
    hpl, input_pair, prop1, prop2, expected
):
    package = object.__new__(hpl.PropertyPackage)
    package.molar_mass = 2.0
    package.TC_to_K = 273.15

    assert package._convert_props_from_SI_with_C(input_pair, prop1, prop2) == (
        pytest.approx(expected)
    )


def test_property_package_rejects_unknown_si_input_pair(hpl):
    package = object.__new__(hpl.PropertyPackage)
    package.molar_mass = 2.0
    package.TC_to_K = 273.15

    with pytest.raises(Exception, match="Unrecognised"):
        package._convert_props_from_SI_with_C("BAD", 1.0, 2.0)


def test_property_package_selects_named_refprop_instance(hpl, monkeypatch):
    package = object.__new__(hpl.PropertyPackage)
    package.RP_ref = object()
    package.RP_src = object()
    package.RP_snk = object()
    selected = []

    monkeypatch.setattr(package, "_set_database_fluid", lambda: selected.append(package.RP))

    package._set_RP_lib("src")
    assert package.RP is package.RP_src
    package._set_RP_lib("sink")
    assert package.RP is package.RP_snk
    package._set_RP_lib("default")
    assert package.RP is package.RP_ref
    assert selected == [package.RP_src, package.RP_snk, package.RP_ref]


def test_refprop_platform_paths_are_explicit_and_testable():
    darwin_ref, darwin_src, darwin_snk, darwin_root = (
        properties_package._get_refprop_paths("darwin")
    )
    windows_ref, windows_src, windows_snk, windows_root = (
        properties_package._get_refprop_paths("win32")
    )

    assert darwin_ref.endswith("librefprop.dylib")
    assert darwin_src.endswith("librefprop-source.dylib")
    assert darwin_snk.endswith("librefprop-sink.dylib")
    assert darwin_root == r"/Library/REFPROP/"
    assert windows_ref.endswith("REFPROP.DLL")
    assert windows_src.endswith("REFPROP-source.DLL")
    assert windows_snk.endswith("REFPROP-sink.DLL")
    assert windows_root == "C:/Program Files (x86)/REFPROP"


def test_refprop_loader_sets_the_same_root_for_each_property_context():
    class FakeLibrary:
        def __init__(self, path):
            self.path = path
            self.roots = []

        def SETPATHdll(self, root):
            self.roots.append(root)

    rp_ref, rp_src, rp_snk = properties_package._load_refprop_libraries(
        "win32",
        library_factory=FakeLibrary,
    )

    assert rp_ref.path.endswith("REFPROP.DLL")
    assert rp_src.path.endswith("REFPROP-source.DLL")
    assert rp_snk.path.endswith("REFPROP-sink.DLL")
    assert rp_ref.roots == ["C:/Program Files (x86)/REFPROP"]
    assert rp_src.roots == ["C:/Program Files (x86)/REFPROP"]
    assert rp_snk.roots == ["C:/Program Files (x86)/REFPROP"]


def test_property_package_checks_loaded_fluid_identity(hpl, monkeypatch):
    package = object.__new__(hpl.PropertyPackage)
    package.PP_fld_names = ["PROPANE", "BUTANE"]

    monkeypatch.setattr(
        package,
        "_get_database_fluid_name",
        lambda: ["PROPANE", "BUTANE"],
    )
    assert package._check_correct_fluid() is True

    monkeypatch.setattr(package, "_get_database_fluid_name", lambda: ["BUTANE"])
    assert package._check_correct_fluid() is False


def test_pure_working_fluid_case_solves_process_targets(
    solve_recuperative_case, assert_successful_cycle
):
    hp = solve_recuperative_case("Butane")

    assert_successful_cycle(hp)
    assert hp.Tref_evap_outlet == pytest.approx(hp.Tref_evap_inlet, abs=1e-6)


def test_binary_working_fluid_case_solves_process_targets_with_glide(
    solve_recuperative_case, assert_successful_cycle
):
    hp = solve_recuperative_case("Propane;Butane|0.5;0.5")

    assert_successful_cycle(hp)
    assert hp.Tref_evap_outlet > hp.Tref_evap_inlet + 5.0


def test_recuperative_model_respects_evaporator_outlet_initial_guess(
    hpl, base_rhp_args
):
    hp = hpl.RecuperativeHP(
        **(
            base_rhp_args
            | {
                "refrigerant": "Propane;Butane|0.5;0.5|mass",
                "Tref_evap_outlet": 21.5,
            }
        )
    )

    assert hp.Tref_evap_outlet == pytest.approx(21.5)
    assert hp.evap_sup.sink_in.temperature == pytest.approx(21.5)


def test_mass_basis_binary_case_matches_equivalent_mole_basis_case(
    solve_recuperative_case, assert_successful_cycle
):
    mass_basis = solve_recuperative_case("Propane;Butane|0.5;0.5|mass")
    equivalent_mole_basis = solve_recuperative_case(
        "Propane;Butane|0.5681983211;0.4318016789|mol"
    )

    assert_successful_cycle(mass_basis)
    assert_successful_cycle(equivalent_mole_basis)
    assert mass_basis.COP == pytest.approx(equivalent_mole_basis.COP, abs=5e-3)
    assert mass_basis.comp_discharge.pressure == pytest.approx(
        equivalent_mole_basis.comp_discharge.pressure, abs=5e-2
    )


def test_high_recuperator_lift_mass_basis_binary_case_converges(
    solve_recuperative_case, assert_successful_cycle
):
    hp = solve_recuperative_case(
        "Propane;Butane|0.5;0.5|mass",
        ihx_on=False,
        dT_rhx1_cold=0.0,
        dT_rhx2_cold=60.0,
        dT_rhx3_cold=0.0,
        comp2_power_frac=0.65,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=0.0,
    )

    assert_successful_cycle(hp)
    assert hp.rhx2_on == 1
    assert hp.ihx_on == 0


def test_low_lift_mass_basis_binary_case_from_unsolved_results_converges(
    solve_recuperative_case, assert_successful_cycle
):
    hp = solve_recuperative_case(
        "Propane;Butane|0.5;0.5|mass",
        ihx_on=False,
        dT_rhx1_cold=2.5,
        dT_rhx2_cold=5.0,
        dT_rhx3_cold=0.0,
        comp2_power_frac=0.3,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=0.0,
    )

    assert_successful_cycle(hp)
    assert hp.rhx1_on == 1
    assert hp.rhx2_on == 1
    assert hp.ihx_on == 0


@pytest.mark.parametrize(
    (
        "dT_rhx2_cold",
        "comp2_power_frac",
        "dT_gc_out_min_offset",
    ),
    [
        (5.0, 0.3, 0.259133195283001),
        (7.5, 0.3, 0.050046581209243775),
        (27.5, 0.35, 0.0),
        (32.5, 0.65, 0.2),
        (35.0, 0.45, 0.6585786437626905),
    ],
)
def test_mass_basis_binary_cases_from_unsolved_results_converge(
    solve_recuperative_case,
    assert_successful_cycle,
    dT_rhx2_cold,
    comp2_power_frac,
    dT_gc_out_min_offset,
):
    hp = solve_recuperative_case(
        "Propane;Butane|0.5;0.5|mass",
        ihx_on=False,
        dT_rhx1_cold=2.5,
        dT_rhx2_cold=dT_rhx2_cold,
        dT_rhx3_cold=0.0,
        comp2_power_frac=comp2_power_frac,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=dT_gc_out_min_offset,
    )

    assert_successful_cycle(hp)
    assert hp.rhx1_on == 1
    assert hp.rhx2_on == 1
    assert hp.ihx_on == 0


@pytest.mark.parametrize(
    ("refrigerant", "dT_rhx1_cold", "dT_rhx2_cold"),
    [
        ("Butane;Neopentane|0.5;0.5|mass", 0.0, 66.66666666666667),
    ],
)
def test_solver_resolves_feasible_binary_cases_from_unsolved_csv(
    solve_recuperative_case,
    assert_successful_cycle,
    refrigerant,
    dT_rhx1_cold,
    dT_rhx2_cold,
):
    hp = solve_recuperative_case(
        refrigerant,
        ihx_on=False,
        dT_rhx1_cold=dT_rhx1_cold,
        dT_rhx2_cold=dT_rhx2_cold,
        dT_rhx3_cold=0.0,
        comp2_power_frac=0.5,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=0.0,
    )

    assert_successful_cycle(hp)


@pytest.mark.parametrize(
    "refrigerant",
    [
        "Propane;Isopentane|0.5;0.5|mass",
        "Propane;Pentane|0.5;0.5|mass",
    ],
)
def test_binary_cases_accept_colder_expansion_outlet_at_superheater_limit(
    solve_recuperative_case,
    assert_successful_cycle,
    refrigerant,
):
    hp = solve_recuperative_case(
        refrigerant,
        ihx_on=False,
        dT_rhx1_cold=0.0,
        dT_rhx2_cold=0.0,
        dT_rhx3_cold=0.0,
        comp2_power_frac=0.5,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=0.0,
        record_failed_solve=False,
    )

    assert_successful_cycle(hp)
    assert hp.evaporator.sink_in.temperature < hp.Tref_evap_inlet
    assert hp.evap_sup.duty >= 0.0
    assert hp.Tref_evap_outlet == pytest.approx(
        hp.T_source_in - hp.dTmin_rhx,
        abs=2e-4,
    )


def test_mixture_phase_envelope_removes_false_internal_temperature_cross(
    solve_recuperative_case,
    monkeypatch,
):
    monkeypatch.setattr(cycle_solver, "print_result_to_csv", lambda *args, **kwargs: None)
    hp = solve_recuperative_case(
        "Isopentane;Neopentane|0.5;0.5|mass",
        ihx_on=False,
        dT_rhx1_cold=0.0,
        dT_rhx2_cold=61.80339887498948,
        dT_rhx3_cold=0.0,
        comp2_power_frac=0.5,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=0.0,
    )

    assert hp.feasible is True
    assert hp.converged is True
    # The old -53 K spike used a PH error-248 output whose enthalpy missed
    # its input by 334 kJ/kg. Check the whole profile, not just solver status.
    import numpy as np
    from HeatPumpLab.cycles.RecuperativeHeatPump.solver._gas_cooler_pinch import (
        _pinch_temperature_difference,
    )

    for exchanger, minimum in [(gc, 19.95) for gc in hp.gas_coolers] + [(hp.rhx2, 5.0)]:
        exchanger.source_pinch.pressure = exchanger.source_out.pressure
        exchanger.sink_pinch.pressure = exchanger.sink_in.pressure
        approaches = [
            _pinch_temperature_difference(float(duty), exchanger)
            for duty in np.linspace(0, exchanger.duty, 501)
        ]
        assert min(approaches) >= minimum - 1e-6


def test_recuperated_binary_case_uses_superheater_limited_suction_pressure(
    solve_recuperative_case,
    assert_successful_cycle,
    monkeypatch,
):
    monkeypatch.setattr(cycle_solver, "print_result_to_csv", lambda *args, **kwargs: None)
    hp = solve_recuperative_case(
        "Propane;Isopentane|0.5;0.5|mass",
        ihx_on=False,
        dT_rhx1_cold=33.333333333333336,
        dT_rhx2_cold=0.0,
        dT_rhx3_cold=0.0,
        comp2_power_frac=0.5,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=0.0,
    )

    assert_successful_cycle(hp)
    assert hp.evaporator.sink_in.temperature < hp.Tref_evap_inlet
    assert hp.evap_sup.duty >= 0.0


@pytest.mark.parametrize(
    ("refrigerant", "dT_rhx1_cold", "dT_rhx2_cold"),
    [
        (
            "Propane;Pentane|0.5;0.5|mass",
            41.1111111111111,
            82.22222222222221,
        ),
        (
            "Propane;Pentane|0.5;0.5|mass",
            47.10914138751296,
            76.22419194582037,
        ),
        (
            "Propane;Isopentane|0.5;0.5|mass",
            56.23058987490536,
            34.75241575014722,
        ),
    ],
)
def test_outer_design_cases_solve_at_superheater_limited_suction_pressure(
    solve_recuperative_case,
    assert_successful_cycle,
    refrigerant,
    dT_rhx1_cold,
    dT_rhx2_cold,
):
    hp = solve_recuperative_case(
        refrigerant,
        ihx_on=False,
        dT_rhx1_cold=dT_rhx1_cold,
        dT_rhx2_cold=dT_rhx2_cold,
        dT_rhx3_cold=0.0,
        comp2_power_frac=0.5,
        comp3_power_frac=0.0,
        dT_gc_out_min_offset=0.0,
        record_failed_solve=False,
    )

    assert_successful_cycle(hp)
    assert hp.evaporator.sink_in.temperature <= hp.Tref_evap_inlet
    assert hp.evap_sup.duty >= 0.0
    assert hp.comp_discharge.pressure < hp.comp_discharge.P_max


@pytest.mark.parametrize(
    "refrigerant",
    [
        "Propane;Cyclopentane|0.5;0.5|mass",
        "Butane;Cyclopentane|0.5;0.5|mass",
        "Isobutane;Cyclopentane|0.25;0.75|mass",
    ],
)
def test_large_glide_mixture_initializes_water_before_evaporator_alignment(
    hpl, base_rhp_args, assert_successful_cycle, refrigerant
):
    hp = hpl.RecuperativeHP(
        **(
            base_rhp_args
            | {
                "refrigerant": refrigerant,
                "comp2_power_frac": 0.1,
                "record_failed_solve": False,
                "report_failed_solve": False,
            }
        )
    )

    # The unaligned refrigerant state would give a subzero water seed.
    assert hp.evaporator.sink_in.temperature + hp.dTmin_rhx < 0.0
    assert hp.evaporator.source_out.temperature == pytest.approx(hp.T_source_in)

    hp.solve()

    assert_successful_cycle(hp)
    source_in = hp.evaporator.source_in
    source_out = hp.evaporator.source_out
    assert 0.0 < source_out.temperature < source_in.temperature
    assert source_out.mass_flow * (
        source_in.enthalpy - source_out.enthalpy
    ) == pytest.approx(hp.evaporator.duty, rel=1e-7)


def test_three_stage_recuperative_model_initialises_optional_heat_recovery_path(
    hpl, base_rhp_args
):
    hp = hpl.RecuperativeHP(
        **(
            base_rhp_args
            | {
                "refrigerant": "Butane",
                "n_comp": 3,
                "comp2_power_frac": 0.3,
                "comp3_power_frac": 0.2,
                "dT_rhx1_cold": 5.0,
                "dT_rhx2_cold": 6.0,
                "dT_rhx3_cold": 7.0,
                "dT_gc_out_min_offset": 20.0,
                "P_H_est": 5.0,
            }
        )
    )

    assert hp.feasible is True
    assert hp.n_comp == 3
    assert len(hp.compressors) == 3
    assert len(hp.gas_coolers) == 4
    assert hp.rhx1_on == 1
    assert hp.rhx2_on == 1
    assert hp.rhx3_on == 1
    assert hp.ihx_on == 1
    assert hp.comp_discharge is hp.comp3.flow_out
    assert hp.gas_cooler4.source_in is hp.comp_discharge
    assert hp.rhx3.source_in is hp.gas_cooler4.source_out
    assert hp.rhx2.source_in is hp.gas_cooler3.source_out
    assert hp.rhx1.source_in is hp.gas_cooler2.source_out
    assert hp.ihx.source_in is hp.gas_cooler1.source_out
    assert hp.exp_valve.flow_in is hp.ihx.source_out
    assert hp.comp_discharge.pressure == pytest.approx(5.0)


def test_discharge_pressure_estimate_uses_pressure_cap_when_scan_has_no_root(
    hpl, base_rhp_args, monkeypatch
):
    from HeatPumpLab.cycles.RecuperativeHeatPump import rhp_helpers

    monkeypatch.setattr(
        rhp_helpers,
        "root_scalar",
        lambda **kwargs: (_ for _ in ()).throw(RuntimeError("no bracket")),
    )

    hp = hpl.RecuperativeHP(
        **(
            base_rhp_args
            | {
                "refrigerant": "Butane",
                "n_comp": 1,
                "comp2_power_frac": 0.0,
                "comp3_power_frac": 0.0,
                "ihx_on": False,
                "P_H_est": 0.0,
            }
        )
    )

    assert hp.feasible is True
    assert hp.comp_discharge.pressure == pytest.approx(
        min(20.0, hp.comp1.flow_in.pressure * 36.0)
    )


@pytest.mark.parametrize(
    ("fallback_result", "expected_pressure"),
    [
        (SimpleNamespace(converged=True, root=7.0), 7.0),
        (SimpleNamespace(converged=False, root=7.0), 20.0),
        (RuntimeError("fallback failed"), 20.0),
    ],
)
def test_discharge_pressure_estimate_uses_critical_pressure_fallback(
    hpl, base_rhp_args, monkeypatch, fallback_result, expected_pressure
):
    from HeatPumpLab.cycles.RecuperativeHeatPump import rhp_helpers

    monkeypatch.setattr(
        rhp_helpers.np,
        "linspace",
        lambda *args, **kwargs: (_ for _ in ()).throw(RuntimeError("scan failed")),
    )

    def root_scalar_fallback(**kwargs):
        if isinstance(fallback_result, Exception):
            raise fallback_result
        kwargs["f"](kwargs["x0"])
        return fallback_result

    monkeypatch.setattr(rhp_helpers, "root_scalar", root_scalar_fallback)

    hp = hpl.RecuperativeHP(
        **(
            base_rhp_args
            | {
                "refrigerant": "Butane",
                "n_comp": 1,
                "comp2_power_frac": 0.0,
                "comp3_power_frac": 0.0,
                "ihx_on": False,
                "P_H_est": 0.0,
            }
        )
    )

    assert hp.feasible is True
    assert hp.comp_discharge.pressure == pytest.approx(expected_pressure)


def test_infeasible_design_is_reported_without_exception(hpl, base_rhp_args):
    hp = hpl.RecuperativeHP(
        **(
            base_rhp_args
            | {
                "dT_rhx1_cold": 90.0,
                "dT_rhx2_cold": 90.0,
                "refrigerant": "Butane",
            }
        )
    )

    hp.solve()

    assert hp.feasible is False
    assert hp.converged is False
