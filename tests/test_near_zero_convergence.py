"""Regression coverage for delayed suction-state updates and stale pinches."""

from types import SimpleNamespace

import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.solver import (
    _cycle_analysis as cycle_analysis,
    _discharge_pressure as pressure_solver,
    _gas_cooler_pinch as pinch_solver,
)


@pytest.mark.parametrize("old_version", ["v11", "v12"])
def test_completed_cases_from_old_objective_are_not_reused(monkeypatch, tmp_path, old_version):
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import _cases

    options = dict(
        file_name="completed",
        results_dir=tmp_path,
        param_kwargs={"refrigerant": "Propane;Pentane|0.3;0.7|mass"},
    )
    with monkeypatch.context() as old_model:
        old_model.setattr(
            _cases, "OPTIMISER_VERSION", f"isolated-refprop-pattern-search-{old_version}"
        )
        [old_case] = _cases._get_optimisation_case_list_to_solve(**options)
    (tmp_path / "completed.csv").write_text(
        f"case_signature,batch_success\n{old_case.signature},True\n"
    )

    [new_case] = _cases._get_optimisation_case_list_to_solve(**options)

    assert new_case.signature != old_case.signature


@pytest.mark.parametrize("property_name", ["pressure", "enthalpy", "mass_flow"])
def test_late_compressor_inlet_update_requires_another_sweep(
    monkeypatch, property_name
):
    inlet = SimpleNamespace(pressure=0.06, enthalpy=440.0, mass_flow=1.0)
    seen = []
    compressor = SimpleNamespace(flow_in=inlet, power=10.0, converged=True)
    compressor.solve = lambda: seen.append(getattr(inlet, property_name))
    original = getattr(inlet, property_name)
    downstream = SimpleNamespace(
        converged=True,
        solve=lambda: setattr(inlet, property_name, 2 * original),
    )
    hp = SimpleNamespace(
        tol=1e-6,
        debug=False,
        compressors=[compressor],
        comp_discharge=SimpleNamespace(pressure=4.3),
        gas_cooler1=SimpleNamespace(source_out=SimpleNamespace(pressure=4.3)),
        exp_valve=SimpleNamespace(flow_in=SimpleNamespace(enthalpy=100.0)),
        gas_coolers=[SimpleNamespace(duty=30.0)],
        operations=[compressor, downstream],
    )
    monkeypatch.setattr(
        cycle_analysis, "set_gas_cooler_sink_temperatures", lambda hp: hp
    )

    cycle_analysis.run_sequential_cycle_analysis(4.3, hp)

    assert hp.converged
    assert seen == pytest.approx([original, 2 * original])


@pytest.mark.parametrize("minimum", [2.618325, 19.94, float("nan")])
def test_final_pressure_rejects_invalid_pinch_despite_component_convergence(
    monkeypatch, minimum
):
    hp = SimpleNamespace(
        feasible=True,
        iteration_converged=True,
        tol=1e-6,
        dT_GC_min_tar=20.0,
        compressors=[SimpleNamespace(power=1.0)],
        operations=[SimpleNamespace(converged=True)],
    )
    monkeypatch.setattr(
        pressure_solver, "run_sequential_cycle_analysis", lambda p, hp: hp
    )
    monkeypatch.setattr(
        pressure_solver,
        "update_gc_pinch_properties",
        lambda p, hp: setattr(hp, "gc_minimum_approach", minimum),
    )

    failure = pressure_solver.solve_at_discharge_pressure(4.3, hp)

    assert failure is not None
    assert "pinch" in failure.message.lower()
    if minimum == minimum:
        assert failure.constraint_errors["gas_cooler_pinch_shortfall"] == pytest.approx(
            20 - minimum - 0.05
        )


def test_pinch_refresh_converts_relative_location_to_duty(monkeypatch):
    gc = SimpleNamespace(duty=200.0)
    hp = SimpleNamespace(gas_coolers=[gc], gas_cooler1=gc)
    monkeypatch.setattr(
        pinch_solver,
        "get_all_gc_pinch_candidates",
        lambda p, hp: [(30.0, 1.0, True), (19.97, 0.25, True), (20.0, 0.0, True)],
    )
    duties = []
    monkeypatch.setattr(
        pinch_solver, "_set_pinch_properties", lambda duty, gc: duties.append(duty)
    )

    pinch_solver.update_gc_pinch_properties(4.3, hp)

    assert duties == [50.0]
    assert hp.gc_minimum_approach == pytest.approx(19.97)
    assert hp.gc_pinch_active == "gc1_internal"


def test_active_pinch_label_uses_the_gas_cooler_number(monkeypatch):
    # hp.gas_coolers runs hot to cold: position 0 here is gas_cooler2
    gc1, gc2 = SimpleNamespace(duty=100.0), SimpleNamespace(duty=100.0)
    hp = SimpleNamespace(gas_coolers=[gc2, gc1], gas_cooler1=gc1, gas_cooler2=gc2)
    monkeypatch.setattr(
        pinch_solver,
        "get_all_gc_pinch_candidates",
        lambda p, hp: [(19.0, 1.0, True), (25.0, 0.5, True), (30.0, 0.0, True),
                       (40.0, 1.0, True), (35.0, 0.5, True), (30.0, 0.0, True)],
    )
    monkeypatch.setattr(pinch_solver, "_set_pinch_properties", lambda duty, gc: None)

    pinch_solver.update_gc_pinch_properties(4.3, hp)

    assert hp.gc_pinch_active == "gc2_hot"


@pytest.mark.parametrize("minimum", [19.95, 20.0, 25.0])
def test_final_pressure_accepts_feasible_pinch_band(monkeypatch, minimum):
    hp = SimpleNamespace(
        feasible=True,
        iteration_converged=True,
        tol=1e-6,
        dT_GC_min_tar=20.0,
        compressors=[SimpleNamespace(power=1.0)],
        operations=[SimpleNamespace(converged=True)],
    )
    monkeypatch.setattr(
        pressure_solver, "run_sequential_cycle_analysis", lambda p, hp: hp
    )
    monkeypatch.setattr(
        pressure_solver,
        "update_gc_pinch_properties",
        lambda p, hp: setattr(hp, "gc_minimum_approach", minimum),
    )

    assert pressure_solver.solve_at_discharge_pressure(6.0, hp) is None


def _assert_final_pinch(hp):
    # Recompute the profile from final stream states, rather than trusting the
    # cached residual or the new validation attribute.
    for gc in hp.gas_coolers:
        for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
            approach = pinch_solver._pinch_temperature_difference(
                fraction * gc.duty, gc
            )
            assert approach >= hp.dT_GC_min_tar - 0.05 - 1e-8


def test_cop_remains_continuous_across_near_zero_rhx1_lifts(
    solve_recuperative_case,
    assert_successful_cycle,
):
    cops = []
    for lift in (0, 1e-6, 0.001, 0.003, 0.005809540466003726, 0.01, 0.015, 0.1):
        hp = solve_recuperative_case(
            "Propane;Pentane|0.3;0.7|mass",
            dT_rhx1_cold=lift,
            dT_rhx2_cold=31.029212817770762,
            comp2_power_frac=0.2947185847601135,
            record_failed_solve=False,
            report_failed_solve=False,
        )
        assert_successful_cycle(hp)
        _assert_final_pinch(hp)
        cops.append(hp.COP)
    assert max(cops) - min(cops) < 0.003


@pytest.mark.parametrize("seed", [0, 3, 4.3, 5.5, 7])
def test_final_near_zero_candidate_is_feasible_across_pressure_guesses(
    solve_recuperative_case,
    assert_successful_cycle,
    seed,
):
    hp = solve_recuperative_case(
        "Propane;Pentane|0.3;0.7|mass",
        dT_rhx1_cold=0.016142113039702075,
        dT_rhx2_cold=16.268979769802385,
        comp2_power_frac=0.16987922757884472,
        P_H_est=seed,
        record_failed_solve=False,
        report_failed_solve=False,
    )
    assert_successful_cycle(hp)
    _assert_final_pinch(hp)
    assert 2.66 < hp.COP < 2.67
