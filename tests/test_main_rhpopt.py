import main_RHPopt


def test_main_rhpopt_run_delegates_to_public_batch_optimiser(monkeypatch, tmp_path):
    calls = []
    result = [
        {
            "case_signature": "abc",
            "best_objective": 0.42,
            "best_phase": "split",
        }
    ]

    def fake_optimise_recuperative_heat_pump_cases(**kwargs):
        calls.append(kwargs)
        return result

    monkeypatch.setattr(
        main_RHPopt,
        "optimise_recuperative_heat_pump_cases",
        fake_optimise_recuperative_heat_pump_cases,
    )
    config = main_RHPopt.OuterDesignOptimisationConfig(max_exact_evaluations=3)

    returned = main_RHPopt.run(config=config, n_cpu=1, results_dir=tmp_path)

    assert returned is result
    assert calls == [
        {
            "filename": main_RHPopt.DEFAULT_FILENAME,
            "param_kwargs": main_RHPopt.DEFAULT_PARAM_KWARGS,
            "config": config,
            "n_cpu": 1,
            "results_dir": tmp_path,
        }
    ]


def test_main_rhpopt_print_summary_reports_best_candidate(capsys):
    result = [
        {
            "best_objective": 0.42,
            "best_phase": "split",
            "exact_evaluations": 7,
            "failed_evaluations": 1,
            "best_dT_rhx1_cold": 1.0,
            "best_dT_rhx2_cold": 2.0,
            "best_comp2_power_frac": 0.4,
            "best_first_component_mass_fraction": 0.6,
            "batch_success": True,
        },
    ]

    main_RHPopt.print_summary(result)

    output = capsys.readouterr().out
    assert "RHP optimisation batch complete" in output
    assert "Cases completed: 1" in output
    assert "Best SEC: 0.42" in output
    assert "Best phase: split" in output
    assert "Exact evaluations: 7" in output
    assert "dT_rhx1_cold: 1.0" in output


def test_main_rhpopt_print_summary_handles_empty_batch(capsys):
    main_RHPopt.print_summary([])

    output = capsys.readouterr().out
    assert "RHP optimisation batch complete" in output
    assert "Cases completed: 0" in output
