from types import SimpleNamespace

import main_RHPopt


def test_main_rhpopt_run_delegates_to_public_optimiser(monkeypatch, tmp_path):
    calls = []
    result = SimpleNamespace(
        best=SimpleNamespace(x=(1.0, 2.0, 0.4, 0.6)),
        exact_evaluations=7,
        failed_evaluations=1,
        to_summary_dict=lambda: {
            "best_objective": 0.42,
            "best_phase": "split",
        },
    )

    def fake_optimise_recuperative_heat_pump(**kwargs):
        calls.append(kwargs)
        return result

    monkeypatch.setattr(
        main_RHPopt,
        "optimise_recuperative_heat_pump",
        fake_optimise_recuperative_heat_pump,
    )
    config = main_RHPopt.OuterDesignOptimisationConfig(max_exact_evaluations=3)

    returned = main_RHPopt.run(config=config, results_dir=tmp_path)

    assert returned is result
    assert calls == [
        {
            "fluids": main_RHPopt.DEFAULT_FLUIDS,
            "config": config,
            "results_dir": tmp_path,
        }
    ]


def test_main_rhpopt_print_summary_reports_best_candidate(capsys):
    result = SimpleNamespace(
        best=SimpleNamespace(x=(1.0, 2.0, 0.4, 0.6)),
        exact_evaluations=7,
        failed_evaluations=1,
        to_summary_dict=lambda: {
            "best_objective": 0.42,
            "best_phase": "split",
        },
    )

    main_RHPopt.print_summary(result)

    output = capsys.readouterr().out
    assert "RHP optimisation complete" in output
    assert "Best SEC: 0.42" in output
    assert "Best phase: split" in output
    assert "Exact evaluations: 7" in output
    assert "dT_rhx1_cold: 1.0" in output
