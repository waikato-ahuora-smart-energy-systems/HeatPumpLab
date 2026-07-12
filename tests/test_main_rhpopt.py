import main_RHPopt


def test_main_rhpopt_run_many_delegates_to_public_batch_optimiser(
    monkeypatch,
    tmp_path,
):
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

    returned = main_RHPopt.run_many(config=config, n_cpu=1, results_dir=tmp_path)

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


def test_main_rhpopt_run_one_delegates_to_single_case_optimiser(
    monkeypatch,
    tmp_path,
):
    calls = []
    result = object()

    def fake_optimise_recuperative_heat_pump(**kwargs):
        calls.append(kwargs)
        return result

    monkeypatch.setattr(
        main_RHPopt,
        "optimise_recuperative_heat_pump",
        fake_optimise_recuperative_heat_pump,
    )
    base_case = {"T_air_tar": 373.15}
    fluids = ("Propane", "Isobutane")
    config = main_RHPopt.OuterDesignOptimisationConfig(max_exact_evaluations=3)

    returned = main_RHPopt.run_one(
        base_case,
        fluids=fluids,
        config=config,
        results_dir=tmp_path,
    )

    assert returned is result
    assert calls == [
        {
            "base_case": base_case,
            "fluids": fluids,
            "config": config,
            "results_dir": tmp_path,
        }
    ]


def test_main_rhpopt_run_one_uses_default_config(monkeypatch):
    calls = []

    def fake_optimise_recuperative_heat_pump(**kwargs):
        calls.append(kwargs)
        return object()

    monkeypatch.setattr(
        main_RHPopt,
        "optimise_recuperative_heat_pump",
        fake_optimise_recuperative_heat_pump,
    )

    main_RHPopt.run_one()

    assert calls[0]["base_case"] is None
    assert calls[0]["fluids"] == main_RHPopt.DEFAULT_FLUIDS
    assert calls[0]["config"] == main_RHPopt.OuterDesignOptimisationConfig()
    assert calls[0]["results_dir"] == main_RHPopt.DEFAULT_RESULTS_DIR
