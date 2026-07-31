import main_RHPopt


def test_main_rhpopt_default_param_kwargs_include_pure_binary_and_ternary_options():
    refrigerants = main_RHPopt.DEFAULT_PARAM_KWARGS["refrigerant_ls"]

    assert refrigerants[:9] == [
        "Propane",
        "Isobutane",
        "Butane",
        "Isopentane",
        "Neopentane",
        "Pentane",
        "CO2",
        "Water",
        "Ammonia",
    ]
    assert ("Propane", "Butane") in refrigerants
    assert ("Propane", "Isobutane", "Butane") in refrigerants
    assert ("Water", "Ammonia") in refrigerants
    assert refrigerants.index(("Water", "Ammonia")) < refrigerants.index(
        ("Propane", "Isobutane", "Butane")
    )
    assert len(refrigerants) == 66
    assert len(main_RHPopt.DEFAULT_PARAM_KWARGS["T_air_tar_ls"]) == 25


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


def test_main_rhpopt_run_many_can_reverse_default_case_order(monkeypatch, tmp_path):
    calls = []

    def fake_optimise_recuperative_heat_pump_cases(**kwargs):
        calls.append(kwargs)
        return []

    monkeypatch.setattr(
        main_RHPopt,
        "optimise_recuperative_heat_pump_cases",
        fake_optimise_recuperative_heat_pump_cases,
    )
    original_refrigerants = list(main_RHPopt.DEFAULT_PARAM_KWARGS["refrigerant_ls"])
    original_air_targets = list(main_RHPopt.DEFAULT_PARAM_KWARGS["T_air_tar_ls"])

    main_RHPopt.run_many(reverse=True, results_dir=tmp_path)

    assert calls[0]["param_kwargs"]["refrigerant_ls"] == list(
        reversed(original_refrigerants)
    )
    assert calls[0]["param_kwargs"]["T_air_tar_ls"] == list(
        reversed(original_air_targets)
    )
    assert main_RHPopt.DEFAULT_PARAM_KWARGS["refrigerant_ls"] == original_refrigerants
    assert main_RHPopt.DEFAULT_PARAM_KWARGS["T_air_tar_ls"] == original_air_targets


def test_main_rhpopt_run_one_delegates_to_single_case_optimiser(
    monkeypatch,
    tmp_path,
):
    calls = []
    result = {"case_signature": "abc", "best_objective": 0.42}

    def fake_optimise_recuperative_heat_pump_case(**kwargs):
        calls.append(kwargs)
        return result

    monkeypatch.setattr(
        main_RHPopt,
        "optimise_recuperative_heat_pump_case",
        fake_optimise_recuperative_heat_pump_case,
    )
    base_case = {"T_air_tar": 373.15}
    refrigerant = "Propane;Isobutane|0.5;0.5|mass"
    config = main_RHPopt.OuterDesignOptimisationConfig(max_exact_evaluations=3)

    returned = main_RHPopt.run_one(
        base_case,
        refrigerant=refrigerant,
        config=config,
        n_cpu=2,
        results_dir=tmp_path,
    )

    assert returned is result
    assert calls == [
        {
            "filename": main_RHPopt.DEFAULT_FILENAME,
            "base_case": base_case,
            "refrigerant": refrigerant,
            "config": config,
            "n_cpu": 2,
            "results_dir": tmp_path,
        }
    ]


def test_main_rhpopt_run_one_uses_default_config(monkeypatch):
    calls = []

    def fake_optimise_recuperative_heat_pump_case(**kwargs):
        calls.append(kwargs)
        return {"case_signature": "abc"}

    monkeypatch.setattr(
        main_RHPopt,
        "optimise_recuperative_heat_pump_case",
        fake_optimise_recuperative_heat_pump_case,
    )

    main_RHPopt.run_one()

    assert calls[0]["filename"] == main_RHPopt.DEFAULT_FILENAME
    assert calls[0]["base_case"] is None
    assert calls[0]["refrigerant"] == main_RHPopt.DEFAULT_REFRIGERANT
    assert calls[0]["config"] == main_RHPopt.OuterDesignOptimisationConfig()
    assert calls[0]["n_cpu"] >= 1
    assert calls[0]["results_dir"] == main_RHPopt.DEFAULT_RESULTS_DIR
