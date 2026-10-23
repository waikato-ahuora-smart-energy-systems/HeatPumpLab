"""Batch case generation for RHP outer-design optimisation."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, replace
from itertools import product
from pathlib import Path
from typing import Any

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import (
    RefrigerantInput,
    normalise_refrigerant_input,
)

from ._incumbents import attach_historical_incumbents
from ._results_io import _result_row_is_complete, _results_csv_path
from ._types import (
    DEFAULT_BASE_CASE,
    OPTIMISER_VERSION,
    OptimisationCase,
    OuterDesignOptimisationConfig,
)

_DEFAULT_REFRIGERANT = "Propane;Butane|0.5;0.5|mass"

_OPTIMISER_LIST_FIELDS = {
    "dT_rhx1_cold": (
        "dT_rhx1_cold_ls",
        "fixed_dT_rhx1_cold",
        "vary_dT_rhx1_cold",
    ),
    "dT_rhx2_cold": (
        "dT_rhx2_cold_ls",
        "fixed_dT_rhx2_cold",
        "vary_dT_rhx2_cold",
    ),
    "comp2_power_frac": (
        "comp2_power_frac_ls",
        "fixed_comp2_power_frac",
        "vary_comp2_power_frac",
    ),
}

_CASE_INPUT_KEYS = tuple(
    key
    for key in DEFAULT_BASE_CASE
    if key not in {"dT_rhx1_cold", "dT_rhx2_cold", "comp2_power_frac", "refrigerant"}
)


def _get_optimisation_case_list_to_solve(
    *,
    file_name: str,
    results_dir: str | Path,
    param_kwargs: dict[str, Any] | None = None,
    config: OuterDesignOptimisationConfig | None = None,
) -> list[OptimisationCase]:
    config = OuterDesignOptimisationConfig() if config is None else config
    cases = _build_optimisation_cases(param_kwargs=param_kwargs, config=config)
    cases = _attach_incumbents_and_refresh_signatures(
        cases,
        results_dir=results_dir,
    )
    completed = _completed_case_signatures(file_name=file_name, results_dir=results_dir)
    return [case for case in cases if case.signature not in completed]


def _build_optimisation_cases(
    *,
    param_kwargs: dict[str, Any] | None = None,
    config: OuterDesignOptimisationConfig | None = None,
) -> list[OptimisationCase]:
    config = OuterDesignOptimisationConfig() if config is None else config
    kwargs = dict(param_kwargs or {})
    refrigerants = _refrigerant_options(kwargs)
    case_options = _case_input_options(kwargs)
    optimiser_options = _optimiser_variable_options(kwargs)
    if kwargs:
        unknown = ", ".join(sorted(kwargs))
        raise ValueError(f"Unsupported optimisation case parameter(s): {unknown}")

    cases: list[OptimisationCase] = []
    case_keys = tuple(case_options)
    optimiser_keys = tuple(optimiser_options)
    for refrigerant, case_values, optimiser_values in product(
        refrigerants,
        product(*(case_options[key] for key in case_keys)),
        product(*(optimiser_options[key] for key in optimiser_keys)),
    ):
        base_case = dict(DEFAULT_BASE_CASE)
        base_case.update(dict(zip(case_keys, case_values)))
        case_config = _config_for_case(
            config,
            dict(zip(optimiser_keys, optimiser_values)),
        )
        cases.append(
            _make_case(
                base_case=base_case,
                refrigerant=refrigerant,
                config=case_config,
            )
        )
    return cases


def _build_optimisation_case(
    *,
    base_case: dict[str, Any] | None,
    refrigerant: RefrigerantInput,
    config: OuterDesignOptimisationConfig | None = None,
) -> OptimisationCase:
    complete_base_case = dict(DEFAULT_BASE_CASE)
    complete_base_case.update(base_case or {})
    return _make_case(
        base_case=complete_base_case,
        refrigerant=normalise_refrigerant_input(refrigerant).format(),
        config=OuterDesignOptimisationConfig() if config is None else config,
    )


def _make_case(
    *,
    base_case: dict[str, Any],
    refrigerant: str,
    config: OuterDesignOptimisationConfig,
) -> OptimisationCase:
    signature = _case_signature(
        base_case=base_case,
        refrigerant=refrigerant,
        config=config,
    )
    return OptimisationCase(
        base_case=base_case,
        refrigerant=refrigerant,
        config=config,
        signature=signature,
        temp_suffix=f"-{signature}",
    )


def _refrigerant_options(kwargs: dict[str, Any]) -> list[str]:
    supplied = [key for key in ("refrigerant", "refrigerant_ls") if key in kwargs]
    if len(supplied) > 1:
        raise ValueError("Use only one of refrigerant or refrigerant_ls.")
    if not supplied:
        return [_DEFAULT_REFRIGERANT]
    key = supplied[0]
    value = kwargs.pop(key)
    if key == "refrigerant":
        raw_values = [value]
    else:
        if not isinstance(value, (list, tuple)):
            raise ValueError("refrigerant_ls must be a sequence of refrigerant inputs.")
        raw_values = list(value)
    return [normalise_refrigerant_input(item).format() for item in raw_values]


def _case_input_options(kwargs: dict[str, Any]) -> dict[str, list[Any]]:
    options: dict[str, list[Any]] = {}
    for key in _CASE_INPUT_KEYS:
        list_key = f"{key}_ls"
        if list_key in kwargs:
            options[key] = _as_list(kwargs.pop(list_key))
        elif key in kwargs:
            options[key] = [kwargs.pop(key)]
        else:
            options[key] = [DEFAULT_BASE_CASE[key]]
    return options


def _optimiser_variable_options(kwargs: dict[str, Any]) -> dict[str, list[Any | None]]:
    options: dict[str, list[Any | None]] = {}
    for variable, (list_key, _fixed_attr, _vary_attr) in _OPTIMISER_LIST_FIELDS.items():
        if list_key in kwargs:
            options[variable] = _as_list(kwargs.pop(list_key))
        elif variable in kwargs:
            options[variable] = [kwargs.pop(variable)]
        else:
            options[variable] = [None]
    return options


def _config_for_case(
    base_config: OuterDesignOptimisationConfig,
    optimiser_values: dict[str, Any | None],
) -> OuterDesignOptimisationConfig:
    updates: dict[str, Any] = {}
    for variable, value in optimiser_values.items():
        if value is None:
            continue
        _list_key, fixed_attr, vary_attr = _OPTIMISER_LIST_FIELDS[variable]
        updates[fixed_attr] = float(value)
        updates[vary_attr] = False
    return replace(base_config, **updates) if updates else base_config


def _case_signature(
    *,
    base_case: dict[str, Any],
    refrigerant: str,
    config: OuterDesignOptimisationConfig,
    incumbents=(),
) -> str:
    payload = {
        "optimiser_version": OPTIMISER_VERSION,
        "base_case": base_case,
        "refrigerant": refrigerant,
        "config": asdict(config),
        "incumbents": [asdict(incumbent) for incumbent in incumbents],
    }
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


def _attach_incumbents_and_refresh_signatures(
    cases: list[OptimisationCase],
    *,
    results_dir: str | Path,
) -> list[OptimisationCase]:
    attached = attach_historical_incumbents(cases, results_dir=results_dir)
    refreshed = []
    for case in attached:
        signature = _case_signature(
            base_case=case.base_case,
            refrigerant=case.refrigerant,
            config=case.config,
            incumbents=case.incumbents,
        )
        refreshed.append(
            replace(
                case,
                signature=signature,
                temp_suffix=f"-{signature}",
            )
        )
    return refreshed


def _completed_case_signatures(
    *,
    file_name: str,
    results_dir: str | Path,
) -> set[str]:
    path = _results_csv_path(file_name=file_name, results_dir=results_dir)
    if not path.exists():
        return set()
    with path.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)
        return {
            signature
            for row in reader
            if (signature := row.get("case_signature"))
            and signature != "case_signature"
            and _result_row_is_complete(row)
        }


def _as_list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else list(value) if isinstance(value, tuple) else [value]
