"""Batch case generation for RHP outer-design optimisation."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict, replace
from itertools import product
from pathlib import Path
from typing import Any

from ._results_io import _parse_binary_refrigerant
from ._types import (
    DEFAULT_BASE_CASE,
    OPTIMISER_VERSION,
    OptimisationCase,
    OuterDesignOptimisationConfig,
)

_DEFAULT_FLUIDS = ("Propane", "Butane")

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
    "first_component_mass_fraction": (
        "first_component_mass_fraction_ls",
        "fixed_first_component_mass_fraction",
        "vary_first_component_mass_fraction",
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
    """Return remaining batch optimisation cases after completed-result filtering."""
    config = OuterDesignOptimisationConfig() if config is None else config
    cases = _build_optimisation_cases(param_kwargs=param_kwargs, config=config)
    completed = _completed_case_signatures(file_name=file_name, results_dir=results_dir)
    return [case for case in cases if case.signature not in completed]


def _build_optimisation_cases(
    *,
    param_kwargs: dict[str, Any] | None = None,
    config: OuterDesignOptimisationConfig | None = None,
) -> list[OptimisationCase]:
    """Build all batch optimisation case payloads."""
    config = OuterDesignOptimisationConfig() if config is None else config
    kwargs = dict(param_kwargs or {})
    fluid_options = _fluid_options(kwargs)
    case_options = _case_input_options(kwargs)
    optimiser_options = _optimiser_variable_options(kwargs)
    if kwargs:
        unknown = ", ".join(sorted(kwargs))
        raise ValueError(f"Unsupported optimisation case parameter(s): {unknown}")

    cases: list[OptimisationCase] = []
    case_keys = tuple(case_options)
    optimiser_keys = tuple(optimiser_options)
    for fluid_option, case_values, optimiser_values in product(
        fluid_options,
        product(*(case_options[key] for key in case_keys)),
        product(*(optimiser_options[key] for key in optimiser_keys)),
    ):
        base_case = dict(DEFAULT_BASE_CASE)
        base_case.update(dict(zip(case_keys, case_values)))
        case_config = _config_for_case(
            base_config=config,
            fluid_option=fluid_option,
            optimiser_values=dict(zip(optimiser_keys, optimiser_values)),
        )
        signature = _case_signature(
            base_case=base_case,
            fluids=fluid_option["fluids"],
            config=case_config,
        )
        cases.append(
            OptimisationCase(
                base_case=base_case,
                fluids=fluid_option["fluids"],
                config=case_config,
                signature=signature,
                temp_suffix=f"-{signature}",
            )
        )
    return cases


def _fluid_options(kwargs: dict[str, Any]) -> list[dict[str, Any]]:
    supplied = [
        key
        for key in (
            "fluids",
            "fluids_ls",
            "fluid_pairs_ls",
            "refrigerant",
            "refrigerant_ls",
        )
        if key in kwargs
    ]
    if len(supplied) > 1:
        raise ValueError(f"Use only one fluid input, got: {', '.join(supplied)}")
    if not supplied:
        return [{"fluids": _DEFAULT_FLUIDS, "first_component_mass_fraction": None}]

    key = supplied[0]
    value = kwargs.pop(key)
    if key in {"refrigerant", "refrigerant_ls"}:
        refrigerants = (
            list(value)
            if key == "refrigerant_ls" and isinstance(value, tuple)
            else _as_list(value)
        )
        options = []
        for refrigerant in refrigerants:
            fluids, composition = _parse_binary_refrigerant(str(refrigerant))
            options.append(
                {
                    "fluids": fluids,
                    "first_component_mass_fraction": composition[0],
                }
            )
        return options

    fluid_values = _as_list(value)
    return [
        {"fluids": _normalise_fluids(fluids), "first_component_mass_fraction": None}
        for fluids in fluid_values
    ]


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
    *,
    base_config: OuterDesignOptimisationConfig,
    fluid_option: dict[str, Any],
    optimiser_values: dict[str, Any | None],
) -> OuterDesignOptimisationConfig:
    updates: dict[str, Any] = {}
    fluid_fraction = fluid_option.get("first_component_mass_fraction")
    if (
        fluid_fraction is not None
        and optimiser_values["first_component_mass_fraction"] is None
    ):
        updates["fixed_first_component_mass_fraction"] = float(fluid_fraction)

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
    fluids: tuple[str, str],
    config: OuterDesignOptimisationConfig,
) -> str:
    payload = {
        "optimiser_version": OPTIMISER_VERSION,
        "base_case": base_case,
        "fluids": fluids,
        "config": asdict(config),
    }
    encoded = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()[:16]


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
        }


def _results_csv_path(*, file_name: str, results_dir: str | Path) -> Path:
    relative_path = Path(file_name + ".csv")
    if relative_path.parts and relative_path.parts[0] == "results":
        relative_path = Path(*relative_path.parts[1:])
    return Path(results_dir) / relative_path


def _normalise_fluids(fluids: Any) -> tuple[str, str]:
    if isinstance(fluids, str):
        parts = tuple(part.strip() for part in fluids.split(";"))
    else:
        parts = tuple(fluids)
    if len(parts) != 2:
        raise ValueError(f"Expected two fluids, got {fluids!r}.")
    return str(parts[0]), str(parts[1])


def _as_list(value: Any) -> list[Any]:
    if isinstance(value, list):
        return value
    if isinstance(value, tuple) and not _looks_like_fluid_pair(value):
        return list(value)
    return [value]


def _looks_like_fluid_pair(value: tuple[Any, ...]) -> bool:
    return len(value) == 2 and all(isinstance(part, str) for part in value)
