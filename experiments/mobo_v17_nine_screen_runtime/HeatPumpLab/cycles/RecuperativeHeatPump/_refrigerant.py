"""Refrigerant input normalisation for recuperative heat-pump APIs."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import TypeAlias

import numpy as np

RefrigerantInput: TypeAlias = (
    str
    | list[str]
    | tuple[str, ...]
    | list[tuple[str, float]]
    | tuple[tuple[str, float], ...]
)


@dataclass(frozen=True)
class RefrigerantSpec:
    """Canonical components, composition, and caller-facing basis."""

    components: tuple[str, ...]
    fractions: tuple[float, ...]
    basis: str

    @property
    def component_count(self) -> int:
        return len(self.components)

    def format(self, fractions: tuple[float, ...] | None = None, basis: str | None = None) -> str:
        values = self.fractions if fractions is None else _normalise_fractions(fractions)
        if len(values) != self.component_count:
            raise ValueError("Fluid composition length does not match component count.")
        output_basis = self.basis if basis is None else _validate_basis(basis)
        if self.component_count == 1:
            return self.components[0]
        names = ";".join(self.components)
        composition = ";".join(f"{value:.12g}" for value in values)
        return f"{names}|{composition}|{output_basis}"


def normalise_refrigerant_input(refrigerant: RefrigerantInput) -> RefrigerantSpec:
    """Return a validated RHP refrigerant specification."""
    if isinstance(refrigerant, str):
        return _parse_refrigerant_string(refrigerant)
    if not isinstance(refrigerant, (list, tuple)) or not refrigerant:
        raise ValueError("Refrigerant must be a non-empty string, list, or tuple.")

    if all(isinstance(component, str) for component in refrigerant):
        components = _normalise_components(refrigerant)
        fractions = tuple(1.0 / len(components) for _ in components)
        return RefrigerantSpec(components, fractions, "mass")

    if all(
        isinstance(item, tuple)
        and len(item) == 2
        and isinstance(item[0], str)
        for item in refrigerant
    ):
        components = _normalise_components([item[0] for item in refrigerant])
        try:
            fractions = _normalise_fractions(tuple(float(item[1]) for item in refrigerant))
        except (TypeError, ValueError) as exc:
            raise ValueError("Structured refrigerant fractions must be numeric.") from exc
        return RefrigerantSpec(components, fractions, "mass")

    raise ValueError(
        "Refrigerant sequences must contain only names or only (name, fraction) pairs."
    )


def convert_composition(
    fractions: tuple[float, ...] | list[float] | np.ndarray,
    *,
    from_basis: str,
    to_basis: str,
    molar_masses: tuple[float, ...] | list[float] | np.ndarray,
) -> tuple[float, ...]:
    """Convert a normalised composition between mass and mole basis."""
    source_basis = _validate_basis(from_basis)
    target_basis = _validate_basis(to_basis)
    composition = np.asarray(_normalise_fractions(tuple(fractions)), dtype=float)
    masses = np.asarray(molar_masses, dtype=float)
    if masses.shape != composition.shape or np.any(~np.isfinite(masses)) or np.any(masses <= 0.0):
        raise ValueError("Molar masses must be finite, positive, and match the composition.")
    if source_basis == target_basis:
        return tuple(float(value) for value in composition)
    amounts = composition / masses if source_basis == "mass" else composition * masses
    converted = amounts / amounts.sum()
    return tuple(float(value) for value in converted)


def _parse_refrigerant_string(refrigerant: str) -> RefrigerantSpec:
    parts = [part.strip() for part in refrigerant.split("|")]
    if len(parts) > 3 or not parts[0]:
        raise ValueError(f"Invalid refrigerant description {refrigerant!r}.")
    components = _normalise_components(parts[0].split(";"))
    if len(parts) == 1:
        if len(components) == 1:
            return RefrigerantSpec(components, (1.0,), "mass")
        fractions = tuple(1.0 / len(components) for _ in components)
        return RefrigerantSpec(components, fractions, "mass")

    try:
        fractions = _normalise_fractions(tuple(float(value) for value in parts[1].split(";")))
    except ValueError as exc:
        raise ValueError(f"Invalid refrigerant composition in {refrigerant!r}.") from exc
    if len(fractions) != len(components):
        raise ValueError("Fluid composition length does not match component count.")
    basis = "mol" if len(parts) == 2 else _validate_basis(parts[2])
    return RefrigerantSpec(components, fractions, basis)


def _normalise_components(components) -> tuple[str, ...]:
    values = tuple(str(component).strip() for component in components)
    if not values or any(not value for value in values):
        raise ValueError("Refrigerant component names cannot be empty.")
    if len({value.casefold() for value in values}) != len(values):
        raise ValueError("Refrigerant component names must be unique.")
    return values


def _normalise_fractions(fractions) -> tuple[float, ...]:
    values = tuple(float(value) for value in fractions)
    if not values or any(not math.isfinite(value) for value in values):
        raise ValueError("Fluid composition fractions must be finite.")
    if any(value < 0.0 for value in values):
        raise ValueError("Fluid composition cannot contain negative fractions.")
    total = sum(values)
    if total <= 0.0:
        raise ValueError("Fluid composition must contain at least one positive fraction.")
    return tuple(value / total for value in values)


def _validate_basis(basis: str) -> str:
    value = str(basis).strip().lower()
    if value not in {"mass", "mol"}:
        raise ValueError("Fluid composition basis must be 'mol' or 'mass'.")
    return value
