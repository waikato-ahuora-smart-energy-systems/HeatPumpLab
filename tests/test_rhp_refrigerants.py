import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump._refrigerant import (
    convert_composition,
    normalise_refrigerant_input,
)


def test_component_name_list_defaults_to_equal_mass_split():
    spec = normalise_refrigerant_input(["Propane", "Butane", "Pentane"])

    assert spec.components == ("Propane", "Butane", "Pentane")
    assert spec.fractions == pytest.approx((1 / 3, 1 / 3, 1 / 3))
    assert spec.basis == "mass"
    assert spec.format() == (
        "Propane;Butane;Pentane|0.333333333333;0.333333333333;0.333333333333|mass"
    )


def test_component_name_tuple_defaults_to_equal_mass_split():
    spec = normalise_refrigerant_input(("Propane", "Butane", "Pentane"))

    assert spec.components == ("Propane", "Butane", "Pentane")
    assert spec.fractions == pytest.approx((1 / 3, 1 / 3, 1 / 3))
    assert spec.basis == "mass"


def test_structured_component_list_normalises_mass_fractions():
    spec = normalise_refrigerant_input(
        [("Propane", 2.0), ("Butane", 3.0), ("Pentane", 5.0)]
    )

    assert spec.fractions == pytest.approx((0.2, 0.3, 0.5))
    assert spec.basis == "mass"


def test_structured_component_tuple_normalises_mass_fractions():
    spec = normalise_refrigerant_input(
        (("Propane", 2.0), ("Butane", 3.0), ("Pentane", 5.0))
    )

    assert spec.fractions == pytest.approx((0.2, 0.3, 0.5))
    assert spec.basis == "mass"


def test_equivalent_mass_input_forms_share_canonical_description():
    inputs = [
        "Propane;Butane;Pentane|0.2;0.3;0.5|mass",
        [("Propane", 0.2), ("Butane", 0.3), ("Pentane", 0.5)],
    ]

    assert {
        normalise_refrigerant_input(refrigerant).format() for refrigerant in inputs
    } == {"Propane;Butane;Pentane|0.2;0.3;0.5|mass"}


def test_component_only_string_defaults_to_equal_mass_split():
    spec = normalise_refrigerant_input("Propane;Butane;Pentane")

    assert spec.fractions == pytest.approx((1 / 3, 1 / 3, 1 / 3))
    assert spec.basis == "mass"


def test_explicit_composition_without_basis_retains_mole_default():
    spec = normalise_refrigerant_input("Propane;Butane|0.4;0.6")

    assert spec.fractions == pytest.approx((0.4, 0.6))
    assert spec.basis == "mol"


def test_mass_and_mole_composition_conversion_round_trips():
    mole = (0.2, 0.3, 0.5)
    molar_masses = (44.0, 58.0, 72.0)

    mass = convert_composition(
        mole,
        from_basis="mol",
        to_basis="mass",
        molar_masses=molar_masses,
    )
    recreated = convert_composition(
        mass,
        from_basis="mass",
        to_basis="mol",
        molar_masses=molar_masses,
    )

    assert recreated == pytest.approx(mole)


@pytest.mark.parametrize(
    "refrigerant",
    [
        [],
        ["Propane", ("Butane", 0.5)],
        [("Propane", 0.5), ("propane", 0.5)],
        [("Propane", -0.1), ("Butane", 1.1)],
        [("Propane", 0.0), ("Butane", 0.0)],
    ],
)
def test_invalid_structured_refrigerants_are_rejected(refrigerant):
    with pytest.raises(ValueError):
        normalise_refrigerant_input(refrigerant)


def test_recuperative_hp_normalises_list_before_model_construction(monkeypatch):
    from HeatPumpLab.cycles.RecuperativeHeatPump import rhp_base

    captured = []
    monkeypatch.setattr(rhp_base, "set_rhp_model", lambda hp: captured.append(hp.refrigerant))

    hp = rhp_base.RecuperativeHP(refrigerant=["Propane", "Butane", "Pentane"])

    assert hp.refrigerant == (
        "Propane;Butane;Pentane|0.333333333333;0.333333333333;0.333333333333|mass"
    )
    assert captured == [hp.refrigerant]
