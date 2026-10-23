import pandas as pd

from HeatPumpLab.cycles.RecuperativeHeatPump.utils import fluid_properties_table


def test_fluid_property_overview_writes_requested_properties(monkeypatch, tmp_path):
    calls = []

    class FakePropertyPackage:
        def __init__(self, refrigerant):
            self.refrigerant = refrigerant

        def get_specified_properties(self, propsOut, propsIn, prop1, prop2):
            calls.append((self.refrigerant, propsOut, propsIn, prop1, prop2))
            values = {
                ("PC", "PQ"): 4.0,
                ("T", "PQ"): 101.0,
                ("DSDT_P", "TQ"): 0.2,
                ("H", "PQ"): 300.0,
                ("H", "TQ", 0): 100.0,
                ("H", "TQ", 1): 500.0,
            }
            key = (
                (propsOut, propsIn, prop2)
                if propsOut == "H" and propsIn == "TQ"
                else (propsOut, propsIn)
            )
            return [values[key]]

    monkeypatch.setattr(
        fluid_properties_table,
        "PropertyPackage",
        FakePropertyPackage,
    )

    destination = tmp_path / "props.csv"
    fluid_properties_table.get_fluid_prop_overview_table(
        refrigerant_ls=["Butane"],
        dest_filename=destination,
        props_list=["PC", "T", "DSDT_P", "OFF"],
    )

    df = pd.read_csv(destination, index_col=0)
    assert df.loc["Butane", "PC"] == 4.0
    assert df.loc["Butane", "T"] == 101.0
    assert df.loc["Butane", "DSDT_P"] == 0.2
    assert df.loc["Butane", "OFF"] == 0.0
    assert ("Butane", "PC", "PQ", 0.101325, 0) in calls
    assert ("Butane", "DSDT_P", "TQ", 15, 1) in calls
