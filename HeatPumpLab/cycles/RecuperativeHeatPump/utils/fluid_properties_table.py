import pandas as pd
from ....utils.properties_package import PropertyPackage

PROPS_LIST = [
    "T",
    "PC",
    "TC",
    "ODP",
    "GWP",
    "ACF",
    "OFF",
    "Z",
    "DSDT_P",
    "RIEM",
    "GRUN",
]


def get_fluid_prop_overview_table(
    refrigerant_ls: list[str],
    dest_filename: str = "results/Refrigerant props.csv",
    props_list: list[str] = PROPS_LIST,
):
    """Generate a table of key properties for a list of refrigerants."""
    table = {}
    for r in refrigerant_ls:
        props = {}
        PP = PropertyPackage(r)
        for p in props_list:
            if p not in ["OFF", "DSDT_P"]:
                props[p] = PP.get_specified_properties(
                    propsOut=p,
                    propsIn="PQ",
                    prop1=0.101325,
                    prop2=0,
                )[0]
            elif p in ["DSDT_P"]:
                props[p] = PP.get_specified_properties(
                    propsOut=p,
                    propsIn="TQ",
                    prop1=15,
                    prop2=1,
                )[0]
            elif p in ["OFF"]:
                hc = PP.get_specified_properties(
                    propsOut="H",
                    propsIn="PQ",
                    prop1=props["PC"],
                    prop2=0,
                )[0]

                hl = PP.get_specified_properties(
                    propsOut="H",
                    propsIn="TQ",
                    prop1=15,
                    prop2=0,
                )[0]

                hv = PP.get_specified_properties(
                    propsOut="H",
                    propsIn="TQ",
                    prop1=15,
                    prop2=1,
                )[0]

                props[p] = (hc - (hv + hl) / 2) / (hv - hl)
            pass
        table[r] = props

    df = pd.DataFrame(table).T
    print(df)

    df.to_csv(dest_filename)
