import pandas as pd
from RecuperativeHeatPump.streams import PropertyPackage

def get_fluid_prop_overview_table():
    props_ls = [
        'T',
        'PC',
        'TC',
        'ODP',
        'GWP',
        'ACF',
        'OFF',
        'Z',
        'DSDT_P',
        'RIEM',
        'GRUN',
    ]

    refrigerant_ls = [
        'Ammonia',
        'Propane', 
        'Butane', 
        'Isobutane',
        'Pentane',
        'Isopentane',
        'Neopentane',
        'Cyclopentane',
        'Carbon dioxide', 
        'R1233zd(E)', 
        'R1336mzz(E)', 
        'R1336mzz(Z)', 
        'R1234ze(E)', 
        'R1234ze(Z)', 
        'R1234yf', 
        'R1243zf', 
        'R1224yd(Z)', 
        'PROPYLEN',
        'DME',
        'Methanol',
        'propane;butane|.5;.5|mol',
        'propane;Isobutane|.5;.5|mol',
        'propane;Pentane|.5;.5|mol',
        'propane;Isopentane|.5;.5|mol',
        'propane;Neopentane|.5;.5|mol',
        'Isobutane;butane|.5;.5|mol',
        'Isobutane;Pentane|.5;.5|mol',
        'Isobutane;Isopentane|.5;.5|mol',
        'Isobutane;Neopentane|.5;.5|mol',
        'butane;Pentane|.5;.5|mol',
        'butane;Isopentane|.5;.5|mol',
        'butane;Neopentane|.5;.5|mol',
    ]

    table = {}

    for r in refrigerant_ls:
        props = {}
        PP = PropertyPackage(r)
        for p in props_ls:
            if p not in ['OFF', 'DSDT_P']:
                props[p] = PP.get_specified_properties(
                    propsOut=p,
                    propsIn='PQ',
                    prop1=0.101325,
                    prop2=0,
                )[0]
            elif p in ['DSDT_P']:
                props[p] = PP.get_specified_properties(
                    propsOut=p,
                    propsIn='TQ',
                    prop1=15,
                    prop2=1,
                )[0]
            elif p in ['OFF']:
                hc = PP.get_specified_properties(
                    propsOut='H',
                    propsIn='PQ',
                    prop1=props['PC'],
                    prop2=0,
                )[0]

                hl = PP.get_specified_properties(
                    propsOut='H',
                    propsIn='TQ',
                    prop1=15,
                    prop2=0,
                )[0]

                hv = PP.get_specified_properties(
                    propsOut='H',
                    propsIn='TQ',
                    prop1=15,
                    prop2=1,
                )[0]                

                props[p] = (hc - (hv + hl) / 2) / (hv - hl)         
            pass
        table[r] = props

    df = pd.DataFrame(table).T 
    print(df)

    df.to_csv('results/Refrigerant props.csv')