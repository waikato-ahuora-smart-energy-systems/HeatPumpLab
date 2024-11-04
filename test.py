from streams import MassStream, PropertyPackage
from functions.decorator_functions import time_func
import time

@time_func
def test1(flds: list[str]):
    for _ in range(100):
        fld: str
        for fld in flds:
            m = MassStream(name='test', 
                        fluid=fld,
                        pressure=1, #MPa
                        temperature=200, #degC
                        mass_flow=1 #kg/s
                        )

@time_func
def test2(flds:str):
    for _ in range(100):
        fld: str
        for fld in flds:
            m = MassStream(name='test', 
                        fluid=fld,
                        pressure=1, #MPa
                        temperature=200, #degC
                        mass_flow=1 #kg/s
                        )

@time_func
def test3(flds:str):
    for _ in range(100):
        fld: str
        for fld in flds:
            m = MassStream(name='test', 
                        fluid=fld
                        )
 
refrigerant_ls: list[str] = [
        'Butane', 
        'Carbon dioxide', 
        'Propane', 
        'Isobutane', 
        'R1233zd(E)', 
        'R1336mzz(Z)', 
        'R1234ze(E)', 
        'R1234ze(Z)', 
        'R1234yf', 
        'R1243zf', 
        'R1224yd(Z)', 
        'PROPYLEN',
        'DME',
        'Methanol',
        'Ammonia',
        'propane;butane|.5;.5|mol',
        'Methanol;butane|.5;.5|mol',
        'R1234ze(E);R1336mzz(Z)|.5;.5|mol',
        'propane;Isobutane|.5;.5|mol',
        'Propane;DME|.5;.5|mol',
        ]

dt = test3(refrigerant_ls)[1]
print(dt)
print(test2(refrigerant_ls)[1] - dt)
print(test1(refrigerant_ls)[1] - dt)
print(test2(refrigerant_ls)[1] - dt)
print(test1(refrigerant_ls)[1] - dt)
