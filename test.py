from streams import MassStream, PropertyPackage
import time

def test1(flds: list[str]):
    for i in range(100):
        fld: str
        for fld in flds:
            m1 = MassStream(name='test', 
                        fluid=PropertyPackage(fld),
                        temperature=50, #MPa
                        quality=0.5, #kJ/kg 
                        mass_flow=1 #kg/s
                        )
            m1.props_in_mTQ_old()

def test2(flds:str):
    for i in range(100):
        fld: str
        for fld in flds:
            m1 = MassStream(name='test', 
                        fluid=PropertyPackage(fld),
                        pressure=1.4612098550237957, #MPa
                        enthalpy=555.6308939805065, #kJ/kg 
                        mass_flow=1 #kg/s
                        )
            m1.props_in_mPH()
            m2 = MassStream(name='test', 
                        fluid=PropertyPackage(fld),
                        pressure=1.4612098550237957, #MPa
                        enthalpy=555.6308939805065, #kJ/kg 
                        mass_flow=1 #kg/s
                        )            
            m2.props_in_mPH_old()
            None

def test3(flds:str):
    for i in range(100):
        fld: str
        for fld in flds:
            m = MassStream(name='test', 
                        fluid=PropertyPackage(fld)
                        )

refrigerant = [
                # 'propane;butane|.5;.5|mol',
                # 'water',
                'R1336mzz(Z)', 
                'Air',
                'Butane', 
                'Carbon dioxide', 
                'Propane', 
                'Isobutane', 
                'R1233zd(E)', 
                'R1234ze(E)', 
                'R1234ze(Z)', 
                'R1234yf', 
                'R1243zf', 
                'R1224yd(Z)', 
                'PROPYLEN',
                'DME',
                'Methanol',
                'Ammonia',
                ]

fld = refrigerant
test3(fld)

tic = time.perf_counter()
test3(fld)
toc = time.perf_counter()
test3_time = toc - tic
print(test3_time)

tic = time.perf_counter()
test2(fld)
toc = time.perf_counter()
test2_time = toc - tic
print(test2_time - test3_time)

tic = time.perf_counter()
test1(fld)
toc = time.perf_counter()
test1_time = toc - tic
print(test1_time - test3_time)