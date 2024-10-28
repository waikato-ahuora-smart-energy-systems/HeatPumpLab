from streams import MassStream, PropertyPackage
# from timeit import default_timer as timer
import time

def test1(flds: list[str]):
    for i in range(100):
        fld: str
        for fld in flds:
            m = MassStream(name='test', 
                        fluid=PropertyPackage(fld),
                        pressure=2, #MPa
                        quality=1, 
                        mass_flow=1 #kg/s
                        )
            m.props_in_mPX()

def test2(flds:str):
    for i in range(100):
        fld: str
        for fld in flds:
            m = MassStream(name='test', 
                        fluid=PropertyPackage(fld),
                        pressure=2, #MPa
                        quality=1, 
                        mass_flow=1 #kg/s
                        )
            m.get_H_from_PX()

def test3(flds:str):
    for i in range(100):
        fld: str
        for fld in flds:
            m = MassStream(name='test', 
                        fluid=PropertyPackage(fld)
                        )

refrigerant = [
                'propane;butane|.5;.5|mol',
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
                'Methanol'
                'Ammonia',
                ]

fld = refrigerant[7]
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