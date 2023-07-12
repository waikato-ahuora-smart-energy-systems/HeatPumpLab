__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley, Keri Adamson and Alex Geary, University of Waikato, 2021 (tim.walmsley@waikato.ac.nz)
"""

from timeit import default_timer as timer
from functions.optimise import optimiser_function, TRHP

if __name__=='__main__':
    TIC = timer()

    Refrigerant = ['R1233zd(E)', 'R1234yf','R1234ze(Z)', 'R1243zf','R1336mzz(Z)', 'R1224yd(Z)', 'R1234ze(E)', 'Carbon dioxide', 'Isobutane', 'Butane', 'Propane']
    dT_GC_min_tar = [j for j in range(15,35,5)]
    T_air_tar = [k for k in range(100,250,10)]
    T_air_inlet = [l for l in range(10,30,5)]
    T_air_exh = [65.0]
    n_comp = [n + 1 for n in range(3)]

    if 0:
        [optimiser_function(indep_var=[Refrigerant[0], dT_GC_min_tar[2], T_air_tar[12], T_air_inlet[1], T_air_exh[0], n_comp[1]]) ] # for j in range(len(Refrigerant))
    else:
        m = TRHP(
                    Refrigerant = Refrigerant[0],
                    dT_GC_min_tar = dT_GC_min_tar[2],
                    T_air_tar = T_air_tar[12],
                    T_air_inlet = T_air_inlet[1],
                    T_air_exh = T_air_exh[0],
                    n_comp = n_comp[1],
                    dT_phx1_max = 56.90464836084587,
                    dT_phx2_max = 46.15299827209145,
                    dT_phx3_max = 10.0, 
                    PR_2_init = 6.250023469515739, 
                    PR_3_init = 1.2,
                    print_output=True,
                    graph_output=True
                )
        m.solve()

    TOC = timer()
    print('Total time:',TOC - TIC) 