__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley, Keri Adamson and Alex Geary, University of Waikato, 2021 (tim.walmsley@waikato.ac.nz)
"""
from timeit import default_timer as timer
#from classes.HTTHP import HTTHP
from functions.optimiseHTTHP import optimiser_function, HTTHP

if __name__=='__main__':
    TIC = timer()

    Refrigerant = ['R1233zd(E)', 'R1234yf','R1234ze(Z)', 'R1243zf','R1336mzz(Z)', 'R1224yd(Z)', 'R1234ze(E)', 'Carbon dioxide', 'Ammonia', 'Isobutane', 'Butane']
    dT_GC_min_tar = [i for i in range(15,30,5)]
    T_air_tar = [i for i in range(180,220,10)]
    T_air_inlet = [i for i in range(10,30,5)]
    T_air_exh = 65.0
    n_comp = [i + 1 for i in range(3)]

    if 0:
        [optimiser_function(indep_var=[Refrigerant[j], dT_GC_min_tar[1], T_air_tar[2], T_air_inlet[1], T_air_exh, n_comp[1]]) for j in range(len(Refrigerant))]
    else:
        m = HTTHP(Refrigerant = Refrigerant[2],
                dT_GC_min_tar = dT_GC_min_tar[1],
                T_air_tar = T_air_tar[2],
                T_air_inlet = T_air_inlet[1],
                T_air_exh = T_air_exh,
                n_comp = n_comp[1],
                dT_phx1_max = 0.0,
                dT_phx2_max = 0.0,
                dT_phx3_max = 0.0, 
                PR_2_init = 1.0, 
                PR_3_init = 1.0,
                print_output=True)
        m.solve()

    TOC = timer()
    print('Total time:',TOC - TIC) 