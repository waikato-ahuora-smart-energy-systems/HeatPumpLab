__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from functions import single_case_analysis, massive_parallel_case_analysis, unsolved_case_analysis

def get_analysis_parameters():
        # Refrigerants analysed
    refrigerant = [
                    # 'propane;butane|.5;.5|mol',
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

    # Cycle parameters and their ranges tested
    dT_GC_min_tar = [20.0]
    T_air_tar = [float(k) for k in range(250, 99, -5)] 
    T_air_inlet = [15.0] #, 35.0] 
    T_source_in = [40.0]
    n_comp = [2, 1]
    Tsat_evap = [10.0]
    comp_PR_uniform = [True] #[True, False]
    w = [float(o/100) for o in range(100, 40, -10)]

    return {
            'refrigerant' : refrigerant,
            'dT_GC_min_tar' : dT_GC_min_tar,
            'T_air_tar' : T_air_tar,
            'T_air_inlet' : T_air_inlet,
            'T_source_in' : T_source_in,
            'n_comp' : n_comp,
            'Tsat_evap' : Tsat_evap,
            'comp_PR_uniform' : comp_PR_uniform,
            'w' : w,
           }

if __name__=='__main__':
    # single_case_analysis()
    # unsolved_case_analysis()
    parameters = get_analysis_parameters()
    massive_parallel_case_analysis(parameters)
    None   