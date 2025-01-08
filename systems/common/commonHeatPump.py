__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

from functions.graph import plot_cycle, plot_driving_force
from functions.readWriteOutput import print_result_to_terminal
import numpy as np

class HeatPump():
    gas_coolers: list
    evaporators: list
    compressors: list
    operations: list

    feasible: bool
    cycle_type: str = None
    converged: bool = False
    COP: np.float32 = -np.inf
    SEC: np.float32 = np.inf
    VHC: np.float32 = 0
    eDTL: np.float32 = 0
    SLE: np.float32 = 0
    LMTD_gcs: np.float32 = 0

    refrigerant: float
    source: float ='Water'
    sink: float ='Air'
    tol: float =1e-8
    debug: float =False
    output: dict

    T_air_tar: float =200.0 
    T_air_inlet: float =15.0 
    dT_GC_min_tar: float =20.0 
    dT_gc_out_min_offset: float =0.0
    T_gc_out: float 
    T_source_in: float =40.0 
    Tref_evap_inlet: float =10.0
    Tref_evap_outlet: float =10.0    
    eff_motor: float =0.96 
    ihx_on: int = 1
    P_H_est: float =0.0
    min_dT_ihx : float = 5.0
    comp_model_index: int =1     

    compute_jac: bool =False
    additional_data: float =None                

    def __init__(self) -> None:
        pass
   
    def print_currrent_case(self):
        print_result_to_terminal(
            self.get_summary()
            )

    def get_cycle_plot(self, x_axis: str ='h', y_axis: str ='T') -> any:
        cycle_obj, sink_obj = self.get_components_sets_to_graph()
        return plot_cycle(cycle_obj=cycle_obj,
                sink_obj=sink_obj,
                x_axis=x_axis,
                y_axis=y_axis)

    def get_driving_force_plot(self, hx_obj_ls: list, dtmin: float =0.0, x_axis: str ='h') -> any:
        return plot_driving_force(hx_obj=hx_obj_ls,
                x_axis=x_axis,
                dtmin=dtmin)
