"""Shared functionality for heat pump cycle classes."""

import numpy as np
from ...utils.graph import plot_cycle, plot_driving_force
from ...utils.readWriteOutput import print_result_to_terminal

class BaseHeatPump():
    """Base class exposing visualisation and reporting helpers for cycles."""
    gas_coolers: list
    evaporators: list
    compressors: list
    operations: list

    feasible: bool
    args: dict
    cycle_type: str = None
    converged: bool = False
    COP: np.float32 = -np.inf
    SEC: np.float32 = np.inf
    VHC: np.float32 = 0
    eDTL: np.float32 = 0
    SLE: np.float32 = 0
    dQdW: np.float32 = 0
    LMTD_gcs: np.float32 = 0

    refrigerant: float
    source: float ='Water'
    sink: float ='Air'
    tol: float =1e-6
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
    motor_eff: float =1.0 
    ihx_on: int = 1
    P_H_est: float =0.0
    dTmin_rhx : float = 5.0
    comp_model_index: int =1     

    compute_jac: bool =False
    additional_data: float =None                

    def __init__(self) -> None:
        """Initialise default attributes for derived heat pump implementations."""
        pass
   
    def print_currrent_case(self):
        """Print the current case summary to the terminal."""
        print_result_to_terminal(
            self.get_report()
            )

    def get_cycle_plot(self, x_axis: str ='h', y_axis: str ='T', refrigerant='default', include_sink=True) -> any:
        """Generate a cycle plot using the configured component graph."""
        cycle_obj, sink_obj = self.get_components_sets_to_graph()
        return plot_cycle(
            cycle_obj=cycle_obj,
            sink_obj=sink_obj,
            x_axis=x_axis,
            y_axis=y_axis,
            refrigerant=refrigerant,
            include_sink=include_sink,
        )

    def get_driving_force_plot(self, hx_obj_ls: list, dtmin: float =0.0, x_axis: str ='h', refrigerant='default') -> any:
        """Generate a driving force plot for the supplied heat exchangers."""
        return plot_driving_force(
            hx_obj=hx_obj_ls,
            x_axis=x_axis,
            dtmin=dtmin,
            refrigerant=refrigerant,
        )
