__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

from functions.graph import Graph
from functions.readWriteOutput import print_result_to_terminal
import numpy as np


class HeatPump():
    gas_coolers: list
    compressors: list
    feasible: bool
    cycle_type: str = None
    converged: bool = False
    COP: np.float32 = -np.inf
    refrigerant: float
    source: float ='Water'
    sink: float ='Air'
    tol: float =1e-8
    debug: float =False
    output: dict

    def __init__(self) -> None:
        pass
   
    def print_currrent_case(self):
        print_result_to_terminal(
            self.get_summary()
            )

    def get_cycle_plot(self, x_axis: str ='h', y_axis: str ='T') -> Graph:
        cycle_obj, sink_obj = self.get_components_sets_to_graph()
        g = Graph()
        g.plot_cycle(cycle_obj=cycle_obj,
                     sink_obj=sink_obj,
                     x_axis=x_axis,
                     y_axis=y_axis)
        return g

    def get_driving_force_plot(self, hx_obj_ls: list, dtmin: float =0.0, x_axis: str ='h') -> Graph:
        g: Graph = Graph()
        g.plot_driving_force(hx_obj=hx_obj_ls,
                             x_axis=x_axis,
                             dtmin=dtmin)
        return g
