"""Generic plate heat exchanger component."""

from .common.base_heat_exchanger import BaseHeatExchanger


class HeatExchanger(BaseHeatExchanger):
    """Configurable heat exchanger wrapper around the base implementation."""

    def __init__(self, source_fluid: str, sink_fluid: str, name: str='Heat exchanger', 
                 HX_Eff: float =0.9, dTmin: float =5, dP_h: float =0.0, dP_c: float =0.0, 
                 model: str ='dTmin_terminal', bypassed: bool =False, converged: bool =False, 
                 tol: float =1e-6) -> None:
        """Initialise the exchanger with effectiveness and pressure data."""
        self.name = name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.HX_Eff = HX_Eff
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.set_port_mass_streams(source_fluid, sink_fluid)
