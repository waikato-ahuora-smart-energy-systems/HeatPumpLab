"""Isentropic turbine expander model."""

from .stream import Stream


class Turbine:
    """Simple turbine with specified isentropic and generator efficiency."""
    name: str
    flow_in: Stream
    flow_out: Stream
    flow_out_isen: Stream
    eff_gen: float
    eff_isen: float
    pressure_ratio: float =1.0
    power: float =0.0 
    dh: float =0.0
    converged: bool = False

    def __init__(self, name: str ='Turbine', fluid: str =None, eff_isen: float =1.0, eff_gen: float =1.0) -> None:
        """Initialise the turbine streams and efficiency parameters."""
        self.name = name
        self.flow_in = Stream(name='flow_in', fluid=fluid)
        self.flow_out = Stream(name='flow_out', fluid=fluid)
        self.eff_isen = eff_isen
        self.eff_gen = eff_gen

    def solve(self) -> None:
        """Solve the turbine outlet state using the isentropic relation."""
        self.flow_out.entropy = self.flow_in.entropy 
        self.flow_out.props_in_mPS()
        self.dh_isen = self.flow_in.enthalpy - self.flow_out.enthalpy
        self.dh = self.dh_isen / self.eff_isen
        self.flow_out.enthalpy = self.flow_in.enthalpy - self.dh 
        self.flow_out.props_in_mPH()
        self.power = self.dh * self.flow_out.mass_flow * self.eff_gen
        self.converged = True
        
