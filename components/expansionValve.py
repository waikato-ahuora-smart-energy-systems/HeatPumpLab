__author__ = 'Tim Walmsley'
__credits__ = ['Keri Adamson']

from streams import MassStream

class ExpansionValve:
    def __init__(self, name='Expansion valve', flow_in = MassStream(), flow_out = MassStream(), specific_heat_loss=0.0, tol=1e-6):
        self.name = name
        self.flow_in = flow_in
        self.flow_out = flow_out
        self.specific_heat_loss = specific_heat_loss
        self.tol = tol
        self.power = 0.0

    def get_outlet_enthalpy(self):
        if self.flow_in.mass_flow > self.tol:
            return self.flow_in.enthalpy - self.specific_heat_loss
        else:
            print('No mass flow through the expansion valve: ' + self.name)
            return self.flow_in.enthalpy
    
    def solve(self):
        self.flow_out.enthalpy = self.get_outlet_enthalpy()
        self.flow_out.props_in_mPH()