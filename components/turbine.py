__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams import MassStream

class Turbine: 
    def __init__(self, name='Turbine', eff_isen = 1.0, eff_gen = 1.0):
        self.name = name
        self.flow_in = MassStream()
        self.flow_out = MassStream()
        self.eff_isen = eff_isen
        self.eff_gen = eff_gen
        self.dh = 0.0
        self.pressure_ratio = 1.0
        self.power = 0.0

    def solve(self):
        self.flow_out.entropy = self.flow_in.entropy 
        self.flow_out.props_in_mPS()
        self.dh_isen = self.flow_in.enthalpy - self.flow_out.enthalpy
        self.dh = self.dh_isen / self.eff_isen
        self.flow_out.enthalpy = self.flow_in.enthalpy - self.dh 
        self.flow_out.props_in_mPH()
        self.power = self.dh * self.flow_out.mass_flow * self.eff_gen
        