__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams import get_property_value, MassStream

class Separator:
    def __init__(self, name='Separator', dP=0.0, RP=None):
        self.name = name
        self.flow_in = MassStream()
        self.vap_out = MassStream()
        self.liquid_out = MassStream()
        self.dP = dP
        self.RP = RP

    def solve(self):
        self.vap_out.pressure = self.flow_in.pressure - self.dP
        self.vap_out.quality = 1.0
        self.liquid_out.pressure = self.flow_in.pressure - self.dP
        self.liquid_out.quality = 0.0

        vap_frac = get_property_value('QMASS', self.flow_in.fluid, 'PH', self.vap_out.pressure / 1000, self.flow_in.enthalpy, self.RP)
        
        self.vap_out.mass_flow = self.flow_in.mass_flow * vap_frac
        self.vap_out.props_in_mPX()
        
        self.liquid_out.mass_flow = self.flow_in.mass_flow - self.vap_out.mass_flow
        self.liquid_out.props_in_mPX()
        