__author__ = 'Tim Walmsley'
__credits__ = ['Keri Adamson']

from streams import MassStream, get_property_value
from math import log, exp

class ExpansionValve:
    def __init__(self, name='Expansion valve', flow_in = MassStream(), flow_out = MassStream(), specific_heat_loss=0.0, tol=1e-6, RP=None):
        self.name = name
        self.flow_in = flow_in
        self.flow_out = flow_out
        self.specific_heat_loss = specific_heat_loss
        self.tol = tol
        self.power = 0.0
        self.RP = RP

    def get_outlet_enthalpy(self):
        if self.flow_in.mass_flow > self.tol:
            return self.flow_in.enthalpy - self.specific_heat_loss
        else:
            print('No mass flow through the expansion valve: ' + self.name)
            return self.flow_in.enthalpy
    
    def solve(self):
        self.flow_out.enthalpy = self.get_outlet_enthalpy()
        self.flow_out.props_in_mPH()

    def get_graph_data(self, divisions=50):
        data = {'siso': 
                    {'P': [], 'h': [], 'T': [], 's': [], 'v': [], 'H': []}   
                }

        for i in range(divisions + 1):
            dlog_P = (log(self.flow_in.pressure) - log(self.flow_out.pressure)) * i / divisions
            data['siso']['P'].append(exp(log(self.flow_in.pressure) - dlog_P))
            data['siso']['h'].append(self.flow_in.enthalpy)
            data['siso']['T'].append(    get_property_value('T', self.flow_out.fluid, 'PH', data['siso']['P'][-1] / 1000, data['siso']['h'][-1], self.RP))
            data['siso']['s'].append(    get_property_value('S', self.flow_out.fluid, 'PH', data['siso']['P'][-1] / 1000, data['siso']['h'][-1], self.RP))
            data['siso']['v'].append(1 / get_property_value('D', self.flow_out.fluid, 'PH', data['siso']['P'][-1] / 1000, data['siso']['h'][-1], self.RP))
            data['siso']['H'].append(0.0)
        
        return data