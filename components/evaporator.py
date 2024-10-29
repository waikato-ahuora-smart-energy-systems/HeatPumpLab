__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from components.common.commonHeatExchanger import BaseHeatExchanger

class Evaporator(BaseHeatExchanger):
    def __init__(self, source_fluid: str, sink_fluid: str, name='Evaporator', dP_h=0.0, dP_c=0.0, model='saturated gas outlet evaporator', bypassed=False, converged=False, tol=1e-6):
        self.name = name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.variable_dP_model = False
        self.set_port_mass_streams(source_fluid, sink_fluid)
        
    def get_detailed_performance(self):
        raise ValueError('No detailed evaporator model defined: ' + self.name)
