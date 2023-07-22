__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes.components.common.commonHeatExchanger import BaseHeatExchanger

class Evaporator(BaseHeatExchanger):

    def __init__(self, Name='Evaporator', dP_h=0.0, dP_c=0.0, model='fixed cold duty', bypassed=False, converged=False, tol=1e-6):
        self.Name = Name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.variable_dP_model = False
        self.init_stream_connections()
        self.init_performance_metrics()        
        
    
    def get_detailed_performance(self):
        raise ValueError('No detailed evaporator model defined: ' + self.Name)
