__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes.unit_operations.common.commonHeatExc import BaseHeatExc

class Evaporator(BaseHeatExc):

    def __init__(self, Name='Evaporator', dP_h=0.0, dP_c=0.0, model='fixed cold duty', bypassed=False, tol=1e-8):
        self.Name = Name
        self.model = model
        self.bypassed = bypassed
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.variable_dP_model = False
        self.init_stream_connections()
        self.init_performance_metrics()        
        
    
    def solve(self):       
        try:
            self.get_base_performance(model=self.model)
        
        except:
            raise ValueError('Evaporator model failed: ' + self.Name)
