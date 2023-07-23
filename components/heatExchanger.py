__author__ = 'Tim Walmsley'
__credits__ = ['']

from common.commonHeatExchanger import BaseHeatExchanger

class HeatExchanger(BaseHeatExchanger):
    def __init__(self, Name='Heat exchanger', HX_Eff=0.9, dTmin=5, dP_h=0.0, dP_c=0.0, model='dTmin_terminal', bypassed=False, converged=False, tol=1e-6):
        self.Name = Name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.HX_Eff = HX_Eff
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.variable_dP_model = False
        self.init_stream_connections()
        self.init_performance_metrics()
