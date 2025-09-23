from .common.baseHeatExchanger import BaseHeatExchanger

class GasCooler(BaseHeatExchanger):

    def __init__(self, source_fluid: str, sink_fluid: str, name='Gas cooler', dTmin=15.0, dP_h=0.0, dP_c=0.0, model='fixed hot temperature outlet', bypassed=False, converged=False, tol=1e-6):
        self.name = name
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.tol = tol
        self.set_port_mass_streams(source_fluid, sink_fluid)
