"""Gas cooler / condenser implementation."""

from .common.base_heat_exchanger import BaseHeatExchanger


class GasCooler(BaseHeatExchanger):
    """Gas cooler with configurable temperature approach and pressure drops."""

    def __init__(self, source_fluid: str, sink_fluid: str, name='Gas cooler', dTmin=15.0, dP_h=0.0, dP_c=0.0, model='fixed hot temperature outlet', bypassed=False, converged=False, tol=1e-6):
        """Initialise the gas cooler configuration and streams."""
        self.name = name
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.tol = tol
        self.set_port_mass_streams(source_fluid, sink_fluid)
