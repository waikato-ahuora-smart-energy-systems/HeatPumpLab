"""Evaporator component built on the base heat exchanger implementation."""

from .common.base_heat_exchanger import BaseHeatExchanger


class Evaporator(BaseHeatExchanger):
    """Evaporator with configurable pressure drops and operating model."""

    def __init__(
        self,
        source_fluid: str,
        sink_fluid: str,
        name="Evaporator",
        dP_h=0.0,
        dP_c=0.0,
        model="saturated gas outlet evaporator",
        bypassed=False,
        converged=False,
        tol=1e-6,
    ):
        """Initialise the evaporator streams and configuration."""
        self.name = name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.variable_dP_model = False
        self.set_port_mass_streams(source_fluid, sink_fluid)
