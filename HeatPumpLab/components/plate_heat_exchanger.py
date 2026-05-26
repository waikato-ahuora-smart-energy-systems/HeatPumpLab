"""Plate heat exchanger specialised settings."""

from .common.base_heat_exchanger import BaseHeatExchanger


class PlateHeatExchanger(BaseHeatExchanger):
    """Plate heat exchanger used for recuperators and internal HX."""

    def __init__(
        self,
        source_fluid: str,
        sink_fluid: str,
        name="Plate heat exchanger",
        hx_eff=0.0,
        dTmin=0.0,
        dT_h=0.0,
        dT_c=0.0,
        dP_h=0.0,
        dP_c=0.0,
        model="dTmin_terminal",
        bypassed=False,
        converged=False,
        tol=1e-6,
    ):
        """Initialise plate heat exchanger configuration."""
        self.name = name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.hx_eff = hx_eff
        self.dTmin = dTmin
        self.dT_h = dT_h
        self.dT_c = dT_c
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.variable_dP_model = False
        self.set_port_mass_streams(source_fluid, sink_fluid)
