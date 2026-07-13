"""Abstract base support for heat exchanger component models."""

from math import log
from ..stream import Stream


class BaseHeatExchanger:
    """Shared heat-exchanger behaviour extended by concrete component models."""

    # Set the stream connection points
    name: str = ""
    model: str
    source_fluid: str
    source_in: Stream
    source_pinch: Stream
    source_out: Stream
    source_out_max: Stream
    sink_fluid: str
    sink_in: Stream
    sink_pinch: Stream
    sink_out: Stream
    sink_out_max: Stream
    pinch_length = None
    bypassed: bool = False
    converged: bool = False
    HX_Eff: float = 0.0
    NTU: float = 0.0
    c: float = 0.0
    U: float = 0.0
    Area: float = 0.0
    duty: float = 0.0
    delta_entropy_src: float = 0.0
    delta_entropy_snk: float = 0.0
    dP_h: float = 0.0
    dP_c: float = 0.0
    variable_dP_model: bool = False
    tol: float = 1e-6

    def __init__(self) -> None:
        """Initialise placeholder attributes for subclasses."""
        pass

    def set_port_mass_streams(self, source_fluid: str, sink_fluid: str) -> None:
        """Create inlet, pinch, and outlet streams for both sides."""
        self.source_in = Stream(name="source in", fluid=source_fluid)
        self.source_pinch = Stream(name="source pinch", fluid=source_fluid)
        self.source_out = Stream(name="source out", fluid=source_fluid)
        self.source_out_max = Stream(name="source out max", fluid=source_fluid)

        self.sink_in = Stream(name="sink in", fluid=sink_fluid)
        self.sink_pinch = Stream(name="sink pinch", fluid=sink_fluid)
        self.sink_out = Stream(name="sink out", fluid=sink_fluid)
        self.sink_out_max = Stream(name="sink out max", fluid=sink_fluid)

    def solve(self):
        """Evaluate the exchanger according to the selected performance model."""
        if self.bypassed:
            self.source_out = self.source_in
            self.sink_out = self.sink_in
            self.converged = True
        else:
            if (
                self.source_in.mass_flow <= self.tol
                or self.sink_in.mass_flow <= self.tol
            ):
                self.source_out.copy_fld_props(fld=self.source_in)
                self.source_out.set_flows(self.source_in.mass_flow)
                self.sink_out.copy_fld_props(fld=self.sink_in)
                self.sink_out.set_flows(self.sink_in.mass_flow)
                self.converged = True

            else:
                self.get_performance()

        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow
        if self._has_active_heat_transfer():
            self.converged = self.converged and self._has_valid_heat_transfer_direction()
        self.delta_entropy_src = (
            abs(self.source_in.entropy - self.source_out.entropy)
            * self.source_out.mass_flow
        )
        self.delta_entropy_snk = (
            abs(self.sink_in.entropy - self.sink_out.entropy) * self.sink_out.mass_flow
        )

    def _has_active_heat_transfer(self) -> bool:
        """Return whether the exchanger should enforce heat-flow direction."""
        return (
            not self.bypassed
            and self.source_in.mass_flow > self.tol
            and self.sink_in.mass_flow > self.tol
        )

    def _has_valid_heat_transfer_direction(self) -> bool:
        """Return whether heat flows from the source side to the sink side."""
        return (
            self.duty >= -self.tol
            and self.source_out.temperature <= self.source_in.temperature + self.tol
            and self.sink_out.temperature >= self.sink_in.temperature - self.tol
        )

    def get_performance(self):
        """Dispatch to the underlying heat exchanger formulation."""
        if self.model == "dTmin_terminal":
            self.get_dTmin_terminal_model()

        elif self.model == "effectiveness":
            self.get_effectiveness_model()

        elif self.model == "fixed cold duty":
            self.get_fixed_cold_duty_model()

        elif self.model == "fixed hot duty":
            self.get_fixed_hot_duty_model()

        elif self.model == "fixed cold temperature outlet":
            self.get_fixed_cold_T_outlet_model()

        elif self.model == "fixed hot temperature outlet":
            self.get_fixed_hot_T_outlet_model()

        elif self.model == "fixed hot quality outlet":
            self.get_fixed_hot_Q_outlet_model()

        elif self.model == "fixed cold delta temperature":
            self.get_fixed_cold_delta_T_model()

        elif self.model == "fixed hot delta temperature":
            self.get_fixed_hot_delta_T_model()

        elif self.model == "saturated gas outlet evaporator":
            self.get_saturated_gas_sink_out_model()

        elif self.model == "specialised":
            self.get_detailed_performance()

        else:
            raise ValueError("No heat exchanger model selected.")

    def get_dTmin_terminal_model(self):
        """Solve for the duty limited by the specified terminal delta T."""
        try:
            self.source_out.pressure = self.source_in.pressure - self.dP_h
            self.source_out_max.props_in_ab(
                "PT", self.source_out.pressure, self.sink_in.temperature + self.dTmin
            )
            self.dH_h_max = self.source_in.mass_flow * (
                self.source_in.enthalpy - self.source_out_max.enthalpy
            )

            self.sink_out.pressure = self.sink_in.pressure - self.dP_c
            self.sink_out_max.props_in_ab(
                "PT", self.sink_out.pressure, self.source_in.temperature - self.dTmin
            )
            self.dH_c_max = self.sink_in.mass_flow * (
                self.sink_out_max.enthalpy - self.sink_in.enthalpy
            )

            self.duty = min(self.dH_h_max, self.dH_c_max)

            self.source_out.enthalpy = (
                self.source_in.enthalpy - self.duty / self.source_in.mass_flow
            )
            self.source_out.props_in_mPH()
            self.sink_out.enthalpy = (
                self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
            )
            self.sink_out.props_in_mPH()
        except Exception as exc:
            raise ValueError(
                "Warning: Phase change occuring in heat exchanger " + self.name
            ) from exc

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_effectiveness_model(self):
        """Apply an effectiveness definition to compute energy transfer."""
        try:
            self.dT_max = abs(self.source_in.temperature - self.sink_in.temperature)

            self.source_out.pressure = self.source_in.pressure - self.dP_h
            self.source_out_max.pressure = self.source_out.pressure
            self.source_out_max.temperature = (
                self.sink_in.temperature + self.dT_max * self.HX_Eff
            )
            self.source_out_max.props_in_mPT()
            self.dH_h_max = self.source_in.mass_flow * (
                self.source_in.enthalpy - self.source_out_max.enthalpy
            )

            self.sink_out.pressure = self.sink_in.pressure - self.dP_c
            self.sink_out_max.pressure = self.sink_out.pressure
            self.sink_out_max.temperature = (
                self.source_in.temperature - self.dT_max * self.HX_Eff
            )
            self.sink_out_max.props_in_mPT()
            self.dH_c_max = self.sink_in.mass_flow * (
                self.sink_out_max.enthalpy - self.sink_in.enthalpy
            )

            self.duty = min(self.dH_h_max, self.dH_c_max)

            self.source_out.enthalpy = (
                self.source_in.enthalpy - self.duty / self.source_in.mass_flow
            )
            self.source_out.props_in_mPH()
            self.sink_out.enthalpy = (
                self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
            )
            self.sink_out.props_in_mPH()

        except Exception as exc:
            raise ValueError(
                "Warning: Phase change occuring in heat exchanger " + self.name
            ) from exc

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_fixed_cold_duty_model(self):
        """Resolve state changes when cold-side duty is fixed."""
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.props_in_mPH()

        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = (
            self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        )
        self.source_out.props_in_mPH()

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_fixed_hot_duty_model(self):
        """Resolve state changes when hot-side duty is fixed."""
        self.duty = self.source_in.heat_flow - self.source_out.heat_flow
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.props_in_mPH()

        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.enthalpy = (
            self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
        )
        self.sink_out.props_in_mPH()

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_fixed_cold_T_outlet_model(self):
        """Solve with a fixed cold outlet temperature constraint."""
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.props_in_mPT()
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow

        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = (
            self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        )
        self.source_out.props_in_mPH()

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_fixed_hot_T_outlet_model(self):
        """Solve with a fixed hot outlet temperature constraint."""
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.props_in_mPT()
        self.duty = self.source_in.heat_flow - self.source_out.heat_flow

        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.enthalpy = (
            self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
        )
        self.sink_out.props_in_mPH()

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_fixed_hot_Q_outlet_model(self):
        """Solve with a fixed hot outlet vapour quality constraint."""
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.props_in_mPQ()
        self.duty = self.source_in.heat_flow - self.source_out.heat_flow

        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.enthalpy = (
            self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
        )
        self.sink_out.props_in_mPH()

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_fixed_cold_delta_T_model(self):
        """Resolve states when enforcing a cold-side temperature rise."""
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.temperature = self.sink_in.temperature + self.dT_c
        self.sink_out.props_in_mPT()
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow

        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = (
            self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        )
        self.source_out.props_in_mPH()

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_fixed_hot_delta_T_model(self):
        """Resolve states when enforcing a hot-side temperature drop."""
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.temperature = self.source_in.temperature - self.dT_h
        self.source_out.props_in_mPT()
        self.duty = self.source_in.heat_flow - self.source_out.heat_flow

        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.enthalpy = (
            self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
        )
        self.sink_out.props_in_mPH()

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_saturated_gas_sink_out_model(self):
        """Evaluate an evaporator with saturated-vapour outlet constraint."""
        self.sink_out.props_in_mTQ()
        self.sink_in.pressure = self.sink_out.pressure + self.dP_c
        self.sink_in.props_in_mPH()
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow

        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = (
            self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        )
        self.source_out.props_in_mPH()

        # if self.source_out.temperature < self.sink_in.temperature or self.source_in.temperature < self.sink_out.temperature:
        #      print('temperature cross occuring in the heat exchanger: ' + self.name)

        self.converged = (
            True
            if self.source_out.temperature > self.sink_in.temperature
            and self.source_in.temperature > self.sink_out.temperature
            else False
        )

    def get_graph_data(self, divisions=50):
        """Sample temperature profiles along the exchanger length."""
        data = {"source": [], "sink": []}

        def get_inter_data(s_in: Stream, s_out: Stream):
            """Interpolate intermediate stream states for plotting."""
            ls = []
            dP_tot = s_in.pressure - s_out.pressure
            dh_tot = s_in.enthalpy - s_out.enthalpy

            for i in range(divisions + 1):
                if i == 0:
                    m_str = Stream(
                        name="source i",
                        fluid=s_in.fluid,
                        pp_instance=s_in.pp_instance,
                        mass_flow=s_in.mass_flow,
                    )
                    m_str.copy_fld_props(s_in)
                    m_str.set_flows(s_in.mass_flow)
                    ls.append(m_str)
                    continue
                if i == divisions:
                    m_str = Stream(
                        name="source i",
                        fluid=s_out.fluid,
                        pp_instance=s_out.pp_instance,
                        mass_flow=s_out.mass_flow,
                    )
                    m_str.copy_fld_props(s_out)
                    m_str.set_flows(s_out.mass_flow)
                    ls.append(m_str)
                    continue
                m_str = Stream(
                    name="source i",
                    fluid=s_in.fluid,
                    pp_instance=s_in.pp_instance,
                    pressure=s_in.pressure - dP_tot * i / divisions,
                    enthalpy=s_in.enthalpy - dh_tot * i / divisions,
                    mass_flow=s_in.mass_flow,
                )
                m_str.props_in_mPH()
                ls.append(m_str)
            return ls

        data["sink"] = get_inter_data(self.sink_in, self.sink_out)
        data["source"] = get_inter_data(self.source_in, self.source_out)
        return data

    def get_LMTD_ave(self):
        """Average the segment-wise log-mean temperature difference."""
        data = self.get_graph_data()
        div = len(data["source"]) - 1
        LMTD_seg = []
        for i in range(div):
            LMTD_seg.append(
                _get_LMTD(
                    Thin=data["source"][div - i].temperature,
                    Thout=data["source"][div - 1 - i].temperature,
                    Tcin=data["sink"][i].temperature,
                    Tcout=data["sink"][i + 1].temperature,
                )
            )
        return sum(LMTD_seg) / len(LMTD_seg)


def _get_LMTD(Thin, Thout, Tcin, Tcout):
    """Compute log mean temperature difference for a segment."""
    dt_hot = Thin - Tcout
    dt_cold = Thout - Tcin
    if dt_hot <= 0.0 or dt_cold <= 0.0:
        raise ValueError("LMTD is undefined because a temperature cross occurs.")
    if dt_hot == dt_cold:
        return dt_hot
    return (dt_hot - dt_cold) / log(dt_hot / dt_cold)
