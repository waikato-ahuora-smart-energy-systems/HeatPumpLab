"""Compressor component model."""

import math

from .stream import Stream

# Pareto front of published compressor efficiency vs pressure ratio,
# eta(PR) = eta_1 - a * ln(PR)**b, fitted on 1 <= PR <= 10 and projected to PR 20
# (0.692 at 10, 0.613 at 14, 0.513 at 20; no data above 10). Source data, fit
# and limitations: paper/compressor_research/efficiency-pr-front-20261004/.
PARETO_FRONT_COEFFS = {
    "isentropic": (0.87096, 0.01991, 2.63312),
    "overall": (0.81982, 0.02940, 2.29041),
}
PARETO_FRONT_PR_RANGE = (1.0, 20.0)  # 1-10 fitted, 10-20 projected

# v2 design front (Oct 2026). One smooth formula,
#   eta(PR) = eta_1 - a * ln(PR)**b - d / (1 + exp(-(ln(PR) - L0)/w)),
# a gradual decline that starts flat at PR = 1 (b > 1) plus a logistic fall-off centred at
# PR = exp(L0) (~19), with a transition width w >= 0.2 in ln(PR) so the curve stays smooth.
# Fitted as a lower envelope (on or below at every PR from 1 to 20) of the design front: the
# highest efficiency of any eligible source (manufacturer software and tables, measured or
# peer-reviewed data, plus Piller single-stage MVR fans at 0.86-0.87; suction superheat
# <= 10 K; no claims), with dips between sources filled (the front never drops and recovers).
# Shaft and electrical data are pooled with a motor efficiency of 0.93
# (overall = 0.93 x isentropic). See
# paper/compressor_research/efficiency-pr-front-20261004/v2-manufacturer-software/.
PARETO_FRONT_V2_COEFFS = {  # (eta_1, a, b, d, L0, w)
    "isentropic": (0.8695, 0.025392, 1.476, 0.29915, 2.9599, 0.2),  # shaft-power basis
    "overall": (0.80859, 0.02366, 1.47491, 0.27829, 2.96035, 0.2),  # electrical input (= 0.93 x isentropic)
}
PARETO_FRONT_V2_PR_RANGE = (1.0, 20.0)


def pareto_front_efficiency(
    pressure_ratio: float, basis: str = "isentropic", version: int = 1
) -> float:
    """Design compressor efficiency at ``pressure_ratio``.

    ``basis`` is "isentropic" (shaft) or "overall" (electrical input). PR below
    1 is treated as 1. ``version=1``: original v1 front, eta_1 - a*ln(PR)**b (fit to PR 10,
    projected to 20). ``version=2``: v2 design front, one smooth formula (see above).
    """
    pr = max(float(pressure_ratio), 1.0)
    L = math.log(pr)
    if version == 2:
        eta_1, a, b, d, L0, w = PARETO_FRONT_V2_COEFFS[basis]
        return eta_1 - a * L**b - d / (1.0 + math.exp(-(L - L0) / w))
    eta_1, a, b = PARETO_FRONT_COEFFS[basis]
    return eta_1 - a * L**b


class Compressor:
    """Simple compressor with configurable efficiency correlations."""

    name: str
    flow_in: Stream = 0.0
    flow_out: Stream = 0.0
    flow_out_isen: Stream = 0.0
    isentropic_eff: float = 0.0
    pressure_ratio: float = 1.0
    overall_eff: float = 0.0
    shaft_work: float = 0.0
    power: float = 0.0
    delta_entropy: float = 0.0
    motor_eff: float = 1.0
    volume_eff: float = 1.0
    model_index: int = 1
    model_coeffs: list[float] = []
    converged: bool = False

    def __init__(
        self,
        name="Compressor",
        fluid: str = None,
        model_index: int = 0,
        motor_eff: float = 1.0,
    ) -> None:
        """Initialise the compressor with fluid, model, and motor efficiency."""
        self.name = name
        self.flow_in = Stream("flow_in", fluid=fluid)
        self.flow_out = Stream("flow_out", fluid=fluid)
        self.flow_out_isen = Stream("flow_out_isen", fluid=fluid)
        self.motor_eff = motor_eff
        self.model_index = model_index
        self.model_coeffs = self.set_model_coeff(model_index)

    def set_model_coeff(self, model_index: int) -> list[float]:
        """Return coefficient sets for the selected efficiency correlation."""
        available_model_coeffs = [
            ["Proportional", 1.0],  # Ideal compression
            ["Proportional", 0.7 / 0.96],  # Typical isentropic efficiency
            ["Polynomial"],  # Compressor model constants (Wang & Cleland, 2018)
            ["Brendel et al"],
            ["Pareto front"],  # Best published eta_is vs PR (2026 survey)
            ["Pareto front v2"],  # Design front, manufacturer software and plant data (2026)
        ]
        return available_model_coeffs[model_index]

    def get_isentropic_efficiency(self) -> float:
        """Calculate the isentropic efficiency based on the chosen model."""
        self.pressure_ratio = self.flow_out.pressure / self.flow_in.pressure
        if self.model_coeffs[0] == "Proportional":
            shaft_eff = self.model_coeffs[1]
        elif self.model_coeffs[0] == "Polynomial":
            shaft_eff = (
                0.65
                + 0.015 * self.pressure_ratio
                + -0.0015 * pow(self.pressure_ratio, 2.0)
            ) / 0.96
        elif self.model_coeffs[0] == "Brendel et al":
            if self.pressure_ratio > 1:
                shaft_eff = (
                    0.66981
                    - 0.6
                    / pow(self.pressure_ratio - 0.01466, 8.38 * self.flow_in.pressure)
                    - 0.00102 * pow(self.pressure_ratio, 1.8)
                )
            else:
                shaft_eff = 0.4
        elif self.model_coeffs[0] == "Pareto front":
            shaft_eff = pareto_front_efficiency(self.pressure_ratio, "isentropic")
        elif self.model_coeffs[0] == "Pareto front v2":
            shaft_eff = pareto_front_efficiency(
                self.pressure_ratio, "isentropic", version=2
            )
        return min(max(shaft_eff, 0.4), 1.0)

    def get_compressor_power(self) -> float:
        """Return the mechanical power requirement for the current operating point."""
        self.flow_out_isen.pp_instance = self.flow_out.pp_instance
        self.flow_out_isen.props_in_ab(
            "PS", self.flow_out.pressure, self.flow_in.entropy
        )
        dh_shaft = (
            self.flow_out_isen.enthalpy - self.flow_in.enthalpy
        ) / self.isentropic_eff
        self.shaft_work = dh_shaft * self.flow_in.mass_flow
        self.flow_out.enthalpy = (
            self.flow_in.enthalpy + self.shaft_work / self.flow_in.mass_flow
        )
        self.flow_out.props_in_mPH()
        return self.shaft_work / self.motor_eff

    def get_volumetric_efficiency(self) -> float:
        """Estimate volumetric efficiency using the default correlation."""
        # Equation from Brendel et al
        return 1 - 0.08244 * pow(self.pressure_ratio - 1, 0.72773)

    def solve(self) -> None:
        """Solve the compressor work and outlet state."""
        self.flow_out.mass_flow = self.flow_in.mass_flow
        self.isentropic_eff = self.get_isentropic_efficiency()
        self.overall_eff = self.isentropic_eff * self.motor_eff
        self.power = self.get_compressor_power()
        self.delta_entropy = (
            self.flow_out.entropy - self.flow_in.entropy
        ) * self.flow_in.mass_flow
        self.volume_eff = self.get_volumetric_efficiency()
        self.converged = True

    def get_graph_data(self, divisions: int = 50) -> dict:
        """Generate intermediate states for plotting the compression path."""
        data = {"siso": []}

        def get_inter_data(s_in: Stream, s_out: Stream) -> list[Stream]:
            """Interpolate states along the compression process."""
            ls = []
            # PV^n = const
            PR_i = (self.flow_out.pressure / self.flow_in.pressure) ** (1 / divisions)
            dh_tot = s_in.enthalpy - s_out.enthalpy
            for i in range(divisions + 1):
                m_str = Stream(
                    name="source i",
                    fluid=s_in.fluid,
                    pp_instance=s_in.pp_instance,
                    pressure=s_in.pressure * (PR_i) ** i,
                    enthalpy=s_in.enthalpy - dh_tot * i / divisions,
                    mass_flow=s_in.mass_flow,
                )
                m_str.props_in_mPH()
                ls.append(m_str)
            return ls

        data["siso"] = get_inter_data(self.flow_in, self.flow_out)
        return data
