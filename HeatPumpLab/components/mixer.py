"""Two-stream mixing junction."""

from ..streams import MassStream


class Mixer:
    """Adiabatic mixer with optional specific heat loss."""
    flow_in_1: MassStream = None
    flow_in_2: MassStream = None
    flow_out: MassStream = None
    converged: bool = True

    def __init__(self, name='Mixer', fluid: str = None, specific_heat_loss: float =0.0, tol: float =1e-6):
        """Initialise inlet/outlet streams and optional heat loss."""
        self.name = name
        self.flow_in_1 = MassStream(name='flow_in_1', fluid=fluid)
        self.flow_in_2 = MassStream(name='flow_in_2', fluid=fluid)
        self.flow_out = MassStream(name='flow_out', fluid=fluid)
        self.specific_heat_loss = specific_heat_loss
        self.tol = tol

    def get_outlet_enthalpy(self):
        """Return the mixed enthalpy weighted by inlet mass fractions."""
        if self.flow_out.mass_flow > self.tol:
            y_1 = self.flow_in_1.mass_flow / self.flow_out.mass_flow
            y_2 = 1 - y_1
            return self.flow_in_1.enthalpy * y_1 + self.flow_in_2.enthalpy * y_2 - self.specific_heat_loss
        else:
            print('No mass flow through the mixer: ' + self.name)
            return self.flow_in_1.enthalpy
    
    def solve(self):
        """Balance mass and energy to determine the outlet state."""
        self.flow_out.mass_flow = self.flow_in_1.mass_flow + self.flow_in_2.mass_flow
        self.flow_out.pressure = min(self.flow_in_1.pressure, self.flow_in_2.pressure)
        self.flow_out.enthalpy = self.get_outlet_enthalpy()
        self.flow_out.props_in_mPH()
        self.converged = True

    # def get_graph_data(self, divisions=50):
    #     data = {'siso': []}
    #     def get_inter_data(s_in, s_out):
    #         ls = []
    #         # PV^n = const
    #         PR_i = (self.flow_out.pressure / self.flow_in.pressure) ** (1 / divisions)
    #         dh_tot = s_in.enthalpy - s_out.enthalpy
    #         for i in range(divisions + 1):
    #             m_str = MassStream(name='source i',
    #                                pp_instance =s_in.pp_instance,
    #                                fluid=s_in.fluid, 
    #                                pressure=s_in.pressure * (PR_i) ** i, 
    #                                enthalpy=s_in.enthalpy - dh_tot * i / divisions, 
    #                                mass_flow=s_in.mass_flow, 
    #                                )
    #             # m_str.props_in_mPH()
    #             ls.append(m_str)
            
    #         return ls

    #     data['siso'] = get_inter_data(self.flow_in, self.flow_out)

    #     return data
