"""Binary splitter component."""

from ..streams import MassStream


class Splitter:
    """Splits a stream into two branches based on a fixed ratio."""
    flow_in: MassStream = None
    flow_out_1: MassStream = None
    flow_out_2: MassStream = None
    r_ratio: float =0.0
    converged: bool = True

    def __init__(self, name='Splitter', fluid: str = None, r_ratio: float =0.0, specific_heat_loss: float =0.0, tol: float =1e-6):
        """Initialise inlet/outlet streams and the flow split ratio."""
        self.name = name
        self.flow_in = MassStream(name='flow_in', fluid=fluid)
        self.flow_out_1 = MassStream(name='flow_out_1', fluid=fluid)
        self.flow_out_2 = MassStream(name='flow_out_2', fluid=fluid)
        self.r_ratio = r_ratio
        self.specific_heat_loss = specific_heat_loss
        self.tol = tol
    
    def solve(self):
        """Split the inlet flow and copy thermodynamic properties."""
        self.flow_out_1.copy_fld_props(self.flow_in)
        self.flow_out_2.copy_fld_props(self.flow_in)
        m_1 = self.flow_in.mass_flow * 1 / (1 + self.r_ratio)
        m_2 = self.flow_in.mass_flow - m_1
        self.flow_out_1.set_flows(mass_flow=m_1)
        self.flow_out_2.set_flows(mass_flow=m_2)
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
