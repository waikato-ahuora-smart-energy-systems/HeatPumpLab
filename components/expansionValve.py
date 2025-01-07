__author__ = 'Tim Walmsley'
__credits__ = ['Keri Adamson']

from streams import MassStream

class ExpansionValve:
    flow_in: MassStream = None
    flow_out: MassStream = None
    converged: bool = True

    def __init__(self, name='Expansion valve', fluid: str = None, specific_heat_loss: float =0.0, tol: float =1e-6):
        self.name = name
        self.flow_in = MassStream(name='flow_in', fluid=fluid)
        self.flow_out = MassStream(name='flow_out', fluid=fluid)
        self.specific_heat_loss = specific_heat_loss
        self.tol = tol
        self.power = 0.0

    def get_outlet_enthalpy(self):
        return self.flow_in.enthalpy - self.specific_heat_loss
    
    def solve(self):
        self.flow_out.mass_flow = self.flow_in.mass_flow
        self.flow_out.enthalpy = self.get_outlet_enthalpy()
        self.flow_out.props_in_mPH()
        self.converged = True

    def get_graph_data(self, divisions=50):
        data = {'siso': []}
        def get_inter_data(s_in, s_out):
            ls = []
            # PV^n = const
            PR_i = (self.flow_out.pressure / self.flow_in.pressure) ** (1 / divisions)
            dh_tot = s_in.enthalpy - s_out.enthalpy
            for i in range(divisions + 1):
                m_str = MassStream(name='source i',
                                   pp_instance =s_in.pp_instance,
                                   fluid=s_in.fluid, 
                                   pressure=s_in.pressure * (PR_i) ** i, 
                                   enthalpy=s_in.enthalpy - dh_tot * i / divisions, 
                                   mass_flow=s_in.mass_flow, 
                                   )
                ls.append(m_str)
            
            return ls

        data['siso'] = get_inter_data(self.flow_in, self.flow_out)

        return data