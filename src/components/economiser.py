from scipy.optimize import root_scalar
from ..streams import MassStream

class Economiser():
    #Set the stream connection points
    name: str = ''
    model: str
    source_fluid: str
    source_in: MassStream
    source_out: MassStream
    sink_fluid: str
    sink_in: MassStream
    sink_out: MassStream
    bypassed: bool = False
    converged: bool = False
    y_econ: float = 0.0
    dTmin: float = 0,0
    HX_Eff: float = 0.0
    NTU: float = 0.0
    c: float = 0.0
    U: float = 0.0
    Area: float = 0.0
    duty: float = 0.0
    dP_h: float = 0.0
    dP_c: float = 0.0
    variable_dP_model: bool = False
    tol: float = 1e-6    

    comp1_flow_out: MassStream
    comp2_flow_out: MassStream
    gc_2_dP_h: float

    def __init__(self, source_fluid: str, sink_fluid: str, name: str='Economiser', 
                 y_econ: float =0.0, dTmin: float =5, dP_h: float =0.0, dP_c: float =0.0, 
                 converged: bool =False, tol: float =1e-6) -> None:
        self.name = name
        self.converged = converged
        self.y_econ = y_econ
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.set_port_mass_streams(source_fluid, sink_fluid)

    def set_port_mass_streams(self, source_fluid: str, sink_fluid: str) -> None:
        self.source_in = MassStream(name='source in', fluid=source_fluid)
        self.source_out = MassStream(name='source out', fluid=source_fluid)
        
        self.sink_in = MassStream(name='sink in', fluid=sink_fluid)
        self.sink_out = MassStream(name='sink out', fluid=sink_fluid)

    def set_add_obj(self, comp1_flow_out: MassStream, comp2_flow_out: MassStream, gc_2_dP_h: float =0.0, T_gc_out: float =0.0):
        self.comp1_flow_out = comp1_flow_out
        self.comp2_flow_out = comp2_flow_out
        self.gc_2_dP_h = gc_2_dP_h
        self.T_gc_out = T_gc_out

    def get_duty(self): 
        if self.sink_in.pressure < self.sink_in.P_crit and self.sink_in.mass_flow > 0.0: 

            def _calc_pinch_dtmin(Q):
                self.source_in.enthalpy = self.source_out.enthalpy + Q / self.source_out.mass_flow
                self.source_in.props_in_mPH()

                self.sink_out.enthalpy = self.sink_in.enthalpy + Q / self.sink_in.mass_flow
                self.sink_out.props_in_mPH()
                return self.source_in.temperature - self.sink_out.temperature - self.dTmin
        
            Q_max = self.source_out.mass_flow * (self.comp2_flow_out.enthalpy - self.source_out.enthalpy)
            r = root_scalar(f=_calc_pinch_dtmin, method='brentq', bracket=[0.01, Q_max])
            # print(r.root)
            return max(r.root, 0.0)
        else:
            return 0.0

    def solve(self):
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.source_in.pressure = self.source_out.pressure + self.dP_h
        self.source_in.mass_flow = self.source_out.mass_flow
        self.sink_out.mass_flow = self.sink_in.mass_flow

        self.duty = self.get_duty()
        
        self.source_in.enthalpy = self.source_out.enthalpy + self.duty / self.source_out.mass_flow if self.source_out.mass_flow > 0.0 else self.source_out.enthalpy
        self.source_in.props_in_mPH()

        self.sink_out.enthalpy = self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow if self.sink_in.mass_flow > 0.0 else self.sink_in.enthalpy
        self.sink_out.props_in_mPH()
        self.converged = True

    def get_graph_data(self, divisions=50):
        data = {'source': [], 'sink': []}
        
        def get_inter_data(s_in: MassStream, s_out: MassStream):
            ls = []
            dP_tot = s_in.pressure - s_out.pressure
            dh_tot = s_in.enthalpy - s_out.enthalpy

            for i in range(divisions + 1):
                m_str = MassStream(name='source i', 
                                   fluid=s_in.fluid,
                                   pp_instance=s_in.pp_instance,
                                   pressure=s_in.pressure - dP_tot * i / divisions, 
                                   enthalpy=s_in.enthalpy - dh_tot * i / divisions, 
                                   mass_flow=s_in.mass_flow, 
                                   )
                m_str.props_in_mPH()
                ls.append(m_str)
            return ls

        data['sink'] = get_inter_data(self.sink_in, self.sink_out)
        data['source'] = get_inter_data(self.source_in, self.source_out)
        return data