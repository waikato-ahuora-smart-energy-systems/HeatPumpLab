"""Economiser heat exchanger with pinch-based duty determination."""

from scipy.optimize import root_scalar
from .stream import Stream

class Economiser():
    """Economiser component acting between two compressor stages."""
    #Set the stream connection points
    name: str = ''
    model: str
    source_fluid: str
    source_in: Stream
    source_out: Stream
    sink_fluid: str
    sink_in: Stream
    sink_out: Stream
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

    comp1_flow_out: Stream
    comp2_flow_out: Stream
    gc_2_dP_h: float

    def __init__(self, source_fluid: str, sink_fluid: str, name: str='Economiser', 
                 y_econ: float =0.0, dTmin: float =5, dP_h: float =0.0, dP_c: float =0.0, 
                 converged: bool =False, tol: float =1e-6) -> None:
        """Initialise exchanger geometry and connect stream placeholders."""
        self.name = name
        self.converged = converged
        self.y_econ = y_econ
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.set_port_mass_streams(source_fluid, sink_fluid)

    def set_port_mass_streams(self, source_fluid: str, sink_fluid: str) -> None:
        """Create inlet and outlet streams for both sides of the economiser."""
        self.source_in = Stream(name='source in', fluid=source_fluid)
        self.source_out = Stream(name='source out', fluid=source_fluid)
        
        self.sink_in = Stream(name='sink in', fluid=sink_fluid)
        self.sink_out = Stream(name='sink out', fluid=sink_fluid)

    def set_add_obj(self, comp1_flow_out: Stream, comp2_flow_out: Stream, gc_2_dP_h: float =0.0, T_gc_out: float =0.0):
        """Attach upstream/downstream compressor states used during solution."""
        self.comp1_flow_out = comp1_flow_out
        self.comp2_flow_out = comp2_flow_out
        self.gc_2_dP_h = gc_2_dP_h
        self.T_gc_out = T_gc_out

    def get_duty(self): 
        """Compute the heat exchange duty respecting the pinch constraint."""
        if self.sink_in.pressure < self.sink_in.P_crit and self.sink_in.mass_flow > 0.0: 

            def _calc_pinch_dtmin(Q):
                """Return pinch delta-T minus target for a trial duty."""
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
        """Solve mass/energy balances for the current duty."""
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
        """Generate temperature/enthalpy profiles for plotting."""
        data = {'source': [], 'sink': []}
        
        def get_inter_data(s_in: Stream, s_out: Stream):
            """Interpolate intermediate stream states."""
            ls = []
            dP_tot = s_in.pressure - s_out.pressure
            dh_tot = s_in.enthalpy - s_out.enthalpy

            for i in range(divisions + 1):
                m_str = Stream(name='source i', 
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
