from ..streams import MassStream

class Compressor:
    name: str
    flow_in: MassStream = 0.0
    flow_out: MassStream = 0.0
    flow_out_isen: MassStream = 0.0
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

    def __init__(self, name='Compressor', fluid: str =None, model_index: int =0, motor_eff: float =1.0) -> None:
        self.name = name
        self.flow_in = MassStream('flow_in', fluid=fluid) 
        self.flow_out = MassStream('flow_out', fluid=fluid)
        self.flow_out_isen = MassStream('flow_out_isen', fluid=fluid) 
        self.motor_eff = motor_eff
        self.model_index = model_index
        self.model_coeffs = self.set_model_coeff(model_index)

    def set_model_coeff(self, model_index: int) -> list[float]:    
        available_model_coeffs = [
            ['Proportional', 1.0],  # Ideal compression
            ['Proportional', 0.7],  # Typical isentropic efficiency
            ['Polynomial'],  # Compressor model constants (Wang & Cleland, 2018)  
            ['Brendel et al'],
        ]     
        return available_model_coeffs[model_index]

    def get_isentropic_efficiency(self) -> float:
        self.pressure_ratio = self.flow_out.pressure / self.flow_in.pressure
        if self.model_coeffs[0] == 'Proportional':
            shaft_eff = self.model_coeffs[1]
        elif self.model_coeffs[0] == 'Polynomial':
            shaft_eff = (0.65 + 0.015 * self.pressure_ratio + -0.0015 * pow( self.pressure_ratio, 2.0 )) 
        elif self.model_coeffs[0] == 'Brendel et al':
            if self.pressure_ratio > 1:
                shaft_eff = (0.66981 - 0.6 / pow( self.pressure_ratio - 0.01466, 8.38 * self.flow_in.pressure ) - 0.00102 * pow( self.pressure_ratio, 1.8 ))
            else:
                shaft_eff = 0.4
        return min(max(shaft_eff, 0.4), 1.0) / self.motor_eff
        
    def get_compressor_power(self) -> float:
        self.flow_out_isen.pp_instance = self.flow_out.pp_instance
        self.flow_out_isen.props_in_ab('PS', self.flow_out.pressure, self.flow_in.entropy)
        dh_shaft = (self.flow_out_isen.enthalpy - self.flow_in.enthalpy)  / self.isentropic_eff
        self.shaft_work = dh_shaft * self.flow_in.mass_flow
        self.flow_out.enthalpy = self.flow_in.enthalpy + self.shaft_work / self.flow_in.mass_flow
        self.flow_out.props_in_mPH()
        return self.shaft_work / self.motor_eff
    
    def get_volumetric_efficiency(self) -> float:
        # Equation from Brendel et al
        return 1 - 0.08244 * pow(self.pressure_ratio - 1, 0.72773)

    def solve(self) -> None:
        self.flow_out.mass_flow = self.flow_in.mass_flow
        self.isentropic_eff = self.get_isentropic_efficiency()
        self.overall_eff = self.isentropic_eff * self.motor_eff
        self.power = self.get_compressor_power()
        self.delta_entropy = (self.flow_out.entropy - self.flow_in.entropy) * self.flow_in.mass_flow
        self.volume_eff = self.get_volumetric_efficiency()
        self.converged = True

    def get_graph_data(self, divisions: int =50) -> dict:
        data = {'siso': []}
        def get_inter_data(s_in: MassStream, s_out: MassStream) -> list[MassStream]:
            ls = []
            # PV^n = const
            PR_i = (self.flow_out.pressure / self.flow_in.pressure) ** (1 / divisions)
            dh_tot = s_in.enthalpy - s_out.enthalpy
            for i in range(divisions + 1):
                m_str = MassStream(name='source i',
                                   fluid=s_in.fluid, 
                                   pp_instance=s_in.pp_instance,
                                   pressure=s_in.pressure * (PR_i) ** i, 
                                   enthalpy=s_in.enthalpy - dh_tot * i / divisions, 
                                   mass_flow=s_in.mass_flow)
                m_str.props_in_mPH()
                ls.append(m_str)
            return ls

        data['siso'] = get_inter_data(self.flow_in, self.flow_out)
        return data