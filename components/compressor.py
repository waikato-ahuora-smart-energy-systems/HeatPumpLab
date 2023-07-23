__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams import get_property, MassStream

class Compressor:
    def __init__(self, name='Compressor', flow_in = MassStream(), flow_out = MassStream(), model_index=0, motor_eff=1.0):
        self.name = name
        self.flow_in = flow_in
        self.flow_out = flow_out
        self.motor_eff = motor_eff
        self.isentropic_eff = 1.0
        self.pressure_ratio = 1.0
        self.overall_eff = 1.0
        self.shaft_work = 0.0
        self.power = 0.0
        self.model_index = model_index
        self.model_coeffs = self.set_model_coeff()

    def set_model_coeff(self):    
        available_model_coeffs = [[       1.0,   0.0,     0.0],
                                  [(0.7/0.96),   0.0,     0.0],         #Typical isentropic efficiency
                                  [      0.65, 0.015, -0.0015]]  #Compressor model constants (Wang & Cleland, 2018)       
        return available_model_coeffs[self.model_index]

    def get_isentropic_efficiency(self):
        self.pressure_ratio = self.flow_out.pressure / self.flow_in.pressure
        self.overall_eff = self.model_coeffs[0] + self.model_coeffs[1] * self.pressure_ratio + self.model_coeffs[2] * pow( self.pressure_ratio, 2 )
        self.overall_eff = min(max(self.overall_eff, 0.0), 1.0)
        return self.overall_eff / self.motor_eff
        
    def get_compressor_power(self):
        h_isen = get_property('H', self.flow_in.fluid, 'PS', self.flow_out.pressure / 1000, self.flow_in.entropy, 0)
        dh_shaft = (h_isen - self.flow_in.enthalpy)  / self.isentropic_eff
        self.shaft_work = dh_shaft * self.flow_in.mass_flow
        self.flow_out.enthalpy = self.flow_in.enthalpy + self.shaft_work / self.flow_in.mass_flow
        self.flow_out.props_in_mPH()
        return self.shaft_work / self.motor_eff
    
    def solve(self):
        self.isentropic_eff = self.get_isentropic_efficiency()
        self.power = self.get_compressor_power()
