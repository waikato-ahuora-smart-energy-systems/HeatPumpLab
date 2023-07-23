__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams import get_property, MassStream

class Compressor:
    def __init__(self, Name='Compressor', FlowIn = MassStream(), FlowOut = MassStream(), model_index=0, Motor_Eff=1.0):
        self.Name = Name
        self.FlowIn = FlowIn
        self.FlowOut = FlowOut
        self.Motor_Eff = Motor_Eff
        self.Isentropic_Eff = 1.0
        self.Pressure_Ratio = 1.0
        self.Overall_Eff = 1.0
        self.Shaft_Work = 0.0
        self.Power = 0.0
        self.model_index = model_index
        self.model_coeffs = self.set_model_coeff()

    def set_model_coeff(self):    
        available_model_coeffs = [[       1.0,   0.0,     0.0],
                                  [(0.7/0.96),   0.0,     0.0],         #Typical isentropic efficiency
                                  [      0.65, 0.015, -0.0015]]  #Compressor model constants (Wang & Cleland, 2018)       
        return available_model_coeffs[self.model_index]

    def get_isentropic_efficiency(self):
        self.Pressure_Ratio = self.FlowOut.Pressure / self.FlowIn.Pressure
        self.Overall_Eff = self.model_coeffs[0] + self.model_coeffs[1] * self.Pressure_Ratio + self.model_coeffs[2] * pow( self.Pressure_Ratio, 2 )
        self.Overall_Eff = min(max(self.Overall_Eff, 0.0), 1.0)
        return self.Overall_Eff / self.Motor_Eff
        
    def get_compressor_power(self):
        h_isen = get_property('H', self.FlowIn.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.FlowIn.Entropy, 0)
        dh_shaft = (h_isen - self.FlowIn.Enthalpy)  / self.Isentropic_Eff
        self.Shaft_Work = dh_shaft * self.FlowIn.Mass_Flow
        self.FlowOut.Enthalpy = self.FlowIn.Enthalpy + self.Shaft_Work / self.FlowIn.Mass_Flow
        self.FlowOut.PropsIn_mPH()
        return self.Shaft_Work / self.Motor_Eff
    
    def solve(self):
        self.Isentropic_Eff = self.get_isentropic_efficiency()
        self.Power = self.get_compressor_power()
