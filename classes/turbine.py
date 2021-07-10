__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

class Turbine:
    
    def __init__(self, eff_isen = 1.0, eff_gen = 1.0):
        self.FlowIn = None
        self.FlowOut = None
        self.dh = 0.0
        self.Pressure_Ratio = 1.0
        self.eff_isen = eff_isen
        self.eff_gen = eff_gen
        self.Power = 0.0
    
    def solve(self):
        self.FlowOut.Entropy = self.FlowIn.Entropy 
        self.FlowOut.PropsIn_mPS()
        self.dh_isen = self.FlowIn.Enthalpy - self.FlowOut.Enthalpy
        self.dh = self.dh_isen / self.eff_isen
        self.FlowOut.Enthalpy = self.FlowIn.Enthalpy - self.dh 
        self.FlowOut.PropsIn_mPH()
        self.Power = self.dh * self.FlowOut.Mass_Flow * self.eff_gen
        