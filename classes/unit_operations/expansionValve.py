__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

class Expansion_Valve:
    
    def __init__(self):
        self.FlowIn = None
        self.FlowOut = None
        #Isenthalpic valve
        self.dH_loss = 0
        self.Pressure_drop = 1
        self.Power = 0
    
    def solve(self):
        self.FlowOut.Enthalpy = self.FlowIn.Enthalpy - self.dH_loss
        self.FlowOut.PropsIn_mPH()
        self.Pressure_drop = self.FlowIn.Pressure - self.FlowOut.Pressure         