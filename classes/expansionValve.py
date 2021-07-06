__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

class Expansion_Valve:
    
    def __init__(self):
        self.FlowIn = None
        self.FlowOut = None
        #Isenthalpic valve
        self.dH_loss = 0
        self.Pressure_Ratio = 1
        self.Power = 0
    
    def solve(self):
        self.FlowOut.Enthalpy = self.FlowIn.Enthalpy - self.dH_loss
        self.FlowOut.PropsIn_mPH()

        #Turbine 
        #self.FlowOut.Entropy = self.FlowIn.Entropy 
        #self.FlowOut.PropsIn_mPS()
        #self.Power = self.FlowIn.Heat_Flow - self.FlowOut.Heat_Flow
        