__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes.streams.properties import REFPROP
from classes.streams.massStream import MassStream

class Compressor:

    def __init__(self, Name='Compressor', Motor_Eff=1.0, model=0):
        self.Name = Name
        self.FlowIn = MassStream()
        self.FlowOut = MassStream()

        self.Motor_Eff = Motor_Eff
        self.Pressure_Ratio = 1.0
        self.dP = 0
        self.Isen_Eff = 0
        self.dh_shaft = 0
        self.Shaft_Work = 0
        self.Power = 0
        self.model = model
        self.set_model_coeff()


    def set_model_coeff(self):    
        compressor_model = [[0.7, 0.0, 0.0],         #Typical isentropic efficiency
                            [0.65, 0.015, -0.0015]]  #Compressor model constants (Wang & Cleland, 2018)       

        self.COMP_COEFF_A = compressor_model[self.model][0]
        self.COMP_COEFF_B = compressor_model[self.model][1]
        self.COMP_COEFF_C = compressor_model[self.model][2]


    def solve(self):
        self.Pressure_Ratio = self.FlowOut.Pressure / self.FlowIn.Pressure
        self.dP = self.FlowOut.Pressure - self.FlowIn.Pressure
        
        #Determines efficiency based on compressor model, shaft work and power use (including motor efficiency)
        self.Overall_Eff = min(max(self.COMP_COEFF_A + self.COMP_COEFF_B * self.Pressure_Ratio + self.COMP_COEFF_C * pow( self.Pressure_Ratio, 2 ), 0.25), 1.0)
        self.Isen_Eff = self.Overall_Eff / self.Motor_Eff
        
        self.dh_comp = ( REFPROP('H', self.FlowOut.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.FlowIn.Entropy, 0) - self.FlowIn.Enthalpy )  / self.Overall_Eff
        self.dh_shaft = self.dh_comp * self.Motor_Eff
        self.FlowOut.Enthalpy = self.FlowIn.Enthalpy + self.dh_shaft
        self.FlowOut.PropsIn_mPH()
        self.Shaft_Work = self.dh_shaft * self.FlowIn.Mass_Flow
        self.Power = self.dh_comp * self.FlowIn.Mass_Flow