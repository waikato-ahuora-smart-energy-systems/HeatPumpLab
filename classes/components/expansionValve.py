__author__ = 'Tim Walmsley'
__credits__ = ['Keri Adamson']

from classes.streams.massStream import MassStream

class ExpansionValve:
    
    def __init__(self, Name='Expansion valve', heat_loss=0.0, tol=1e-6):
        self.Name = Name
        self.FlowIn = MassStream()
        self.FlowOut = MassStream()
        
        self.heat_loss = heat_loss
        self.Pressure_drop = 0.0
        self.Power = 0.0
        self.tol = tol
    
    
    def solve(self):
        if self.FlowIn.Mass_Flow > self.tol:
            self.FlowOut.Enthalpy = self.FlowIn.Enthalpy - self.heat_loss / self.FlowIn.Mass_Flow
            self.FlowOut.PropsIn_mPH()
            self.Pressure_drop = self.FlowIn.Pressure - self.FlowOut.Pressure

        else:
            raise ValueError('No mass flow through the expansion valve: ' + self.Name)         