__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']


from streams import get_property, MassStream


class Separator:

    def __init__(self, Name='Separator', dP=0.0):
        self.Name = Name
        self.FlowIn = MassStream()
        self.VapOut = MassStream()
        self.LiquidOut = MassStream()
        self.dP = dP


    def solve(self):
        self.VapOut.Pressure = self.FlowIn.Pressure - self.dP
        self.VapOut.Quality = 1.0
        self.LiquidOut.Pressure = self.FlowIn.Pressure - self.dP
        self.LiquidOut.Quality = 0.0

        vap_frac = get_property('QMASS', self.FlowIn.Fluid, 'PH', self.VapOut.Pressure / 1000, self.FlowIn.Enthalpy, 0)
        
        self.VapOut.Mass_Flow = self.FlowIn.Mass_Flow * vap_frac
        self.VapOut.PropsIn_mPX()
        
        self.LiquidOut.Mass_Flow = self.FlowIn.Mass_Flow - self.VapOut.Mass_Flow
        self.LiquidOut.PropsIn_mPX()