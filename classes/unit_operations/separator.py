__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

class Separator:
    def __init__(self):
        self.FlowIn = None
        self.VapOut = None
        self.LiquidOut = None

    def solve(self):
        self.VapOut.Pressure = self.FlowIn.Pressure
        self.VapOut.Mass_Flow = self.FlowIn.Mass_Flow - self.LiquidOut.Mass_Flow
        self.VapOut.PropsIn_mPX()
        self.LiquidOut.Pressure = self.FlowIn.Pressure
        self.LiquidOut.Mass_Flow = 1
        self.LiquidOut.PropsIn_mPX()