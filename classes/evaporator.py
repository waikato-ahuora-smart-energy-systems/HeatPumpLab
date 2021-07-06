__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from globalFunctions import REFPROP

class Evaporator:

    def __init__(self):
        self.HotIn = None
        self.HotOut = None
        self.ColdIn = None
        self.ColdOut = None
        self.HX_Eff = 0.8
        self.NTU = 0
        self.c = 0
        self.U = 1
        self.Duty = 0
        self.Area = 0
        self.dP_h = 0
        self.dP_c = 0
        self.Beta = 587 # m^2/m^3
        self.FLength = 0.5 # m
        self.Depth = 1 # m
    
    def solve(self):       
        self.Duty = self.ColdOut.Heat_Flow - self.ColdIn.Heat_Flow
        self.ColdIn.Pressure = self.ColdOut.Pressure + self.dP_c
        self.ColdIn.PropsIn_mPH()