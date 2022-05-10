__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from functions.refprop import REFPROPS

class Stream:
    def __init__(self, Name, Fluid = '', Mass_Flow = 1.0, Pressure = -1000, Temperature = -1000, Enthalpy = -1000, Quality = -1000, underdefined = 0):
        self.Name = Name
        self.Fluid = Fluid
        self.Mass_Flow = Mass_Flow
        self.Pressure = Pressure
        self.Temperature = Temperature
        self.Enthalpy = Enthalpy
        self.Entropy = 0
        self.Quality = Quality

        if not (self.Fluid == '') and underdefined == 0:
            if Temperature == -1000 and Quality == -1000:
                self.PropsIn_mPH()
            elif Enthalpy == -1000 and Quality == -1000:
                self.PropsIn_mPT()
            elif Temperature == -1000 and Enthalpy == -1000:
                self.PropsIn_mPX()
            elif Pressure == -1000 and Enthalpy == -1000:
                self.PropsIn_mTX()

    def PropsIn_mPH(self):
        props = REFPROPS('QMASS;T;S;CP', self.Fluid, 'PH', self.Pressure / 1000, self.Enthalpy)
        self.Quality = props[0];  self.Temperature = props[1];  self.Entropy = props[2];  self.Heat_Capacity = props[3] 
        self.Heat_Flow = self.Enthalpy * self.Mass_Flow
        
    def PropsIn_mPT(self):
        props = REFPROPS('QMASS;H;S;CP', self.Fluid, 'PT', self.Pressure / 1000, self.Temperature)
        self.Quality = props[0];  self.Enthalpy = props[1];  self.Entropy = props[2];  self.Heat_Capacity = props[3]
        self.Heat_Flow = self.Enthalpy * self.Mass_Flow

    def PropsIn_mPX(self):
        props = REFPROPS('T;H;S;CP', self.Fluid, 'PQMASS', self.Pressure / 1000, self.Quality)
        self.Temperature = props[0];  self.Enthalpy = props[1];  self.Entropy = props[2];  self.Heat_Capacity = props[3]
        self.Heat_Flow = self.Enthalpy * self.Mass_Flow

    def PropsIn_mPS(self):
        props = REFPROPS('T;H;QMASS;CP', self.Fluid, 'PS', self.Pressure / 1000, self.Entropy)
        self.Temperature = props[0];  self.Enthalpy = props[1];  self.Quality = props[2];  self.Heat_Capacity = props[3]
        self.Heat_Flow = self.Enthalpy * self.Mass_Flow

    def PropsIn_mTX(self):
        props = REFPROPS('P;H;S;CP', self.Fluid, 'TQMASS', self.Temperature, self.Quality)
        self.Pressure = props[0] * 1000;  self.Enthalpy = props[1];  self.Entropy = props[2];  self.Heat_Capacity = props[3]
        self.Heat_Flow = self.Enthalpy * self.Mass_Flow

    def Combine(self, fld1, fld2):
        #Assumes the fluids have the same composition
        self.Fluid = fld1.Fluid
        self.Mass_Flow = fld1.Mass_Flow + fld2.Mass_Flow
        self.Pressure = min(fld1.Pressure, fld2.Pressure)
        self.Enthalpy = (fld1.Heat_Flow + fld2.Heat_Flow) / self.Mass_Flow
        self.PropsIn_mPH()
    
    def PropsExtend(self): 
        props = REFPROPS('D;PRANDTL;VIS;TCX', self.Fluid, 'PT', self.Pressure / 1000, self.Temperature)
        self.Density = props[0];  self.PRANDTL = props[1];  self.Viscosity = props[2] /1000000000; self.Thermal_conductivity = props[3] / 1000000;  self.Specific_volume = 1 / props[0]