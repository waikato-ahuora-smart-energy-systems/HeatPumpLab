__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams.properties import Properties

class MassStream:
    def __init__(self, name='', fluid=Properties(), mass_flow=1.0, pressure=-1e-6, temperature=-1e-6, enthalpy=-1e-6, entropy=-1e-6, quality=-1e-6):
        self.name = name
        self.RP_fluid = fluid if fluid is not str else Properties(fluid)
        self.mass_flow = mass_flow
        self.pressure = pressure
        self.temperature = temperature
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.quality = quality
        self.specific_volume = 0

        if temperature == -1e-6 and quality == -1e-6 and entropy == -1e-6:
            self.props_in_mPH()
        elif enthalpy == -1e-6 and quality == -1e-6 and entropy == -1e-6:
            self.props_in_mPT()
        elif temperature == -1e-6 and enthalpy == -1e-6 and entropy == -1e-6:
            self.props_in_mPX()
        elif pressure == -1e-6 and enthalpy == -1e-6 and entropy == -1e-6:
            self.props_in_mTX()
        elif quality == -1e-6 and enthalpy == -1e-6 and temperature == -1e-6:
            self.props_in_mPS()

    def props_in_mPH(self):
        self.quality, self.temperature, self.entropy, self.heat_capacity, self.density = self.RP_fluid.get_properties('QMASS;T;S;CP;D', 'PH', self.pressure, self.enthalpy)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPT(self):
        self.quality, self.enthalpy, self.entropy, self.heat_capacity, self.density = self.RP_fluid.get_properties('QMASS;H;S;CP;D', 'PT', self.pressure, self.temperature)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPX(self):
        self.enthalpy, self.temperature, self.entropy, self.heat_capacity, self.density = self.RP_fluid.get_properties('H;T;S;CP;D', 'PQ', self.pressure, self.quality)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPS(self):
        self.quality, self.temperature, self.enthalpy, self.heat_capacity, self.density = self.RP_fluid.get_properties('QMASS;T;H;CP;D', 'PS', self.pressure, self.entropy)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mTX(self):
        self.pressure, self.enthalpy, self.entropy, self.heat_capacity, self.density = self.RP_fluid.get_properties('P;H;S;CP;D', 'TQ', self.temperature, self.quality)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow
        
    def props_extend(self): 
        self.prandtl, self.viscosity, self.thermal_conductivity = self.RP_fluid.get_properties('PRANDTL;VIS;TCX', 'PH', self.pressure, self.enthalpy)[0:3]
        
    # def get_combined_stream(self, fld1, fld2):
    #     #Assumes the fluids have the same composition
    #     self.fluid = fld1.fluid
    #     self.mass_flow = fld1.mass_flow + fld2.mass_flow
    #     self.pressure = min(fld1.pressure, fld2.pressure)
    #     self.enthalpy = (fld1.heat_flow + fld2.heat_flow) / self.mass_flow
    #     self.props_in_mPH()

    def set_equal(self, fld):
        self.mass_flow = fld.mass_flow
        self.pressure = fld.pressure
        self.temperature = fld.temperature
        self.enthalpy = fld.enthalpy
        self.entropy = fld.entropy
        self.heat_capacity = fld.heat_capacity        
        self.quality = fld.quality
        self.density = fld.density
        self.specific_volume = fld.specific_volume
