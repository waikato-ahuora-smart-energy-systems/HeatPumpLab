__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams.propertiesPackage import PropertyPackage

class MassStream:
    name: str
    RP_fluid: PropertyPackage
    mass_flow: float
    pressure: float
    temperature: float
    enthalpy: float
    entropy: float
    quality: float | str
    specific_volume: float
    prandtl: float
    viscosity: float
    thermal_conductivity: float

    def __init__(self, 
                 name:str ='', 
                 fluid:PropertyPackage =PropertyPackage(), 
                 mass_flow:float =1.0, 
                 pressure:float =-1e-6, 
                 temperature:float =-1e-6, 
                 enthalpy:float =-1e-6, 
                 entropy:float =-1e-6, 
                 quality:float =-1e-6,
                 specific_volume:float = -1e-6):
        
        self.name = name
        self.RP_fluid = fluid
        self.mass_flow = mass_flow
        self.pressure = pressure
        self.temperature = temperature
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.quality = quality
        self.specific_volume = specific_volume

        if pressure != -1e-6 and enthalpy != -1e-6:
            self.props_in_mPH()
        elif pressure != -1e-6 and temperature != -1e-6:
            self.props_in_mPT()
        elif pressure != -1e-6 and quality != -1e-6:
            self.props_in_mPX()
        elif temperature != -1e-6 and quality != -1e-6:
            self.props_in_mTX()
        elif pressure != -1e-6 and entropy != -1e-6:
            self.props_in_mPS()
        else:
            None

    def get_H_from_PX(self):
        self.enthalpy = self.RP_fluid.get_properties('H', 'PX', self.pressure, self.quality)[0]

    def props_in_mPH(self):
        self.quality, self.temperature, self.entropy, self.heat_capacity, self.density, self.prandtl, self.viscosity, self.thermal_conductivity = self.RP_fluid.get_properties('QMASS;T;S;CP;D;PRANDTL;VIS;TCX', 'PH', self.pressure, self.enthalpy)[0:8]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPT(self):
        self.quality, self.enthalpy, self.entropy, self.heat_capacity, self.density, self.prandtl, self.viscosity, self.thermal_conductivity = self.RP_fluid.get_properties('QMASS;H;S;CP;D;PRANDTL;VIS;TCX', 'PT', self.pressure, self.temperature)[0:8]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPX(self):
        self.enthalpy, self.temperature, self.entropy, self.heat_capacity, self.density, self.prandtl, self.viscosity, self.thermal_conductivity = self.RP_fluid.get_properties('H;T;S;CP;D;PRANDTL;VIS;TCX', 'PQ', self.pressure, self.quality)[0:8]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPS(self):
        self.quality, self.temperature, self.enthalpy, self.heat_capacity, self.density, self.prandtl, self.viscosity, self.thermal_conductivity = self.RP_fluid.get_properties('QMASS;T;H;CP;D;PRANDTL;VIS;TCX', 'PS', self.pressure, self.entropy)[0:8]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mTX(self):
        self.pressure, self.enthalpy, self.entropy, self.heat_capacity, self.density, self.prandtl, self.viscosity, self.thermal_conductivity = self.RP_fluid.get_properties('P;H;S;CP;D;PRANDTL;VIS;TCX', 'TQ', self.temperature, self.quality)[0:8]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def copy_fld_props(self, fld):
        self.mass_flow = fld.mass_flow
        self.pressure = fld.pressure
        self.temperature = fld.temperature
        self.enthalpy = fld.enthalpy
        self.entropy = fld.entropy
        self.heat_capacity = fld.heat_capacity        
        self.quality = fld.quality
        self.density = fld.density
        self.specific_volume = fld.specific_volume
        self.prandtl = fld.prandtl
        self.viscosity = fld.viscosity
        self.thermal_conductivity = fld.thermal_conductivity
