__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams.properties import get_properties

class MassStream:
    def __init__(self, name='', fluid='', mass_flow=1.0, pressure=-1000, temperature=-1000, enthalpy=-1000, quality=-1000, underdefined=0, RP=None):
        self.name = name
        self.fluid = fluid
        self.mass_flow = mass_flow
        self.pressure = pressure
        self.temperature = temperature
        self.enthalpy = enthalpy
        self.entropy = 0
        self.quality = quality
        self.RP = RP

        if not (self.fluid == '') and underdefined == 0:
            if temperature == -1000 and quality == -1000:
                self.props_in_mPH()
            elif enthalpy == -1000 and quality == -1000:
                self.props_in_mPT()
            elif temperature == -1000 and enthalpy == -1000:
                self.props_in_mPX()
            elif pressure == -1000 and enthalpy == -1000:
                self.props_in_mTX()

    def props_in_mPH(self):
        props = get_properties('QMASS;T;S;CP', self.fluid, 'PH', self.pressure / 1000, self.enthalpy, RP=self.RP)
        self.quality = props[0];  self.temperature = props[1];  self.entropy = props[2];  self.heat_capacity = props[3] 
        self.heat_flow = self.enthalpy * self.mass_flow
        del props

    def props_in_mPT(self):
        props = get_properties('QMASS;H;S;CP', self.fluid, 'PT', self.pressure / 1000, self.temperature, RP=self.RP)
        self.quality = props[0];  self.enthalpy = props[1];  self.entropy = props[2];  self.heat_capacity = props[3]
        self.heat_flow = self.enthalpy * self.mass_flow
        del props

    def props_in_mPX(self):
        props = get_properties('T;H;S;CP', self.fluid, 'PQMASS', self.pressure / 1000, self.quality, RP=self.RP)
        self.temperature = props[0];  self.enthalpy = props[1];  self.entropy = props[2];  self.heat_capacity = props[3]
        self.heat_flow = self.enthalpy * self.mass_flow
        del props

    def props_in_mPS(self):
        props = get_properties('T;H;QMASS;CP', self.fluid, 'PS', self.pressure / 1000, self.entropy, RP=self.RP)
        self.temperature = props[0];  self.enthalpy = props[1];  self.quality = props[2];  self.heat_capacity = props[3]
        self.heat_flow = self.enthalpy * self.mass_flow
        del props

    def props_in_mTX(self):
        props = get_properties('P;H;S;CP', self.fluid, 'TQMASS', self.temperature, self.quality, RP=self.RP)
        self.pressure = props[0] * 1000;  self.enthalpy = props[1];  self.entropy = props[2];  self.heat_capacity = props[3]
        self.heat_flow = self.enthalpy * self.mass_flow
        del props
        
    def props_extend(self): 
        props = get_properties('D;PRANDTL;VIS;TCX', self.fluid, 'PT', self.pressure / 1000, self.temperature, RP=self.RP)
        self.density = props[0];  self.prandtl = props[1];  self.viscosity = props[2] /1000000000; self.thermal_conductivity = props[3] / 1000000;  self.specific_volume = 1 / props[0]
        del props
        
    def get_combined_stream(self, fld1, fld2):
        #Assumes the fluids have the same composition
        self.fluid = fld1.fluid
        self.mass_flow = fld1.mass_flow + fld2.mass_flow
        self.pressure = min(fld1.pressure, fld2.pressure)
        self.enthalpy = (fld1.heat_flow + fld2.heat_flow) / self.mass_flow
        self.props_in_mPH()

    def set_equal(self, fld):
        self.mass_flow = fld.mass_flow
        self.pressure = fld.pressure
        self.temperature = fld.temperature
        self.enthalpy = fld.enthalpy
        self.entropy = fld.entropy
        self.heat_capacity = fld.heat_capacity        
        self.quality = fld.quality
