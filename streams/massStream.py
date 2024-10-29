__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams.propertiesPackage import PropertyPackage


class MassStream:
    name: str
    fluid: str
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
    TC_to_K: float = 273.15
    DEFAULT_VAL: float = -9999990.0

    def __init__(self, 
                 name:str ='', 
                 fluid:str =None,
                 mass_flow:float =1.0, 
                 pressure:float =None, 
                 temperature:float =None, 
                 enthalpy:float =None, 
                 entropy:float =None, 
                 quality:float =None,
                 specific_volume:float = None) -> None:
        
        self.name = name
        self.fluid = fluid
        self.mass_flow = mass_flow
        self.pressure = pressure
        self.temperature = temperature
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.quality = quality
        self.specific_volume = specific_volume
        if fluid is not None:
            self.set_fluid(fluid)
            self.auto_calc_props()

    def set_fluid(self, fluid: str) -> None:
        self.fluid = fluid
        self.PP = PropertyPackage(fluid)
        self.P_crit = self.PP.P_crit
        self.T_crit = self.PP.T_crit
        self.D_crit = self.PP.D_crit
        self.molar_mass = self.PP.molar_mass
        self.fld_composition = self.PP.fld_composition

    def auto_calc_props(self) -> None:
        if self.pressure != None and self.enthalpy != None:
            self.props_in_mPH()
        elif self.pressure != None and self.temperature != None:
            self.props_in_mPT()
        elif self.pressure != None and self.quality != None:
            self.props_in_mPQ()
        elif self.temperature != None and self.quality != None:
            self.props_in_mTQ()
        elif self.pressure != None and self.entropy != None:
            self.props_in_mPS()
        elif self.temperature != None and self.entropy != None:
            self.props_in_mTS()    
        else:
            None        

    def props_in_ab(self, input_str: str, prop1: float, prop2: float, mass_flow: float =None, name: str =None) -> None:
        if mass_flow != None:
            self.mass_flow = mass_flow

        if name != None:
            self.name = name

        if input_str == 'PH':
            self.pressure = prop1
            self.enthalpy = prop2
            self.props_in_mPH()
        elif input_str == 'HP':
            self.pressure = prop2
            self.enthalpy = prop1
            self.props_in_mPH()

        elif input_str == 'PT':
            self.pressure = prop1
            self.temperature = prop2
            self.props_in_mPT()
        elif input_str == 'TP':
            self.pressure = prop2
            self.temperature = prop1
            self.props_in_mPT()

        elif input_str == 'PQ':
            self.pressure = prop1
            self.quality = prop2
            self.props_in_mPQ()
        elif input_str == 'QP':
            self.pressure = prop2
            self.quality = prop1
            self.props_in_mPQ()   

        elif input_str == 'TQ':
            self.temperature = prop1
            self.quality = prop2
            self.props_in_mTQ()
        elif input_str == 'QT':
            self.temperature = prop2
            self.quality = prop1
            self.props_in_mTQ()         
        
        elif input_str == 'PS':
            self.pressure = prop1
            self.entropy = prop2
            self.props_in_mPS()
        elif input_str == 'SP':
            self.pressure = prop2
            self.entropy = prop1
            self.props_in_mPS()              

        elif input_str == 'TS':
            self.temperature = prop1
            self.entropy = prop2
            self.props_in_mTS()
        elif input_str == 'ST':
            self.temperature = prop2
            self.entropy = prop1
            self.props_in_mTS()   

        else:
            raise Exception("Unrecognised fluid property input combination.")         

    def props_in_mPH(self) -> None:
        # self.props_in_mPH_old()
        self._update_props(
            self.PP.get_properties_fast('PH', self.pressure * 1000, self.enthalpy * self.molar_mass, self.fld_composition)
        )

    def props_in_mPQ(self) -> None:
        # self.props_in_mPQ_old()
        self._update_props(
            self.PP.get_properties_fast('PQ', self.pressure * 1000, self.quality, self.fld_composition)
        )   

    def props_in_mPT(self) -> None:
        self.props_in_mPT_old()
        # self._update_props(
        #     self.PP.get_properties_fast('PT', self.pressure * 1000, self.temperature + self.TC_to_K, self.fld_composition)
        # ) 

        # if (abs(1 - self.enthalpy / self.enthalpy1) > 1e-4 or
        #     abs(1 - self.quality / self.quality1) > 1e-4 or
        #     abs(1 - self.entropy / self.entropy1) > 1e-4 or
        #     abs(1 - self.heat_capacity / self.heat_capacity1) > 1e-4 or
        #     abs(1 - self.density / self.density1) > 1e-4 or
        #     abs(1 - self.specific_volume / self.specific_volume) > 1e-4
        #     ):
        #     print('Check') # 550.7620311829611

    def props_in_mPS(self) -> None:
        # self.props_in_mPS_old()
        self._update_props(
            self.PP.get_properties_fast('PS', self.pressure * 1000, self.entropy * self.molar_mass, self.fld_composition)
        )           

    def props_in_mTQ(self) -> None:
        # self.props_in_mTQ_old()
        self._update_props(
            self.PP.get_properties_fast('TQ', self.temperature + self.TC_to_K, self.quality, self.fld_composition)
        )        

    def props_in_mTS(self) -> None:
        self.props_in_mTS_old()
        # self._update_props(
        #     self.PP.get_properties_fast('TQ', self.temperature + self.TC_to_K, self.entropy * self.molar_mass, self.fld_composition)
        # ) 

    def _update_props(self, props) -> None:
        self.pressure = props.P / 1000
        self.enthalpy = props.h / self.molar_mass
        self.quality = props.q if (props.q >= 0.0 and props.q <= 1.0) else self.DEFAULT_VAL
        self.temperature = props.T - self.TC_to_K
        self.entropy = props.s / self.molar_mass
        self.heat_capacity = props.Cp / self.molar_mass if props.Cp > 0.0 else self.DEFAULT_VAL
        self.density = props.D * self.molar_mass
        self.specific_volume = 1 / self.density if self.density > 0.0 else 0.0
        self.heat_flow = self.enthalpy * self.mass_flow  



    def props_in_mPH_old(self):
        self.quality, self.temperature, self.entropy, self.heat_capacity, self.density = self.PP.get_properties('QMASS;T;S;CP;D', 'PH', self.pressure, self.enthalpy)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPT_old(self):
        # test = self.PP.get_properties('QMASS;H;S;CP;D', 'PT', self.pressure, self.temperature)
        self.quality, self.enthalpy, self.entropy, self.heat_capacity, self.density = self.PP.get_properties('QMASS;H;S;CP;D', 'PT', self.pressure, self.temperature)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPQ_old(self):
        self.enthalpy, self.temperature, self.entropy, self.heat_capacity, self.density = self.PP.get_properties('H;T;S;CP;D', 'PQ', self.pressure, self.quality)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mPS_old(self):
        self.quality, self.temperature, self.enthalpy, self.heat_capacity, self.density = self.PP.get_properties('QMASS;T;H;CP;D', 'PS', self.pressure, self.entropy)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mTQ_old(self):
        self.pressure, self.enthalpy, self.entropy, self.heat_capacity, self.density = self.PP.get_properties('P;H;S;CP;D', 'TQ', self.temperature, self.quality)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_in_mTS_old(self):
        self.quality, self.pressure, self.enthalpy, self.heat_capacity, self.density = self.PP.get_properties('QMASS;P;H;CP;D', 'TS', self.temperature, self.entropy)[0:5]
        self.specific_volume = 1 / self.density
        self.heat_flow = self.enthalpy * self.mass_flow



    def props_extend(self):
        self.prandtl, self.viscosity, self.thermal_conductivity = self.PP.get_properties('PRANDTL;VIS;TCX', 'PH', self.pressure, self.enthalpy)[0:3]

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
