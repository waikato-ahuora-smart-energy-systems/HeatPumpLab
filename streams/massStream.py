__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams.propertiesPackage import PropertyPackage


class MassStream():
    name: str
    fluid: str
    pp_instance: str
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
                 pp_instance:str ='Default',
                 mass_flow:float =1.0, 
                 pressure:float =None, 
                 temperature:float =None, 
                 enthalpy:float =None, 
                 entropy:float =None, 
                 quality:float =None,
                 specific_volume:float = None) -> None:
        
        self.name = name
        self.fluid = fluid
        self.pp_instance = pp_instance
        self.mass_flow = mass_flow
        self.pressure = pressure
        self.temperature = temperature
        self.enthalpy = enthalpy
        self.entropy = entropy
        self.quality = quality
        self.specific_volume = specific_volume
        if fluid is not None:
            self.set_fluid(fluid, pp_instance)
            self.auto_calc_props()

    def set_fluid(self, fluid: str, pp_instance: str) -> None:
        self.PP = PropertyPackage(fluid, pp_instance)
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

    def props_in_ab(self, input_str: str, prop1: float, prop2: float, mass_flow: float =None, name: str =None, pp_instance: str =None) -> None:
        if mass_flow is not None:
            self.mass_flow = mass_flow
        if name is not None:
            self.name = name
        if pp_instance is not None:
            self.pp_instance = pp_instance
            self.PP._set_RP_lib(pp_instance)
        self._update_props(
            self.PP.get_std_properties(input_str, prop1, prop2, self.fld_composition)
        )

    def props_in_mPH(self) -> None:
        self._update_props(
            self.PP.get_std_properties('PH', self.pressure, self.enthalpy, self.fld_composition)
        )

    def props_in_mPQ(self) -> None:
        self._update_props(
            self.PP.get_std_properties('PQ', self.pressure, self.quality, self.fld_composition)
        )

    def props_in_mPT(self) -> None:
        self._update_props(
            self.PP.get_std_properties('PT', self.pressure, self.temperature, self.fld_composition)
        )         

    def props_in_mPS(self) -> None:
        self._update_props(
            self.PP.get_std_properties('PS', self.pressure, self.entropy, self.fld_composition)
        )           

    def props_in_mTQ(self) -> None:
        self._update_props(
            self.PP.get_std_properties('TQ', self.temperature, self.quality, self.fld_composition)
        )        

    def props_in_mTS(self) -> None:
        self._update_props(
            self.PP.get_std_properties('TQ', self.temperature, self.entropy, self.fld_composition)
        )      

    def _update_props(self, props) -> None:
        self.pressure = props['P']
        self.enthalpy = props['h']
        self.quality = props['q']
        self.temperature = props['T']
        self.entropy = props['s']
        self.heat_capacity = props['Cp']
        self.density = props['D']
        self.specific_volume = 1 / self.density if self.density > 0.0 else 0.0
        self.heat_flow = self.enthalpy * self.mass_flow

    def props_extend(self):
        self.prandtl, self.viscosity, self.thermal_conductivity = self.PP.get_specified_properties('PRANDTL;VIS;TCX', 'PH', self.pressure, self.enthalpy)[0:3]

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
