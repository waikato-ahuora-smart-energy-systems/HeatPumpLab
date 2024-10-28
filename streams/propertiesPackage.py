__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from os import environ
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from sys import platform
import numpy as np

class PropertyPackage:
    RP: REFPROPFunctionLibrary
    if platform == "darwin":
        environ['RPPREFIX'] = r'/Library/REFPROP'
        RP = REFPROPFunctionLibrary(environ['RPPREFIX'])
        RP.SETPATHdll(environ['RPPREFIX'])
    else:
        try:
            RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/REFPROP.DLL')
        except:
            RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/')
        RP.SETPATHdll('C:/Program Files (x86)/REFPROP')
    
    iFlag: int = 0
    UNITS: int = RP.GETENUMdll(0,'SI with C').iEnum
    iMass: int = 1
    P_crit: float
    T_crit: float
    D_crit: float
    _fld: list[str]
    fld_name: str
    fld_composition: list[float]
    _composition: np.ndarray[float]
    
    def __init__(self, fld_name: str ='water'):
        self.set_fluid(fld_name)
        self.P_crit, self.T_crit, self.D_crit = self.get_properties('PCRIT;TCRIT;DCRIT', 'PH', 1.0, 500)[0:3]

    def set_iMass(self, basis: str ='mass'):
        self.iMass = 0 if basis == 'mol' else 1

    def set_fluid(self, fld: str):
        self._fld = fld.split('|')
        if len(self._fld) == 1:
            self.fld_name: str = self._fld[0]
            self.fld_composition: float = [1.0]
        elif len(self._fld) == 2 or len(self._fld) == 3:
            self.fld_name = '*'.join([f+'.FLD' for f in self._fld[0].split(';')])
            self._composition = np.array(self._fld[1].split(';'), dtype=float)
            self.fld_composition = self._composition / self._composition.sum()
            if len(self._fld) == 3:
                self.set_iMass(self._fld[2])
        else:
            ValueError('Invalid fluid string entry')

    def get_properties(self, propsOut:str, propsIn:str, prop1:float, prop2:float):
        return self.RP.REFPROPdll(self.fld_name, propsIn, propsOut, self.UNITS, self.iMass, self.iFlag, float(prop1), float(prop2), self.fld_composition).Output
    
    def set_critical_properties(self):
        self.P_crit, self.T_crit, self.D_crit = self.get_properties('PCRIT;TCRIT;DCRIT', 'PH', 1.0, 500)[0:3]
