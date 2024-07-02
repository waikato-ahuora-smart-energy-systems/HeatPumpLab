__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from sys import platform
import numpy as np

# Identify the location of REFPROP library
class Properties:
    def __init__(self, fld='water'):
        if platform == "darwin":
            os.environ['RPPREFIX'] = r'/Library/REFPROP'
            self.RP = REFPROPFunctionLibrary(os.environ['RPPREFIX'])
            self.RP.SETPATHdll(os.environ['RPPREFIX'])
        else:
            try:
                self.RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/REFPROP.DLL')
                self.RP.SETPATHdll('C:/Program Files (x86)/REFPROP')
            except:
                self.RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/')
                self.RP.SETPATHdll('C:/Program Files (x86)/REFPROP')
        self.set_units()
        self.set_iFlag()
        self.set_iMass()
        self.set_fluid(fld)
        self.set_critical_properties()
    
    def set_units(self):
        self.UNITS = self.RP.GETENUMdll(0,'SI with C').iEnum

    def set_iMass(self, basis='mass'):
        self.iMass = 0 if basis == 'mol' else 1

    def set_fluid(self, fld):
        self._fld = fld.split('|')
        if len(self._fld) == 1:
            self.fld_name = self._fld[0]
            self.fld_composition = [1.0]
        elif len(self._fld) == 2 or len(self._fld) == 3:
            self.fld_name = '*'.join([f+'.FLD' for f in self._fld[0].split(';')])
            self._composition = np.array(self._fld[1].split(';'), dtype=float)
            self.fld_composition = self._composition / self._composition.sum()
            if len(self._fld) == 3:
                self.set_iMass(self._fld[2])
        else:
            ValueError('Invalid fluid string entry')

    def set_iFlag(self):
        self.iFlag = 0

    def get_properties(self, propsOut, propsIn, prop1, prop2):
        return self.RP.REFPROPdll(self.fld_name, propsIn, propsOut, self.UNITS, self.iMass, self.iFlag, float(prop1), float(prop2), self.fld_composition).Output
    
    def set_critical_properties(self):
        v = self.get_properties('PCRIT;TCRIT;DCRIT', 'PH', 1.0, 500)
        self.pcrit, self.tcrit, self.dcrit = v[0], v[1], v[2]
