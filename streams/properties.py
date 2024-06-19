__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from sys import platform

# Identify the location of REFPROP library
class Properties:
    def __init__(self) -> None:
        pass

    def init_properties(self):
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
        return self.RP

# Retrieve more than one property value from REFPROP
def get_properties_values(propsOut, fluid, propsIn, prop1, prop2, RP):
    UNITS = RP.GETENUMdll(0,'SI with C').iEnum
    return RP.REFPROPdll(fluid, propsIn, propsOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output

# Retrieve one property value from REFPROP
def get_property_value(propOut, fluid, propsIn, prop1, prop2, RP): 
    UNITS = RP.GETENUMdll(0,'SI with C').iEnum
    if propOut == 'P' or propOut == 'PCRIT' or propOut == 'Pcrit':
        var = RP.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0] * 1000
    else:
        var = RP.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0]
    return var

