__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

def get_property(propOut, fluid, propsIn, prop1, prop2, RP): 
    UNITS = RP.GETENUMdll(0,'SI with C').iEnum
    return RP.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0]

def get_properties(propsOut, fluid, propsIn, prop1, prop2, RP):
    UNITS = RP.GETENUMdll(0,'SI with C').iEnum
    return RP.REFPROPdll(fluid, propsIn, propsOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output

class Properties:
    def __init__(self) -> None:
        pass

    def init_properties(self):
        try:
            self.RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/') 
        except:
            self.RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/REFPROP.DLL')
        self.RP.SETPATHdll('C:/Program Files (x86)/REFPROP')