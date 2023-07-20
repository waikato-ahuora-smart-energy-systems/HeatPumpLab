__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

try:
    RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/REFPROP.DLL') 
except:
    RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/')
RP.SETPATHdll('C:/Program Files (x86)/REFPROP')
UNITS = RP.GETENUMdll(0,'SI with C').iEnum

def REFPROP(propOut, fluid, propsIn, prop1, prop2, arg = [1.0]): 
    return RP.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0]

def REFPROPS(propsOut, fluid, propsIn, prop1, prop2, arg = [1.0]):
    return RP.REFPROPdll(fluid, propsIn, propsOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output

def HAPROPS(propOut, propsIn, P, T_or_H, AH, arg = [1.0]):
    pass