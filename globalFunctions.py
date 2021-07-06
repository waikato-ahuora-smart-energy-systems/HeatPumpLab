__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from globalVariables import CONVERGANCE
import math
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/') #/REFPROP.DLL
RP.SETPATHdll('C:/Program Files (x86)/REFPROP')
UNITS = RP.GETENUMdll(0,'SI with C').iEnum

def REFPROP(propOut, fluid, propsIn, prop1, prop2, arg = [1.0]): #Arguments order matches Excel call function
    return RP.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0]

def REFPROPS(propsOut, fluid, propsIn, prop1, prop2, arg = [1.0]):
    return RP.REFPROPdll(fluid, propsIn, propsOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output
