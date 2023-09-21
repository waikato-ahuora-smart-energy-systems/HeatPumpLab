__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
# try:
#     RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/') 
# except:
#     RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/REFPROP.DLL')
# RP.SETPATHdll('C:/Program Files (x86)/REFPROP')
# UNITS = RP.GETENUMdll(0,'SI with C').iEnum


def get_property(propOut, fluid, propsIn, prop1, prop2, arg = [1.0]): 
    """ 
    Fully encapsulated for multiprocessing support 
    Fresh copy of REFPROP will be loaded into each instance, for parallelism
    """
    root = 'C:/Program Files (x86)/REFPROP/'
    R = REFPROPFunctionLibrary(root)
    R.SETPATHdll(root)
    UNITS = R.GETENUMdll(0,'SI with C').iEnum
    return R.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0]
    # return RP.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0]

def get_properties(propsOut, fluid, propsIn, prop1, prop2, arg = [1.0]):
    root = 'C:/Program Files (x86)/REFPROP/'
    R = REFPROPFunctionLibrary(root)
    R.SETPATHdll(root)
    UNITS = R.GETENUMdll(0,'SI with C').iEnum
    return R.REFPROPdll(fluid, propsIn, propsOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output
    # return RP.REFPROPdll(fluid, propsIn, propsOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output
