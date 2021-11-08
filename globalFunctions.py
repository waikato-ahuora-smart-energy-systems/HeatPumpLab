__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

try:
    RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/REFPROP.DLL') #/REFPROP.DLL
except:
    RP = REFPROPFunctionLibrary('C:/Program Files (x86)/REFPROP/')
RP.SETPATHdll('C:/Program Files (x86)/REFPROP')
UNITS = RP.GETENUMdll(0,'SI with C').iEnum

def REFPROP(propOut, fluid, propsIn, prop1, prop2, arg = [1.0]): #Arguments order matches Excel call function
    return RP.REFPROPdll(fluid, propsIn, propOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output[0]

def REFPROPS(propsOut, fluid, propsIn, prop1, prop2, arg = [1.0]):
    return RP.REFPROPdll(fluid, propsIn, propsOut, UNITS, 0, 0, prop1, prop2, [1.0]).Output

def HAPROPS(propOut, propsIn, P, T_or_H, AH, arg = [1.0]):

    fract_h2o = AH / (1 + AH)
    fract_air = 1 - fract_h2o
    P_vap = P * fract_h2o
    
    if (propOut == 'h' or propOut == 'H') and propsIn == '.+T':
        T = T_or_H
        P_sat = RP.REFPROPdll('Water', 'TQ', 'P', UNITS, 0, 0, T, 1.0, [1.0]).Output[0]
        RH = P_vap / P_sat
        
        h_air = RP.REFPROPdll('Air', 'PT', propOut, UNITS, 0, 0, P - P_vap, T, [1.0]).Output[0]
        h_vap = RP.REFPROPdll('Water', 'PT', propOut, UNITS, 0, 0, P_vap, T, [1.0]).Output[0]
        return fract_h2o * h_vap + fract_air * h_air
    elif (propOut == 'hdew' or propOut == 'Hdew') and propsIn == 'P.+':
        T_sat = RP.REFPROPdll('Water', 'PQ', 'T', UNITS, 0, 0, P, 1.0, [1.0]).Output[0]

        h_air = RP.REFPROPdll('Air', 'PT', propOut, UNITS, 0, 0, P - P_vap, T_sat, [1.0]).Output[0]
        h_vap = RP.REFPROPdll('Water', 'PQ', propOut, UNITS, 0, 0, P_vap, 1.0, [1.0]).Output[0]

        return fract_h2o * h_vap + fract_air * h_air

def write_output(filename, to_write):
    f = open(filename, 'a')
    [print(elements, file=f) for elements in to_write]         
    print('\n', file=f)
    f.close()  