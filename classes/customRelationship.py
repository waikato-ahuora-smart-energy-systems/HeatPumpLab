__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from globalVariables import CONVERGANCE
from globalFunctions import REFPROP

class CustomRelationship:
    def __init__(self, source_obj=0, target_obj=0, var1=0, var2=0, relation=0):
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.var1 = var1
        self.var2 = var2
        self.relation = relation
        self.solve()

    def solve(self):
        if self.relation == 0:
            self.target_obj.Temperature = max(min(self.source_obj[0].Temperature, self.source_obj[1].Temperature) + self.var1, self.source_obj[0].Temperature, self.source_obj[1].Temperature)

        elif self.relation == 1:
            self.Determine_Pint(self.source_obj[0].Pressure, self.source_obj[1].Pressure, self.source_obj[2], self.var1)

        elif self.relation == 2:
            self.Air_Mass_Flow()
    
    #Determines the n-1 intermediate pressures between n compressors
    def Determine_Pint(self, P_in, P_out, n_comp=1, PR_offset=[0.0]):
        
        for i in range(1, n_comp):
            self.target_obj[i - 1].Pressure = P_in ** (1 - i / n_comp) * P_out ** (i / n_comp)

            if (PR_offset[i - 1] > 0 and PR_offset[i - 1] <= 1.0) or (PR_offset[i - 1] < 0 and PR_offset[i - 1] >= -1.0): 
                if i == 1:
                    n_exp = math.log(self.target_obj[i - 1].Pressure) - math.log(P_in)
                else:
                    n_exp = math.log(self.target_obj[i - 1].Pressure) - math.log(self.target_obj[i - 2].Pressure)

                self.target_obj[i - 1].Pressure = self.target_obj[i - 1].Pressure * math.exp(n_exp * PR_offset[i - 1])

    def Air_Mass_Flow(self):
        mcp_air_list = []
        #Hot side terminal temperature
        dH = 0
        for gc in self.source_obj:
            gc_duty = gc.HotIn.Heat_Flow - gc.HotOut.Heat_Flow
            dH += gc_duty
            dT = gc.HotIn.Temperature - gc.dT_GC_min_tar - self.source_obj[0].ColdIn.Temperature
            mcp_air_list.append(dH / dT)

            if dH > gc.HotIn.Mass_Flow * gc.HotIn.Heat_Capacity * dT and dT > 0:
                x = gc_duty / gc.HotIn.Mass_Flow / 2
                h = gc.HotIn.Enthalpy
                dH_1 = 0

                while x > CONVERGANCE * gc_duty and gc_duty > 0 and h >= gc.HotOut.Enthalpy:
                    h -= x
                    dT = REFPROP('T', gc.HotIn.Fluid, 'PH', gc.HotIn.Pressure / 1000, h) - gc.dT_GC_min_tar - self.source_obj[0].ColdIn.Temperature
                    dH_1 = dH - gc.HotIn.Mass_Flow * (gc.HotIn.Enthalpy - h)
                    cp_i = REFPROP('CP', gc.HotIn.Fluid, 'PH', gc.HotIn.Pressure / 1000, h)

                    if dH_1 < gc.HotIn.Mass_Flow * cp_i * dT: 
                        h += x
                        x /= 2
                
                if dH_1 > 0.0:
                    mcp_air_list.append(dH_1 / dT)
        
        mcp_air = max(mcp_air_list)

        self.target_obj[0].ColdIn.Mass_Flow = mcp_air / self.source_obj[0].ColdIn.Heat_Capacity
        self.target_obj[0].ColdIn.PropsIn_mPT()
        
        for gc in self.target_obj:
            gc.ColdOut.Mass_Flow = gc.ColdIn.Mass_Flow
            gc.ColdOut.Enthalpy = gc.ColdIn.Enthalpy + (gc.HotIn.Heat_Flow - gc.HotOut.Heat_Flow) / gc.ColdIn.Mass_Flow
            gc.ColdOut.PropsIn_mPH()