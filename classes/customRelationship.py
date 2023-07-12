__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from globalVariables import CONVERGANCE
from functions.refprop import REFPROP
from scipy import optimize

class CustomRelationship:
    def __init__(self, source_obj=0, target_obj=0, relation=0):
        self.CONVERGANCE = CONVERGANCE
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.relation = relation
        self.solve()

    def solve(self):
        if self.relation == 0:
            self.Determine_T()

        elif self.relation == 1:
            self.Determine_Pint()

        elif self.relation == 2:
            self.Air_Mass_Flow()

        elif self.relation == 3:
            self.Airsource_Mass_Flow()

        else:
            pass
    
    def Determine_T(self):
        self.target_obj.Temperature = max(
            min(self.source_obj[0].Temperature, self.source_obj[1].Temperature) + self.var, 
            self.source_obj[0].Temperature, 
            self.source_obj[1].Temperature
            )
    
    def Determine_Pint(self):
        #Determines the n-1 intermediate pressures between n compressors
        P_comp1_in = self.source_obj[0].Pressure 
        P_comp1_out = self.source_obj[1].Pressure
        comp1_power_frac = self.source_obj[2]
        comp2_power_frac = self.source_obj[3]
        comp3_power_frac = self.source_obj[4]

        PR_comp1 = P_comp1_out / P_comp1_in
        self.target_obj[1].Pressure = P_comp1_out * PR_comp1 * comp2_power_frac / comp1_power_frac if comp2_power_frac > 0.0 else P_comp1_out
        self.target_obj[2].Pressure = self.target_obj[1].Pressure * PR_comp1 * comp3_power_frac / comp1_power_frac if comp3_power_frac > 0.0 else self.target_obj[1].Pressure


    def Air_Mass_Flow(self):
        m_air_list = []
        #Hot side terminal temperature
        dQ = 0
        for gc in self.source_obj:
            gc.HotOut.Pressure = gc.HotIn.Pressure - gc.dP_refrigerant
            gc.HotOut.PropsIn_mPT()
            gc.Duty = gc.HotIn.Heat_Flow - gc.HotOut.Heat_Flow
            
            gc_duty = gc.HotIn.Heat_Flow - gc.HotOut.Heat_Flow
            dQ += gc_duty
            T_air_i = gc.HotIn.Temperature - gc.dT_GC_min_tar
            dh_air = REFPROP('H', gc.ColdIn.Fluid, 'PT', gc.ColdIn.Pressure / 1000, T_air_i) - self.source_obj[0].ColdIn.Enthalpy
            
            if dh_air > 0:
                m_air_list.append(dQ / dh_air)
                self.mcp_air = dQ / (T_air_i - self.source_obj[0].ColdIn.Temperature)

                P_crit = REFPROP('PCRIT', gc.HotIn.Fluid, 'PT', gc.HotIn.Pressure / 1000, gc.HotIn.Enthalpy) * 1000

                if self.mcp_air > gc.HotIn.Mass_Flow * gc.HotIn.Heat_Capacity and self.mcp_air < gc.HotOut.Mass_Flow * gc.HotOut.Heat_Capacity:
                    #a = self.mcp_air - gc.HotOut.Mass_Flow * REFPROP('CP', gc.HotOut.Fluid, 'PH', gc.HotOut.Pressure / 1000, gc.HotOut.Enthalpy)
                    #b = self.mcp_air - gc.HotIn.Mass_Flow * REFPROP('CP', gc.HotIn.Fluid, 'PH', gc.HotIn.Pressure / 1000, gc.HotIn.
                    h_ref = optimize.brentq(f=lambda h_ref: self.mcp_air - gc.HotIn.Mass_Flow * REFPROP('CP', gc.HotIn.Fluid, 'PH', (gc.HotIn.Pressure * (h_ref - gc.HotOut.Enthalpy) + gc.HotOut.Pressure * (gc.HotIn.Enthalpy - h_ref)) / (gc.HotIn.Enthalpy - gc.HotOut.Enthalpy) / 1000, h_ref), a=gc.HotIn.Enthalpy, b=gc.HotOut.Enthalpy)
                    dQ_i = dQ - gc.HotIn.Mass_Flow * (gc.HotIn.Enthalpy - h_ref)
                    if dQ_i > self.CONVERGANCE:
                        T_air_i = REFPROP('T', gc.HotIn.Fluid, 'PH', gc.HotIn.Pressure / 1000, h_ref) - gc.dT_GC_min_tar
                        dh_air = REFPROP('H', gc.ColdIn.Fluid, 'PT', gc.ColdIn.Pressure / 1000, T_air_i) - self.source_obj[0].ColdIn.Enthalpy
                        m_air_list.append(dQ_i / dh_air)                

                if gc.HotIn.Pressure < REFPROP('PCRIT', gc.HotIn.Fluid, 'PT', gc.HotIn.Pressure / 1000, gc.HotIn.Enthalpy) * 1000:
                    h_ref = REFPROP('H', gc.HotIn.Fluid, 'PQMASS', gc.HotIn.Pressure / 1000, 1.0)
                    
                    if h_ref < gc.HotIn.Enthalpy and h_ref > gc.HotOut.Enthalpy:
                        dQ_i = dQ - gc.HotIn.Mass_Flow * (gc.HotIn.Enthalpy - h_ref)
                        if dQ_i > self.CONVERGANCE:
                            T_air_i = REFPROP('T', gc.HotIn.Fluid, 'PH', gc.HotIn.Pressure / 1000, h_ref) - gc.dT_GC_min_tar
                            dh_air = REFPROP('H', gc.ColdIn.Fluid, 'PT', gc.ColdIn.Pressure / 1000, T_air_i) - self.source_obj[0].ColdIn.Enthalpy
                            m_air_list.append(dQ_i / dh_air)  

        self.target_obj[0].ColdIn.Mass_Flow = max(m_air_list)
        self.target_obj[0].ColdIn.PropsIn_mPT()
        
        for gc in self.target_obj:
            gc.ColdOut.Mass_Flow = gc.ColdIn.Mass_Flow
            gc.ColdOut.Enthalpy = gc.ColdIn.Enthalpy + (gc.HotIn.Heat_Flow - gc.HotOut.Heat_Flow) / gc.ColdIn.Mass_Flow
            gc.ColdOut.PropsIn_mPH()


    def Airsource_Mass_Flow(self):
        m_air_list = []
        #Hot side terminal temperature
        dQ = 0
        for ev in self.source_obj:
            ev.ColdOut.Pressure = ev.ColdIn.Pressure
            ev.HotOut.PropsIn_mPT()
            ev.Duty = ev.ColdOut.Heat_Flow - ev.ColdIn.Heat_Flow
            
            ev_duty = ev.ColdOut.Heat_Flow - ev.ColdIn.Heat_Flow
            dQ += ev_duty
            T_air_i = ev.ColdOut.Temperature - 25
            dh_air = self.source_obj[0].HotOut.Enthalpy - REFPROP('H', ev.HotOut.Fluid, 'PT', ev.HotOut.Pressure / 1000, T_air_i) 
            
            if dh_air > 0:
                m_air_list.append(dQ / dh_air)
                self.mcp_air = dQ / (T_air_i - self.source_obj[0].HotOut.Temperature)

        self.target_obj[0].HotIn.Mass_Flow = max(m_air_list)
        self.target_obj[0].HotOut.Mass_Flow = self.target_obj[0].HotIn.Mass_Flow
        self.target_obj[0].HotIn.PropsIn_mPT()