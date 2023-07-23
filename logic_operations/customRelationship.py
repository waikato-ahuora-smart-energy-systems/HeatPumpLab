__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from streams import get_property
from scipy.optimize import brentq

class CustomRelationship:
    def __init__(self, source_obj=0, target_obj=0, relation=0, tol=1e-6):
        self.tol = tol
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.relation = relation

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
        self.target_obj.temperature = max(
            min(self.source_obj[0].temperature, self.source_obj[1].temperature) + self.var, 
            self.source_obj[0].temperature, 
            self.source_obj[1].temperature
            )
    
    def Determine_Pint(self):
        #Determines the n-1 intermediate pressures between n compressors
        S0, S3, F2, F3 = self.source_obj
        log_P3_star = math.log(S3.pressure / S0.pressure)
        log_P1_star = (1 - F2 - F3) * log_P3_star
        log_P2_star = (1 - F3) * log_P3_star
        self.target_obj[0].pressure = S0.pressure * math.exp(log_P1_star)
        self.target_obj[1].pressure = S0.pressure * math.exp(log_P2_star)

    def Air_Mass_Flow(self):
        m_air_list = []
        #Hot side terminal temperature
        dQ = 0
        for gc in self.source_obj:
            gc.source_out.pressure = gc.source_in.pressure - gc.dP_h
            gc.source_out.props_in_mPT()
            gc.duty = gc.source_in.heat_flow - gc.source_out.heat_flow
            
            gc_duty = gc.source_in.heat_flow - gc.source_out.heat_flow
            dQ += gc_duty
            T_air_i = gc.source_in.temperature - gc.dTmin
            dh_air = get_property('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i) - self.source_obj[0].sink_in.enthalpy
            
            if dh_air > 0:
                m_air_list.append(dQ / dh_air)
                self.mcp_air = dQ / (T_air_i - self.source_obj[0].sink_in.temperature)

                if self.mcp_air > gc.source_in.mass_flow * gc.source_in.heat_capacity and self.mcp_air < gc.source_out.mass_flow * gc.source_out.heat_capacity:
                    #a = self.mcp_air - gc.source_out.mass_flow * get_property('CP', gc.source_out.fluid, 'PH', gc.source_out.pressure / 1000, gc.source_out.enthalpy)
                    #b = self.mcp_air - gc.source_in.mass_flow * get_property('CP', gc.source_in.fluid, 'PH', gc.source_in.pressure / 1000, gc.source_in.
                    h_ref = brentq(f=lambda h_ref: self.mcp_air - gc.source_in.mass_flow * get_property('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref), a=gc.source_in.enthalpy, b=gc.source_out.enthalpy)
                    dQ_i = dQ - gc.source_in.mass_flow * (gc.source_in.enthalpy - h_ref)
                    if dQ_i > self.tol:
                        T_air_i = get_property('T', gc.source_in.fluid, 'PH', gc.source_in.pressure / 1000, h_ref) - gc.dTmin
                        dh_air = get_property('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i) - self.source_obj[0].sink_in.enthalpy
                        m_air_list.append(dQ_i / dh_air)                

                if gc.source_in.pressure < get_property('PCRIT', gc.source_in.fluid, 'PT', gc.source_in.pressure / 1000, gc.source_in.enthalpy) * 1000:
                    h_ref = get_property('H', gc.source_in.fluid, 'PQMASS', gc.source_in.pressure / 1000, 1.0)
                    
                    if h_ref < gc.source_in.enthalpy and h_ref > gc.source_out.enthalpy:
                        dQ_i = dQ - gc.source_in.mass_flow * (gc.source_in.enthalpy - h_ref)
                        if dQ_i > self.tol:
                            T_air_i = get_property('T', gc.source_in.fluid, 'PH', gc.source_in.pressure / 1000, h_ref) - gc.dTmin
                            dh_air = get_property('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i) - self.source_obj[0].sink_in.enthalpy
                            m_air_list.append(dQ_i / dh_air)  

        self.target_obj[0].sink_in.mass_flow = max(m_air_list)
        self.target_obj[0].sink_in.props_in_mPT()
        
        for gc in self.target_obj:
            gc.sink_out.mass_flow = gc.sink_in.mass_flow
            gc.sink_out.enthalpy = gc.sink_in.enthalpy + (gc.source_in.heat_flow - gc.source_out.heat_flow) / gc.sink_in.mass_flow
            gc.sink_out.props_in_mPH()

    def Airsource_Mass_Flow(self):
        m_air_list = []
        #Hot side terminal temperature
        dQ = 0
        for ev in self.source_obj:
            ev.sink_out.pressure = ev.sink_in.pressure
            ev.source_out.props_in_mPT()
            ev.duty = ev.sink_out.heat_flow - ev.sink_in.heat_flow
            
            ev_duty = ev.sink_out.heat_flow - ev.sink_in.heat_flow
            dQ += ev_duty
            T_air_i = ev.sink_out.temperature - 25
            dh_air = self.source_obj[0].source_out.enthalpy - get_property('H', ev.source_out.fluid, 'PT', ev.source_out.pressure / 1000, T_air_i) 
            
            if dh_air > 0:
                m_air_list.append(dQ / dh_air)
                self.mcp_air = dQ / (T_air_i - self.source_obj[0].source_out.temperature)

        self.target_obj[0].source_in.mass_flow = max(m_air_list)
        self.target_obj[0].source_out.mass_flow = self.target_obj[0].source_in.mass_flow
        self.target_obj[0].source_in.props_in_mPT()
        