__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math, pprint
from streams import get_property_value
from scipy.optimize import brentq, root_scalar

class CustomRelationship:
    def __init__(self, source_obj=0, target_obj=0, additional_obj=0, relation=0, tol=1e-9, RP=None):
        self.tol = tol
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.additional_obj = additional_obj
        self.relation = relation
        self.RP = RP

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

        return 0
    
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

        Q_tot = 0
        for gc in self.source_obj:
            Q_tot += gc.duty
        self.mcp_air_0 = Q_tot / (self.additional_obj[0] - self.source_obj[0].sink_in.temperature)

        for gc in self.source_obj:
            dQ += gc.duty
            T_air_i = gc.source_in.temperature - gc.dTmin
            dh_air = get_property_value('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i, RP=self.RP) - self.source_obj[0].sink_in.enthalpy
            
            if dh_air > 0.0 and gc.duty > 0.0:
                
                def get_mid_pinch_point_error(h_ref, mcp_air):
                    return mcp_air[0] - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP)

                def get_bounds_for_mid_pinch_point(mcp_air):
                    bnds = False
                    if mcp_air - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity or mcp_air + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity:
                        if mcp_air - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity and mcp_air + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity:
                            bnds = (gc.source_in.enthalpy, gc.source_out.enthalpy)
                        else:
                            dh = gc.source_in.enthalpy - gc.source_out.enthalpy
                            x = 0.1
                            if mcp_air - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity:
                                e0 = get_mid_pinch_point_error(h_ref=gc.source_in.enthalpy, mcp_air=[mcp_air])
                                e1 = e0
                                while e0 * e1 > 0.0 and x < 1.0:
                                    h = gc.source_in.enthalpy - x * dh
                                    e1 = get_mid_pinch_point_error(h_ref=h, mcp_air=[mcp_air])
                                    x += 0.1
                                if e0 * e1 < 0.0:
                                    bnds = (gc.source_in.enthalpy, h)

                            else: #mcp_air + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity
                                e0 = get_mid_pinch_point_error(h_ref=gc.source_out.enthalpy, mcp_air=[mcp_air])
                                e1 = e0
                                while e0 * e1 > 0.0 and x < 1.0:
                                    h = gc.source_out.enthalpy + x * dh
                                    e1 = get_mid_pinch_point_error(h_ref=h, mcp_air=[mcp_air])
                                    x += 0.1
                                if e0 * e1 < 0.0:
                                    bnds = (h, gc.source_out.enthalpy)

                    return bnds
                
                def get_mid_pinch_point(mcp_air, method, bnds):
                    result = root_scalar(f=get_mid_pinch_point_error,
                                         args=[mcp_air],
                                         method=method, 
                                         bracket=bnds, 
                                         x0=(bnds[0] + bnds[1]) / 2, 
                                         xtol=self.tol)
                    h_ref = result['root']
                    e = abs(mcp_air - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP)) / mcp_air
                    return h_ref, e
                
                def get_m_air(mcp_air):
                    m_air = 0
                    bnds = get_bounds_for_mid_pinch_point(mcp_air)
                    if bnds != False:
                        h_ref, e = get_mid_pinch_point(mcp_air=mcp_air, method='brentq', bnds=bnds)
                        if abs(e) > 1e-6:
                            h_ref, e = get_mid_pinch_point(mcp_air=mcp_air, method='bisect', bnds=bnds)

                        if abs(e) < 1e-6:
                            dQ_i = dQ - gc.source_in.mass_flow * (gc.source_in.enthalpy - h_ref)
                            if dQ_i > 1e-3:
                                T_ref_i = get_property_value('T', gc.source_in.fluid, 'PH', gc.source_in.pressure / 1000, h_ref, RP=self.RP)
                                T_air_i = T_ref_i - gc.dTmin
                                dh_air = get_property_value('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i, RP=self.RP) - self.source_obj[0].sink_in.enthalpy
                                m_air = dQ_i / dh_air
                    return m_air

                m_air_list.append(dQ / dh_air)
                m_air_list.append(get_m_air(self.mcp_air_0))
                self.mcp_air_i = dQ / (T_air_i - self.source_obj[0].sink_in.temperature)
                m_air_list.append(get_m_air(self.mcp_air_i))
                
                # if self.mcp_air_0 - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity or self.mcp_air_0 + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity:
                #     result = root_scalar(f=lambda h_ref: self.mcp_air_0 - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP), 
                #                         method='brentq', 
                #                         bracket=(gc.source_in.enthalpy, gc.source_out.enthalpy), 
                #                         x0=(gc.source_in.enthalpy + gc.source_out.enthalpy) / 2, 
                #                         rtol=self.tol)
                #     h_ref = result['root']
                #     e = abs(self.mcp_air_0 - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP)) / self.mcp_air_0
                #     if e > self.tol:
                #         result = root_scalar(f=lambda h_ref: self.mcp_air_0 - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP), 
                #                             method='bisect', 
                #                             bracket=(gc.source_in.enthalpy, gc.source_out.enthalpy), 
                #                             x0=gc.source_in.enthalpy, 
                #                             rtol=self.tol)
                #         h_ref = result['root']
                #         e = abs(self.mcp_air_0 - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP)) / self.mcp_air_0

                #     dQ_i = dQ - gc.source_in.mass_flow * (gc.source_in.enthalpy - h_ref)
                #     if dQ_i > 1e-3:
                #         T_air_i = get_property_value('T', gc.source_in.fluid, 'PH', gc.source_in.pressure / 1000, h_ref, RP=self.RP) - gc.dTmin
                #         dh_air = get_property_value('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i, RP=self.RP) - self.source_obj[0].sink_in.enthalpy
                #         m_air_list.append(dQ_i / dh_air)                

                # self.mcp_air_i = dQ / (T_air_i - self.source_obj[0].sink_in.temperature)
                # if self.mcp_air_i - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity or self.mcp_air_i + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity:
                #     result = root_scalar(f=lambda h_ref: self.mcp_air_i - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP), 
                #                         method='brentq', 
                #                         bracket=(gc.source_in.enthalpy, gc.source_out.enthalpy), 
                #                         x0=(gc.source_in.enthalpy + gc.source_out.enthalpy) / 2, 
                #                         rtol=self.tol)
                #     h_ref = result['root']
                #     e = abs(self.mcp_air_i - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP)) / self.mcp_air_i
                #     if e > self.tol:
                #         result = root_scalar(f=lambda h_ref: self.mcp_air_i - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP), 
                #                             method='bisect', 
                #                             bracket=(gc.source_in.enthalpy, gc.source_out.enthalpy), 
                #                             x0=gc.source_in.enthalpy, 
                #                             rtol=self.tol)
                #         h_ref = result['root']
                #         e = abs(self.mcp_air_i - gc.source_in.mass_flow * get_property_value('CP', gc.source_in.fluid, 'PH', (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy) / 1000, h_ref, RP=self.RP)) / self.mcp_air_i

                #     dQ_i = dQ - gc.source_in.mass_flow * (gc.source_in.enthalpy - h_ref)
                #     if dQ_i > 1e-3:
                #         T_air_i = get_property_value('T', gc.source_in.fluid, 'PH', gc.source_in.pressure / 1000, h_ref, RP=self.RP) - gc.dTmin
                #         dh_air = get_property_value('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i, RP=self.RP) - self.source_obj[0].sink_in.enthalpy
                #         m_air_list.append(dQ_i / dh_air)                

                P_crit = get_property_value('PCRIT', gc.source_in.fluid, 'PT', gc.source_in.pressure / 1000, gc.source_in.enthalpy, RP=self.RP)
                if gc.source_in.pressure < P_crit:
                    h_ref = get_property_value('H', gc.source_in.fluid, 'PQMASS', gc.source_in.pressure / 1000, 1.0, RP=self.RP)
                    
                    if h_ref < gc.source_in.enthalpy and h_ref > gc.source_out.enthalpy:
                        dQ_i = dQ - gc.source_in.mass_flow * (gc.source_in.enthalpy - h_ref)
                        
                        if dQ_i > self.tol:
                            T_air_i = get_property_value('T', gc.source_in.fluid, 'PH', gc.source_in.pressure / 1000, h_ref, RP=self.RP) - gc.dTmin
                            dh_air = get_property_value('H', gc.sink_in.fluid, 'PT', gc.sink_in.pressure / 1000, T_air_i, RP=self.RP) - self.source_obj[0].sink_in.enthalpy
                            m_air_list.append(dQ_i / dh_air)  
        
        if not m_air_list == []:
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
            dh_air = self.source_obj[0].source_out.enthalpy - get_property_value('H', ev.source_out.fluid, 'PT', ev.source_out.pressure / 1000, T_air_i, RP=self.RP) 
            
            if dh_air > 0:
                m_air_list.append(dQ / dh_air)
                self.mcp_air = dQ / (T_air_i - self.source_obj[0].source_out.temperature)

        self.target_obj[0].source_in.mass_flow = max(m_air_list)
        self.target_obj[0].source_out.mass_flow = self.target_obj[0].source_in.mass_flow
        self.target_obj[0].source_in.props_in_mPT()
        