__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from streams import MassStream
from components import GasCooler
from scipy.optimize import root_scalar

class CustomRelationship:
    converged: bool = False
    r: bool = False

    def __init__(self, source_obj=0, target_obj=0, relation=0, tol=1e-9):
        self.tol = tol
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.relation = relation

    def solve(self):
        if self.relation == 0:
            r = self.Determine_T()
        elif self.relation == 1:
            r = self.Determine_Pint()
        elif self.relation == 2:
            r = self.Air_Mass_Flow()
        elif self.relation == 3:
            r = self.Airsource_Mass_Flow()
        elif self.relation == 4:
            r = self.Determine_EV_sm_outlet()
        self.converged = True
        return r
    
    def Determine_T(self): 
        self.target_obj.temperature = max(
            min(self.source_obj[0].temperature, self.source_obj[1].temperature) + self.var, 
            self.source_obj[0].temperature, 
            self.source_obj[1].temperature,
            )
        return True
    
    def Determine_Pint(self):
        #Determines the n-1 intermediate pressures between n compressors
        S0, S3, F2, F3 = self.source_obj
        log_P3_star = math.log(S3.pressure / S0.pressure)
        log_P1_star = (1 - F2 - F3) * log_P3_star
        log_P2_star = (1 - F3) * log_P3_star
        self.target_obj[0].pressure = S0.pressure * math.exp(log_P1_star)
        self.target_obj[1].pressure = S0.pressure * math.exp(log_P2_star)
        return True

    def Air_Mass_Flow(self):
        gc: GasCooler
        m_air_list = []
        #Hot side terminal temperature
        dQ = 0
        
        Q_tot = 0
        for gc in self.source_obj:
            Q_tot += gc.duty

        for gc in self.source_obj:
            dQ += gc.duty
            T_air_i = gc.source_in.temperature - gc.dTmin
            sink_in_i = MassStream(fluid=gc.sink_in.fluid, pp_instance=gc.sink_in.pp_instance, pressure=gc.sink_in.pressure, temperature=T_air_i)
            dh_air = sink_in_i.enthalpy - self.source_obj[0].sink_in.enthalpy
            
            if dh_air > 0.0 and gc.duty > 0.0:
                
                def get_mid_pinch_point_error(h_ref, mcp_air):
                    P_ave = (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy)
                    source_in_i = MassStream(fluid=gc.source_in.fluid, pp_instance=gc.source_in.pp_instance, pressure=P_ave, temperature=h_ref)
                    return mcp_air[0] - gc.source_in.mass_flow * source_in_i.heat_capacity

                def get_bounds_for_mid_pinch_point(mcp_air):
                    bnds = False
                    if mcp_air - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity or mcp_air + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity:
                        if mcp_air - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity and mcp_air + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity:
                            e0 = get_mid_pinch_point_error(h_ref=gc.source_in.enthalpy, mcp_air=[mcp_air])
                            e1 = get_mid_pinch_point_error(h_ref=gc.source_out.enthalpy, mcp_air=[mcp_air])
                            if e0 * e1 < 0.0:
                                bnds = (gc.source_in.enthalpy, gc.source_out.enthalpy)
                        else:
                            dh = gc.source_in.enthalpy - gc.source_out.enthalpy
                            if mcp_air - 1e-6 > gc.source_in.mass_flow * gc.source_in.heat_capacity:
                                h0 = gc.source_in.enthalpy
                                dh = gc.source_in.enthalpy - gc.source_out.enthalpy

                            else: #mcp_air + 1e-6 < gc.source_out.mass_flow * gc.source_out.heat_capacity
                                h0 = gc.source_out.enthalpy
                                dh = gc.source_out.enthalpy - gc.source_in.enthalpy

                            x = 0.1
                            e0 = get_mid_pinch_point_error(h_ref=h0, mcp_air=[mcp_air])
                            e1 = e0
                            while (e0 * e1 > 0.0 or abs(e1) < 1e-3) and (x < 1.0):
                                h = h0 - x * dh
                                e1 = get_mid_pinch_point_error(h_ref=h, mcp_air=[mcp_air])
                                x += 0.1
                            if e0 * e1 < 0.0:
                                bnds = (h0, h)

                    return bnds
                
                def get_mid_pinch_point(mcp_air, method, bnds):
                    result = root_scalar(f=get_mid_pinch_point_error,
                                         args=[mcp_air],
                                         method=method, 
                                         bracket=bnds, 
                                         x0=(bnds[0] + bnds[1]) / 2, 
                                         xtol=self.tol)
                    h_ref = result['root']
                    P_ref = (gc.source_in.pressure * (h_ref - gc.source_out.enthalpy) + gc.source_out.pressure * (gc.source_in.enthalpy - h_ref)) / (gc.source_in.enthalpy - gc.source_out.enthalpy)
                    source_in_i = MassStream(fluid=gc.source_in.fluid, pp_instance=gc.source_in.pp_instance, pressure=P_ref, enthalpy=h_ref)
                    e = abs(mcp_air - gc.source_in.mass_flow * source_in_i.heat_capacity) / mcp_air
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
                                source_in_i = MassStream(fluid=gc.source_in.fluid, pp_instance=gc.source_in.pp_instance, pressure=gc.source_in.pressure, enthalpy=h_ref)
                                T_ref_i = source_in_i.temperature
                                T_air_i = T_ref_i - gc.dTmin
                                sink_in_i = MassStream(fluid=gc.sink_in.fluid, pp_instance=gc.sink_in.pp_instance, pressure=gc.sink_in.pressure, temperature=T_air_i)
                                dh_air = sink_in_i.enthalpy - self.source_obj[0].sink_in.enthalpy
                                m_air = dQ_i / dh_air
                    return m_air

                m_air_list.append(dQ / dh_air)
                self.mcp_air_i = dQ / (T_air_i - self.source_obj[0].sink_in.temperature)
                m_air_list.append(get_m_air(self.mcp_air_i))              

                if gc.source_in.pressure < gc.source_in.P_crit:
                    source_in_i = MassStream(fluid=gc.source_in.fluid, pp_instance=gc.source_in.pp_instance, pressure=gc.source_in.pressure, quality=1.0)
                    h_ref = source_in_i.enthalpy
                    
                    if h_ref < gc.source_in.enthalpy and h_ref > gc.source_out.enthalpy:
                        dQ_i = dQ - gc.source_in.mass_flow * (gc.source_in.enthalpy - h_ref)
                        
                        if dQ_i > self.tol:
                            source_in_i = MassStream(fluid=gc.source_in.fluid, pp_instance=gc.source_in.pp_instance, pressure=gc.source_in.pressure, enthalpy=h_ref)
                            T_air_i = source_in_i.temperature - gc.dTmin
                            sink_in_i = MassStream(fluid=gc.sink_in.fluid, pp_instance=gc.sink_in.pp_instance, pressure=gc.sink_in.pressure, temperature=T_air_i)
                            dh_air = sink_in_i.enthalpy - self.source_obj[0].sink_in.enthalpy
                            m_air_list.append(dQ_i / dh_air)  
        
        if not m_air_list == []:
            self.target_obj[0].sink_in.mass_flow = max(m_air_list)
            self.target_obj[0].sink_in.props_in_mPT()
            
            for gc in self.target_obj:
                gc.sink_out.mass_flow = gc.sink_in.mass_flow
                gc.sink_out.enthalpy = gc.sink_in.enthalpy + (gc.source_in.heat_flow - gc.source_out.heat_flow) / gc.sink_in.mass_flow
                gc.sink_out.props_in_mPH()
        
        return m_air_list

    def Airsource_Mass_Flow(self):
        ev: GasCooler
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
            source_out_i = MassStream(fluid=ev.source_out.fluid, pp_instance=ev.source_out.pp_instance, pressure=ev.source_out.pressure, temperature=T_air_i)
            dh_air = self.source_obj[0].source_out.enthalpy - source_out_i.enthalpy
            
            if dh_air > 0:
                m_air_list.append(dQ / dh_air)
                self.mcp_air = dQ / (T_air_i - self.source_obj[0].source_out.temperature)

        self.target_obj[0].source_in.mass_flow = max(m_air_list)
        self.target_obj[0].source_out.mass_flow = self.target_obj[0].source_in.mass_flow
        self.target_obj[0].source_in.props_in_mPT()
        return True

    def Determine_EV_sm_outlet(self):
        EV, SP = self.target_obj
        Comp1, Econ, Comp_discharge, gc_2_dP_h = self.source_obj

        EV.flow_in.pressure = Comp_discharge.pressure - gc_2_dP_h - Econ.dP_h
        EV.flow_out.pressure = Comp1.flow_out.pressure + Econ.dP_c

        # First estimate
        if EV.flow_out.pressure < EV.flow_out.P_crit:
            EV.flow_in.mass_flow = Comp1.flow_out.mass_flow * SP.r_ratio
            Econ.source_out.mass_flow = Comp1.flow_out.mass_flow + EV.flow_in.mass_flow

            EV.flow_out.quality = 0.001
            EV.flow_out.props_in_mPQ()

            h_EV = EV.flow_out.enthalpy

            def _calc_dt_drop(h, F_in, F_out, dTmin):
                F_in.enthalpy = h
                F_in.props_in_mPH()
                F_out.enthalpy = h
                F_out.props_in_mPH()
                return F_in.temperature - F_out.temperature - dTmin

            h_min, h_max = EV.flow_out.enthalpy, Comp_discharge.enthalpy
            if _calc_dt_drop(h_min, EV.flow_in, EV.flow_out, EV.dTmin) < 0.0:
                r = root_scalar(f=lambda x: _calc_dt_drop(x, EV.flow_in, EV.flow_out, EV.dTmin), method='brentq', x0=h_min, bracket=[h_min, h_max])
                h_EV = r.root
    
        else:
            EV.flow_in.mass_flow = 0.0
            Econ.source_out.mass_flow = Comp1.flow_out.mass_flow
            h_EV = Comp_discharge.enthalpy


        EV.flow_out.enthalpy = h_EV
        EV.flow_out.props_in_mPH()

        EV.flow_in.enthalpy = h_EV
        EV.flow_in.props_in_mPH()

        SP.flow_in.pressure = Comp_discharge.pressure - gc_2_dP_h
        SP.flow_in.enthalpy = EV.flow_in.enthalpy
        SP.flow_in.props_in_mPH()
        
        pass
