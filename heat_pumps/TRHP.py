__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator
from streams import MassStream, Properties
from logic_operations import SetRelationship, CustomRelationship
from scipy.optimize import root_scalar 
from scipy.interpolate import CubicSpline, interp1d
from functions import Graph, print_result_to_csv
import numpy as np

class TRHP():
    def __init__(
            self, 
            refrigerant='R1336mzz(Z)', 
            source='Water',
            sink='Air',
            T_air_tar=200.0, 
            T_air_inlet=15.0, 
            dT_GC_min_tar=20.0, 
            T_gc_out_min=0.0, 
            T_source_in=40.0, 
            Tsat_evap=10.0,
            dTmin_phx=5.0,
            min_dT_ihx = 5.0,
            n_comp=3, 
            comp_PR_uniform=False, 
            eff_motor=0.96, 
            dT_phx1_cold=0.0, 
            dT_phx2_cold=0.0, 
            dT_phx3_cold=0.0, 
            comp2_power_frac=0.0, 
            comp3_power_frac=0.0,
            P_H_est=0.0,
            tol=1e-12,
            display=False,
            result_parameters=None
        ):
        
        self.refrigerant = refrigerant
        self.cycle_type = None
        self.T_air_tar = float(T_air_tar)
        self.dT_GC_min_tar = float(dT_GC_min_tar)
        self.T_air_inlet = float(T_air_inlet)
        self.T_source_in = float(T_source_in)
        self.T_gc_out = max(self.T_air_inlet + self.dT_GC_min_tar, T_gc_out_min)
        self.Tsat_evap_tar = float(Tsat_evap)
        self.Tsat_evap = float(Tsat_evap)
        self.n_comp = int(n_comp)
        self.comp2_power_frac = float(comp2_power_frac)
        self.comp3_power_frac = float(comp3_power_frac)
        self.comp_PR_uniform = bool(comp_PR_uniform)
        self.dTmin = float(dTmin_phx)
        self.eff_motor = float(eff_motor)
        self.min_dT_ihx = float(min_dT_ihx)
        self.dT_phx1_cold = float(dT_phx1_cold) #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.dT_phx2_cold = float(dT_phx2_cold) #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.dT_phx3_cold = float(dT_phx3_cold) #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.P_H_est = float(P_H_est)
        self.tol = float(tol)
        self.converged = False
        self.COP = -np.inf
        self.disp = display

        self.RP_ref = Properties(fld=refrigerant)
        self.RP_sink = Properties(fld=sink)
        self.RP_src = Properties(fld=source)

        self.init_components()
        self.init_structure()
        self.init_streams()
        self.init_logic_blocks()
        self.init_solver_pathway()
   
    def init_components(self):
        #Initialisation of unit operations
        self.comp1 = Compressor(name='Compressor', motor_eff=self.eff_motor, model_index=1)
        self.comp2 = Compressor(name='Compressor', motor_eff=self.eff_motor, model_index=1) 
        self.comp3 = Compressor(name='Compressor', motor_eff=self.eff_motor, model_index=1) 
        self.gas_cooler1 = GasCooler(name='Gas cooler', dTmin=self.dT_GC_min_tar)
        self.gas_cooler2 = GasCooler(name='Gas cooler', dTmin=self.dT_GC_min_tar)
        self.gas_cooler3 = GasCooler(name='Gas cooler', dTmin=self.dT_GC_min_tar)
        self.gas_cooler4 = GasCooler(name='Gas cooler', dTmin=self.dT_GC_min_tar)
        self.ihx = PlateHeatExchanger(name='Internal heat exchanger', dTmin=self.dTmin, model='dTmin_terminal')
        self.phx1 = PlateHeatExchanger(name='Recuperator', dT_c=self.dT_phx1_cold, model='fixed cold delta temperature')
        self.phx2 = PlateHeatExchanger(name='Recuperator', dT_c=self.dT_phx2_cold, model='fixed cold delta temperature')        
        self.phx3 = PlateHeatExchanger(name='Recuperator', dT_c=self.dT_phx3_cold, model='fixed cold delta temperature')
        self.exp_valve = ExpansionValve(name='Expansion valve')
        self.evaporator = Evaporator(name='Evaporator', model='saturated gas outlet evaporator')
        self.evap_sup = PlateHeatExchanger(name='Evaporator', model='fixed cold temperature outlet')

    def init_structure(self):
        #Check the PHE combo
        self.dT_phx1_cold = self.dT_phx1_cold if self.dT_phx1_cold > 0.0 else 0.0
        self.dT_phx2_cold = self.dT_phx2_cold if self.n_comp > 1 and self.dT_phx2_cold > 0.0 else 0.0
        self.dT_phx3_cold = self.dT_phx3_cold if self.n_comp > 2 and self.dT_phx3_cold > 0.0 else 0.0
        
        #Check the compressor combo
        if self.comp_PR_uniform or self.comp2_power_frac + self.comp3_power_frac >= 1.0:
            self.comp2_power_frac = 1 / self.n_comp if self.n_comp > 1 else 0.0
            self.comp3_power_frac = 1 / self.n_comp if self.n_comp > 2 else 0.0
        else:
            self.comp2_power_frac = self.comp2_power_frac if self.n_comp > 1 else 0.0
            self.comp3_power_frac = self.comp3_power_frac if self.n_comp > 2 else 0.0
        
        if self.comp2_power_frac == 0.0:
            self.comp3_power_frac = 0.0
        
        if self.n_comp == 3 and self.comp3_power_frac == 0.0:
            self.n_comp = 2
        
        if self.n_comp == 2 and self.comp2_power_frac == 0.0:
            self.n_comp = 1

        self.comp1_power_frac = 1 - (self.comp2_power_frac + self.comp3_power_frac)

        #Binary variables for the heat pump structure
        self.comp2_on = 1  if self.n_comp > 1 else 0
        self.comp3_on = 1  if self.n_comp > 2 else 0
        self.ihx_on = 1    if self.T_source_in < self.T_gc_out else 0
                
        self.phx1_on = 1   if self.dT_phx1_cold > 0.0 else 0
        self.phx2_on = 1   if self.dT_phx2_cold > 0.0 else 0
        self.phx3_on = 1   if self.dT_phx3_cold > 0.0 else 0

    def init_streams(self):
        #Evaporator superheat section
        self.evap_sup.source_in = MassStream( name='source In', fluid=self.RP_src, pressure=0.101, temperature=self.T_source_in, mass_flow=10^6)
        self.evap_sup.source_out = MassStream( name='source Out', fluid=self.RP_src, pressure=0.101, temperature=self.T_source_in, mass_flow=10^6)
        self.evap_sup.sink_in = MassStream( name='LP Ref: Evap to Super', fluid=self.RP_ref, temperature=self.Tsat_evap, quality=1.0)
        P_L = self.evap_sup.sink_in.pressure
        self.evap_sup.sink_out = MassStream( name='LP Ref: Super to IHX', fluid=self.RP_ref, pressure=P_L, temperature=(self.T_source_in - self.dTmin))
        
        #Internal heat exchanger
        self.ihx.sink_in = self.evap_sup.sink_out
        if self.ihx_on == 1:
            self.ihx.sink_out = MassStream( name='LP Ref: IHX to PHX1', fluid=self.RP_ref, pressure=P_L, temperature=(self.T_gc_out - self.dTmin))
        else:   
            self.ihx.sink_out = self.ihx.sink_in

        #Plate heat exchanger 1 - cold-side
        self.phx1.sink_in = self.ihx.sink_out
        if self.phx1_on == 1: 
            self.phx1.sink_out = MassStream( name='LP Ref: PHX1 to Comp1', fluid=self.RP_ref, pressure=P_L, temperature=(self.phx1.sink_in.temperature + self.dT_phx1_cold))
        else:
            self.phx1.sink_out = self.phx1.sink_in

        #Final discharge pressure
        if self.P_H_est > 0.0:
            P_H = self.P_H_est
        else:
            P_H = self.RP_ref.get_properties('P', 'TS', self.T_air_tar + self.dT_GC_min_tar * 2, self.phx1.sink_out.entropy)[0]
            if P_H < P_L:
                P_H = P_L * 16
            
        PR = (P_H / P_L) ** (1 / self.n_comp)
        
        #Compressor 1
        self.comp1.flow_in = self.phx1.sink_out
        self.comp1.flow_out = MassStream( name='Comp1_out', fluid=self.RP_ref, pressure=(P_L * PR), entropy=self.phx1.sink_out.entropy)

        #Plate heat exchanger 2 - cold-side
        self.phx2.sink_in = self.comp1.flow_out
        if self.phx2_on == 1:
            self.phx2.sink_out = MassStream( name='Comp2_in', fluid=self.RP_ref, pressure=self.phx2.sink_in.pressure, temperature=self.phx2.sink_in.temperature + self.dT_phx2_cold)
        else:
            self.phx2.sink_out = self.comp1.flow_out
        
        #Compressor 2
        self.comp2.flow_in = self.phx2.sink_out
        if self.comp2_on == 1:
            self.comp2.flow_out = MassStream( name='Comp2_out', fluid=self.RP_ref, pressure=(P_L * PR * PR), entropy=self.phx1.sink_out.entropy)  
        else:
            self.comp2.flow_out = self.comp2.flow_in

        #Plate heat exchanger 3 - cold-side
        self.phx3.sink_in = self.comp2.flow_out
        if self.phx3_on == 1:
            self.phx3.sink_out = MassStream( name='Comp3_in', fluid=self.RP_ref, pressure=self.phx3.sink_in.pressure, temperature=self.phx3.sink_in.temperature + self.dT_phx3_cold)
        else:
            self.phx3.sink_out = self.phx3.sink_in

        #Compressor 3
        self.comp3.flow_in = self.phx3.sink_out
        if self.comp3_on == 1:
            self.comp3.flow_out = MassStream( name='Ref: Comp3 to GC4', fluid=self.RP_ref, pressure=(P_L * PR * PR * PR), entropy=self.phx1.sink_out.entropy)
        else:
            self.comp3.flow_out = self.comp3.flow_in

        #Gas cooler 4 and Plate heat exchanger 3 - hot-sides
        self.gas_cooler4.sink_out = MassStream( name='Air: GC4 to Dryer', fluid=self.RP_sink, pressure=0.101, temperature=self.T_air_tar )
        self.gas_cooler4.sink_in = MassStream( name='Air: GC3 to GC4', fluid=self.RP_sink, pressure=0.101, temperature=self.T_air_tar )
        
        self.gas_cooler4.source_in = self.comp3.flow_out
        if self.comp3_on == 1 and self.phx3_on == 1:
            self.gas_cooler4.source_out = MassStream( name='Ref: GC4 to PHX3', fluid=self.RP_ref )
            self.phx3.source_in = self.gas_cooler4.source_out
            self.phx3.source_out = MassStream( name='Ref: PHX3 to GC3', fluid=self.RP_ref )       
        else:
            self.gas_cooler4.source_out = self.gas_cooler4.source_in
            self.phx3.source_in = self.gas_cooler4.source_out
            self.phx3.source_out = self.phx3.source_in       
        
        #Gas cooler 3 and Plate heat exchanger 2 - hot-sides
        self.gas_cooler3.sink_out = self.gas_cooler4.sink_in
        self.gas_cooler3.sink_in = MassStream( name='Air: GC2 to GC3', fluid=self.RP_sink, pressure=0.101, temperature=self.T_air_tar )
        
        self.gas_cooler3.source_in = self.phx3.source_out
        if self.comp2_on == 1 and self.phx2_on == 1:
            self.gas_cooler3.source_out = MassStream( name='Ref: GC3 to PHX2', fluid=self.RP_ref )
            self.phx2.source_in = self.gas_cooler3.source_out
            self.phx2.source_out = MassStream( name='Ref: PHX2 to GC2', fluid=self.RP_ref )
        else:
            self.gas_cooler3.source_out = self.gas_cooler3.source_in
            self.phx2.source_in = self.gas_cooler3.source_out
            self.phx2.source_out = self.phx2.source_in

        #Gas cooler 2 and Plate heat exchanger 1 - hot-sides
        self.gas_cooler2.sink_out = self.gas_cooler3.sink_in
        self.gas_cooler2.sink_in = MassStream( name='Air: GC1 to GC2', fluid=self.RP_sink, pressure=0.101, temperature=self.T_air_tar )

        self.gas_cooler2.source_in = self.phx2.source_out
        if self.phx1_on == 1:
            self.gas_cooler2.source_out = MassStream( name='Ref: GC2 to PHX1', fluid=self.RP_ref )
            self.phx1.source_in = self.gas_cooler2.source_out
            self.phx1.source_out = MassStream( name='Ref: PHX1 to GC1', fluid=self.RP_ref )
        else:
            self.gas_cooler2.source_out = self.gas_cooler2.source_in
            self.phx1.source_in = self.gas_cooler2.source_out
            self.phx1.source_out = self.phx1.source_in

        #Gas cooler 1
        self.gas_cooler1.sink_in = MassStream( name='Air: Inlet to GC1', fluid=self.RP_sink, pressure=0.101, temperature=self.T_air_inlet )
        self.gas_cooler1.sink_out = self.gas_cooler2.sink_in

        self.gas_cooler1.source_out = MassStream( name='Ref: GC1 to IHX', fluid=self.RP_ref, pressure=P_H, temperature=self.T_gc_out )
        self.gas_cooler1.source_in = self.phx1.source_out

        #Internal heat exchanger
        self.ihx.source_in = self.gas_cooler1.source_out
        if self.ihx_on == 1:
            self.ihx.source_out = MassStream( name='Ref: IHX to Exp V', fluid=self.RP_ref, pressure=P_H, temperature=self.T_gc_out )
        else:   
            self.ihx.source_out = self.ihx.source_in
        
        #Expansion valve
        self.exp_valve.flow_in = self.ihx.source_out
        self.exp_valve.flow_out = MassStream( name='Evap_in', fluid=self.RP_ref, pressure=P_L, enthalpy=self.exp_valve.flow_in.enthalpy )

        #Evaporator
        self.evaporator.sink_in = self.exp_valve.flow_out
        self.evaporator.sink_out = self.evap_sup.sink_in
        self.evaporator.source_in = MassStream( name='Exhaust air in', fluid=self.RP_src, pressure=0.101, temperature=self.T_source_in, mass_flow=1000)
        self.evaporator.source_out = MassStream( name='Exhaust air out', fluid=self.RP_src, pressure=0.101, temperature=min(self.evaporator.sink_in.temperature + self.dTmin, self.T_source_in), mass_flow=1000)

        self.streams = [self.evaporator.sink_out, 
                        self.comp1.flow_out, 
                        self.comp2.flow_out, 
                        self.comp3.flow_out, 
                        self.gas_cooler1.sink_out, 
                        self.gas_cooler2.sink_out, 
                        self.gas_cooler3.sink_out, 
                        self.gas_cooler4.sink_out, 
                        self.ihx.source_out, 
                        self.exp_valve.flow_out]

    def init_logic_blocks(self):
        #Set special relationships between points in the cycle 
        if self.phx3_on == 1:
            self.set_phx3_dT_cold = SetRelationship( source_obj=self.phx3.sink_in, target_obj=[self.phx3.sink_out], c=self.dT_phx3_cold, prop='temperature', prop2='pressure' ) 
            self.set_gc4_hot_out = SetRelationship( source_obj=self.phx3.sink_out, target_obj=[self.phx3.source_in], c=self.dTmin, prop='temperature', prop2='pressure' ) 

        if self.phx2_on == 1:    
            self.set_phx2_dT_cold = SetRelationship( source_obj=self.phx2.sink_in, target_obj=[self.phx2.sink_out], c=self.dT_phx2_cold, prop='temperature', prop2='pressure' ) 
            self.set_gc3_hot_out = SetRelationship( source_obj=self.phx2.sink_out, target_obj=[self.phx2.source_in], c=self.dTmin, prop='temperature', prop2='pressure' ) 
        
        if self.phx1_on == 1:
            self.set_phx1_dT_cold = SetRelationship( source_obj=self.phx1.sink_in, target_obj=[self.phx1.sink_out], c=self.dT_phx1_cold, prop='temperature', prop2='pressure' ) 
            self.set_gc2_hot_out = SetRelationship( source_obj=self.phx1.sink_out, target_obj=[self.phx1.source_in], c=self.dTmin, prop='temperature', prop2='pressure' ) 

        self.set_P_ints = CustomRelationship( source_obj=[self.comp1.flow_in, self.comp3.flow_out, self.comp2_power_frac, self.comp3_power_frac], target_obj=[self.comp1.flow_out, self.comp2.flow_out], relation=1, tol=self.tol )
        self.set_gc_air_massflow = CustomRelationship( source_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], target_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], relation=2, tol=self.tol )
        #self.set_evap_air_massflow = CustomRelationship( source_obj=[self.evaporator, self.evap_sup], target_obj=[self.evaporator, self.evap_sup], relation=3)

    def init_solver_pathway(self):
        #Solution pathway set-up
        self.operations = [self.evap_sup, self.set_P_ints]

        if self.ihx_on == 1:
            self.operations.append(self.ihx)

        if self.phx1_on == 1:
            self.operations.append(self.set_phx1_dT_cold)
            self.operations.append(self.set_gc2_hot_out)
        self.operations.append(self.comp1)
        if self.phx2_on == 1:
            self.operations.append(self.set_phx2_dT_cold)
            self.operations.append(self.set_gc3_hot_out)
        if self.comp2_on == 1:
            self.operations.append(self.comp2)
        if self.phx3_on == 1:
            self.operations.append(self.set_phx3_dT_cold)
            self.operations.append(self.set_gc4_hot_out)
        if self.comp3_on == 1:
            self.operations.append(self.comp3)
        
        if self.phx3_on == 1:
            self.operations.append(self.gas_cooler4)
            self.operations.append(self.phx3)
        if self.phx2_on == 1:
            self.operations.append(self.gas_cooler3)
            self.operations.append(self.phx2)
        if self.phx1_on == 1:
            self.operations.append(self.gas_cooler2)        
            self.operations.append(self.phx1)
        self.operations.append(self.gas_cooler1)

        if self.ihx_on == 1:
            self.operations.append(self.ihx)

        self.operations.append(self.exp_valve)
        self.operations.append(self.evaporator)

    def get_pressure_bounds(self):
        feasible = 0
        P_lb = self.comp3.flow_out.pressure
        P_ub = P_lb
        P_i = P_lb

        dt_plus = self.get_solve_state(P_lb)
        dt_neg = dt_plus
        P = np.array(P_lb)
        DT = np.array(dt_plus)

        i = 0
        m_step = 0.75
        if dt_neg > 0.0:
            #Find a lower bound pressure resulting in an air temperature less than the target
            while dt_neg > 0.0 and i < 100:
                dt_plus = dt_neg
                P_ub = P_lb
                P_lb = self.comp1.flow_in.pressure + (P_lb - self.comp1.flow_in.pressure) * m_step
                dt_neg = self.get_solve_state(P_lb)
                P = np.append(P, P_lb)
                DT = np.append(DT, dt_neg)
                i += 1

        else:
            #Find a upper bound pressure resulting in an air temperature greater than the target
            while dt_plus < 0.0 and i < 100:
                dt_neg = dt_plus
                P_lb = P_ub
                P_ub = self.comp1.flow_in.pressure + (P_ub - self.comp1.flow_in.pressure) / m_step
                dt_plus = self.get_solve_state(P_ub)
                P = np.append(P, P_ub)
                DT = np.append(DT, dt_plus)
                i += 1

        if P_lb < self.comp1.flow_in.pressure:
            #Check the low pressure bound remains above the evaporator pressure
            P_lb = self.comp1.flow_in.pressure
            dt_neg = self.get_solve_state(P_lb)
            if P_ub < P_lb:
                #Check the upper pressure bound stays above the low pressure bound
                dt_neg = dt_plus
        
        if abs(dt_plus * dt_neg) == -dt_plus * dt_neg and abs(dt_plus * dt_neg) > 0.0:
            feasible = 1
            P_i = P_lb - dt_neg / (dt_plus - dt_neg) * (P_ub - P_lb)
            if P_i * 1.01 < P_ub:
                dt = self.get_solve_state(P_i * 1.01)
                if dt > 0.0:
                    P_ub = P_i * 1.01
                    dt_plus = dt
            if P_i * 0.99 > P_lb:
                dt = self.get_solve_state(P_i * 0.99)
                if dt < 0.0:
                    P_lb = max(P_i * 0.99, P_lb)
                    dt_neg = dt
            
            P_i = P_lb - dt_neg / (dt_plus - dt_neg) * (P_ub - P_lb)

        return P_lb, P_ub, P_i, feasible

    def get_solve_state(self, P_comp_out):
        self.comp3.flow_out.pressure = P_comp_out
        i = 0; i_max = 10
        h = 0.0; m = 0.0; P = 0.0; W = 0.0
        dh = 0; dm = 0; dP = 0; dW = 0

        while ((i == 0 or min(dP, dW, dQ, dh, dm) < 0.0) and i < i_max): #Check for model convergance: dh, dm, 
            # Saves pressure and temperature properties to check for solution convergance
            h = self.exp_valve.flow_in.enthalpy
            m = self.gas_cooler4.sink_out.mass_flow
            P = self.gas_cooler4.source_in.pressure - self.gas_cooler1.source_out.pressure
            W = self.comp1.power + self.comp2.power + self.comp3.power
            Q = self.gas_cooler4.duty + self.gas_cooler3.duty + self.gas_cooler2.duty + self.gas_cooler1.duty
            i += 1

            #Solve the performance of each component and update each special relationship
            for obj in self.operations:
                obj.solve()

            dh = self.tol * h - abs(self.exp_valve.flow_in.enthalpy - h)
            dm = self.tol * m - abs(self.gas_cooler4.sink_out.mass_flow - m)
            dP = self.tol * P - abs(self.gas_cooler4.source_in.pressure - self.gas_cooler1.source_out.pressure - P)
            dW = self.tol * W - abs(self.comp1.power + self.comp2.power + self.comp3.power - W)
            dQ = self.tol * Q - abs(self.gas_cooler4.duty + self.gas_cooler3.duty + self.gas_cooler2.duty + self.gas_cooler1.duty - Q)

        self.set_gc_air_massflow.solve()

        if 0:
            self.get_cycle_plot(x_axis='h', y_axis='T')
            None

        if i < i_max:
            self.converged = True
            e = self.gas_cooler4.sink_out.temperature - self.T_air_tar
        else:
            self.converged = False
            e = -self.T_air_tar

        return e

    def get_solve_gas_pressure(self, method, bnds, x0):
        r = root_scalar(f=self.get_solve_state, 
                        method=method, 
                        bracket=bnds, 
                        x0=x0, 
                        xtol=self.tol,
                        )
        r['error'] = 1e-4 - abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature) / self.T_air_tar
        if self.disp:
            # self.get_cycle_plot()
            print(r['error'])
            
        return r

    def set_quality_relationship(self):
        e = []
        gcs = [self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4]
        for g in gcs:
            T_sat = self.RP_ref.get_properties('T', 'PQ', g.source_out.pressure, 0.0)[0]
            e.append(abs(T_sat - g.source_out.temperature))
        
        self.q_gc_out = [1.0]

        n = e.index(min(e))
        new_rel = SetRelationship( source_obj=self.q_gc_out, target_obj=[gcs[n].source_out], prop='quality', prop2='pressure' )        
        if n == 1:
            old_rel = self.set_gc2_hot_out
            self.set_gc2_hot_out = new_rel
            self.phx_fixed_src_quality = self.phx1
        elif n == 2:
            old_rel = self.set_gc3_hot_out
            self.set_gc3_hot_out = new_rel
            self.phx_fixed_src_quality = self.phx2
        elif n == 3:
            old_rel = self.set_gc4_hot_out
            self.set_gc4_hot_out = new_rel
            self.phx_fixed_src_quality = self.phx3
        else:
            ValueError('Error in the model.')

        for i in range(len(self.operations)):
            if self.operations[i] == old_rel:
                self.operations[i] = new_rel

        gcs[n].model = 'fixed hot quality outlet'
        return gcs[n]

    def get_phx_min_dT_for_solve_quality_method(self, QMASS_gc_out):
        self.q_gc_out[0] = QMASS_gc_out
        r = self.get_solve_gas_pressure(method='secant', 
                                        x0=self.P_i, 
                                        bnds=(self.P_lb, self.P_ub))
        hx = self.phx_fixed_src_quality
        e = hx.source_in.temperature - hx.sink_out.temperature - self.dTmin
        return e
    
    def get_solve_quality(self):
        e_one = self.get_phx_min_dT_for_solve_quality_method(1.0)
        P_one = self.comp3.flow_out.pressure
        e_zero = self.get_phx_min_dT_for_solve_quality_method(0.0)
        P_zero = self.comp3.flow_out.pressure
        if e_zero * e_one < 0:
            r = root_scalar(f=self.get_phx_min_dT_for_solve_quality_method,
                            method='brentq',
                            bracket=(0.0, 1.0),
                            x0=0.5,
                            xtol=1e-3
                            )
            r['error'] = 1e-4 - abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature)
        else:
            if abs(P_zero - P_one) > 1:
                None
            r = {'error' : 1.0}
        return r

    def get_p_test(self):
        P = np.linspace(start=self.P_lb / 1.5, stop=self.P_ub * 3, num=1000)
        dT = np.array([self.T_air_tar + self.get_solve_state(P_i) for P_i in P])
        Graph().simple(x=P, y=dT)
        ddT = np.diff(dT)
        i = ddT.argmax()
        P0, P1 = P[i], P[i+1]
        # P0, P1 = 16.181339509277915, 16.221350320439466
        # self.get_solve_state(P0)
        # print(self.set_gc_air_massflow.solve())
        # self.get_cycle_plot()
        # self.get_solve_state(P1)
        # print(self.set_gc_air_massflow.solve())
        # self.get_cycle_plot()
        None

    def get_solve_T_evap(self, T_evap_sat_out):
        #Solve heat pump cycle
        self.evaporator.sink_out.temperature = T_evap_sat_out
        self.evaporator.sink_out.props_in_mTX()
        self.converged = False
        self.COP = -np.inf
        
        if (self.T_air_tar - self.T_air_inlet) * 2/3 + 2 >= self.dT_phx1_cold + self.dT_phx2_cold + self.dT_phx3_cold:
            self.P_lb, self.P_ub, self.P_i, self.feasible = self.get_pressure_bounds()

            if self.feasible == 1:
                # self.get_p_test()
                try:
                    r = self.get_solve_gas_pressure(method='secant', x0=self.P_i, bnds=(self.P_lb, self.P_ub))
                except:
                    self.print_failed_to_solve_report()
                    r = {'error' : -1,
                         'root' : self.P_i}

                # if r['error'] < 0.0:
                #     # Try bisect method instead
                #     r = self.get_solve_gas_pressure(method='bisect', x0=r['root'], bnds=(self.P_lb, self.P_ub))
                #     if self.disp:
                #         self.get_cycle_plot()
                #     if r['error'] > 0.0:
                #         None

                if r['error'] < 0.0:
                    # Possibly unconverged due to a fixed hot outlet temperature of a gas cooler 
                    # being near the saturation temperature, causing instability
                    if r['root'] < self.RP_ref.pcrit:
                        self.set_quality_relationship()
                        self.P_i = r['root']
                        r = self.get_solve_quality()
                        if self.disp:
                            self.get_cycle_plot()
                    else:
                        self.P_lb = self.RP_ref.pcrit + 10e-4
                        dt_neg = self.get_solve_state(self.P_lb)
                        
                        if dt_neg < 0.0:
                            r = self.get_solve_gas_pressure(method='brentq', x0=r['root'], bnds=(self.P_lb, self.P_ub))
                        if self.disp:
                            self.get_cycle_plot()

                    if r['error'] < 0.0:
                        self.print_failed_to_solve_report()
                             
                if r['error'] > 0.0 and self.converged:
                    self.COP = (self.gas_cooler4.duty + self.gas_cooler3.duty + self.gas_cooler2.duty + self.gas_cooler1.duty) / (self.comp1.power + self.comp2.power + self.comp3.power)                   
                else:
                    self.converged = False
        
                self.cycle_type = 'Transcritical' if self.comp3.flow_out.pressure > self.RP_ref.pcrit else 'Subcritical'

        return self.evaporator.sink_in.temperature - self.Tsat_evap_tar

    def print_failed_to_solve_report(self, filename='results/Cases-Unsolved'):
        print('Heat pump model failed to solve.')
        self.get_summary()
        print_result_to_csv(self.output, file_name=filename)        

    def solve(self):
        #Solve heat pump cycle
        self.converged = False
        self.COP = -np.inf
        # try:
        r = self.get_solve_T_evap(self.Tsat_evap_tar)
        if abs(r) > 1e-6 and abs(self.evaporator.sink_in.temperature - self.evaporator.sink_out.temperature) > 0.1:

            r = root_scalar(
                            f=self.get_solve_T_evap,
                            method='secant',
                            x0=self.Tsat_evap_tar,
                            bracket=(self.Tsat_evap_tar, self.T_source_in),
                            xtol=self.tol
                        )
        # except:
        #     self.print_failed_to_solve_report()
                            
        return self.COP, self.comp3.flow_out.pressure

    def print_currrent_case(self):
        print({
            'refrigerant' : self.RP_ref_str,
            'dT_GC_min_tar' : self.dT_GC_min_tar,
            'T_air_tar' : self.T_air_tar,
            'T_air_inlet' : self.T_air_inlet,
            'T_source_in' : self.T_source_in,
            'n_comp' : self.n_comp,
            'dT_phx1_cold' : self.dT_phx1_cold,
            'dT_phx2_cold' : self.dT_phx2_cold,
            'dT_phx3_cold' : self.dT_phx3_cold,
            'comp2_power_frac' : self.comp2_power_frac,
            'comp3_power_frac' : self.comp3_power_frac,
            'comp_PR_uniform' : self.comp_PR_uniform
            })

    def get_summary(self):
        self.output = {
            #Parameters
            'Refrigerant:' : self.refrigerant,
            'Type:' : self.cycle_type,
            'n_comp:' : self.n_comp,
            'T_source_in:' : self.evaporator.source_in.temperature,
            'dT_min Gas Cooler:' : self.dT_GC_min_tar,
            'T_air_inlet:' : self.gas_cooler1.sink_in.temperature,
            'T_air_tar:' : self.T_air_tar,
            'Comp_PR_uniform:' : self.comp_PR_uniform,
            #Key results
            'COP:' : self.COP,
            'Tol:' : abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature) / self.T_air_tar,
            'Converged:' : self.converged,
            #Additional results
            'Actual air T:' : self.gas_cooler4.sink_out.temperature,            
            'Evap T:' : self.Tsat_evap,
            'Evap P:' : self.evaporator.sink_in.pressure * 1000,
            'Comp 1 P out:' : self.comp1.flow_out.pressure * 1000,
            'Comp 2 P out:' : self.comp2.flow_out.pressure * 1000,
            'Comp 3 P out:' : self.comp3.flow_out.pressure * 1000,
            'Comp 1 power frac:' : self.comp1_power_frac, 
            'Comp 2 power frac:' : self.comp2_power_frac, 
            'Comp 3 power frac:' : self.comp3_power_frac,            
            'dT phx max 1:' : self.dT_phx1_cold,
            'dT phx max 2:' : self.dT_phx2_cold,
            'dT phx max 3:' : self.dT_phx3_cold,
            'PR 1:' : self.comp1.pressure_ratio,
            'PR 2:' : self.comp2.pressure_ratio,
            'PR 3:' : self.comp3.pressure_ratio,
            'Overall P ratio:' : self.gas_cooler4.source_in.pressure / self.evaporator.sink_in.pressure,
            'Gas cooler(s) inlet T:' : self.gas_cooler4.source_in.temperature,
            'Gas cooler(s) outlet T:' : self.gas_cooler1.source_out.temperature
        }

    def get_cycle_plot(self, x_axis='h', y_axis='T'):
        g = Graph()
        g.cycle(cycle_obj=[
            [self.ihx, 'source'],
            [self.exp_valve, 'siso'],
            [self.evaporator, 'sink'],
            [self.evap_sup, 'sink'],
            [self.ihx, 'sink'],
            [self.phx1, 'sink'],
            [self.comp1, 'siso'],
            [self.phx2, 'sink'],
            [self.comp2, 'siso'],
            [self.phx3, 'sink'],
            [self.comp3, 'siso'],
            [self.gas_cooler4, 'source'],
            [self.phx3, 'source'],
            [self.gas_cooler3, 'source'],
            [self.phx2, 'source'],
            [self.gas_cooler2, 'source'],
            [self.phx1, 'source'],
            [self.gas_cooler1, 'source'], 
            ], 
            sink_obj=[
                [self.gas_cooler4, 'sink'],
                [self.gas_cooler3, 'sink'],
                [self.gas_cooler2, 'sink'],
                [self.gas_cooler1, 'sink'],
            ],
            x_axis=x_axis, y_axis=y_axis)
