__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2023'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021 (tim.walmsley@waikato.ac.nz)
"""

from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator
from streams import MassStream, Properties
from streams.properties import get_property_value
from logic_operations import SetRelationship, CustomRelationship, FlowSheetSolver
from scipy.optimize import root_scalar 
from scipy.interpolate import CubicSpline
from outputs import Graph, print_result_to_csv
import numpy as np

class TRHP(FlowSheetSolver, Properties):
    def __init__(
            self, 
            refrigerant='R1336mzz(Z)', 
            source='Water',
            Sink='Air',
            T_air_tar=200.0, 
            T_air_inlet=15.0, 
            dT_GC_min_tar=20.0, 
            T_gc_out_min=0.0, 
            T_source_in=40.0, 
            Tsat_evap=10.0,
            dTmin_phx=5.0,
            n_comp=3, 
            comp_PR_uniform=False, 
            eff_motor=0.96, 
            dT_phx1_cold=0.0, 
            dT_phx2_cold=0.0, 
            dT_phx3_cold=0.0, 
            comp2_power_frac=0.0, 
            comp3_power_frac=0.0,
            P_H_est=0.0,
            tol=1e-9,
            result_parameters=None
        ):
        
        self.refrigerant = refrigerant
        self.cycle_type = None
        self.source = source
        self.Sink = Sink
        self.T_air_tar = T_air_tar
        self.dT_GC_min_tar = dT_GC_min_tar
        self.T_air_inlet = T_air_inlet
        self.T_source_in = T_source_in
        self.T_gc_out = max(self.T_air_inlet + self.dT_GC_min_tar, T_gc_out_min)
        self.Tsat_evap = Tsat_evap
        self.n_comp = n_comp
        self.comp2_power_frac = comp2_power_frac
        self.comp3_power_frac = comp3_power_frac
        self.comp_PR_uniform = comp_PR_uniform
        self.dTmin = dTmin_phx
        self.eff_motor = eff_motor
        self.min_dT_ihx = 5.0
        self.dT_phx1_cold = dT_phx1_cold #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.dT_phx2_cold = dT_phx2_cold #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.dT_phx3_cold = dT_phx3_cold #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.P_H_est = P_H_est
        self.tol = tol
        self.converged = False
        self.COP = -np.inf
        self.disp = False
        
        self.init_properties()
        self.init_components()
        self.init_structure()
        self.init_streams()
        self.init_logic_blocks()
        self.init_solver_pathway()

    def init_components(self):
        #Initialisation of unit operations
        self.comp1 = Compressor(motor_eff=self.eff_motor, model_index=1, RP=self.RP)
        self.comp2 = Compressor(motor_eff=self.eff_motor, model_index=1, RP=self.RP) 
        self.comp3 = Compressor(motor_eff=self.eff_motor, model_index=1, RP=self.RP) 
        self.gas_cooler1 = GasCooler(dTmin=self.dT_GC_min_tar, RP=self.RP)
        self.gas_cooler2 = GasCooler(dTmin=self.dT_GC_min_tar, RP=self.RP)
        self.gas_cooler3 = GasCooler(dTmin=self.dT_GC_min_tar, RP=self.RP)
        self.gas_cooler4 = GasCooler(dTmin=self.dT_GC_min_tar, RP=self.RP)
        self.ihx = PlateHeatExchanger(dTmin=self.dTmin, model='dTmin_terminal', RP=self.RP)
        self.phx1 = PlateHeatExchanger(dT_c=self.dT_phx1_cold, model='fixed cold delta temperature', RP=self.RP)
        self.phx2 = PlateHeatExchanger(dT_c=self.dT_phx2_cold, model='fixed cold delta temperature', RP=self.RP)        
        self.phx3 = PlateHeatExchanger(dT_c=self.dT_phx3_cold, model='fixed cold delta temperature', RP=self.RP)
        self.exp_valve = ExpansionValve(RP=self.RP)
        self.evaporator = Evaporator(model='saturated gas outlet evaporator', RP=self.RP)
        self.evap_sup = PlateHeatExchanger(model='fixed cold temperature outlet', RP=self.RP)

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
        self.evap_sup.source_in = MassStream( name='source In', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6, RP=self.RP)
        self.evap_sup.source_out = MassStream( name='source Out', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6, RP=self.RP)
        self.evap_sup.sink_in = MassStream( name='LP Ref: Evap to Super', fluid=self.refrigerant, temperature=self.Tsat_evap, quality=1.0, RP=self.RP)
        P_L = self.evap_sup.sink_in.pressure
        self.evap_sup.sink_out = MassStream( name='LP Ref: Super to IHX', fluid=self.refrigerant, pressure=P_L, temperature=(self.T_source_in - self.dTmin), RP=self.RP)
        
        #Internal heat exchanger
        self.ihx.sink_in = self.evap_sup.sink_out
        if self.ihx_on == 1:
            self.ihx.sink_out = MassStream( name='LP Ref: IHX to PHX1', fluid=self.refrigerant, pressure=P_L, temperature=(self.T_gc_out - self.dTmin), RP=self.RP)
        else:   
            self.ihx.sink_out = self.ihx.sink_in

        #Plate heat exchanger 1 - cold-side
        self.phx1.sink_in = self.ihx.sink_out
        if self.phx1_on == 1: 
            self.phx1.sink_out = MassStream( name='LP Ref: PHX1 to Comp1', fluid=self.refrigerant, pressure=P_L, temperature=(self.phx1.sink_in.temperature + self.dT_phx1_cold), RP=self.RP)
        else:
            self.phx1.sink_out = self.phx1.sink_in

        #Final discharge pressure
        if self.P_H_est > 0.0:
            P_H = self.P_H_est
        else:
            P_H = get_property_value('P', self.refrigerant, 'TS', self.T_air_tar, self.phx1.sink_out.entropy, RP=self.RP)
            
        PR = (P_H / P_L) ** (1 / self.n_comp)
        
        #Compressor 1
        self.comp1.flow_in = self.phx1.sink_out
        self.comp1.flow_out = MassStream( name='Comp1_out', fluid=self.refrigerant, pressure=(P_L * PR), temperature=get_property_value('T', self.refrigerant, 'PS', (P_L * PR), self.phx1.sink_out.entropy, RP=self.RP), RP=self.RP)

        #Plate heat exchanger 2 - cold-side
        self.phx2.sink_in = self.comp1.flow_out
        if self.phx2_on == 1:
            self.phx2.sink_out = MassStream( name='Comp2_in', fluid=self.refrigerant, pressure=self.phx2.sink_in.pressure, temperature=self.phx2.sink_in.temperature + self.dT_phx2_cold, RP=self.RP)
        else:
            self.phx2.sink_out = self.comp1.flow_out
        
        #Compressor 2
        self.comp2.flow_in = self.phx2.sink_out
        if self.comp2_on == 1:
            self.comp2.flow_out = MassStream( name='Comp2_out', fluid=self.refrigerant, pressure=(P_L * PR * PR), temperature=get_property_value('T', self.refrigerant, 'PS', (P_L * PR * PR), self.phx1.sink_out.entropy, RP=self.RP), RP=self.RP)       
        else:
            self.comp2.flow_out = self.comp2.flow_in

        #Plate heat exchanger 3 - cold-side
        self.phx3.sink_in = self.comp2.flow_out
        if self.phx3_on == 1:
            self.phx3.sink_out = MassStream( name='Comp3_in', fluid=self.refrigerant, pressure=self.phx3.sink_in.pressure, temperature=self.phx3.sink_in.temperature + self.dT_phx3_cold, RP=self.RP)
        else:
            self.phx3.sink_out = self.phx3.sink_in

        #Compressor 3
        self.comp3.flow_in = self.phx3.sink_out
        if self.comp3_on == 1:
            self.comp3.flow_out = MassStream( name='Ref: Comp3 to GC4', fluid=self.refrigerant, pressure=(P_L * PR * PR * PR), temperature=get_property_value('T', self.refrigerant, 'PS', P_H, self.phx1.sink_out.entropy, RP=self.RP), RP=self.RP)       
        else:
            self.comp3.flow_out = self.comp3.flow_in

        #Gas cooler 4 and Plate heat exchanger 3 - hot-sides
        self.gas_cooler4.sink_out = MassStream( name='Air: GC4 to Dryer', fluid=self.Sink, pressure=101, temperature=self.T_air_tar, RP=self.RP )
        self.gas_cooler4.sink_in = MassStream( name='Air: GC3 to GC4', fluid=self.Sink, pressure=101, temperature=self.T_air_tar, RP=self.RP )
        
        self.gas_cooler4.source_in = self.comp3.flow_out
        if self.comp3_on == 1 and self.phx3_on == 1:
            self.gas_cooler4.source_out = MassStream( name='Ref: GC4 to PHX3', fluid=self.refrigerant, RP=self.RP )
            self.phx3.source_in = self.gas_cooler4.source_out
            self.phx3.source_out = MassStream( name='Ref: PHX3 to GC3', fluid=self.refrigerant, RP=self.RP )       
        else:
            self.gas_cooler4.source_out = self.gas_cooler4.source_in
            self.phx3.source_in = self.gas_cooler4.source_out
            self.phx3.source_out = self.phx3.source_in       
        
        #Gas cooler 3 and Plate heat exchanger 2 - hot-sides
        self.gas_cooler3.sink_out = self.gas_cooler4.sink_in
        self.gas_cooler3.sink_in = MassStream( name='Air: GC2 to GC3', fluid=self.Sink, pressure=101, temperature=self.T_air_tar, RP=self.RP )
        
        self.gas_cooler3.source_in = self.phx3.source_out
        if self.comp2_on == 1 and self.phx2_on == 1:
            self.gas_cooler3.source_out = MassStream( name='Ref: GC3 to PHX2', fluid=self.refrigerant, RP=self.RP )
            self.phx2.source_in = self.gas_cooler3.source_out
            self.phx2.source_out = MassStream( name='Ref: PHX2 to GC2', fluid=self.refrigerant, RP=self.RP )
        else:
            self.gas_cooler3.source_out = self.gas_cooler3.source_in
            self.phx2.source_in = self.gas_cooler3.source_out
            self.phx2.source_out = self.phx2.source_in

        #Gas cooler 2 and Plate heat exchanger 1 - hot-sides
        self.gas_cooler2.sink_out = self.gas_cooler3.sink_in
        self.gas_cooler2.sink_in = MassStream( name='Air: GC1 to GC2', fluid=self.Sink, pressure=101, temperature=self.T_air_tar, RP=self.RP )

        self.gas_cooler2.source_in = self.phx2.source_out
        if self.phx1_on == 1:
            self.gas_cooler2.source_out = MassStream( name='Ref: GC2 to PHX1', fluid=self.refrigerant, RP=self.RP )
            self.phx1.source_in = self.gas_cooler2.source_out
            self.phx1.source_out = MassStream( name='Ref: PHX1 to GC1', fluid=self.refrigerant, RP=self.RP )
        else:
            self.gas_cooler2.source_out = self.gas_cooler2.source_in
            self.phx1.source_in = self.gas_cooler2.source_out
            self.phx1.source_out = self.phx1.source_in

        #Gas cooler 1
        self.gas_cooler1.sink_in = MassStream( name='Air: Inlet to GC1', fluid=self.Sink, pressure=101, temperature=self.T_air_inlet, RP=self.RP )
        self.gas_cooler1.sink_out = self.gas_cooler2.sink_in

        self.gas_cooler1.source_out = MassStream( name='Ref: GC1 to IHX', fluid=self.refrigerant, pressure=P_H, temperature=self.T_gc_out, RP=self.RP )
        self.gas_cooler1.source_in = self.phx1.source_out

        #Internal heat exchanger
        self.ihx.source_in = self.gas_cooler1.source_out
        if self.ihx_on == 1:
            self.ihx.source_out = MassStream( name='Ref: IHX to Exp V', fluid=self.refrigerant, pressure=P_H, temperature=self.T_gc_out, RP=self.RP )
        else:   
            self.ihx.source_out = self.ihx.source_in
        
        #Expansion valve
        self.exp_valve.flow_in = self.ihx.source_out
        self.exp_valve.flow_out = MassStream( name='Evap_in', fluid=self.refrigerant, pressure=P_L, enthalpy=self.exp_valve.flow_in.enthalpy, RP=self.RP )

        #Evaporator
        self.evaporator.sink_in = self.exp_valve.flow_out
        self.evaporator.sink_out = self.evap_sup.sink_in
        self.evaporator.source_in = MassStream( name='Exhaust air in', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=1000, RP=self.RP)
        self.evaporator.source_out = MassStream( name='Exhaust air out', fluid=self.source, pressure=101, temperature=min(self.evaporator.sink_in.temperature + self.dTmin, self.T_source_in), mass_flow=1000, RP=self.RP)

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

        self.set_P_ints = CustomRelationship( source_obj=[self.comp1.flow_in, self.comp3.flow_out, self.comp2_power_frac, self.comp3_power_frac], target_obj=[self.comp1.flow_out, self.comp2.flow_out], relation=1, RP=self.RP, tol=self.tol )
        self.set_gc_air_massflow = CustomRelationship( source_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], target_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], relation=2, RP=self.RP, tol=self.tol )
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
        
        feasible = int((-dt_plus * dt_neg) / abs(dt_plus * dt_neg))
        
        if feasible == 1:
            try:
                cs = CubicSpline(y=P, x=DT) if DT[0] < DT[1] else CubicSpline(y=P[::-1], x=DT[::-1])
                P_i = float(cs(0.0))
            except:
                P_i = P_lb + dt_plus / (dt_plus - dt_neg) * (P_ub - P_lb)

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

        if self.disp:
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
                        rtol=self.tol,
                        )
        r['error'] = self.tol * 100 - abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature)
        return r

    def set_quality_relationship(self):
        e = []
        gcs = [self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4]
        for g in gcs:
            T_sat = get_property_value('T', self.refrigerant, 'PQ', g.source_out.pressure / 1000, 0.0, RP=self.RP)
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
                i0 = i

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
        r = root_scalar(f=self.get_phx_min_dT_for_solve_quality_method,
                        method='brentq',
                        bracket=(0.0, 1.0),
                        x0=0.5,
                        rtol=self.tol
                        )
        r['error'] = self.tol * 100 - abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature)
        return r

    def print_currrent_case(self):
        print({
            'refrigerant' : self.refrigerant,
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

    def solve(self):
        #Solve heat pump cycle
        self.converged = False
        self.COP = -np.inf
        
        if (self.T_air_tar - self.T_air_inlet) * 2/3 + 2 >= self.dT_phx1_cold + self.dT_phx2_cold + self.dT_phx3_cold:
            self.P_lb, self.P_ub, self.P_i, self.feasible = self.get_pressure_bounds()     
            if self.feasible == 1:
                r = self.get_solve_gas_pressure(method='secant', x0=self.P_i, bnds=(self.P_lb, self.P_ub))
                if r['error'] < 0.0:
                    # Try brentq method instead
                    r = self.get_solve_gas_pressure(method='brentq', x0=r['root'], bnds=(self.P_lb, self.P_ub))
                    if r['error'] < 0.0:
                        # Possibly unconverged due to a fixed hot outlet temperature of a gas cooler 
                        # being near the saturation temperature, causing instability
                        self.set_quality_relationship()
                        self.P_i = r['root']
                        r = self.get_solve_quality()
                        
                        if r['error'] < 0.0:
                            print('Heat pump model failed to solve.')
                            self.get_summary()
                            print_result_to_csv(self.output,file_name='results/Cases', suffix= '-Unsolved')
                            None
                
                if r['error'] > 0.0 and self.converged:
                    self.COP = (self.gas_cooler4.duty + self.gas_cooler3.duty + self.gas_cooler2.duty + self.gas_cooler1.duty) / (self.comp1.power + self.comp2.power + self.comp3.power)                   
                else:
                    self.converged = False
        
        P_crit = get_property_value('PCRIT', self.refrigerant, 'PT', 1, 100, RP=self.RP)
        self.cycle_type = 'Transcritical' if self.comp3.flow_out.pressure > P_crit else 'Subcritical'

        return self.COP, self.comp3.flow_out.pressure

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
            #Key results
            'COP:' : self.COP,
            'Tol:' : abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature) / self.T_air_tar,
            'Converged:' : self.converged,
            #Additional results
            'Actual air T:' : self.gas_cooler4.sink_out.temperature,            
            'Evap T:' : self.Tsat_evap,
            'Evap P:' : self.evaporator.sink_in.pressure,
            'Comp 1 P out:' : self.comp1.flow_out.pressure,
            'Comp 2 P out:' : self.comp2.flow_out.pressure,
            'Comp 3 P out:' : self.comp3.flow_out.pressure,
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
