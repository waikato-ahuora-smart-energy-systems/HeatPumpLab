__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes.components.compressor import Compressor
from classes.components.plateHeatExc import PlateHeatExchanger
from classes.components.gasCooler import GasCooler
from classes.components.expansionValve import ExpansionValve
from classes.components.evaporator import Evaporator
from classes.streams.massStream import MassStream
from classes.logic_operations.setRelationship import SetRelationship
from classes.logic_operations.customRelationship import CustomRelationship
from classes.streams.properties import REFPROP
from scipy.optimize import brentq

class TRHP:
    def __init__(
            self, 
            Refrigerant='R1336mzz(Z)', 
            Source='Water',
            Sink='Air',
            T_air_tar=200.0, 
            T_air_inlet=15.0, 
            dT_GC_min_tar=20.0, 
            T_gc_out_min=0.0, 
            T_source_in=45.0, 
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
            tol=1e-6
        ):

        self.Refrigerant = Refrigerant
        self.Source = Source
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
        self.dT_phx1_cold = dT_phx1_cold
        self.dT_phx2_cold = dT_phx2_cold 
        self.dT_phx3_cold = dT_phx3_cold 
        self.tol = tol
        self.converged = False
        self.COP = 0.0
        
        self.init_model()


    def init_model(self):
        self.check_logic()
        self.init_components()
        self.init_streams()
        self.set_relationships()
        self.define_solver_pathway()


    def check_logic(self):
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


    def init_components(self):
        #Initialisation of unit operations
        self.comp1 = Compressor(Motor_Eff=self.eff_motor, model=0)
        self.comp2 = Compressor(Motor_Eff=self.eff_motor, model=0) 
        self.comp3 = Compressor(Motor_Eff=self.eff_motor, model=0) 
        self.gas_cooler1 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.gas_cooler2 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.gas_cooler4 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.gas_cooler3 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.ihx = PlateHeatExchanger(dTmin=self.dTmin, model='dTmin_terminal')
        self.phx1 = PlateHeatExchanger(dT_c=self.dT_phx1_cold, model='fixed cold delta temperature')
        self.phx2 = PlateHeatExchanger(dT_c=self.dT_phx2_cold, model='fixed cold delta temperature')        
        self.phx3 = PlateHeatExchanger(dT_c=self.dT_phx3_cold, model='fixed cold delta temperature')
        self.exp_valve = ExpansionValve()
        self.evaporator = Evaporator(model='fixed cold duty')
        self.evap_sup = PlateHeatExchanger(model='fixed cold temperature outlet')


    def init_streams(self):

        #Evaporator superheat section
        self.evap_sup.HotIn = MassStream( Name='Source In', Fluid=self.Source, Pressure=101, Temperature=self.T_source_in, Mass_Flow=10^6)
        self.evap_sup.HotOut = MassStream( Name='Source Out', Fluid=self.Source, Pressure=101, Temperature=self.T_source_in, Mass_Flow=10^6)
        self.evap_sup.ColdIn = MassStream( Name='LP Ref: Evap to Super', Fluid=self.Refrigerant, Temperature=self.Tsat_evap, Quality=1.0 )
        P_L = self.evap_sup.ColdIn.Pressure
        self.evap_sup.ColdOut = MassStream( Name='LP Ref: Super to IHX', Fluid=self.Refrigerant, Pressure=P_L, Temperature=(self.T_source_in - self.dTmin) )
        
        #Internal heat exchanger
        self.ihx.ColdIn = self.evap_sup.ColdOut
        if self.ihx_on == 1:
            self.ihx.ColdOut = MassStream( Name='LP Ref: IHX to PHX1', Fluid=self.Refrigerant, Pressure=P_L, Temperature=(self.T_gc_out - self.dTmin) )
        else:   
            self.ihx.ColdOut = self.ihx.ColdIn

        #Plate heat exchanger 1 - cold-side
        self.phx1.ColdIn = self.ihx.ColdOut
        if self.phx1_on == 1: 
            self.phx1.ColdOut = MassStream( Name='LP Ref: PHX1 to Comp1', Fluid=self.Refrigerant, Pressure=P_L, Temperature=(self.phx1.ColdIn.Temperature + self.dT_phx1_cold) )
        else:
            self.phx1.ColdOut = self.phx1.ColdIn

        #Final discharge pressure
        P_H = REFPROP('P', self.comp1.FlowIn.Fluid, 'TS', self.T_air_tar + (self.dT_GC_min_tar - self.dT_phx1_cold + self.dT_phx2_cold + self.dT_phx3_cold) / 2, self.phx1.ColdOut.Entropy, 0) * 1000
        PR = (P_H / P_L)**(1 / self.n_comp) 

        #Compressor 1
        self.comp1.FlowIn = self.phx1.ColdOut
        self.comp1.FlowOut = MassStream( Name='Comp1_out', Fluid=self.Refrigerant, Pressure=(P_L * PR), Temperature=REFPROP('T', self.Refrigerant, 'PS', (P_L * PR), self.phx1.ColdOut.Entropy, 0))

        #Plate heat exchanger 2 - cold-side
        self.phx2.ColdIn = self.comp1.FlowOut
        if self.phx2_on == 1:
            self.phx2.ColdOut = MassStream( Name='Comp2_in', Fluid=self.Refrigerant, Pressure=self.phx2.ColdIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_cold)
        else:
            self.phx2.ColdOut = self.comp1.FlowOut
        
        #Compressor 2
        self.comp2.FlowIn = self.phx2.ColdOut
        if self.comp2_on == 1:
            self.comp2.FlowOut = MassStream( Name='Comp2_out', Fluid=self.Refrigerant, Pressure=(P_L * PR * PR), Temperature=REFPROP('T', self.Refrigerant, 'PS', (P_L * PR * PR), self.phx1.ColdOut.Entropy, 0))       
        else:
            self.comp2.FlowOut = self.comp2.FlowIn

        #Plate heat exchanger 3 - cold-side
        self.phx3.ColdIn = self.comp2.FlowOut
        if self.phx3_on == 1:
            self.phx3.ColdOut = MassStream( Name='Comp3_in', Fluid=self.Refrigerant, Pressure=self.phx3.ColdIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_cold )
        else:
            self.phx3.ColdOut = self.phx3.ColdIn

        #Compressor 3
        self.comp3.FlowIn = self.phx3.ColdOut
        if self.comp3_on == 1:
            self.comp3.FlowOut = MassStream( Name='Ref: Comp3 to GC4', Fluid=self.Refrigerant, Pressure=P_H, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_H, self.phx1.ColdOut.Entropy, 0))       
        else:
            self.comp3.FlowOut = self.comp3.FlowIn

        #Gas cooler 4 and Plate heat exchanger 3 - hot-sides
        self.gas_cooler4.ColdOut = MassStream( Name='Air: GC4 to Dryer', Fluid=self.Sink, Pressure=101, Temperature=self.T_air_tar )
        self.gas_cooler4.ColdIn = MassStream( Name='Air: GC3 to GC4', Fluid=self.Sink, Pressure=101, Temperature=self.T_air_tar )
        
        self.gas_cooler4.HotIn = self.comp3.FlowOut
        if self.comp3_on == 1 and self.phx3_on == 1:
            self.gas_cooler4.HotOut = MassStream( Name='Ref: GC4 to PHX3', Fluid=self.Refrigerant )
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = MassStream( Name='Ref: PHX3 to GC3', Fluid=self.Refrigerant )       
        else:
            self.gas_cooler4.HotOut = self.gas_cooler4.HotIn
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = self.phx3.HotIn       
        
        #Gas cooler 3 and Plate heat exchanger 2 - hot-sides
        self.gas_cooler3.ColdOut = self.gas_cooler4.ColdIn
        self.gas_cooler3.ColdIn = MassStream( Name='Air: GC2 to GC3', Fluid=self.Sink, Pressure=101, Temperature=self.T_air_tar )
        
        self.gas_cooler3.HotIn = self.phx3.HotOut
        if self.comp2_on == 1 and self.phx2_on == 1:
            self.gas_cooler3.HotOut = MassStream( Name='Ref: GC3 to PHX2', Fluid=self.Refrigerant )
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = MassStream( Name='Ref: PHX2 to GC2', Fluid=self.Refrigerant )
        else:
            self.gas_cooler3.HotOut = self.gas_cooler3.HotIn
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = self.phx2.HotIn

        #Gas cooler 2 and Plate heat exchanger 1 - hot-sides
        self.gas_cooler2.ColdOut = self.gas_cooler3.ColdIn
        self.gas_cooler2.ColdIn = MassStream( Name='Air: GC1 to GC2', Fluid=self.Sink, Pressure=101, Temperature=self.T_air_tar )

        self.gas_cooler2.HotIn = self.phx2.HotOut
        if self.phx1_on == 1:
            self.gas_cooler2.HotOut = MassStream( Name='Ref: GC2 to PHX1', Fluid=self.Refrigerant )
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = MassStream( Name='Ref: PHX1 to GC1', Fluid=self.Refrigerant )
        else:
            self.gas_cooler2.HotOut = self.gas_cooler2.HotIn
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = self.phx1.HotIn

        #Gas cooler 1
        self.gas_cooler1.ColdIn = MassStream( Name='Air: Inlet to GC1', Fluid=self.Sink, Pressure=101, Temperature=self.T_air_inlet )
        self.gas_cooler1.ColdOut = self.gas_cooler2.ColdIn

        self.gas_cooler1.HotOut = MassStream( Name='Ref: GC1 to IHX', Fluid=self.Refrigerant, Pressure=P_H, Temperature=self.T_gc_out )
        self.gas_cooler1.HotIn = self.phx1.HotOut

        #Internal heat exchanger
        self.ihx.HotIn = self.gas_cooler1.HotOut
        if self.ihx_on == 1:
            self.ihx.HotOut = MassStream( Name='Ref: IHX to Exp V', Fluid=self.Refrigerant, Pressure=P_H, Temperature=self.T_gc_out )
        else:   
            self.ihx.HotOut = self.ihx.HotIn
        
        #Expansion valve
        self.exp_valve.FlowIn = self.ihx.HotOut
        self.exp_valve.FlowOut = MassStream( Name='Evap_in', Fluid=self.Refrigerant, Pressure=P_L, Enthalpy=self.exp_valve.FlowIn.Enthalpy )

        #Evaporator
        self.evaporator.ColdIn = self.exp_valve.FlowOut
        self.evaporator.ColdOut = self.evap_sup.ColdIn
        self.evaporator.HotIn = MassStream( Name='Exhaust air in', Fluid=self.Source, Pressure=101, Temperature=45.0, Mass_Flow=1000)
        self.evaporator.HotOut = MassStream( Name='Exhaust air out', Fluid=self.Source, Pressure=101, Temperature=45.0, Mass_Flow=1000)


    def set_relationships(self):
        #Set special relationships between points in the cycle 
        if self.phx3_on == 1:
            self.set_phx3_dT_cold = SetRelationship( source_obj=self.phx3.ColdIn, target_obj=[self.phx3.ColdOut], c=self.dT_phx3_cold, prop='Temperature', prop2='Pressure' ) 
            self.set_phx3_hot_T_in = SetRelationship( source_obj=self.phx3.ColdOut, target_obj=[self.phx3.HotIn], c=self.dTmin, prop='Temperature', prop2='Pressure' ) 

        if self.phx2_on == 1:    
            self.set_phx2_dT_cold = SetRelationship( source_obj=self.phx2.ColdIn, target_obj=[self.phx2.ColdOut], c=self.dT_phx2_cold, prop='Temperature', prop2='Pressure' ) 
            self.set_phx2_hot_T_in = SetRelationship( source_obj=self.phx2.ColdOut, target_obj=[self.phx2.HotIn], c=self.dTmin, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx1_on == 1:
            self.set_phx1_dT_cold = SetRelationship( source_obj=self.phx1.ColdIn, target_obj=[self.phx1.ColdOut], c=self.dT_phx1_cold, prop='Temperature', prop2='Pressure' ) 
            self.set_phx1_hot_T_in = SetRelationship( source_obj=self.phx1.ColdOut, target_obj=[self.phx1.HotIn], c=self.dTmin, prop='Temperature', prop2='Pressure' ) 

        self.set_P_ints = CustomRelationship( source_obj=[self.comp1.FlowIn, self.comp3.FlowOut, self.comp2_power_frac, self.comp3_power_frac], target_obj=[self.comp1.FlowOut, self.comp2.FlowOut], relation=1 )
        self.set_gc_air_massflow = CustomRelationship( source_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], target_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], relation=2 )
        #self.set_evap_air_massflow = CustomRelationship( source_obj=[self.evaporator, self.evap_sup], target_obj=[self.evaporator, self.evap_sup], relation=3)


    def define_solver_pathway(self):
        #Solution pathway set-up
        self.OP = [self.evap_sup, self.set_P_ints]

        if self.ihx_on == 1:
            self.OP.append(self.ihx)

        if self.phx1_on == 1:
            self.OP.append(self.set_phx1_dT_cold)
            self.OP.append(self.set_phx1_hot_T_in)
        self.OP.append(self.comp1)
        if self.phx2_on == 1:
            self.OP.append(self.set_phx2_dT_cold)
            self.OP.append(self.set_phx2_hot_T_in)
        if self.comp2_on == 1:
            self.OP.append(self.comp2)
        if self.phx3_on == 1:
            self.OP.append(self.set_phx3_dT_cold)
            self.OP.append(self.set_phx3_hot_T_in)
        if self.comp3_on == 1:
            self.OP.append(self.comp3)
        
        if self.phx3_on == 1:
            self.OP.append(self.gas_cooler4)
            self.OP.append(self.phx3)
        if self.phx2_on == 1:
            self.OP.append(self.gas_cooler3)
            self.OP.append(self.phx2)
        if self.phx1_on == 1:
            self.OP.append(self.gas_cooler2)        
            self.OP.append(self.phx1)
        self.OP.append(self.gas_cooler1)

        self.OP.append(self.exp_valve)
        self.OP.append(self.evaporator)
        self.OP.append(self.set_gc_air_massflow)


    def solve(self):
        #Solve heat pump cycle                        
        if self.set_pressure_bounds() > 0.0:
            #Solve the heat pump thermodynamics
            r = brentq(f=self.state_solve, a=self.P_low_bound, b=self.P_upp_bound, full_output=True, rtol=self.tol, xtol=self.tol)

            if r[1].converged and self.converged and abs(self.gas_cooler4.ColdOut.Temperature - self.T_air_tar) < 0.05:
                self.COP = (self.gas_cooler4.Duty + self.gas_cooler3.Duty + self.gas_cooler2.Duty + self.gas_cooler1.Duty) / (self.comp1.Power + self.comp2.Power + self.comp3.Power)
                        
        if self.COP == 0:
            self.converged = False          
        
        self.get_summary()


    def set_pressure_bounds(self):
        P_gc_est = REFPROP('P', self.comp1.FlowIn.Fluid, 'TS', self.T_air_tar + (self.dT_GC_min_tar - self.dT_phx1_cold + self.dT_phx2_cold + self.dT_phx3_cold) / 2, self.comp1.FlowIn.Entropy, 0) * 1000

        self.P_low_bound = P_gc_est
        self.P_upp_bound = P_gc_est

        dt_plus = self.state_solve(P_gc_est)
        dt_neg = dt_plus
        i = 0
        if dt_neg > 0.0:
            #Find a lower bound pressure resulting in an air temperature less than the target
            while dt_neg > 0.0 and i < 20:
                self.P_upp_bound = self.P_low_bound
                self.P_low_bound *= 0.9
                dt_neg = self.state_solve(self.P_low_bound)
                i += 1
            
        else:
            #Find a upper bound pressure resulting in an air temperature greater than the target
            while dt_plus < 0.0 and i < 20:
                self.P_low_bound = self.P_upp_bound
                self.P_upp_bound *= 1.1
                dt_plus = self.state_solve(self.P_upp_bound)
                i += 1

        if self.P_low_bound < self.comp1.FlowIn.Pressure:
            #Check the low pressure bound remains above the evaporator pressure
            self.P_low_bound = self.comp1.FlowIn.Pressure
            dt_neg = self.state_solve(self.P_low_bound)
            if self.P_upp_bound < self.P_low_bound:
                #Check the upper pressure bound stays above the low pressure bound
                dt_neg = dt_plus

        return -dt_plus * dt_neg    
   

    def state_solve(self, P_comp_out):
        self.converged = False
        self.comp3.FlowOut.Pressure = P_comp_out
        i = 0
        h = 0.0; m = 0.0; P = 0.0; W = 0.0
        dh = 0; dm = 0; dP = 0; dW = 0

        while ((i < 1 or min(dh, dm, dP, dW) < 0.0) and i < 20): #Check for model convergance
            #Saves critical pressure and temperature properties to check for solution convergance
            h = self.exp_valve.FlowIn.Enthalpy
            m = self.gas_cooler4.ColdOut.Mass_Flow
            P = self.gas_cooler4.HotIn.Pressure - self.gas_cooler1.HotOut.Pressure
            W = self.comp1.Power + self.comp2.Power + self.comp3.Power
            i += 1

            #Solve the performance of each component and update each special relationship
            for obj in self.OP:
                obj.solve()

            dh = self.tol * h - abs(self.exp_valve.FlowIn.Enthalpy - h)
            dm = self.tol * m - abs(self.gas_cooler4.ColdOut.Mass_Flow - m)
            dP = self.tol * P - abs(self.gas_cooler4.HotIn.Pressure - self.gas_cooler1.HotOut.Pressure - P)
            dW = self.tol * W - abs(self.comp1.Power + self.comp2.Power + self.comp3.Power - W)

        self.converged = True if i < 20 else False
            
        return self.gas_cooler4.ColdOut.Temperature - self.T_air_tar


    def get_summary(self):
        self.output = {
            #Key results
            'COP:' : self.COP,
            'Tol:' : abs(self.T_air_tar - self.gas_cooler4.ColdOut.Temperature) / self.T_air_tar,
            'Converged:' : self.converged,
            #Parameters
            'Refrigerant:' : self.Refrigerant,
            'n_comp:' : self.n_comp,
            'Dryer exh air T:' : self.evaporator.HotIn.Temperature,
            'dT_min Gas Cooler:' : self.dT_GC_min_tar,
            'T_air_inlet:' : self.gas_cooler1.ColdIn.Temperature,
            'T_air_tar:' : self.T_air_tar,
            #Additional results
            'Actual air T:' : self.gas_cooler4.ColdOut.Temperature,            
            'Evap T:' : self.Tsat_evap,
            'Evap P:' : self.evaporator.ColdIn.Pressure,
            'Comp 1 P out:' : self.comp1.FlowOut.Pressure,
            'Comp 2 P out:' : self.comp2.FlowOut.Pressure,
            'Comp 3 P out:' : self.comp3.FlowOut.Pressure,
            'dT phx max 1:' : self.dT_phx1_cold,
            'dT phx max 2:' : self.dT_phx2_cold,
            'dT phx max 3:' : self.dT_phx3_cold,
            'PR 1:' : self.comp1.Pressure_Ratio,
            'PR 2:' : self.comp2.Pressure_Ratio,
            'PR 3:' : self.comp3.Pressure_Ratio,
            'Overall P ratio:' : self.gas_cooler4.HotIn.Pressure / self.evaporator.ColdIn.Pressure,
            'Gas cooler(s) inlet T:' : self.gas_cooler4.HotIn.Temperature,
            'Gas cooler(s) outlet T:' : self.gas_cooler1.HotOut.Temperature
        }
    