__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes.unit_operations.compressor import Compressor
from classes.unit_operations.plateHeatExc import PlateHeatExc
from classes.unit_operations.gasCooler import GasCooler
from classes.unit_operations.expansionValve import Expansion_Valve
from classes.unit_operations.evaporator import Evaporator
from classes.streams.massStream import MassStream
from classes.logic_operations.setRelationship import SetRelationship
from classes.logic_operations.customRelationship import CustomRelationship
from classes.streams.properties import REFPROP
from scipy.optimize import brentq

class TRHP:
    def __init__(
            self, 
            Refrigerant='R1336mzz(Z)', 
            T_air_tar=200.0, 
            dT_GC_min_tar=20.0, 
            n_comp=3, 
            dT_phx1_cold=0.0, 
            dT_phx2_cold=0.0, 
            dT_phx3_cold=0.0, 
            comp2_power_frac=0.0, 
            comp3_power_frac=0.0,
            comp_PR_uniform=False, 
            T_air_inlet=15.0, 
            T_air_exh=70.0, 
            T_gc_out_min=0.0, 
            T_ihx_in=40.0,
            Tsat_evap=10.0, 
            eff_motor=0.96, 
            dTmin_phx=5.0,
            tol=1e-6
        ):

        self.Refrigerant = Refrigerant
        self.T_air_tar = T_air_tar
        self.dT_GC_min_tar = dT_GC_min_tar
        self.T_air_inlet = T_air_inlet
        self.T_air_exh = T_air_exh
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
        self.T_ihx_in = T_ihx_in
        self.AIR = 'Air'
        self.tol = tol

        self.converged = False
        self.COP = 0.0

        
        self.init_model()


    def init_model(self):
        self.init_components()
        self.init_compressors()
        self.init_structure()
        self.init_streams()
        self.set_relationships()
        self.define_solver_pathway()


    def init_components(self):
        #Initialisation of unit operations
        self.comp1 = Compressor(Motor_Eff=self.eff_motor, model=0)
        self.comp2 = Compressor(Motor_Eff=self.eff_motor, model=0) 
        self.comp3 = Compressor(Motor_Eff=self.eff_motor, model=0) 
        self.gas_cooler1 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.gas_cooler2 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.gas_cooler4 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.gas_cooler3 = GasCooler(dTmin=self.dT_GC_min_tar)
        self.ihx = PlateHeatExc(dTmin=self.dTmin, model='dTmin_terminal')
        self.phx1 = PlateHeatExc(dTmin=self.dTmin, model='dTmin_terminal')
        self.phx2 = PlateHeatExc(dTmin=self.dTmin, model='dTmin_terminal')        
        self.phx3 = PlateHeatExc(dTmin=self.dTmin, model='dTmin_terminal')
        self.exp_valve = Expansion_Valve()
        self.evaporator = Evaporator(model='fixed cold duty')
        self.evap_sup = PlateHeatExc(model='fixed cold temperature outlet')


    def init_compressors(self):
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


    def init_structure(self):
        #Binary variables for the heat pump structure
        self.comp2_on = 1  if self.n_comp > 1 else 0
        self.comp3_on = 1  if self.n_comp > 2 else 0
        self.ihx_on = 1    if self.T_ihx_in < self.T_gc_out else 0
                
        self.dT_phx1_cold = self.dT_phx1_cold if self.dT_phx1_cold > 0.0 else 0.0
        self.dT_phx2_cold = self.dT_phx2_cold if self.n_comp > 1 and self.dT_phx2_cold > 0.0 else 0.0
        self.dT_phx3_cold = self.dT_phx3_cold if self.n_comp > 2 and self.dT_phx3_cold > 0.0 else 0.0
        self.phx1_on = 1   if self.dT_phx1_cold > 0.0 else 0
        self.phx2_on = 1   if self.dT_phx2_cold > 0.0 else 0
        self.phx3_on = 1   if self.dT_phx3_cold > 0.0 else 0
    

    def init_streams(self):
        #Initialise values in the evaporator streams
        self.evaporator.ColdOut = MassStream( 'Evap_out', self.Refrigerant, Temperature=self.Tsat_evap, Quality = 1 )
        P_low_init = self.evaporator.ColdOut.Pressure
        self.evap_sup.ColdIn = self.evaporator.ColdOut
        self.evap_sup.ColdOut = MassStream( 'IHX_in', self.Refrigerant, Pressure=P_low_init, Temperature=self.T_ihx_in )

        #Initialise values in the compression streams
        self.ihx.ColdIn = self.evap_sup.ColdOut
        if self.ihx_on == 1:
            self.ihx.ColdOut = MassStream( 'PHX1_in', self.Refrigerant, Pressure=P_low_init, Temperature=self.T_ihx_in + max(( self.T_gc_out - self.T_ihx_in ) * self.eff_phx, 0) )
        else:   
            self.ihx.ColdOut = self.ihx.ColdIn

        self.phx1.ColdIn = self.ihx.ColdOut
        if self.phx1_on == 1: 
            self.phx1.ColdOut = MassStream( 'Comp1_in', self.Refrigerant, Pressure=P_low_init, Temperature=self.phx1.ColdIn.Temperature + self.dT_phx1_cold * self.eff_phx )
        else:
            self.phx1.ColdOut = self.phx1.ColdIn

        self.comp1.FlowIn = self.phx1.ColdOut

        P_comp_out_est = REFPROP('P', self.Refrigerant, 'TS', self.T_air_tar + self.dT_GC_min_tar - (self.dT_phx1_cold + self.dT_phx2_cold + self.dT_phx3_cold) / 2, self.comp1.FlowIn.Entropy, 0) * 1000
        
        self.PR_tot = P_comp_out_est / P_low_init
        PR1 = 1
        PR2 = 1
        PR3 = 1
        if self.n_comp == 1:
            PR1 = self.PR_tot
        elif self.n_comp == 2:
            PR1 = (self.PR_tot)**(1/2) 
            PR2 = PR1
        else:
            PR1 = (self.PR_tot)**(1/3) 
            PR2 = PR1; PR3 = PR1
                
        P_comp1_out = P_low_init * PR1
        P_comp2_out = P_comp1_out * PR2
        P_comp3_out = P_comp2_out * PR3

        self.comp1.FlowOut = MassStream( 'Comp1_out', self.Refrigerant, Pressure=P_comp1_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp1_out / 1000, self.comp1.FlowIn.Entropy, 0))

        self.phx2.ColdIn = self.comp1.FlowOut
        if self.phx2_on == 1:
            self.phx2.ColdOut = MassStream( 'Comp2_in', self.Refrigerant, Pressure=self.phx2.ColdIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_cold * self.eff_phx)
        else:
            self.phx2.ColdOut = self.comp1.FlowOut
        
        self.comp2.FlowIn = self.phx2.ColdOut
        if self.comp2_on == 1:
            self.comp2.FlowOut = MassStream( 'Comp2_out', self.Refrigerant, Pressure=P_comp2_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp2_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp2.FlowOut = self.comp2.FlowIn

        self.phx3.ColdIn = self.comp2.FlowOut
        if self.phx3_on == 1:
            self.phx3.ColdOut = MassStream( 'Comp3_in', self.Refrigerant, Pressure=self.phx3.ColdIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_cold * self.eff_phx)
        else:
            self.phx3.ColdOut = self.phx3.ColdIn

        self.comp3.FlowIn = self.phx3.ColdOut
        if self.comp3_on == 1:
            self.comp3.FlowOut = MassStream( 'Comp3_out', self.Refrigerant, Pressure=P_comp3_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp3_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp3.FlowOut = self.comp3.FlowIn


        #Initialisation of gas cooler side of the cycle
        self.gas_cooler4.HotIn = self.comp3.FlowOut
        if self.comp3_on == 1 and self.phx3_on == 1:
            self.gas_cooler4.HotOut = MassStream( 'GC4_out', self.Refrigerant, Pressure= self.gas_cooler4.HotIn.Pressure, Temperature=self.comp2.FlowOut.Temperature + self.dT_phx3_cold )
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = MassStream( 'GC3_in', self.Refrigerant, Pressure=self.phx3.HotIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_cold * self.eff_phx)       
        else:
            self.gas_cooler4.HotOut = self.gas_cooler4.HotIn
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = self.phx3.HotIn  

        self.gas_cooler3.HotIn = self.phx3.HotOut
        if self.comp2_on == 1 and self.phx2_on == 1:
            self.gas_cooler3.HotOut = MassStream( 'GC3_out', self.Refrigerant, Pressure= self.gas_cooler3.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx2_cold )
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = MassStream( 'GC2_in', self.Refrigerant, Pressure=self.phx2.HotIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_cold * self.eff_phx)
        else:
            self.gas_cooler3.HotOut = self.gas_cooler3.HotIn
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = self.phx2.HotIn

        self.gas_cooler2.HotIn = self.phx2.HotOut
        if self.phx1_on == 1:
            self.gas_cooler2.HotOut = MassStream( 'GC2_out', self.Refrigerant, Pressure= self.gas_cooler2.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx1_cold )
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = MassStream( 'GC1_in', self.Refrigerant, Pressure=self.phx1.HotIn.Pressure, Temperature=self.phx1.ColdIn.Temperature + self.dT_phx1_cold * self.eff_phx)
        else:
            self.gas_cooler2.HotOut = self.gas_cooler2.HotIn
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = self.phx1.HotIn 

        self.gas_cooler1.HotIn = self.phx1.HotOut
        self.gas_cooler1.HotOut = MassStream( 'GC1_out', self.Refrigerant, Pressure=self.comp2.FlowOut.Pressure, Temperature=self.T_gc_out )


        #Initialisation of internal exchanger and evaporator system
        self.ihx.HotIn = self.gas_cooler1.HotOut
        if self.ihx_on == 1:
            self.ihx.HotOut = MassStream( 'Ej_pri_in', self.Refrigerant, Pressure=self.gas_cooler1.HotOut.Pressure, Enthalpy=self.gas_cooler1.HotOut.Enthalpy -  ( self.phx1.ColdIn.Enthalpy - self.ihx.ColdIn.Enthalpy ) )
        else:   
            self.ihx.HotOut = self.ihx.HotIn
                
        self.exp_valve.FlowIn = self.ihx.HotOut
        self.exp_valve.FlowOut = MassStream( 'Evap_in', self.Refrigerant, Pressure=P_low_init, Enthalpy=self.exp_valve.FlowIn.Enthalpy )
        self.evaporator.ColdIn = self.exp_valve.FlowOut
        

        #Initialisation of inlet and exhaust air streams
        self.gas_cooler1.ColdIn = MassStream( 'AirIn1', self.AIR, Pressure=101, Temperature=self.T_air_inlet )
        self.gas_cooler1.ColdOut = MassStream( 'AirOut1', self.AIR, Pressure=101, Temperature=self.T_air_tar )
        self.gas_cooler2.ColdIn = self.gas_cooler1.ColdOut
        self.gas_cooler2.ColdOut = MassStream( 'AirOut2', self.AIR, Pressure=101, Temperature=self.T_air_tar )
        self.gas_cooler3.ColdIn = self.gas_cooler2.ColdOut
        self.gas_cooler3.ColdOut = MassStream( 'AirOut3', self.AIR, Pressure=101, Temperature=self.T_air_tar )
        self.gas_cooler4.ColdIn = self.gas_cooler3.ColdOut
        self.gas_cooler4.ColdOut = MassStream( 'AirOut4', self.AIR, Pressure=101, Temperature=self.T_air_tar )

        self.evaporator.HotIn = MassStream( 'exhIn', self.AIR, Pressure=101, Temperature=self.T_air_exh)
        self.evaporator.HotOut = MassStream( 'exhout', self.AIR, Pressure=101, Temperature=self.T_air_exh)
        self.evap_sup.HotIn = MassStream( 'supIn', self.AIR, Pressure=101, Temperature=self.T_air_exh)
        self.evap_sup.HotOut = MassStream( 'suphout', self.AIR, Pressure=101, Temperature=self.T_air_exh)


    def set_relationships(self):
        #Set special relationships between points in the cycle 
        if self.phx3_on == 1:
            self.set_phx3_dTmax = SetRelationship( source_obj=self.phx3.ColdIn, target_obj=[self.phx3.HotIn], c=self.dT_phx3_cold, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx2_on == 1:    
            self.set_phx2_dTmax = SetRelationship( source_obj=self.phx2.ColdIn, target_obj=[self.phx2.HotIn], c=self.dT_phx2_cold, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx1_on == 1:
            self.set_phx1_dTmax = SetRelationship( source_obj=self.phx1.ColdIn, target_obj=[self.phx1.HotIn], c=self.dT_phx1_cold, prop='Temperature', prop2='Pressure' ) 

        self.set_P_ints = CustomRelationship( source_obj=[self.n_comp, self.comp1.FlowIn, self.comp3.FlowOut, self.comp2_power_frac, self.comp3_power_frac], target_obj=[self.comp1.FlowOut, self.comp2.FlowOut], relation=1 )
        self.set_gc_air_massflow = CustomRelationship( source_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], target_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], relation=2 )
        #self.set_evap_air_massflow = CustomRelationship( source_obj=[self.evaporator, self.evap_sup], target_obj=[self.evaporator, self.evap_sup], relation=3)


    def define_solver_pathway(self):
        #Solution pathway set-up
        self.OP = [self.exp_valve, self.evaporator, self.evap_sup, self.set_P_ints]

        if self.ihx_on == 1:
            self.OP.append(self.ihx)
        if self.phx1_on == 1:
            self.OP.append(self.phx1)

        self.OP.append(self.comp1)
        if self.phx2_on == 1:
            self.OP.append(self.phx2)      
        if self.comp2_on == 1:
            self.OP.append(self.comp2)
        if self.phx3_on == 1:
            self.OP.append(self.phx3)       
        if self.comp3_on == 1:
            self.OP.append(self.comp3)

        if self.phx3_on == 1:
            self.OP.append(self.set_phx3_dTmax)
        if self.phx2_on == 1:
            self.OP.append(self.set_phx2_dTmax)
        if self.phx1_on == 1:
            self.OP.append(self.set_phx1_dTmax)  

        self.OP.append(self.set_gc_air_massflow)
        if self.phx3_on == 1:
            self.OP.append(self.gas_cooler4)
        if self.phx2_on == 1:
            self.OP.append(self.gas_cooler3)
        if self.phx1_on == 1:
            self.OP.append(self.gas_cooler2)        
        
        self.OP.append(self.gas_cooler1)

        self.AP = [self.gas_cooler1]
        
        if self.ihx_on == 1:
            self.AP.append(self.ihx)  
        if self.phx1_on == 1:
            self.AP.append(self.phx1)
        if self.phx2_on == 1:
            self.AP.append(self.phx2)        
        if self.phx3_on == 1:
            self.AP.append(self.phx3)  
        if self.phx3_on == 1:
            self.AP.append(self.gas_cooler4)
        if self.phx2_on == 1:
            self.AP.append(self.gas_cooler3)
        if self.phx1_on == 1:
            self.AP.append(self.gas_cooler2)


    def solve(self):
        #Solve heat pump cycle                        
        if self.set_pressure_bounds() > 0.0:
            #Solve the heat pump thermodynamics
            r = brentq(f=self.state_solve, a=self.P_low_bound, b=self.P_upp_bound, full_output=True, rtol=self.tol, xtol=self.tol)

            if r[1].converged and self.converged and abs(self.gas_cooler4.ColdOut.Temperature - self.T_air_tar) < 0.05:
                self.COP = (self.gas_cooler4.Duty + self.gas_cooler3.Duty + self.gas_cooler2.Duty + self.gas_cooler1.Duty) / (self.comp1.Power + self.comp2.Power + self.comp3.Power)
                        
        if self.COP == 0:
            self.converged = False          
        
        self.summarise_results()


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
        h_exp_in = 0.0
        m_air_out = 0.0
        dP_gc = 0.0
        duty_comp = 0.0

        while ( 
                (
                    i < 5 
                    or abs(self.exp_valve.FlowIn.Enthalpy - h_exp_in) > self.tol * h_exp_in 
                    or abs(self.gas_cooler4.ColdOut.Mass_Flow - m_air_out) > self.tol * m_air_out 
                    or abs(self.gas_cooler4.HotIn.Pressure - self.gas_cooler1.HotOut.Pressure - dP_gc) > self.tol * 10 * dP_gc
                    or abs(self.comp1.Power + self.comp2.Power + self.comp3.Power - duty_comp) > self.tol * 10 * duty_comp
                ) 
                and i < 100
              ): #Check for model convergance
            
            #Saves critical pressure and temperature properties to check for solution convergance
            h_exp_in = self.exp_valve.FlowIn.Enthalpy
            m_air_out = self.gas_cooler4.ColdOut.Mass_Flow
            dP_gc = self.gas_cooler4.HotIn.Pressure - self.gas_cooler1.HotOut.Pressure
            duty_comp = self.comp1.Power + self.comp2.Power + self.comp3.Power
            i += 1

            #Solve the performance of each component and update each special relationship
            for obj in self.OP:
                obj.solve()

            if i < 100:            
                self.converged = True
            
        return self.gas_cooler4.ColdOut.Temperature - self.T_air_tar


    def summarise_results(self):
        self.output = {
            #Key results
            'COP:' : self.COP,
            'Tol:' : abs(self.T_air_tar - self.gas_cooler4.ColdOut.Temperature) / self.T_air_tar,
            'Converged:' : self.converged,
            #Parameters
            'Refrigerant:' : self.Refrigerant,
            'dT_min Gas Cooler:' : self.dT_GC_min_tar,
            'T_air_tar:' : self.T_air_tar,
            'T_air_inlet:' : self.gas_cooler1.ColdIn.Temperature,
            'Dryer exh air T:' : self.evaporator.HotIn.Temperature,
            'n_comp:' : self.n_comp,
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
    