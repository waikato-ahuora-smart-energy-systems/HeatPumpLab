__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes.compressor import Compressor
from classes.plateHeatExc import PlateHeatExc
from classes.gasCooler import GasCooler
from classes.expansionValve import Expansion_Valve
from classes.evaporator import Evaporator
from classes.stream import Stream
from classes.set_relationship import Set_Relationship
from classes.customRelationship import CustomRelationship
from globalVariables import CONVERGANCE
from functions.refprop import REFPROP
from scipy.optimize import brentq
from functions.graph import graph_heat_exchangers
from functions.output import write_screen


class TRHP:
    def __init__(self, Refrigerant='R1336mzz(Z)', T_air_tar=200.0, dT_GC_min_tar=20.0, n_comp = 3, dT_phx1_max=0.0, dT_phx2_max=0.0, dT_phx3_max=0.0, 
                 PR_2_init=0.0, PR_3_init=0.0, T_air_inlet=15.0, T_air_exh=70.0, RH_air_exh = 0.15, T_gc_out_min=0.0, T_ihx_in_max=50.0, eff_motor=0.95, 
                 eff_phx = 0.9, print_output=False, graph_output=False):
        
        self.CONVERGANCE = CONVERGANCE
        self.print_output = print_output
        self.graph_output = graph_output

        self.Refrigerant = Refrigerant
        self.T_air_tar = T_air_tar
        self.dT_GC_min_tar = dT_GC_min_tar
        self.T_air_inlet = T_air_inlet
        self.T_air_exh = T_air_exh
        
        if PR_2_init < 1.0 and PR_3_init < 1.0:
            self.comp_PR_type = 'Constant' 
        else:
            self.comp_PR_type = 'Variable'
            if PR_2_init < 1.0:
                PR_2_init = PR_3_init
            if PR_3_init < 1.0:
                PR_3_init = PR_2_init

        self.PR_2_init = max(PR_2_init, 1.0)
        self.PR_3_init = max(PR_3_init, 1.0)
        self.n_comp = n_comp

        self.eff_phx = eff_phx
        self.eff_motor = eff_motor

        self.dT_phx1_max = dT_phx1_max if dT_phx1_max > 0.0 else 0.0
        self.dT_phx2_max = dT_phx2_max if self.n_comp > 1 and dT_phx2_max > 0.0 else 0.0
        self.dT_phx3_max = dT_phx3_max if self.n_comp > 2 and dT_phx3_max > 0.0 else 0.0

        self.T_gc_out = max(T_air_inlet + dT_GC_min_tar, T_gc_out_min)
        self.T_ihx_in = min(T_air_exh - dT_GC_min_tar, T_ihx_in_max)
        
        self.converged = False
        self.COP = 0        
        self.AIR = 'Air'
        self.T_evap_init = 15.0

        self.determine_HP_structure()
        self.create_components()
        self.stream_initialisation()
        self.set_dTmin_relationships()
        self.set_gas_cooler_analysis()
        self.define_solver_pathway()


    def determine_HP_structure(self):
        #Binary variables for the heat pump structure
        self.comp2_on = 1  if self.n_comp > 1 else 0
        self.comp3_on = 1  if self.n_comp > 2 else 0
        self.ihx_on = 1    if self.T_ihx_in < self.T_gc_out else 0
        self.phx1_on = 1   if self.dT_phx1_max > 0.0 else 0
        self.phx2_on = 1   if self.dT_phx2_max > 0.0 else 0
        self.phx3_on = 1   if self.dT_phx3_max > 0.0 else 0

        self.n_comp = 1 + self.comp2_on + self.comp3_on 
    

    def create_components(self):
        #Initialisation of unit operations
        self.comp1 = Compressor(self.eff_motor)
        self.comp2 = Compressor(self.eff_motor) 
        self.comp3 = Compressor(self.eff_motor) 
        self.gas_cooler1 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=self.dT_GC_min_tar)
        self.gas_cooler2 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=self.dT_GC_min_tar)
        self.gas_cooler4 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=self.dT_GC_min_tar)
        self.gas_cooler3 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=self.dT_GC_min_tar)
        self.ihx = PlateHeatExc(self.eff_phx)
        self.phx1 = PlateHeatExc(self.eff_phx)
        self.phx2 = PlateHeatExc(self.eff_phx)        
        self.phx3 = PlateHeatExc(self.eff_phx)
        self.exp_valve = Expansion_Valve()
        self.evaporator = Evaporator()
        self.evap_sup = Evaporator()
        

    def stream_initialisation(self):
        #Initialise values in the evaporator streams
        self.evaporator.ColdOut = Stream( 'Evap_out', self.Refrigerant, Temperature=self.T_evap_init, Quality=1 )
        P_low_init = self.evaporator.ColdOut.Pressure
        self.evap_sup.ColdIn = self.evaporator.ColdOut
        self.evap_sup.ColdOut = Stream( 'IHX_in', self.Refrigerant, Pressure=P_low_init, Temperature=self.T_ihx_in )

        #Initialise values in the compression streams
        self.ihx.ColdIn = self.evap_sup.ColdOut
        if self.ihx_on == 1:
            self.ihx.ColdOut = Stream( 'PHX1_in', self.Refrigerant, Pressure=P_low_init, Temperature=self.T_ihx_in + max(( self.T_gc_out - self.T_ihx_in ) * self.eff_phx, 0) )
        else:   
            self.ihx.ColdOut = self.ihx.ColdIn

        self.phx1.ColdIn = self.ihx.ColdOut
        if self.phx1_on == 1: 
            self.phx1.ColdOut = Stream( 'Comp1_in', self.Refrigerant, Pressure=P_low_init, Temperature=self.phx1.ColdIn.Temperature + self.dT_phx1_max * self.eff_phx )
        else:
            self.phx1.ColdOut = self.phx1.ColdIn

        self.comp1.FlowIn = self.phx1.ColdOut

        P_comp_out_est = REFPROP('P', self.Refrigerant, 'TS', self.T_air_tar + self.dT_GC_min_tar - (self.dT_phx1_max + self.dT_phx2_max + self.dT_phx3_max) / 2, self.comp1.FlowIn.Entropy, 0) * 1000
        
        if self.n_comp == 1:
            P_comp1_out = P_comp_out_est
        elif self.n_comp == 2:
            if self.comp_PR_type == 'Constant':
                P_comp1_out = P_low_init ** (1 / 2) * P_comp_out_est ** (1 / 2)    
            else:
                PR_1 = P_comp_out_est / P_low_init * 1 / self.PR_2_init
                P_comp1_out = P_low_init * PR_1
            P_comp2_out = P_comp_out_est
        elif self.n_comp == 3:
            if self.comp_PR_type == 'Constant':
                P_comp1_out = P_low_init ** (2 / 3) * P_comp_out_est ** (1 / 3)
                P_comp2_out = P_low_init ** (1 / 3) * P_comp_out_est ** (2 / 3)
            else:
                PR_1 = P_comp_out_est / P_low_init * 1 / (self.PR_2_init * self.PR_3_init)
                P_comp1_out = P_low_init * PR_1
                P_comp2_out = P_comp1_out * self.PR_2_init
            P_comp3_out = P_comp_out_est
        
        self.comp1.FlowOut = Stream( 'Comp1_out', self.Refrigerant, Pressure=P_comp1_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp1_out / 1000, self.comp1.FlowIn.Entropy, 0))

        self.phx2.ColdIn = self.comp1.FlowOut
        if self.phx2_on == 1:
            self.phx2.ColdOut = Stream( 'Comp2_in', self.Refrigerant, Pressure=self.phx2.ColdIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_max * self.eff_phx)
        else:
            self.phx2.ColdOut = self.comp1.FlowOut
        
        self.comp2.FlowIn = self.phx2.ColdOut
        if self.comp2_on == 1:
            self.comp2.FlowOut = Stream( 'Comp2_out', self.Refrigerant, Pressure=P_comp2_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp2_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp2.FlowOut = self.comp2.FlowIn

        self.phx3.ColdIn = self.comp2.FlowOut
        if self.phx3_on == 1:
            self.phx3.ColdOut = Stream( 'Comp3_in', self.Refrigerant, Pressure=self.phx3.ColdIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_max * self.eff_phx)
        else:
            self.phx3.ColdOut = self.phx3.ColdIn

        self.comp3.FlowIn = self.phx3.ColdOut
        if self.comp3_on == 1:
            self.comp3.FlowOut = Stream( 'Comp3_out', self.Refrigerant, Pressure=P_comp3_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp3_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp3.FlowOut = self.comp3.FlowIn


        #Initialisation of gas cooler side of the cycle
        self.gas_cooler4.HotIn = self.comp3.FlowOut
        if self.comp3_on == 1 and self.phx3_on == 1:
            self.gas_cooler4.HotOut = Stream( 'GC4_out', self.Refrigerant, Pressure= self.gas_cooler4.HotIn.Pressure, Temperature=self.comp2.FlowOut.Temperature + self.dT_phx3_max )
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = Stream( 'GC3_in', self.Refrigerant, Pressure=self.phx3.HotIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_max * self.eff_phx)       
        else:
            self.gas_cooler4.HotOut = self.gas_cooler4.HotIn
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = self.phx3.HotIn  

        self.gas_cooler3.HotIn = self.phx3.HotOut
        if self.comp2_on == 1 and self.phx2_on == 1:
            self.gas_cooler3.HotOut = Stream( 'GC3_out', self.Refrigerant, Pressure= self.gas_cooler3.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx2_max )
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = Stream( 'GC2_in', self.Refrigerant, Pressure=self.phx2.HotIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_max * eff_phx)
        else:
            self.gas_cooler3.HotOut = self.gas_cooler3.HotIn
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = self.phx2.HotIn

        self.gas_cooler2.HotIn = self.phx2.HotOut
        if self.phx1_on == 1:
            self.gas_cooler2.HotOut = Stream( 'GC2_out', self.Refrigerant, Pressure= self.gas_cooler2.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx1_max )
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = Stream( 'GC1_in', self.Refrigerant, Pressure=self.phx1.HotIn.Pressure, Temperature=self.phx1.ColdIn.Temperature + self.dT_phx1_max * eff_phx)
        else:
            self.gas_cooler2.HotOut = self.gas_cooler2.HotIn
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = self.phx1.HotIn 

        self.gas_cooler1.HotIn = self.phx1.HotOut
        self.gas_cooler1.HotOut = Stream( 'GC1_out', self.Refrigerant, Pressure=self.comp2.FlowOut.Pressure, Temperature=self.T_gc_out )


        #Initialisation of internal exchanger and evaporator system
        self.ihx.HotIn = self.gas_cooler1.HotOut
        if self.ihx_on == 1:
            self.ihx.HotOut = Stream( 'Ej_pri_in', self.Refrigerant, Pressure=self.gas_cooler1.HotOut.Pressure, Enthalpy=self.gas_cooler1.HotOut.Enthalpy -  ( self.phx1.ColdIn.Enthalpy - self.ihx.ColdIn.Enthalpy ) )
        else:   
            self.ihx.HotOut = self.ihx.HotIn
                
        self.exp_valve.FlowIn = self.ihx.HotOut
        self.exp_valve.FlowOut = Stream( 'Evap_in', self.Refrigerant, Pressure=P_low_init, Enthalpy=self.exp_valve.FlowIn.Enthalpy )
        self.evaporator.ColdIn = self.exp_valve.FlowOut
        

        #Initialisation of inlet and exhaust air streams
        self.gas_cooler1.ColdIn = Stream( 'AirIn1', self.AIR, Pressure=101, Temperature=self.T_air_inlet )
        self.gas_cooler1.ColdOut = Stream( 'AirOut1', self.AIR, Pressure=101, Temperature=self.T_air_tar )
        self.gas_cooler2.ColdIn = self.gas_cooler1.ColdOut
        self.gas_cooler2.ColdOut = Stream( 'AirOut2', self.AIR, Pressure=101, Temperature=self.T_air_tar )
        self.gas_cooler3.ColdIn = self.gas_cooler2.ColdOut
        self.gas_cooler3.ColdOut = Stream( 'AirOut3', self.AIR, Pressure=101, Temperature=self.T_air_tar )
        self.gas_cooler4.ColdIn = self.gas_cooler3.ColdOut
        self.gas_cooler4.ColdOut = Stream( 'AirOut4', self.AIR, Pressure=101, Temperature=self.T_air_tar )

        self.evaporator.HotIn = Stream( 'exhIn', self.AIR, Pressure=101, Temperature=self.T_air_exh)
        self.evaporator.HotOut = Stream( 'exhout', self.AIR, Pressure=101, Temperature=self.T_air_exh)
        self.evap_sup.HotIn = Stream( 'supIn', self.AIR, Pressure=101, Temperature=self.T_air_exh)
        self.evap_sup.HotOut = Stream( 'suphout', self.AIR, Pressure=101, Temperature=self.T_air_exh)


    def set_dTmin_relationships(self):
        #Set special relationships between points in the cycle 
        if self.phx3_on == 1:
            self.set_phx3_dTmax = Set_Relationship( source_obj=self.phx3.ColdIn, target_obj=[self.phx3.HotIn], c=self.dT_phx3_max, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx2_on == 1:    
            self.set_phx2_dTmax = Set_Relationship( source_obj=self.phx2.ColdIn, target_obj=[self.phx2.HotIn], c=self.dT_phx2_max, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx1_on == 1:
            self.set_phx1_dTmax = Set_Relationship( source_obj=self.phx1.ColdIn, target_obj=[self.phx1.HotIn], c=self.dT_phx1_max, prop='Temperature', prop2='Pressure' ) 

        self.set_P_ints = CustomRelationship( source_obj=[self.comp1.FlowIn, self.comp1.FlowOut, self.n_comp], target_obj=[self.comp1.FlowOut, self.comp2.FlowOut, self.comp3.FlowOut], var1=[self.comp_PR_type, self.PR_2_init, self.PR_3_init], relation=1 )
        self.set_gc_air_massflow = CustomRelationship( source_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], target_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], relation=2 )
        self.set_evap_air_massflow = CustomRelationship( source_obj=[self.evaporator, self.evap_sup], target_obj=[self.evaporator, self.evap_sup], relation=3)
    

    def set_gas_cooler_analysis(self):
        #Setup heat exchanger elements
        if self.phx3_on == 1: #therefore gas cooler 4,3,2,1 is on 
            self.N_elements_per_pass = 2
        if self.phx2_on == 1 and self.phx3_on == 0: #therefore gas cooler 3,2,1 is on
            self.N_elements_per_pass = 3
        if self.phx1_on == 1 and self.phx2_on == 0 and self.phx3_on == 0: #therefore gas cooler 2,1 is on
            self.N_elements_per_pass = 3
        else:
            self.N_elements_per_pass = 4  

        self.detailed_model = False
        if self.detailed_model: 
            self.gas_cooler1.init_heat_exchanger_elements(self.N_elements_per_pass)
            self.gas_cooler2.init_heat_exchanger_elements(self.N_elements_per_pass)
            self.gas_cooler3.init_heat_exchanger_elements(self.N_elements_per_pass)
            self.gas_cooler4.init_heat_exchanger_elements(self.N_elements_per_pass)


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
            r = brentq(f=self.state_solve, a=self.P_low_bound, b=self.P_upp_bound, full_output=True, rtol=self.CONVERGANCE, xtol=self.CONVERGANCE)

            if r[1].converged and self.converged and abs(self.gas_cooler4.ColdOut.Temperature - self.T_air_tar) < 0.05:
                self.COP = (self.gas_cooler4.Duty + self.gas_cooler3.Duty + self.gas_cooler2.Duty + self.gas_cooler1.Duty) / (self.comp1.Power + self.comp2.Power + self.comp3.Power)
                        
        if self.COP == 0:
            self.converged = False          
        
        if self.print_output:
            write_screen(heat_pump=self)

        if self.graph_output:
            graph_heat_exchangers([self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4, self.phx1, self.phx2, self.phx3])
            graph_heat_exchangers([self.gas_cooler1, self.phx1, self.gas_cooler2, self.phx2, self.gas_cooler3, self.phx3, self.gas_cooler4])

        return self.COP, self    
    

    def set_pressure_bounds(self):
        P_estimate_high = REFPROP('P', self.comp1.FlowIn.Fluid, 'TS', self.T_air_tar + (self.dT_GC_min_tar - self.dT_phx1_max + self.dT_phx2_max + self.dT_phx3_max) / 2, self.comp1.FlowIn.Entropy, 0) * 1000
        if self.comp_PR_type == 'Constant':
            P_comp1_est = ( P_estimate_high / self.comp1.FlowIn.Pressure ) ** ( 1 / self.n_comp) * self.comp1.FlowIn.Pressure
        else:   
            P_comp1_est = P_estimate_high / ( self.PR_2_init * self.PR_3_init ) 
        
        self.P_low_bound = P_comp1_est
        self.P_upp_bound = P_comp1_est

        dt_plus = self.state_solve(P_comp1_est)
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
        self.comp1.FlowOut.Pressure = P_comp_out
        i = 0
        h_exp_in = 0.0
        m_air_out = 0.0
        dP_gc = 0.0
        duty_comp = 0.0

        while ( (
                    i < 5 
                    or abs(self.exp_valve.FlowIn.Enthalpy - h_exp_in) > self.CONVERGANCE * h_exp_in 
                    or abs(self.gas_cooler4.ColdOut.Mass_Flow - m_air_out) > self.CONVERGANCE * m_air_out 
                    or abs(self.gas_cooler4.HotIn.Pressure - self.gas_cooler1.HotOut.Pressure - dP_gc) > self.CONVERGANCE * 10 * dP_gc
                    or abs(self.comp1.Power + self.comp2.Power + self.comp3.Power - duty_comp) > self.CONVERGANCE * 10 * duty_comp
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


    def Stream_outputs(self):
        self.duty_GC4 = self.gas_cooler4.HotIn.Mass_Flow * (self.gas_cooler4.HotIn.Enthalpy - self.gas_cooler4.HotOut.Enthalpy)
        self.duty_GC3 = self.gas_cooler3.HotIn.Mass_Flow * (self.gas_cooler3.HotIn.Enthalpy - self.gas_cooler3.HotOut.Enthalpy)
        self.duty_GC2 = self.gas_cooler2.HotIn.Mass_Flow * (self.gas_cooler2.HotIn.Enthalpy - self.gas_cooler2.HotOut.Enthalpy)
        self.duty_GC1 = self.gas_cooler1.HotIn.Mass_Flow * (self.gas_cooler1.HotIn.Enthalpy - self.gas_cooler1.HotOut.Enthalpy)
        self.duty_GC_total = self.duty_GC4 + self.duty_GC3 + self.duty_GC2 + self.duty_GC1
        
        return [self.comp1.FlowIn.Pressure, self.comp1.FlowOut.Pressure, self.comp2.FlowOut.Pressure, self.comp3.FlowOut.Pressure, self.gas_cooler4.ColdOut.Temperature, self.gas_cooler3.ColdOut.Temperature, self.duty_GC_total]