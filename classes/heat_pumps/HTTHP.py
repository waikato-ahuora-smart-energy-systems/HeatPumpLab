__author__ = 'Tim Walmsley and Lana Kong'
__credits__ = ['']


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


class HTTHP:
    def __init__(
            self, 
            Refrigerant = ['CO2','CO2'], 
            T_air_tar = 200.0, 
            dT_GC_min_tar = 20.0, 
            dT_ihx_top = 0.0,       #Optimisation variable
            dT_ihx_btm = 0.0,       #Optimisation variable
            dT_ihx_tb = 0.0,        #Optimisation variable
            T_gc_out = 0.0,         #Optimisation variable
            comp_PR_uniform = False, 
            T_air_inlet = 15.0, 
            T_air_exh = 70.0, 
            T_ihx_btm_in = 40.0,
            Tsat_evap = [10.0, 10.0], 
            eff_motor = 0.96, 
            dT_min_ihx = 0.9,
            tolerance=1e-6
        ):

        self.refrigerant_btm = Refrigerant[0]
        self.refrigerant_top = Refrigerant[1]

        self.T_air_tar = T_air_tar
        self.dT_GC_min_tar = dT_GC_min_tar
        self.T_air_inlet = T_air_inlet
        self.T_air_exh = T_air_exh
        self.T_gc_out = T_gc_out
        self.Tsat_evap = Tsat_evap
        self.comp_PR_uniform = comp_PR_uniform
        self.dT_min_ihx = dT_min_ihx
        self.eff_motor = eff_motor
        self.dT_ihx_top = dT_ihx_top
        self.dT_ihx_btm = dT_ihx_btm 
        self.dT_ihx_tb = dT_ihx_tb
        self.T_ihx_btm_in = T_ihx_btm_in
        self.AIR = 'Air'
        self.WATER = 'Water'
        self.tolerance = tolerance

        self.create_components()
        self.model_initialisation()


    def create_components(self):
        #Initialisation of unit operations
        self.exp_valve_btm = ExpansionValve()
        self.evaporator_btm = Evaporator()
        self.evap_sup_btm = Evaporator()
        self.ihx_btm = PlateHeatExchanger(self.dT_min_ihx)
        self.comp_btm = Compressor(self.eff_motor)
        self.gas_cooler_cas = GasCooler(HotFluid=self.refrigerant_btm, ColdFluid=self.refrigerant_top, dT_GC_min_tar=self.dT_GC_min_tar)
        self.comp_top = Compressor(self.eff_motor) 
        self.gas_cooler_air = GasCooler(HotFluid=self.refrigerant_top, ColdFluid=self.AIR, dT_GC_min_tar=self.dT_GC_min_tar)
        self.ihx_top = PlateHeatExchanger(self.dT_min_ihx) #Need to update class
        self.ihx_tb = PlateHeatExchanger(self.dT_min_ihx)
        self.exp_valve_top = ExpansionValve()
        self.evaporator_top = Evaporator()
        self.evap_sup_top = Evaporator()
        

    def model_initialisation(self):
        self.converged = False
        self.COP = 0.0        
        
        self.dT_ihx_top = self.dT_ihx_top if self.dT_ihx_top > 0.0 else 0.0
        self.dT_ihx_btm = self.dT_ihx_btm if self.dT_ihx_btm > 0.0 else 0.0
        self.dT_ihx_tb = self.dT_ihx_tb if self.dT_ihx_tb > 0.0 else 0.0

        self.T_gc_out = max(self.T_air_inlet + self.dT_GC_min_tar, self.T_gc_out)
        
        self.determine_HP_structure()
        self.stream_initialisation()
        self.set_relationships()
        self.define_solver_pathway()


    def determine_HP_structure(self):
        #Binary variables for the heat pump structure
        self.ihx_top_on = 1    if self.dT_ihx_top < self.T_gc_out else 0
        self.ihx_btm_on = 1   if self.dT_ihx_btm > 0.0 else 0
        self.ihx_tb_on = 1   if self.dT_ihx_tb > 0.0 else 0


    def stream_initialisation(self):
        #Initialise values in the evaporator streams
        self.evaporator_btm.ColdOut = MassStream( 'Evap_btm_out', self.refrigerant_btm, Temperature=self.Tsat_evap, Quality = 1 )
        self.evap_sup_btm.ColdIn = self.evaporator_btm.ColdOut
        self.evap_sup_btm.ColdOut = MassStream( 'IHX_in', self.refrigerant_btm, Pressure=self.evaporator_btm.ColdOut.Pressure, Temperature=self.T_ihx_btm_in )

        #Initialise values in the compression streams
        self.ihx_btm.ColdIn = self.evap_sup_btm.ColdOut
        if self.ihx_btm_on == 1:
            self.ihx_btm.ColdOut = MassStream( 'PHX1_in', self.refrigerant_btm, Pressure=self.evaporator_btm.ColdOut.Pressure, Temperature=self.T_ihx_btm_in + self.dT_ihx_btm)
        else:   
            self.ihx_btm.ColdOut = self.ihx_btm.ColdIn

        self.ihx_tb.ColdIn = self.ihx_btm.ColdOut
        if self.ihx_tb_on == 1: 
            self.ihx_tb.ColdOut = MassStream( 'Comp_btm_in', self.refrigerant_btm, Pressure=self.evaporator_btm.ColdOut.Pressure, Temperature=self.ihx_tb.ColdIn.Temperature + self.dT_ihx_tb )
        else:
            self.ihx_tb.ColdOut = self.ihx_tb.ColdIn

        self.comp_btm.FlowIn = self.ihx_tb.ColdOut

        P_comp_out_est = REFPROP('P', self.refrigerant_btm, 'TS', self.T_air_tar + self.dT_GC_min_tar - (self.dT_ihx_btm + self.dT_ihx_top + self.dT_ihx_tb) / 3, self.comp_btm.FlowIn.Entropy, 0) * 1000        

        self.comp_btm.FlowOut = MassStream( 'Comp_btm_out', self.refrigerant_btm, Pressure=P_comp_out_est, Temperature=REFPROP('T', self.refrigerant_btm, 'PS', P_comp_btm_out / 1000, self.comp_btm.FlowIn.Entropy, 0))
        self.gas_cooler_cas.HotIn = self.comp_btm.FlowOut
        self.gas_cooler_cas.HotOut = MassStream( 'Gas_cooler_out', self.refrigerant_btm, Pressure=self.gas_cooler_cas.HotIn.Pressure, Temperature=self.T_gc_out)

        if self.ihx_btm_on == 1:
            self.ihx_btm.HotOut = MassStream( 'Comp2_in', self.refrigerant_btm, Pressure=self.ihx_btm.ColdIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_max * self.eff_phx)
        else:
            self.ihx_btm.HotOut = self.gas_cooler_cas.HotOut
            
        self.exp_valve_btm.FlowIn = self.ihx_btm.HotOut

        self.exp_valve_btm.FlowOut = MassStream( 'Evap_btm_in', self.refrigerant_btm, Pressure=self.evaporator_btm.ColdOut.Pressure, Enthalpy=self.exp_valve_btm.FlowIn.Enthalpy )





        self.phx2.ColdIn = self.comp1.FlowOutMassStream( 'Evap_in', self.Refrigerant, Pressure=self.evaporator_btm.ColdOut.Pressure, Enthalpy=self.exp_valve.FlowIn.Enthalpy )

        if self.phx2_on == 1:
            self.phx2.ColdOut = MassStream( 'Comp2_in', self.Refrigerant, Pressure=self.phx2.ColdIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_max * self.eff_phx)
        else:
            self.phx2.ColdOut = self.comp1.FlowOut
        

        


        self.comp2.FlowIn = self.phx2.ColdOut
        if self.comp2_on == 1:
            self.comp2.FlowOut = MassStream( 'Comp2_out', self.Refrigerant, Pressure=P_comp2_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp2_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp2.FlowOut = self.comp2.FlowIn

        self.phx3.ColdIn = self.comp2.FlowOut
        if self.phx3_on == 1:
            self.phx3.ColdOut = MassStream( 'Comp3_in', self.Refrigerant, Pressure=self.phx3.ColdIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_max * self.eff_phx)
        else:
            self.phx3.ColdOut = self.phx3.ColdIn

        self.comp3.FlowIn = self.phx3.ColdOut
        if self.comp3_on == 1:
            self.comp3.FlowOut = MassStream( 'Comp3_out', self.Refrigerant, Pressure=P_comp3_out, Temperature=REFPROP('T', self.Refrigerant, 'PS', P_comp3_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp3.FlowOut = self.comp3.FlowIn


        #P_comp_top_out = P_comp_out_est

        #Initialisation of gas cooler side of the cycle
        self.gas_cooler4.HotIn = self.comp3.FlowOut
        if self.comp3_on == 1 and self.phx3_on == 1:
            self.gas_cooler4.HotOut = MassStream( 'GC4_out', self.Refrigerant, Pressure= self.gas_cooler4.HotIn.Pressure, Temperature=self.comp2.FlowOut.Temperature + self.dT_phx3_max )
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = MassStream( 'GC3_in', self.Refrigerant, Pressure=self.phx3.HotIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_max * self.eff_phx)       
        else:
            self.gas_cooler4.HotOut = self.gas_cooler4.HotIn
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = self.phx3.HotIn  

        self.gas_cooler3.HotIn = self.phx3.HotOut
        if self.comp2_on == 1 and self.phx2_on == 1:
            self.gas_cooler3.HotOut = MassStream( 'GC3_out', self.Refrigerant, Pressure= self.gas_cooler3.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx2_max )
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = MassStream( 'GC2_in', self.Refrigerant, Pressure=self.phx2.HotIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_max * self.eff_phx)
        else:
            self.gas_cooler3.HotOut = self.gas_cooler3.HotIn
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = self.phx2.HotIn

        self.gas_cooler2.HotIn = self.phx2.HotOut
        if self.phx1_on == 1:
            self.gas_cooler2.HotOut = MassStream( 'GC2_out', self.Refrigerant, Pressure= self.gas_cooler2.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx1_max )
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = MassStream( 'GC1_in', self.Refrigerant, Pressure=self.phx1.HotIn.Pressure, Temperature=self.phx1.ColdIn.Temperature + self.dT_phx1_max * self.eff_phx)
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
        self.exp_valve.FlowOut = MassStream( 'Evap_in', self.Refrigerant, Pressure=self.evaporator_btm.ColdOut.Pressure, Enthalpy=self.exp_valve.FlowIn.Enthalpy )
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
            self.set_phx3_dTmax = SetRelationship( source_obj=self.phx3.ColdIn, target_obj=[self.phx3.HotIn], c=self.dT_phx3_max, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx2_on == 1:    
            self.set_phx2_dTmax = SetRelationship( source_obj=self.phx2.ColdIn, target_obj=[self.phx2.HotIn], c=self.dT_phx2_max, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx1_on == 1:
            self.set_phx1_dTmax = SetRelationship( source_obj=self.phx1.ColdIn, target_obj=[self.phx1.HotIn], c=self.dT_phx1_max, prop='Temperature', prop2='Pressure' ) 

        self.set_P_ints = CustomRelationship( source_obj=[self.n_comp, self.comp1.FlowIn, self.comp3.FlowOut, self.comp2_power_frac, self.comp3_power_frac], target_obj=[self.comp1.FlowOut, self.comp2.FlowOut], relation=1 )
        self.set_gc_air_massflow = CustomRelationship( source_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], target_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], relation=2 )
        #self.set_evap_air_massflow = CustomRelationship( source_obj=[self.evaporator, self.evap_sup], target_obj=[self.evaporator, self.evap_sup], relation=3)
    

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
            r = brentq(f=self.state_solve, a=self.P_low_bound, b=self.P_upp_bound, full_output=True, rtol=self.tolerance, xtol=self.tolerance)

            if r[1].converged and self.converged and abs(self.gas_cooler4.ColdOut.Temperature - self.T_air_tar) < 0.05:
                self.COP = (self.gas_cooler4.Duty + self.gas_cooler3.Duty + self.gas_cooler2.Duty + self.gas_cooler1.Duty) / (self.comp1.Power + self.comp2.Power + self.comp3.Power)
                        
        if self.COP == 0:
            self.converged = False          
        
        self.summarise_results()


    def set_pressure_bounds(self):
        P_gc_est = REFPROP('P', self.comp1.FlowIn.Fluid, 'TS', self.T_air_tar + (self.dT_GC_min_tar - self.dT_phx1_max + self.dT_phx2_max + self.dT_phx3_max) / 2, self.comp1.FlowIn.Entropy, 0) * 1000

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
                    or abs(self.exp_valve.FlowIn.Enthalpy - h_exp_in) > self.tolerance * h_exp_in 
                    or abs(self.gas_cooler4.ColdOut.Mass_Flow - m_air_out) > self.tolerance * m_air_out 
                    or abs(self.gas_cooler4.HotIn.Pressure - self.gas_cooler1.HotOut.Pressure - dP_gc) > self.tolerance * 10 * dP_gc
                    or abs(self.comp1.Power + self.comp2.Power + self.comp3.Power - duty_comp) > self.tolerance * 10 * duty_comp
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
            'dT phx max 1:' : self.dT_phx1_max,
            'dT phx max 2:' : self.dT_phx2_max,
            'dT phx max 3:' : self.dT_phx3_max,
            'PR 1:' : self.comp1.Pressure_Ratio,
            'PR 2:' : self.comp2.Pressure_Ratio,
            'PR 3:' : self.comp3.Pressure_Ratio,
            'Overall P ratio:' : self.gas_cooler4.HotIn.Pressure / self.evaporator.ColdIn.Pressure,
            'Gas cooler(s) inlet T:' : self.gas_cooler4.HotIn.Temperature,
            'Gas cooler(s) outlet T:' : self.gas_cooler1.HotOut.Temperature
        }
    