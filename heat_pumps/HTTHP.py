__author__ = 'Tim Walmsley and Lana Kong'
__credits__ = ['']

from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator
from streams import get_property, MassStream
from logic_operations import SetRelationship, CustomRelationship
from scipy.optimize import brentq

class HTTHP:
    def __init__(
            self, 
            refrigerant = ['CO2','CO2'], 
            T_air_tar = 200.0, 
            dT_GC_min_tar = 20.0, 
            dT_ihx_top = 0.0,       #Optimisation variable
            dT_ihx_btm = 0.0,       #Optimisation variable
            dT_ihx_tb = 0.0,        #Optimisation variable
            T_gc_top_out = 0.0,         #Optimisation variable
            comp_PR_uniform = False, 
            T_air_inlet = 15.0, 
            T_air_exh = 70.0, 
            dT_ihx_btm_in = 20.0,   #CHECK
            dT_ihx_top_in = 5.0,    #CHECK
            Tsat_evap = [10.0, 10.0], 
            eff_motor = 0.96, 
            dT_min_ihx = 0.9,
            tolerance=1e-6
        ):

        fluid=self.refrigerant_btm = refrigerant[0]
        fluid=self.refrigerant_top = refrigerant[1]

        self.T_air_tar = T_air_tar
        self.dT_GC_min_tar = dT_GC_min_tar
        self.T_air_inlet = T_air_inlet
        self.T_air_exh = T_air_exh
        self.T_gc_top_out = T_gc_top_out
        self.Tsat_evap = Tsat_evap
        self.comp_PR_uniform = comp_PR_uniform
        self.dT_min_ihx = dT_min_ihx
        self.eff_motor = eff_motor
        self.dT_ihx_top = dT_ihx_top
        self.dT_ihx_btm = dT_ihx_btm 
        self.dT_ihx_tb = dT_ihx_tb
        self.dT_ihx_btm_in = dT_ihx_btm_in
        self.dT_ihx_top_in = dT_ihx_top_in
        self.AIR = 'Air'
        self.WATER = 'Water'
        self.tolerance = tolerance

        self.init_components()
        self.model_initialisation()

    def init_components(self):
        #Initialisation of unit operations
        self.exp_valve_btm = ExpansionValve()
        self.evaporator_btm = Evaporator(model='fixed cold duty')
        self.evap_sup_btm = Evaporator(model='fixed cold duty')
        self.ihx_btm = PlateHeatExchanger(self.dT_min_ihx)
        self.comp_btm = Compressor(motor_eff=self.eff_motor, model_index=1)
        self.gas_cooler_cas = GasCooler(HotFluid=fluid=self.refrigerant_btm, ColdFluid=fluid=self.refrigerant_top, dT_GC_min_tar=self.dT_GC_min_tar)
        self.comp_top = Compressor(motor_eff=self.eff_motor, model_index=1) 
        self.gas_cooler_air = GasCooler(HotFluid=fluid=self.refrigerant_top, ColdFluid=self.AIR, dT_GC_min_tar=self.dT_GC_min_tar)
        self.ihx_top = PlateHeatExchanger(self.dT_min_ihx, model='dTmin_terminal')
        self.ihx_tb = PlateHeatExchanger(self.dT_min_ihx, model='dTmin_terminal')
        self.exp_valve_top = ExpansionValve()
        self.evaporator_top = Evaporator(model='fixed cold duty')
        self.evap_sup_top = Evaporator(model='fixed cold duty')
        
    def model_initialisation(self):
        self.converged = False
        self.COP = 0.0        
        
        self.dT_ihx_top = self.dT_ihx_top if self.dT_ihx_top > 0.0 else 0.0
        self.dT_ihx_btm = self.dT_ihx_btm if self.dT_ihx_btm > 0.0 else 0.0
        self.dT_ihx_tb = self.dT_ihx_tb if self.dT_ihx_tb > 0.0 else 0.0

        self.T_gc_top_out = max(self.T_air_inlet + self.dT_GC_min_tar, self.T_gc_top_out)
        
        self.determine_HP_structure()
        self.init_streams()
        self.init_logic_blocks()
        self.init_solver_pathway()

    def determine_HP_structure(self):
        #Binary variables for the heat pump structure
        self.ihx_top_on = 1    if self.dT_ihx_top < self.T_gc_top_out else 0
        self.ihx_btm_on = 1   if self.dT_ihx_btm > 0.0 else 0
        self.ihx_tb_on = 1   if self.dT_ihx_tb > 0.0 else 0

    def init_streams(self):
        ###Bottom cycle###
        
        #Bottom evaporator
        self.evaporator_btm.source_in = MassStream( name='source In', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6)
        self.evaporator_btm.source_out = MassStream( name='source Out', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6)
        self.evaporator_btm.sink_out = MassStream( 'Evap_btm_out', fluid=self.refrigerant_btm, temperature=self.Tsat_evap[0], quality = 1 )
        
        #Bottom evaporator superheater
        self.evap_sup_btm.source_in = MassStream( name='source In', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6)
        self.evap_sup_btm.source_out = MassStream( name='source Out', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6)
        self.evap_sup_btm.sink_in = self.evaporator_btm.sink_out
        self.evap_sup_btm.sink_out = MassStream( 'IHX_in', fluid=self.refrigerant_btm, pressure=self.evaporator_btm.sink_out.pressure, temperature=self.T_ihx_btm_in )


        #Bottom internal heat exchanger - sink side
        self.ihx_btm.sink_in = self.evap_sup_btm.sink_out
        if self.ihx_btm_on == 1:
            self.ihx_btm.sink_out = MassStream( 'PHX1_in', fluid=self.refrigerant_btm, pressure=self.evaporator_btm.sink_out.pressure, temperature=self.T_ihx_btm_in + self.dT_ihx_btm)
        else:   
            self.ihx_btm.sink_out = self.ihx_btm.sink_in

        #Top to bottom internal heat exchanger - sink side
        self.ihx_tb.sink_in = self.ihx_btm.sink_out
        if self.ihx_tb_on == 1: 
            self.ihx_tb.sink_out = MassStream( 'Comp_btm_in', fluid=self.refrigerant_btm, pressure=self.evaporator_btm.sink_out.pressure, temperature=self.ihx_tb.sink_in.temperature + self.dT_ihx_tb )
        else:
            self.ihx_tb.sink_out = self.ihx_tb.sink_in

        #Bottom compressor
        self.comp_btm.flow_in = self.ihx_tb.sink_out
        P_comp_btm_out_est = get_property('P', self.refrigerant_btm, 'TS', self.T_air_tar - (self.dT_ihx_btm + self.dT_ihx_top + self.dT_ihx_tb) / 3, self.comp_btm.flow_in.entropy, 0) * 1000        
        self.comp_btm.flow_out = MassStream( 'Comp_btm_out', fluid=self.refrigerant_btm, pressure=P_comp_btm_out_est, temperature=(self.T_air_tar - (self.dT_ihx_btm + self.dT_ihx_top + self.dT_ihx_tb) / 3))
        
        #Cascade gas cooler - source side
        self.gas_cooler_cas.source_in = self.comp_btm.flow_out
        self.gas_cooler_cas.source_out = MassStream( 'Gas_cooler_out', fluid=self.refrigerant_btm, pressure=self.gas_cooler_cas.source_in.pressure, temperature=(self.evaporator_btm + 20))

        #Bottom internal heat exchanger - source side
        self.ihx_btm.source_in = self.gas_cooler_cas.source_out
        if self.ihx_btm_on == 1:
            self.ihx_btm.source_out = MassStream( 'Comp2_in', fluid=self.refrigerant_btm, pressure=self.ihx_btm.sink_in.pressure, temperature=(self.evaporator_btm + 20))
        else:
            self.ihx_btm.source_out = self.gas_cooler_cas.source_out
            
        #Bottom expansion valve
        self.exp_valve_btm.flow_in = self.ihx_btm.source_out
        self.exp_valve_btm.flow_out = MassStream( 'Evap_btm_in', fluid=self.refrigerant_btm, pressure=self.evaporator_btm.sink_out.pressure, enthalpy=self.exp_valve_btm.flow_in.enthalpy )
        self.evaporator_btm.sink_in = self.exp_valve_btm.flow_out

        ###Top cycle###

        #Top evaporator
        self.evaporator_top.source_in = MassStream( name='source In', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6)
        self.evaporator_top.source_out = MassStream( name='source Out', fluid=self.source, pressure=101, temperature=self.T_source_in, mass_flow=10^6)
        self.evaporator_top.sink_out = MassStream( 'Evap_top_out', fluid=self.refrigerant_top, temperature=self.Tsat_evap[1], quality = 1 )
        self.evap_sup_top.sink_in = self.evaporator_top.sink_out
        self.evap_sup_top.sink_out = MassStream( 'IHX_in', fluid=self.refrigerant_top, pressure=self.evaporator_top.sink_out.pressure, temperature=(self.Tsat_evap[1] + self.dT_ihx_top_in) )

        #Top internal heat exchanger - sink side
        self.ihx_top.sink_in = self.evap_sup_top.sink_out
        if self.ihx_top_on == 1:
            self.ihx_top.sink_out = MassStream( 'PHX1_in', fluid=self.refrigerant_top, pressure=self.evaporator_top.sink_out.pressure, temperature=(self.evap_sup_top.sink_out.temperature + self.dT_ihx_top) )
        else:   
            self.ihx_top.sink_out = self.ihx_top.sink_in

        #Cascade gas cooler - sink side
        self.gas_cooler_cas.sink_in = self.ihx_top.sink_out
        self.gas_cooler_cas.sink_out = MassStream( 'Gas_cooler_out', fluid=self.refrigerant_top, pressure=self.gas_cooler_cas.source_in.pressure, temperature=(self.gas_cooler_cas.source_in.temperature - 20))

        #Top compressor
        self.comp_top.flow_in = self.gas_cooler_cas.sink_out
        P_comp_top_out_est = get_property('P', self.refrigerant_top, 'TS', self.T_air_tar + self.dT_GC_min_tar, self.comp_top.flow_in.entropy, 0) * 1000        
        self.comp_top.flow_out = MassStream( 'Comp_top_out', fluid=self.refrigerant_top, pressure=P_comp_top_out_est, temperature=(self.T_air_tar + self.dT_GC_min_tar * 2))
        
        #Top gas cooler for air heating
        self.gas_cooler_air.source_in = self.comp_top.flow_out
        self.gas_cooler_air.source_out = MassStream( 'Gas_cooler_out', fluid=self.refrigerant_top, pressure=self.gas_cooler_cas.source_in.pressure, temperature=self.T_gc_top_out )
        self.gas_cooler_air.sink_in = MassStream( 'AirIn1', self.AIR, pressure=101, temperature=self.T_air_inlet )
        self.gas_cooler_air.sink_out = MassStream( 'AirOut1', self.AIR, pressure=101, temperature=self.T_air_tar )

        #Top internal heat exchanger - source side
        self.ihx_top.source_in = self.gas_cooler_air.source_out
        if self.ihx_top_on == 1:
            self.ihx_top.source_out = MassStream( 'Comp2_in', fluid=self.refrigerant_top, pressure=self.ihx_top.sink_in.pressure, temperature=(self.evaporator_top + 20))
        else:
            self.ihx_top.source_out = self.ihx_top.source_in

        #Top to bottom internal heat exchanger - source side
        self.ihx_tb.source_in = self.ihx_top.source_out
        if self.ihx_tb_on == 1: 
            self.ihx_tb.source_out = MassStream( 'Comp_top_in', fluid=self.refrigerant_top, pressure=self.evaporator_top.sink_out.pressure, temperature=(self.evaporator_top + 20) )
        else:
            self.ihx_tb.source_out = self.ihx_tb.source_in  
        
        #Top - expansion valve
        self.exp_valve_top.flow_in = self.ihx_tb.source_out
        self.exp_valve_top.flow_out = MassStream( 'Evap_top_in', fluid=self.refrigerant_top, pressure=self.evaporator_top.sink_out.pressure, enthalpy=self.exp_valve_top.flow_in.enthalpy )
        self.evaporator_top.sink_in = self.exp_valve_top.flow_out

    def init_logic_blocks(self):
        #Set special relationships between points in the cycle 
        if self.phx3_on == 1:
            self.set_phx3_dTmax = SetRelationship( source_obj=self.phx3.sink_in, target_obj=[self.phx3.source_in], c=self.dT_phx3_max, prop='temperature', prop2='pressure' ) 
        
        if self.phx2_on == 1:    
            self.set_phx2_dTmax = SetRelationship( source_obj=self.phx2.sink_in, target_obj=[self.phx2.source_in], c=self.dT_phx2_max, prop='temperature', prop2='pressure' ) 
        
        if self.phx1_on == 1:
            self.set_phx1_dTmax = SetRelationship( source_obj=self.phx1.sink_in, target_obj=[self.phx1.source_in], c=self.dT_phx1_max, prop='temperature', prop2='pressure' ) 

        self.set_P_ints = CustomRelationship( source_obj=[self.n_comp, self.comp1.flow_in, self.comp3.flow_out, self.comp2_power_frac, self.comp3_power_frac], target_obj=[self.comp1.flow_out, self.comp2.flow_out], relation=1 )
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


    def init_solver_pathway(self):
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

            if r[1].converged and self.converged and abs(self.gas_cooler4.sink_out.temperature - self.T_air_tar) < 0.05:
                self.COP = (self.gas_cooler4.duty + self.gas_cooler3.duty + self.gas_cooler2.duty + self.gas_cooler1.duty) / (self.comp1.power + self.comp2.power + self.comp3.power)
                        
        if self.COP == 0:
            self.converged = False          
        
        self.summarise_results()


    def set_pressure_bounds(self):
        P_gc_est = get_property('P', self.comp1.flow_in.fluid, 'TS', self.T_air_tar + (self.dT_GC_min_tar - self.dT_phx1_max + self.dT_phx2_max + self.dT_phx3_max) / 2, self.comp1.flow_in.entropy, 0) * 1000

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

        if self.P_low_bound < self.comp1.flow_in.pressure:
            #Check the low pressure bound remains above the evaporator pressure
            self.P_low_bound = self.comp1.flow_in.pressure
            dt_neg = self.state_solve(self.P_low_bound)
            if self.P_upp_bound < self.P_low_bound:
                #Check the upper pressure bound stays above the low pressure bound
                dt_neg = dt_plus

        return -dt_plus * dt_neg    
   

    def state_solve(self, P_comp_out):
        self.converged = False
        self.comp3.flow_out.pressure = P_comp_out
        i = 0
        h_exp_in = 0.0
        m_air_out = 0.0
        dP_gc = 0.0
        duty_comp = 0.0

        while ( 
                (
                    i < 5 
                    or abs(self.exp_valve.flow_in.enthalpy - h_exp_in) > self.tolerance * h_exp_in 
                    or abs(self.gas_cooler4.sink_out.mass_flow - m_air_out) > self.tolerance * m_air_out 
                    or abs(self.gas_cooler4.source_in.pressure - self.gas_cooler1.source_out.pressure - dP_gc) > self.tolerance * 10 * dP_gc
                    or abs(self.comp1.power + self.comp2.power + self.comp3.power - duty_comp) > self.tolerance * 10 * duty_comp
                ) 
                and i < 100
              ): #Check for model convergance
            
            #Saves critical pressure and temperature properties to check for solution convergance
            h_exp_in = self.exp_valve.flow_in.enthalpy
            m_air_out = self.gas_cooler4.sink_out.mass_flow
            dP_gc = self.gas_cooler4.source_in.pressure - self.gas_cooler1.source_out.pressure
            duty_comp = self.comp1.power + self.comp2.power + self.comp3.power
            i += 1

            #Solve the performance of each component and update each special relationship
            for obj in self.OP:
                obj.solve()

            if i < 100:            
                self.converged = True
            
        return self.gas_cooler4.sink_out.temperature - self.T_air_tar


    def summarise_results(self):
        self.output = {
            #Key results
            'COP:' : self.COP,
            'Tol:' : abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature) / self.T_air_tar,
            'Converged:' : self.converged,
            #Parameters
            'refrigerant:' : self.refrigerant,
            'dT_min Gas Cooler:' : self.dT_GC_min_tar,
            'T_air_tar:' : self.T_air_tar,
            'T_air_inlet:' : self.gas_cooler1.sink_in.temperature,
            'Dryer exh air T:' : self.evaporator.source_in.temperature,
            'n_comp:' : self.n_comp,
            #Additional results
            'Actual air T:' : self.gas_cooler4.sink_out.temperature,            
            'Evap T:' : self.Tsat_evap,
            'Evap P:' : self.evaporator.sink_in.pressure,
            'Comp 1 P out:' : self.comp1.flow_out.pressure,
            'Comp 2 P out:' : self.comp2.flow_out.pressure,
            'Comp 3 P out:' : self.comp3.flow_out.pressure,
            'dT phx max 1:' : self.dT_phx1_max,
            'dT phx max 2:' : self.dT_phx2_max,
            'dT phx max 3:' : self.dT_phx3_max,
            'PR 1:' : self.comp1.pressure_ratio,
            'PR 2:' : self.comp2.pressure_ratio,
            'PR 3:' : self.comp3.pressure_ratio,
            'Overall P ratio:' : self.gas_cooler4.source_in.pressure / self.evaporator.sink_in.pressure,
            'Gas cooler(s) inlet T:' : self.gas_cooler4.source_in.temperature,
            'Gas cooler(s) outlet T:' : self.gas_cooler1.source_out.temperature
        }
    