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
from globalFunctions import REFPROP
from scipy import interpolate
from scipy import optimize
import numpy as np


class HeatPump:

    def __init__(self, Refrigerant='Carbon Dioxide', T_air_tar=200.0, dT_GC_min_tar=20.0, n_comp = 3, dT_phx1_max=20.0, dT_phx2_max=0.0, dT_phx3_max=0.0, 
                 PR_offset1=0.0, PR_offset2=0.0, T_air_inlet=15.0, T_air_exh=70.0, RH_air_exh = 0.15, T_gc_out_min=0.0, T_ihx_in_max=46.0, eff_motor=0.95, eff_phx = 0.9):

        #Continuous variables for the heat pump operation
        self.Refrigerant = Refrigerant
        self.T_air_tar = T_air_tar
        self.dT_GC_min_tar = dT_GC_min_tar
        self.PR_offset1 = PR_offset1 
        self.PR_offset2 = PR_offset2
        self.n_comp = n_comp

        self.eff_phx = eff_phx
        self.eff_motor = eff_motor

        self.dT_phx1_max = dT_phx1_max if dT_phx1_max > 0.0 else 0.0
        self.dT_phx2_max = dT_phx2_max if self.n_comp > 1 and dT_phx2_max > 0.0 else 0.0
        self.dT_phx3_max = dT_phx3_max if self.n_comp > 2 and dT_phx3_max > 0.0 else 0.0

        self.T_gc_out = max(T_air_inlet + dT_GC_min_tar, T_gc_out_min)
        self.T_ihx_in = min(T_air_exh - dT_GC_min_tar, T_ihx_in_max)
        
        self.converged = False
        self.COP = 1        
        self.AIR = 'Air'
        self.T_evap_init = 15.0

        #Binary variables for the heat pump structure
        self.comp2_on = 1  if self.n_comp > 1 else 0
        self.comp3_on = 1  if self.n_comp > 2 else 0
        self.ihx_on = 1    if self.T_ihx_in < self.T_gc_out else 0
        self.phx1_on = 1   if self.dT_phx1_max > 0.0 else 0
        self.phx2_on = 1   if self.dT_phx2_max > 0.0 else 0
        self.phx3_on = 1   if self.dT_phx3_max > 0.0 else 0

        self.n_comp = 1 + self.comp2_on + self.comp3_on 
    
        #Initialisation of unit operations
        self.comp1 = Compressor(eff_motor)
        self.comp2 = Compressor(eff_motor) 
        self.comp3 = Compressor(eff_motor) 
        self.gas_cooler1 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=dT_GC_min_tar)
        self.gas_cooler2 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=dT_GC_min_tar)
        self.gas_cooler3 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=dT_GC_min_tar)
        self.gas_cooler4 = GasCooler(HotFluid=self.Refrigerant, ColdFluid=self.AIR, dT_GC_min_tar=dT_GC_min_tar)
        self.ihx = PlateHeatExc(eff_phx)
        self.phx1 = PlateHeatExc(eff_phx)
        self.phx2 = PlateHeatExc(eff_phx)        
        self.phx3 = PlateHeatExc(eff_phx)
        self.exp_valve = Expansion_Valve()
        self.evaporator = Evaporator()
        self.evap_sup = Evaporator()
        
        #Initialise values in the evaporator streams
        self.evaporator.ColdOut = Stream( 'Evap_out', Refrigerant, Temperature=self.T_evap_init, Quality=1 )
        P_low_init = self.evaporator.ColdOut.Pressure
        self.evap_sup.ColdIn = self.evaporator.ColdOut
        self.evap_sup.ColdOut = Stream( 'IHX_in', Refrigerant, Pressure=P_low_init, Temperature=self.T_ihx_in )

        #Initialise values in the compression streams
        self.ihx.ColdIn = self.evap_sup.ColdOut
        if self.ihx_on == 1:
            self.ihx.ColdOut = Stream( 'PHX1_in', Refrigerant, Pressure=P_low_init, Temperature=self.T_ihx_in + max(( self.T_gc_out - self.T_ihx_in ) * eff_phx, 0) )
        else:   
            self.ihx.ColdOut = self.ihx.ColdIn

        self.phx1.ColdIn = self.ihx.ColdOut
        if self.phx1_on == 1: #self.ihx_on == 1:
            self.phx1.ColdOut = Stream( 'Comp1_in', Refrigerant, Pressure=P_low_init, Temperature=self.phx1.ColdIn.Temperature + self.dT_phx2_max * eff_phx )
        else:
            self.phx1.ColdOut = self.phx1.ColdIn

        self.comp1.FlowIn = self.phx1.ColdOut

        P_comp_out_est = REFPROP('P', Refrigerant, 'TS', self.T_air_tar + self.dT_GC_min_tar - (dT_phx1_max + dT_phx2_max + dT_phx3_max) / 2, self.comp1.FlowIn.Entropy, 0) * 1000
        
        if self.n_comp == 1:
            P_comp1_out = P_comp_out_est
        elif self.n_comp == 2:
            P_comp1_out = P_low_init ** (1 / 2) * P_comp_out_est ** (1 / 2)
            P_comp2_out = P_comp_out_est
        elif self.n_comp == 3:
            P_comp1_out = P_low_init ** (2 / 3) * P_comp_out_est ** (1 / 3)
            P_comp2_out = P_low_init ** (1 / 3) * P_comp_out_est ** (2 / 3)
            P_comp3_out = P_comp_out_est

        self.comp1.FlowOut = Stream( 'Comp1_out', Refrigerant, Pressure=P_comp1_out, Temperature=REFPROP('T', Refrigerant, 'PS', P_comp1_out / 1000, self.comp1.FlowIn.Entropy, 0))
        
        self.phx2.ColdIn = self.comp1.FlowOut
        if self.phx2_on == 1:
            self.phx2.ColdOut = Stream( 'Comp2_in', Refrigerant, Pressure=self.phx2.ColdIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_max * eff_phx)
        else:
            self.phx2.ColdOut = self.comp1.FlowOut
        
        self.comp2.FlowIn = self.phx2.ColdOut
        if self.comp2_on == 1:
            self.comp2.FlowOut = Stream( 'Comp2_out', Refrigerant, Pressure=P_comp2_out, Temperature=REFPROP('T', Refrigerant, 'PS', P_comp2_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp2.FlowOut = self.comp2.FlowIn

        self.phx3.ColdIn = self.comp2.FlowOut
        if self.phx3_on == 1:
            self.phx3.ColdOut = Stream( 'Comp3_in', Refrigerant, Pressure=self.phx3.ColdIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_max * eff_phx)
        else:
            self.phx3.ColdOut = self.phx3.ColdIn

        self.comp3.FlowIn = self.phx3.ColdOut
        if self.comp3_on == 1:
            self.comp3.FlowOut = Stream( 'Comp3_out', Refrigerant, Pressure=P_comp3_out, Temperature=REFPROP('T', Refrigerant, 'PS', P_comp3_out / 1000, self.comp1.FlowIn.Entropy, 0))       
        else:
            self.comp3.FlowOut = self.comp3.FlowIn

        #Initialisation of gas cooler side of the cycle
        self.gas_cooler4.HotIn = self.comp3.FlowOut
        if self.comp3_on == 1 and self.phx3_on == 1:
            self.gas_cooler4.HotOut = Stream( 'GC4_out', Refrigerant, Pressure= self.gas_cooler4.HotIn.Pressure, Temperature=self.comp2.FlowOut.Temperature + self.dT_phx3_max )
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = Stream( 'GC3_in', Refrigerant, Pressure=self.phx3.HotIn.Pressure, Temperature=self.phx3.ColdIn.Temperature + self.dT_phx3_max * eff_phx)       
        else:
            self.gas_cooler4.HotOut = self.gas_cooler4.HotIn
            self.phx3.HotIn = self.gas_cooler4.HotOut
            self.phx3.HotOut = self.phx3.HotIn  

        self.gas_cooler3.HotIn = self.phx3.HotOut
        if self.comp2_on == 1 and self.phx2_on == 1:
            self.gas_cooler3.HotOut = Stream( 'GC3_out', Refrigerant, Pressure= self.gas_cooler3.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx2_max )
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = Stream( 'GC2_in', Refrigerant, Pressure=self.phx2.HotIn.Pressure, Temperature=self.phx2.ColdIn.Temperature + self.dT_phx2_max * eff_phx)
        else:
            self.gas_cooler3.HotOut = self.gas_cooler3.HotIn
            self.phx2.HotIn = self.gas_cooler3.HotOut
            self.phx2.HotOut = self.phx2.HotIn

        self.gas_cooler2.HotIn = self.phx2.HotOut
        if self.phx1_on == 1:
            self.gas_cooler2.HotOut = Stream( 'GC2_out', Refrigerant, Pressure= self.gas_cooler2.HotIn.Pressure, Temperature=self.comp1.FlowOut.Temperature + self.dT_phx1_max )
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = Stream( 'GC1_in', Refrigerant, Pressure=self.phx1.HotIn.Pressure, Temperature=self.phx1.ColdIn.Temperature + self.dT_phx1_max * eff_phx)
        else:
            self.gas_cooler2.HotOut = self.gas_cooler2.HotIn
            self.phx1.HotIn = self.gas_cooler2.HotOut
            self.phx1.HotOut = self.phx1.HotIn 

        self.gas_cooler1.HotIn = self.phx1.HotOut
        self.gas_cooler1.HotOut = Stream( 'GC1_out', Refrigerant, Pressure=self.comp2.FlowOut.Pressure, Temperature=self.T_gc_out )


        #Initialisation of internal exchanger and evaporator system
        self.ihx.HotIn = self.gas_cooler1.HotOut
        if self.ihx_on == 1:
            self.ihx.HotOut = Stream( 'Ej_pri_in', Refrigerant, Pressure=self.gas_cooler1.HotOut.Pressure, Enthalpy=self.gas_cooler1.HotOut.Enthalpy -  ( self.phx1.ColdIn.Enthalpy - self.ihx.ColdIn.Enthalpy ) )
        else:   
            self.ihx.HotOut = self.ihx.HotIn
                
        self.exp_valve.FlowIn = self.ihx.HotOut
        self.exp_valve.FlowOut = Stream( 'Evap_in', Refrigerant, Pressure=P_low_init, Enthalpy=self.exp_valve.FlowIn.Enthalpy )
        self.evaporator.ColdIn = self.exp_valve.FlowOut

        #Initialisation of inlet and exhaust air streams
        self.gas_cooler1.ColdIn = Stream( 'AirIn1', self.AIR, Pressure=101, Temperature=T_air_inlet )
        self.gas_cooler1.ColdOut = Stream( 'AirOut1', self.AIR, Pressure=101, Temperature=T_air_tar )
        self.gas_cooler2.ColdIn = self.gas_cooler1.ColdOut
        self.gas_cooler2.ColdOut = Stream( 'AirOut2', self.AIR, Pressure=101, Temperature=T_air_tar )
        self.gas_cooler3.ColdIn = self.gas_cooler2.ColdOut
        self.gas_cooler3.ColdOut = Stream( 'AirOut3', self.AIR, Pressure=101, Temperature=T_air_tar )
        self.gas_cooler4.ColdIn = self.gas_cooler3.ColdOut
        self.gas_cooler4.ColdOut = Stream( 'AirOut4', self.AIR, Pressure=101, Temperature=T_air_tar )


        #Set special relationships between points in the cycle 
        if self.phx3_on == 1:
            self.set_phx3_dTmax = Set_Relationship( source_obj=self.phx3.ColdIn, target_obj=[self.phx3.HotIn], c=self.dT_phx3_max, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx2_on == 1:    
            self.set_phx2_dTmax = Set_Relationship( source_obj=self.phx2.ColdIn, target_obj=[self.phx2.HotIn], c=self.dT_phx2_max, prop='Temperature', prop2='Pressure' ) 
        
        if self.phx1_on == 1:
            self.set_phx1_dTmax = Set_Relationship( source_obj=self.phx1.ColdIn, target_obj=[self.phx1.HotIn], c=self.dT_phx1_max, prop='Temperature', prop2='Pressure' ) 
            #self.set_phx1_dTmax = CustomRelationship( source_obj=[self.phx1.ColdIn, self.gas_cooler1.HotOut], target_obj=self.phx1.HotIn, var1=self.dT_phx2_max, relation=0 )
        
        self.set_P_ints = CustomRelationship( source_obj=[self.comp1.FlowIn, self.comp3.FlowOut, n_comp], target_obj=[self.comp1.FlowOut, self.comp2.FlowOut], var1=[self.PR_offset1, self.PR_offset2], relation=1 )
        self.set_air_massflow = CustomRelationship( source_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], target_obj=[self.gas_cooler1, self.gas_cooler2, self.gas_cooler3, self.gas_cooler4], relation=2 )
        
        #Setup heat exchanger elements
        self.gas_cooler1.init_heat_exchanger_elements()
        self.gas_cooler2.init_heat_exchanger_elements()
        self.gas_cooler3.init_heat_exchanger_elements()
        self.gas_cooler4.init_heat_exchanger_elements()

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

        self.OP.append(self.set_air_massflow)
        if self.phx3_on == 1:
            self.OP.append(self.gas_cooler4)
        if self.phx2_on == 1:
            self.OP.append(self.gas_cooler3)
        if self.phx1_on == 1:
            self.OP.append(self.gas_cooler2)        
        
        self.OP.append(self.gas_cooler1)


    def state_solve(self, P_comp_out):
        self.converged = False
        self.comp3.FlowOut.Pressure = P_comp_out
        i = 0
        h_exp_in = 0
        m_air_out = 0

        while (i < 5 or abs(self.exp_valve.FlowIn.Enthalpy - h_exp_in) > CONVERGANCE * h_exp_in 
                or abs(self.gas_cooler4.ColdOut.Mass_Flow - m_air_out) > CONVERGANCE * m_air_out) and i < 100: #Check for model convergance
            
            #Saves critical pressure and temperature properties to check for solution convergance
            h_exp_in = self.exp_valve.FlowIn.Enthalpy
            m_air_out = self.gas_cooler4.ColdOut.Mass_Flow
            i += 1

            #Solve the performance of each component and update each special relationship
            for obj in self.OP:
                obj.solve()
        
        if i < 100:            
            self.converged = True
        
        return self.gas_cooler4.ColdOut.Temperature - self.T_air_tar
            

    #Solve heat pump cycle
    def solve(self):
        #Find a lower bound pressure resulting in an air temperature less than the target
        P_low_bound = REFPROP('P', self.comp1.FlowIn.Fluid, 'TS', self.T_air_tar - (self.dT_phx1_max + self.dT_phx2_max + self.dT_phx3_max), self.comp1.FlowIn.Entropy, 0) * 1000
        i = 0 
        while self.state_solve(P_low_bound) > 0.0 and i < 20:
            P_low_bound *= 0.90
            i += 1
        
        #Find a upper bound pressure resulting in an air temperature greater than the target
        P_upp_bound = REFPROP('P', self.comp1.FlowIn.Fluid, 'TS', self.T_air_tar + self.dT_GC_min_tar, self.comp1.FlowIn.Entropy, 0) * 1000
        i = 0
        while self.state_solve(P_upp_bound) < 0.0 and i < 20:
            P_upp_bound *= 1.1
            i += 1
         
        a = self.state_solve(P_low_bound)
        b = self.state_solve(P_upp_bound)

        if a < 0.0 and b > 0.0:
            
            r = optimize.brentq(f=self.state_solve, a=P_low_bound, b=P_upp_bound, full_output=True)
            
            if r[1].converged and self.converged:
                self.COP = (self.gas_cooler4.Duty + self.gas_cooler3.Duty + self.gas_cooler2.Duty + self.gas_cooler1.Duty) / (self.comp1.Power + self.comp2.Power + self.comp3.Power )  # - self.exp_valve.Power           
            
        if self.COP == 1:
            self.converged = False          
        
        return self.COP    

    
    #Writes the outputs of the model to a nested list (2-D)
    def Write_Outputs(self, dTime = 0):
        pass
    '''
        #Summarises the properties of the refrigerant streams
        output_col1 = ['m']
        output_col2 = ['P']
        output_col3 = ['h']

        for stream in [self.evaporator.ColdOut, self.evaporator.ColdIn, self.separator.LiquidOut, self.separator.VapOut, 
                self.ihx.ColdIn, self.phx1.ColdIn, self.comp1.FlowIn, self.comp1.FlowOut, self.comp2.FlowIn, self.comp2.FlowOut, self.gas_cooler3.HotOut, 
                self.gas_cooler2.HotIn, self.gas_cooler2.HotOut, self.gas_cooler1.HotIn, self.gas_cooler1.HotOut]:
            output_col1.append(stream.Mass_Flow)
            output_col2.append(stream.Pressure)
            output_col3.append(stream.Enthalpy)

        #Summarises the performance of the heat transfer operations
        output_col1.append('U')
        output_col2.append('Area')
        output_col3.append('Duty')
        for unit_op in [self.ihx, self.phx1, self.phx2, self.gas_cooler3, self.gas_cooler2, self.gas_cooler1, self.evaporator, self.evap_sup]:
            output_col1.append(unit_op.U)
            output_col2.append(unit_op.Area)
            output_col3.append(unit_op.Duty)

        #Summarises the performance of the work transfer operations
        output_col1.append('')
        output_col2.append('Pressure Ratio')
        output_col3.append('Power')
        for unit_op in [self.comp1, self.comp2]:
            output_col1.append('')
            output_col2.append(unit_op.Pressure_Ratio)
            output_col3.append(unit_op.Power)

        #Write the time taken to solve the cycle
        output_col1.append('')
        output_col2.append('Solve timer')
        output_col3.append(dTime)

        #Combine outputs and transpose to ready export to Excel
        output = [output_col1, output_col2, output_col3]
        transposed_output = [[row[i] for row in output] for i in range(len(output[0]))]

        return transposed_output
    '''