__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from streams import MassStream
from components.common.crossFlowElement import CrossFlowElement
from components.common.commonHeatExchanger import BaseHeatExchanger
from scipy.optimize import brentq

class GasCooler(BaseHeatExchanger):
    def __init__(self, source_fluid: str, sink_fluid: str, name='Gas cooler', dTmin=15.0, dP_h=0.0, dP_c=0.0, model='fixed hot temperature outlet', bypassed=False, converged=False, tol=1e-6):
        self.name = name
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.tol = tol
        self.set_port_mass_streams(source_fluid, sink_fluid)

        if self.model == 'custom':
            self.init_parameters()
            self.init_nodes()

    def init_parameters(self):
        self.NTU = 0
        self.c = 0
        self.duty = 0
        self.area = 0
        self.converged = True

        #Initial guesses
        self.U = 0.1 # kW/[m^2 C]
        self.eff = 0.8
        self.f_length = 0.5 # m
        self.depth = 1 # m

        #Model constants
        self.beta = 587 # m^2/m^3
        self.N_passes = 6
        self.face_velocity = 4 # m/
        self.gc = 9.81 #m/s
        self.length_longitudinal = 0.022 #m
        self.length_transverse = 0.0254 #m

        #Constants air side 
        self.sigma = 0.534
        self.hydralic_diameter = 0.003632 #m
        self.fin_conductivity = 237 #W/m°C - Al fin 
        self.fin_thickness = 0.0003 #m
        self.fin_diameter = 0.0102 #m
        self.fin_area_per_total_area = 0.913 #m^2/m^2
        
        #Constants refrigerant side
        self.tube_O_Diameter = 0.0102 #m
        self.tube_I_Diameter = 0.00855 #m
        self.tube_Spacing = 0.022 #m
        self.tube_Conductivity = 16 #W/m°C - SS tube
        #self.tube_thickness = 0.00165 #m
        #self.tubes_per_pass = 4 #assumption
        #self.gap_spacing = 0.0005 #m 
    
    def init_nodes(self, n_elements_per_pass=1):
        self.n_elements_per_pass = n_elements_per_pass

        #Initialisation of heat exchanger elements
        self.matrix_elements = [[0 for j in range(self.n_elements_per_pass)] for i in range(self.N_passes)]
        self.linear_elements = []

        step = 1
        start_val = 0
        end_val = self.n_elements_per_pass
        
        for i in range(self.N_passes):
            for j in range(start_val, end_val, step): 
                      
                self.matrix_elements[i][j] = CrossFlowElement(face_velocity=self.face_velocity, 
                                                              sigma=self.sigma, 
                                                              hydralic_diameter=self.hydralic_diameter, 
                                                              fin_conductivity=self.fin_conductivity, 
                                                              fin_thickness=self.fin_thickness, 
                                                              fin_diameter=self.fin_diameter, 
                                                              fin_area_per_total_area=self.fin_area_per_total_area, 
                                                              tube_O_Diameter=self.tube_O_Diameter, 
                                                              tube_I_Diameter=self.tube_I_Diameter, 
                                                              tube_Spacing=self.tube_Spacing, 
                                                              tube_Conductivity=self.tube_Conductivity)
                
                if i == 0:
                    self.matrix_elements[i][j].sink_in = self.sink_in # Air inlet 

                    if j == start_val:
                        self.matrix_elements[i][j].source_out = self.source_out # refrigerant outlet
                    else:
                        self.matrix_elements[i][j].source_out = self.matrix_elements[i][j - step].source_in # refrigerant outlet - inlet relationship
                else:
                    self.matrix_elements[i][j].sink_in = self.matrix_elements[i - 1][j].sink_out # Air inlet - outlet relationship

                    if j == start_val:
                        self.matrix_elements[i][j].source_out = self.matrix_elements[i - 1][j].source_in # refrigerant outlet - inlet relationship
                    else:
                        self.matrix_elements[i][j].source_out = self.matrix_elements[i][j - step].source_in # refrigerant outlet - inlet relationship

                self.matrix_elements[i][j].source_in = MassStream( name=str(i) + str(j), fluid=self.source_in.fluid, underdefined=1 ) # refrigerant inlet
                self.matrix_elements[i][j].sink_out = MassStream( name=str(i) + str(j), fluid=self.sink_in.fluid, underdefined=1 ) # Air outlet

                self.linear_elements.append(self.matrix_elements[i][j])

            step *= -1
            k = start_val
            start_val = end_val + step
            end_val = k + step
                   
    def get_detailed_performance(self):
        raise ValueError('No detailed gas cooler model defined: ' + self.name)

    def area_solve(self):
        if self.duty > 0 and self.source_in.temperature > self.sink_out.temperature and self.source_out.temperature > self.sink_in.temperature:
            self.dT_max = ( self.source_in.temperature - self.sink_in.temperature )
            self.LMTD_linear = (( self.source_in.temperature - self.sink_out.temperature ) - ( self.source_out.temperature - self.sink_in.temperature )) / math.log(( self.source_in.temperature - self.sink_out.temperature ) / ( self.source_out.temperature - self.sink_in.temperature ))
            self.C_ratio_linear = ( self.sink_out.temperature - self.sink_in.temperature ) / ( self.source_in.temperature - self.source_out.temperature )
            self.max_duty = (self.source_in.temperature - self.source_out.temperature ) * (( self.source_in.temperature - self.sink_in.temperature ))

            if self.area == 0:
                area_upp = self.duty / ( self.U * self.LMTD_linear * 0.5) #First approximation for area 
                area_low = self.duty / ( self.U * self.LMTD_linear * 1.5) #First approximation for area 
            else:
                area_upp = self.area * 1.1
                area_low = self.area * 0.9                 
            
            a = self.get_area(area_low)
            b = self.get_area(area_upp)

            if a < 0.0 and b > 0.0:
                
                r = brentq(f=self.get_area, a=area_low, b=area_upp, full_output=True, rtol=self.tol, xtol=self.tol)
                
                if r[1].converged:  
                    self.converged = True          
                    
                    #refrigerant side pressure drop
                    
                    #self.dP_h = self.linear_elements[len(self.linear_elements) - 1].source_in.pressure - self.source_out.pressure 
                    element_pressure = self.linear_elements[len(self.linear_elements) - 1].source_in.pressure
                    self.dP_h = element_pressure  - self.source_out.pressure
                    #self.source_out.pressure = self.source_in.pressure - 0#self.dP_h
                    #Air side pressure drop 
                    
                    
                    sum_P = 0
                    sum_h = 0 
                    for element in self.matrix_elements[self.N_passes - 1]:              
                        sum_P += element.sink_out.pressure
                        sum_h += element.sink_out.enthalpy 
                    
                    self.sink_out.pressure = sum_P / self.n_elements_per_pass
                    self.sink_out.enthalpy = sum_h / self.n_elements_per_pass
                    self.dP_c = self.sink_in.pressure - self.sink_out.pressure
                    
                    '''
                    #Accounting for pressure drop refrigerant
                    self.source_out.props_in_mPT()
                    #Accounting for pressure drop air
                    self.sink_out.props_in_mPH()'''
                    
    def get_area(self, area_0):
        self.volume = area_0 / self.beta
        self.face_area = self.sink_in.mass_flow / (self.sink_in.density * self.face_velocity)

        self.depth = self.volume / self.face_area
        self.tube_Length = (self.face_area) ** 0.5
        self.ht_area_per_element = area_0 / (self.N_passes * self.n_elements_per_pass)           
        
        self.duct_height = self.face_area ** 0.5
        self.tube_cross_area = ((math.pi * (self.tube_O_Diameter ** 2) / 4)/(self.length_longitudinal * self.length_transverse)) * (self.duct_height * (self.depth / self.N_passes)) #calculated tube cross sec area using depth and Ll and Lt

        self.area = area_0

        duty_total = 0
        for element in self.linear_elements:
            element.solve(self.ht_area_per_element, self.tube_cross_area)
            duty_total += element.duty

        return self.linear_elements[len(self.linear_elements) - 1].source_in.temperature - self.source_in.temperature
