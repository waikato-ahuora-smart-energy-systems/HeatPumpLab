__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']


import math
from streams import get_property, MassStream
from components.common.crossFlowElement import CrossFlowElement
from components.common.commonHeatExchanger import BaseHeatExchanger
from scipy.optimize import brentq


class GasCooler(BaseHeatExchanger):

    def __init__(self, Name='Gas cooler', dTmin=15.0, dP_h=0.0, dP_c=0.0, model='fixed hot temperature outlet', bypassed=False, converged=False, tol=1e-6):
        self.Name = Name
        self.dTmin = dTmin
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.tol = tol
        
        self.init_stream_connections()
        self.init_performance_metrics()

        if self.model == 'custom':
            self.init_parameters()
            self.init_nodes()


    def init_parameters(self):
        self.NTU = 0
        self.c = 0
        self.Duty = 0
        self.Area = 0
        self.converged = True

        #Initial guesses
        self.U = 0.1 # kW/[m^2 C]
        self.HX_Eff = 0.8
        self.FLength = 0.5 # m
        self.Depth = 1 # m

        #Model constants
        self.Beta = 587 # m^2/m^3
        self.N_passes = 6
        self.Face_Velocity = 4 # m/
        self.gc = 9.81 #m/s
        self.Length_Longitudinal = 0.022 #m
        self.Length_Transverse = 0.0254 #m

        #Constants air side 
        self.sigma = 0.534
        self.Hydralic_diameter = 0.003632 #m
        self.Fin_conductivity = 237 #W/m°C - Al fin 
        self.Fin_thickness = 0.0003 #m
        self.Fin_diameter = 0.0102 #m
        self.Fin_area_per_total_area = 0.913 #m^2/m^2
        
        #Constants refrigerant side
        self.Tube_O_Diameter = 0.0102 #m
        self.Tube_I_Diameter = 0.00855 #m
        self.Tube_Spacing = 0.022 #m
        self.Tube_Conductivity = 16 #W/m°C - SS tube
        #self.Tube_thickness = 0.00165 #m
        #self.Tubes_per_pass = 4 #assumption
        #self.Gap_spacing = 0.0005 #m
        
    
    def init_nodes(self, N_elements_per_pass=1):
        self.N_elements_per_pass = N_elements_per_pass

        #Initialisation of heat exchanger elements
        self.matrix_elements = [[0 for j in range(self.N_elements_per_pass)] for i in range(self.N_passes)]
        self.linear_elements = []

        step = 1
        start_val = 0
        end_val = self.N_elements_per_pass
        
        for i in range(self.N_passes):
            for j in range(start_val, end_val, step): 
                      
                self.matrix_elements[i][j] = CrossFlowElement(Face_Velocity=self.Face_Velocity, 
                                                              sigma=self.sigma, 
                                                              Hydralic_diameter=self.Hydralic_diameter, 
                                                              Fin_conductivity=self.Fin_conductivity, 
                                                              Fin_thickness=self.Fin_thickness, 
                                                              Fin_diameter=self.Fin_diameter, 
                                                              Fin_area_per_total_area=self.Fin_area_per_total_area, 
                                                              Tube_O_Diameter=self.Tube_O_Diameter, 
                                                              Tube_I_Diameter=self.Tube_I_Diameter, 
                                                              Tube_Spacing=self.Tube_Spacing, 
                                                              Tube_Conductivity=self.Tube_Conductivity)
                
                if i == 0:
                    self.matrix_elements[i][j].ColdIn = self.ColdIn # Air inlet 

                    if j == start_val:
                        self.matrix_elements[i][j].HotOut = self.HotOut # Refrigerant outlet
                    else:
                        self.matrix_elements[i][j].HotOut = self.matrix_elements[i][j - step].HotIn # Refrigerant outlet - inlet relationship
                else:
                    self.matrix_elements[i][j].ColdIn = self.matrix_elements[i - 1][j].ColdOut # Air inlet - outlet relationship

                    if j == start_val:
                        self.matrix_elements[i][j].HotOut = self.matrix_elements[i - 1][j].HotIn # Refrigerant outlet - inlet relationship
                    else:
                        self.matrix_elements[i][j].HotOut = self.matrix_elements[i][j - step].HotIn # Refrigerant outlet - inlet relationship

                self.matrix_elements[i][j].HotIn = MassStream( Name=str(i) + str(j), Fluid=self.HotIn.Fluid, underdefined=1 ) # Refrigerant inlet
                self.matrix_elements[i][j].ColdOut = MassStream( Name=str(i) + str(j), Fluid=self.ColdIn.Fluid, underdefined=1 ) # Air outlet

                self.linear_elements.append(self.matrix_elements[i][j])

            step *= -1
            k = start_val
            start_val = end_val + step
            end_val = k + step
                   

    def get_detailed_performance(self):
        raise ValueError('No detailed gas cooler model defined: ' + self.Name)


    def area_solve(self):
        if self.Duty > 0 and self.HotIn.Temperature > self.ColdOut.Temperature and self.HotOut.Temperature > self.ColdIn.Temperature:
            self.ColdIn.PropsExtend()
            self.ColdOut.PropsExtend()
            self.dT_max = ( self.HotIn.Temperature - self.ColdIn.Temperature )
            self.LMTD_linear = (( self.HotIn.Temperature - self.ColdOut.Temperature ) - ( self.HotOut.Temperature - self.ColdIn.Temperature )) / math.log(( self.HotIn.Temperature - self.ColdOut.Temperature ) / ( self.HotOut.Temperature - self.ColdIn.Temperature ))
            self.C_ratio_linear = ( self.ColdOut.Temperature - self.ColdIn.Temperature ) / ( self.HotIn.Temperature - self.HotOut.Temperature )
            self.Max_Duty = (self.HotIn.Temperature - self.HotOut.Temperature ) * (( self.HotIn.Temperature - self.ColdIn.Temperature ))

            if self.Area == 0:
                Area_upp = self.Duty / ( self.U * self.LMTD_linear * 0.5) #First approximation for area 
                Area_low = self.Duty / ( self.U * self.LMTD_linear * 1.5) #First approximation for area 
            else:
                Area_upp = self.Area * 1.1
                Area_low = self.Area * 0.9                 
            
            a = self.get_area(Area_low)
            b = self.get_area(Area_upp)

            if a < 0.0 and b > 0.0:
                
                r = brentq(f=self.get_area, a=Area_low, b=Area_upp, full_output=True, rtol=self.tol, xtol=self.tol)
                
                if r[1].converged:  
                    self.converged = True          
                    
                    #Refrigerant side pressure drop
                    
                    #self.dP_h = self.linear_elements[len(self.linear_elements) - 1].HotIn.Pressure - self.HotOut.Pressure 
                    element_pressure = self.linear_elements[len(self.linear_elements) - 1].HotIn.Pressure
                    self.dP_h = element_pressure  - self.HotOut.Pressure
                    #self.HotOut.Pressure = self.HotIn.Pressure - 0#self.dP_h
                    #Air side pressure drop 
                    
                    
                    sum_P = 0
                    sum_h = 0 
                    for element in self.matrix_elements[self.N_passes - 1]:              
                        sum_P += element.ColdOut.Pressure
                        sum_h += element.ColdOut.Enthalpy 
                    
                    self.ColdOut.Pressure = sum_P / self.N_elements_per_pass
                    self.ColdOut.Enthalpy = sum_h / self.N_elements_per_pass
                    self.dP_c = self.ColdIn.Pressure - self.ColdOut.Pressure
                    
                    '''
                    #Accounting for pressure drop refrigerant
                    self.HotOut.PropsIn_mPT()
                    #Accounting for pressure drop air
                    self.ColdOut.PropsIn_mPH()'''
                    

    def get_area(self, Area_0):
        self.Volume = Area_0 / self.Beta
        self.Face_Area = self.ColdIn.Mass_Flow / (self.ColdIn.Density * self.Face_Velocity)

        self.Depth = self.Volume / self.Face_Area
        self.Tube_Length = (self.Face_Area) ** 0.5
        self.HT_Area_per_Element = Area_0 / (self.N_passes * self.N_elements_per_pass)           
        
        self.Duct_Height = self.Face_Area ** 0.5
        self.Tube_Cross_Area = ((math.pi * (self.Tube_O_Diameter ** 2) / 4)/(self.Length_Longitudinal * self.Length_Transverse)) * (self.Duct_Height * (self.Depth / self.N_passes)) #calculated tube cross sec area using depth and Ll and Lt

        self.Area = Area_0

        duty_total = 0
        for element in self.linear_elements:
            element.solve(self.HT_Area_per_Element, self.Tube_Cross_Area)
            duty_total += element.Duty

        return self.linear_elements[len(self.linear_elements) - 1].HotIn.Temperature - self.HotIn.Temperature
