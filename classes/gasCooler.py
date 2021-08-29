__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from globalFunctions import REFPROP
from globalVariables import CONVERGANCE, AREA_CALC
from classes.crossFlowElement import CrossFlow_Element
from classes.stream import Stream
import matplotlib.pyplot as plt
from scipy import optimize


class GasCooler():

    def __init__(self, HotFluid, ColdFluid, dT_GC_min_tar = 15):
        self.HotIn = Stream( Name='HotIn', Fluid=HotFluid, underdefined=1 )
        self.HotOut = Stream( Name='HotOut', Fluid=HotFluid, underdefined=1 )
        self.ColdIn = Stream( Name='ColdIn', Fluid=ColdFluid, underdefined=1 )
        self.ColdOut = Stream( Name='ColdOut', Fluid=ColdFluid, underdefined=1 )
        self.dT_GC_min_tar = dT_GC_min_tar

        self.NTU = 0
        self.c = 0
        self.Duty = 0
        self.Area = 0
        self.dP_refrigerant = 0
        self.dP_air = 0
        self.converged = True
        #Initial guesses
        self.U = 0.1 # kW/[m^2 C]
        self.HX_Eff = 0.8
        self.FLength = 0.5 # m
        self.Depth = 1 # m
        #Model constants
        self.Beta = 587 # m^2/m^3
        self.N_passes = 6
        #self.N_elements_per_pass = 4
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
        
    
    def init_heat_exchanger_elements(self, N_elements_per_pass = 4):
        self.N_elements_per_pass = N_elements_per_pass
        #Initialisation of heat exchanger elements
        self.matrix_elements = [[0 for j in range(self.N_elements_per_pass)] for i in range(self.N_passes)]
        self.linear_elements = []

        step = 1
        start_val = 0
        end_val = self.N_elements_per_pass
        
        for i in range(self.N_passes):
            for j in range(start_val, end_val, step): 
                      
                self.matrix_elements[i][j] = CrossFlow_Element(self.Face_Velocity, self.sigma, self.Hydralic_diameter, self.Fin_conductivity, self.Fin_thickness, self.Fin_diameter, self.Fin_area_per_total_area, self.Tube_O_Diameter, self.Tube_I_Diameter, self.Tube_Spacing, self.Tube_Conductivity)
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

                self.matrix_elements[i][j].HotIn = Stream( Name=str(i) + str(j), Fluid=self.HotIn.Fluid, underdefined=1 ) # Refrigerant inlet
                self.matrix_elements[i][j].ColdOut = Stream( Name=str(i) + str(j), Fluid=self.ColdIn.Fluid, underdefined=1 ) # Air outlet

                self.linear_elements.append(self.matrix_elements[i][j])

            step *= -1
            k = start_val
            start_val = end_val + step
            end_val = k + step
                   

    def solve(self):
        self.HotOut.Pressure = self.HotIn.Pressure -self.dP_refrigerant
        self.HotOut.PropsIn_mPT()
        self.Duty = self.HotIn.Heat_Flow - self.HotOut.Heat_Flow
        
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_air
        #self.ColdOut.PropsIn_mPH()
        self.ColdIn.Enthalpy = self.ColdOut.Enthalpy - self.Duty / self.ColdIn.Mass_Flow
        self.ColdIn.PropsIn_mPH()


    def area_solve(self):    
        if AREA_CALC and self.Duty > 0 and self.HotIn.Temperature > self.ColdOut.Temperature and self.HotOut.Temperature > self.ColdIn.Temperature:
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
            
            a = self.Calc_Area(Area_low)
            b = self.Calc_Area(Area_upp)

            if a < 0.0 and b > 0.0:
                
                r = optimize.brentq(f=self.Calc_Area, a=Area_low, b=Area_upp, full_output=True)
                
                if r[1].converged:  
                    self.converged = True          
                    
                    #Refrigerant side pressure drop
                    
                    #self.dP_refrigerant = self.linear_elements[len(self.linear_elements) - 1].HotIn.Pressure - self.HotOut.Pressure 
                    element_pressure = self.linear_elements[len(self.linear_elements) - 1].HotIn.Pressure
                    self.dP_refrigerant = element_pressure  - self.HotOut.Pressure
                    #self.HotOut.Pressure = self.HotIn.Pressure - 0#self.dP_refrigerant
                    #Air side pressure drop 
                    
                    
                    sum_P = 0
                    sum_h = 0 
                    for element in self.matrix_elements[self.N_passes - 1]:              
                        sum_P += element.ColdOut.Pressure
                        sum_h += element.ColdOut.Enthalpy 
                    
                    self.ColdOut.Pressure = sum_P / self.N_elements_per_pass
                    self.ColdOut.Enthalpy = sum_h / self.N_elements_per_pass
                    self.dP_air = self.ColdIn.Pressure - self.ColdOut.Pressure
                    
                    '''
                    #Accounting for pressure drop refrigerant
                    self.HotOut.PropsIn_mPT()
                    #Accounting for pressure drop air
                    self.ColdOut.PropsIn_mPH()'''
                    

    def Calc_Area(self, Area_0):
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


    def graph_data(self):
        H = []; T_hot = []; T_cold = []

        for i in range(51):
            dh_hot = (self.HotIn.Enthalpy - self.HotOut.Enthalpy) * i / 50
            dP_hot = (self.HotIn.Pressure - self.HotOut.Pressure) * i / 50
            T_hot.append(REFPROP('T', self.HotOut.Fluid, 'PH', (self.HotOut.Pressure + dP_hot) / 1000, self.HotOut.Enthalpy + dh_hot))
                        
            dh_cold = (self.ColdOut.Enthalpy - self.ColdIn.Enthalpy) * i / 50
            dP_cold = (self.ColdIn.Pressure - self.ColdOut.Pressure) * i / 50
            T_cold.append(REFPROP('T', self.ColdIn.Fluid, 'PH', (self.ColdIn.Pressure - dP_cold) / 1000, self.ColdIn.Enthalpy + dh_cold))

            H.append(dh_hot * self.HotOut.Mass_Flow)

        return H, T_hot, T_cold         