from __future__ import division
__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math


class CrossFlow_Element:

    def __init__(self, Face_Velocity, sigma, Hydralic_diameter, Fin_conductivity, Fin_thickness, Fin_diameter, Fin_area_per_total_area, Tube_O_Diameter, Tube_I_Diameter, Tube_Spacing, Tube_Conductivity): 
        self.HotIn = None
        self.HotOut = None
        self.ColdIn = None 
        self.ColdOut = None
        self.Face_Velocity = Face_Velocity
        
        #Constants air side 
        self.sigma = sigma
        self.Hydralic_diameter = Hydralic_diameter
        self.Fin_conductivity = Fin_conductivity
        self.Fin_thickness = Fin_thickness
        self.Fin_diameter = Fin_diameter
        self.Fin_area_per_total_area = Fin_area_per_total_area
        
        #Constants refrigerant side
        self.Tube_O_Diameter = Tube_O_Diameter
        self.Tube_I_Diameter = Tube_I_Diameter
        self.Tube_Spacing = Tube_Spacing
        self.Tube_Conductivity = Tube_Conductivity


    def solve(self, HT_Area_per_Element = 0, Tube_Cross_Area = 0):
        self.HT_Area_per_Element = HT_Area_per_Element
        self.Tube_Cross_Area = Tube_Cross_Area
        #self.ColdIn.PropsIn_mPH() #test
        self.ColdIn.PropsExtend()
        self.HotOut.PropsExtend()
        self.HotOut.PropsIn_mPH
  
        '''HEAT TRANSFER''' 
        #Calculate h for the air-side        
        self.G_c = self.Face_Velocity * self.ColdIn.Density / self.sigma
        self.RE_c = (self.G_c * self.Hydralic_diameter) / self.ColdIn.Viscosity
        self.j_c = 0.10024987 * self.RE_c ** (-0.31046686) #Correlation for j
        self.h_c = 1000 * self.j_c * ((self.G_c * self.ColdIn.Heat_Capacity)/(self.ColdIn.PRANDTL ** (2/3)))  #W/m^2C        
        self.Fin_efficiency = math.tanh((math.sqrt(2 * self.h_c / self.Fin_conductivity / self.Fin_thickness) * (self.Fin_diameter - self.Fin_thickness) / 2) / (math.sqrt (2 * self.h_c / self.Fin_conductivity / self.Fin_diameter) * (self.Fin_diameter - self.Fin_thickness) / 2))
        self.n0 = 1 - self.Fin_area_per_total_area * (1 - self.Fin_efficiency)
        self.h_c_ave = self.h_c * self.n0
        
        #Calculate h for the refrigerant-side    
        self.G_h = self.HotOut.Mass_Flow / self.Tube_Cross_Area
        self.RE_h = (self.G_h * self.Tube_I_Diameter * 2) / self.HotOut.Viscosity # = 4rG/u
        ''' supercritical correlation?'''
        self.HotOut.PropsIn_mPH
        self.ColdOut.PropsIn_mPH
        if self.HotOut.PRANDTL < 0:
            pass
        else:
            self.Nu_h = self.Nu_Bishop(Re = self.RE_h, Pr = self.HotOut.PRANDTL, D=self.Tube_I_Diameter) #correlation for NU 
            #self.Nu_h = 0.021 * (self.RE_h ** 0.8) * self.HotOut.PRANDTL ** 0.33 #correlation for NU 
            self.h_h = (self.HotOut.Thermal_conductivity * self.Nu_h) / (2 * self.Tube_I_Diameter) #W/m2C = Nuk/4r
            self.h_h_ave = self.h_h * (1 - self.Fin_area_per_total_area)
            
            #Calculate tube heat transfer resistance
            self.R_tube = (self.Tube_O_Diameter * math.log(self.Tube_O_Diameter / self.Tube_I_Diameter)) / (2 * self.Tube_Conductivity * (1 - self.Fin_area_per_total_area))

            #Calculate U value
            self.U = ((1 / self.h_c_ave + self.R_tube + 1 / self.h_h_ave) ** -1)/ 1000 #kW/m2C
            
            #Calculate effectiveness-NTU method for a crossflow heat exchanger element
            self.mcp_h = self.HotOut.Heat_Capacity * self.HotOut.Mass_Flow
            self.mcp_c = self.ColdIn.Heat_Capacity * self.ColdIn.Mass_Flow

            self.NTU_h = (self.U * self.HT_Area_per_Element) / self.mcp_h
            self.R_h = self.mcp_h / self.mcp_c
            self.P_h = self.CrossflowUnmixedEff(self.NTU_h, self.R_h)
            RIn_T_est = (self.P_h * self.ColdIn.Temperature - self.HotOut.Temperature) / (self.P_h - 1)
            
            #Calculate duty and unknown enthalpies of the element 
            self.Duty = self.P_h * self.mcp_h * (RIn_T_est - self.ColdIn.Temperature)
            self.HotIn.Enthalpy = self.HotOut.Enthalpy + self.Duty / self.HotOut.Mass_Flow
            self.HotIn.Pressure = self.HotOut.Pressure
            self.HotIn.PropsIn_mPH()
            #self.HotIn.PropsIn_mPT()
            self.HotIn.PropsExtend()
            self.ColdOut.Enthalpy = self.ColdIn.Enthalpy + self.Duty / self.ColdIn.Mass_Flow
            self.ColdOut.Pressure = self.ColdIn.Pressure
            self.ColdOut.PropsIn_mPH()
            self.ColdOut.PropsExtend()
        
            '''PRESSURE DROP per element'''
            #Calculate pressure drop refrigerant side for each segment 
            self.f_h = (0.790 * math.log( self.RE_h ) - 1.64) ** -2
            ref_i = ((self.G_h ** 2) /(2 * self.HotIn.Density)) 
            ref_o = ((self.G_h ** 2) / (2 * self.HotOut.Density + 1))
            self.dp_h = ((self.f_h * (self.G_h ** 2) * self.Tube_I_Diameter) / (2 * self.HotIn.Density * self.Tube_I_Diameter)) + ref_i  - ref_o 

            #Calculate pressure drop air side for each segment 
            self.f_c = 0.39301916 * self.RE_c **(-0.24867009) #Correlation for f
            air_i = ((self.G_c ** 2) /(2 * self.ColdIn.Density))
            air_o = ((self.G_c ** 2) / (2 * self.ColdOut.Density))
            self.dp_c = ((self.f_c * (self.G_c ** 2) * self.Hydralic_diameter) / (2 * self.ColdIn.Density * self.Hydralic_diameter)) + air_i - air_o 

            #Calculate pressures in/out of element
            self.HotIn.Pressure = self.HotOut.Pressure + self.dp_h
            self.HotIn.PropsIn_mPH()
            self.ColdOut.Pressure = self.ColdIn.Pressure - self.dp_c
            self.ColdOut.PropsIn_mPH()
        

    def CrossflowUnmixedEff(self, NTU, R, i = 1, j = 1, Pn = 0, Pn_0 = 0):
        if j < i + 1:
            dPn = R ** i / math.factorial(i + 1) * ( i - j + 1 ) / math.factorial(j) * NTU ** ( float( i + j ) )
            Pn = self.CrossflowUnmixedEff(NTU, R, i, j + 1, Pn + dPn, Pn_0) if dPn > 10 ** -6 else self.CrossflowUnmixedEff(NTU, R, i, i + 1, Pn + dPn, Pn_0) 
        elif i < 20 and Pn - Pn_0 > 10 ** -6:
            Pn = self.CrossflowUnmixedEff(NTU, R, i + 1, 1, Pn, Pn)
        else:
            Pn = 1 - math.exp(0-NTU) - math.exp(-(1 + R) * NTU) * Pn
        return Pn

    
    def Nu_Bishop(self,Re = 0, Pr = 0, rho_w=None, rho_b=None, D=None, x=None):
        Nu = 0.0069*Re**0.9*Pr**0.66
        '''
        if rho_w is not None and rho_b is not None:
            Nu *= (rho_w/rho_b)**0.43
        if D is not None and x is not None:
            Nu *= (1 + 2.4*D/x)'''
        if Pr < 0:
            print(Pr, Re, self.HotIn.Fluid, self.HotIn.Pressure, self.HotIn.Temperature)
        return Nu
        