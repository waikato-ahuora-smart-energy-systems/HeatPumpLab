__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from globalFunctions import REFPROP
from globalVariables import AREA_CALC

class PlateHeatExc:
    
    def __init__(self, HX_Eff):
        self.HotIn = None
        self.HotOut = None
        self.ColdIn = None
        self.ColdOut = None
        self.HX_Eff = HX_Eff
        self.NTU = 0
        self.c = 0
        self.U = 0
        self.Number_of_Plates = 0
        self.Area = 0
        self.Duty = 0
        self.dP_h = 0
        self.dP_c = 0
        #Plate heat exchanger model constants
        self.PHX_thick = 0.005
        self.PHX_thermalcond = 0.0249
        self.PHX_flowgap = 0.00225
        self.PHX_D_hyd = 0.003913
        self.PHX_LtoW = 2
        self.PHX_Length = 0.5
        self.PHX_Factor = 1.15


    #Internal heat exchanger analysis to determine duty
    def solve(self):
        self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.HotOut.PropsIn_mPH()
        self.ColdOut.PropsIn_mPH() 
        
        self.dT_max = abs( self.HotIn.Temperature - self.ColdIn.Temperature )

        if self.dT_max > 5:
            self.dh_h_max = self.HotIn.Enthalpy - REFPROP('H', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, self.HotIn.Temperature - self.dT_max * self.HX_Eff, 0)
            self.dh_c_max = REFPROP('H', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, self.ColdIn.Temperature + self.dT_max * self.HX_Eff, 0) - self.ColdIn.Enthalpy
        else:
            self.dh_h_max = 0
            self.dh_c_max = 0

        self.dh_PHX = min(self.dh_h_max, self.dh_c_max)
        self.HotOut.Enthalpy = self.HotIn.Enthalpy - self.dh_PHX
        self.HotOut.PropsIn_mPH()
        self.ColdOut.Enthalpy = self.ColdIn.Enthalpy + self.dh_PHX
        self.ColdOut.PropsIn_mPH()
        self.Duty = self.ColdOut.Heat_Flow - self.ColdIn.Heat_Flow
        '''put in the area function'''
        #if AREA_CALC and self.Duty > 0 and self.HotIn.Temperature > self.ColdOut.Temperature and self.HotOut.Temperature > self.ColdIn.Temperature:          
            #self.Calc_Area()

        #Adjust the stream pressures
        self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
        self.HotOut.PropsIn_mPH()
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.ColdOut.PropsIn_mPH()       
    
    #Determines the design of the plate heat exchanger including thermal and hydralic performance
    def area_solve(self):
        if AREA_CALC and self.Duty > 0 and self.HotIn.Temperature > self.ColdOut.Temperature and self.HotOut.Temperature > self.ColdIn.Temperature:
            Ch = self.Duty / abs( self.HotIn.Temperature - self.HotOut.Temperature )
            Cc = self.Duty / abs( self.ColdOut.Temperature - self.ColdIn.Temperature )
            self.c = min( Ch, Cc ) / max( Ch, Cc )
            
            if not ( self.c == 1 ):
                self.NTU = 1 / ( 1 - self.c )  * math.log(( 1 - self.HX_Eff * self.c )  /  ( 1 - self.HX_Eff ))
            else:
                self.NTU = self.HX_Eff /  ( 1 - self.HX_Eff )
            
            self.UA = self.NTU * min( Ch, Cc )
            
            self.ColdIn.PropsExtend()
            self.ColdOut.PropsExtend()
            self.HotIn.PropsExtend()
            self.HotOut.PropsExtend()
            '''
            k_h = ( REFPROP('TCX', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, ( self.HotIn.Temperature + self.HotOut.Temperature ) / 2, 0) ) / 1000
            Pr_h = REFPROP('PRANDTL', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, ( ( self.HotIn.Temperature + self.HotOut.Temperature ) / 2 ), 0)
            Vis_h = ( REFPROP('VIS', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, ( self.HotIn.Temperature + self.HotOut.Temperature ) / 2, 0) ) * ( 10 ** - 6 )
            
            k_c = ( REFPROP('TCX', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, ( self.ColdIn.Temperature + self.ColdOut.Temperature ) / 2, 0) ) / 1000
            Pr_c = REFPROP('PRANDTL', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, ( ( self.ColdIn.Temperature + self.ColdOut.Temperature ) / 2 ), 0)
            Vis_c = ( REFPROP('VIS', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, ( self.ColdIn.Temperature + self.ColdOut.Temperature ) / 2, 0) ) * ( 10 ** - 6 )
            '''
            Wp = 0.5 * self.PHX_Length
            
            self.U = 1.5 #Intitial estimate
            self.Area_per_plate = self.PHX_Length * self.PHX_Factor * Wp
            Nt = self.UA / ( self.U * self.Area_per_plate )
            Nt = 10 #initial guess for # plates
            Nt_0 = 0
            i = 0
            while Nt != Nt_0 and i < 100:
                Nt_0 = Nt
                G = self.HotIn.Mass_Flow /  ( Wp * (Nt - 1) * self.PHX_flowgap )

                #RE_h = ( G * self.PHX_D_hyd )  /  ( Vis_h )
                RE_h = (G * self.PHX_D_hyd )  /  ( self.HotIn.Viscosity)
                #Nu_h = 0.02 * ( Pr_h ** 0.67 ) * RE_h ** 0.4
                Nu_h = 0.02 * ( self.HotIn.PRANDTL ** 0.67 ) * RE_h ** 0.4
                #HTC_h = ( k_h * Nu_h ) / self.PHX_D_hyd
                HTC_h = ( self.HotIn.Thermal_conductivity * Nu_h ) / self.PHX_D_hyd

                #RE_c = ( G * self.PHX_D_hyd )  /  ( Vis_c )
                RE_c = ( G * self.PHX_D_hyd )  /  ( self.ColdIn.Viscosity)
                #Nu_c = 0.02 * ( Pr_c ** 0.67 ) * RE_c ** 0.4
                Nu_c = 0.02 * ( self.ColdIn.PRANDTL ** 0.67 ) * RE_c ** 0.4
                #HTC_c = ( k_c * Nu_c ) / self.PHX_D_hyd
                HTC_c = ( self.ColdIn.Thermal_conductivity * Nu_c ) / self.PHX_D_hyd

                #Assuming no fouling for now (new HX)
                self.U = 1 / ( ( 1 / HTC_h ) + ( 1 / HTC_c ) + ( self.PHX_thick / self.PHX_thermalcond ) )
                Area_min_req = self.UA / self.U

                Nt = math.trunc(Area_min_req / self.Area_per_plate) + 1                                
                i += 1 

            self.Number_of_Plates = Nt
            self.Area = self.Area_per_plate * Nt
            
            #Calculates hotside pressure drop
            self.HotIn.PropsExtend()
            #Density_h = REFPROP('D', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, ( ( self.HotIn.Temperature + self.HotOut.Temperature )  / 2 ), 0)
            den_avg_h = (self.HotIn.Density + self.HotOut.Density) /2
            #f_h = 1.22 /  ( RE_h ** 0.252 )
            f_h = 38 /  ( RE_h )
            #self.dP_h = ( ( 2 * f_h * self.PHX_Length *  ( G ** 2 ) )  /  ( Density_h * self.PHX_D_hyd ) )  / 1000
            self.dP_h = ( ( 2 * f_h * self.PHX_Length *  ( G ** 2 ) )  /  ( den_avg_h* self.PHX_D_hyd ) )  / 1000
            
            #Calculates coldside pressure drop
            #Density_c = REFPROP('D', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, ( ( self.ColdIn.Temperature + self.ColdOut.Temperature )  / 2 ), 0)
            den_avg_c = (self.ColdIn.Density + self.ColdOut.Density) /2
            #f_c = 1.22 /  ( RE_c ** 0.252 )
            f_c = 38 /  ( RE_c)
            #self.dP_c = ( ( 2 * f_c * self.PHX_Length *  ( G ** 2 ) )  /  ( Density_c * self.PHX_D_hyd ) )  / 1000
            self.dP_c = ( ( 2 * f_c * self.PHX_Length *  ( G ** 2 ) )  /  ( den_avg_c * self.PHX_D_hyd ) )  / 1000         
        

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