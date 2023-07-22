__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from classes.streams.properties import REFPROP
from classes.components.common.commonHeatExchanger import BaseHeatExchanger

class PlateHeatExchanger(BaseHeatExchanger):
    
    def __init__(self, Name='Plate heat exchanger', HX_Eff=0.0, dTmin=0.0, dT_h=0.0, dT_c=0.0, dP_h=0.0, dP_c=0.0, model='dTmin_terminal', bypassed=False, converged=False, tol=1e-6):
        self.Name = Name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.HX_Eff = HX_Eff
        self.dTmin = dTmin
        self.dT_h = dT_h
        self.dT_c = dT_c
        self.dP_h = dP_h
        self.dP_c = dP_c
        self.tol = tol
        self.variable_dP_model = False
        self.init_stream_connections()
        self.init_performance_metrics()

        if self.model == 'custom':
            self.init_detailed_model()
      

    def get_detailed_performance(self):
        raise ValueError('No detailed plate heat exchanger model defined: ' + self.Name)


    def init_detailed_model(self):
        #Plate heat exchanger model constants
        self.PHX_thick = 0.0005
        self.PHX_thermalcond = 0.0249
        self.PHX_flowgap = 0.00225
        self.PHX_D_hyd = 0.003913
        self.PHX_LtoW = 4
        self.PHX_Length = 0.5
        self.PHX_Factor = 1.15


    '''
    def area_solve(self):
        if self.Duty > 0 and self.HotIn.Temperature > self.ColdOut.Temperature and self.HotOut.Temperature > self.ColdIn.Temperature:          
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

            Wp = self.PHX_Length / self.PHX_LtoW           
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
            den_avg_h = (self.HotIn.Density + self.HotOut.Density) /2
            f_h = 1.22 /  ( RE_h ** 0.252 )
            self.dP_h = ( ( 2 * f_h * self.PHX_Length *  ( G ** 2 ) )  /  ( den_avg_h* self.PHX_D_hyd ) )  / 1000
            
            #Calculates coldside pressure drop
            den_avg_c = (self.ColdIn.Density + self.ColdOut.Density) /2
            f_c = 1.22 /  ( RE_c ** 0.252 )
            self.dP_c = ( ( 2 * f_c * self.PHX_Length *  ( G ** 2 ) )  /  ( den_avg_c * self.PHX_D_hyd ) )  / 1000  
            test = 0           

    
    def area_solve(self):
        #Determines the design of the plate heat exchanger including thermal and hydralic performance
        if self.area_calc_on and self.Duty > 0 and self.HotIn.Temperature > self.ColdOut.Temperature and self.HotOut.Temperature > self.ColdIn.Temperature:
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
            
            k_h = ( REFPROP('TCX', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, ( self.HotIn.Temperature + self.HotOut.Temperature ) / 2, 0) ) / 1000
            Pr_h = REFPROP('PRANDTL', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, ( ( self.HotIn.Temperature + self.HotOut.Temperature ) / 2 ), 0)
            Vis_h = ( REFPROP('VIS', self.HotIn.Fluid, 'PT', self.HotIn.Pressure / 1000, ( self.HotIn.Temperature + self.HotOut.Temperature ) / 2, 0) ) * ( 10 ** - 6 )
            
            k_c = ( REFPROP('TCX', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, ( self.ColdIn.Temperature + self.ColdOut.Temperature ) / 2, 0) ) / 1000
            Pr_c = REFPROP('PRANDTL', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, ( ( self.ColdIn.Temperature + self.ColdOut.Temperature ) / 2 ), 0)
            Vis_c = ( REFPROP('VIS', self.ColdIn.Fluid, 'PT', self.ColdIn.Pressure / 1000, ( self.ColdIn.Temperature + self.ColdOut.Temperature ) / 2, 0) ) * ( 10 ** - 6 )
            
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
        '''
       