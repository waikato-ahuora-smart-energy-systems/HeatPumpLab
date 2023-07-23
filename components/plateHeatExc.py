__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import math
from streams import get_property
from components.common.commonHeatExchanger import BaseHeatExchanger

class PlateHeatExchanger(BaseHeatExchanger):   
    def __init__(self, name='Plate heat exchanger', hx_eff=0.0, dTmin=0.0, dT_h=0.0, dT_c=0.0, dP_h=0.0, dP_c=0.0, model='dTmin_terminal', bypassed=False, converged=False, tol=1e-6):
        self.name = name
        self.model = model
        self.bypassed = bypassed
        self.converged = converged
        self.hx_eff = hx_eff
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
        raise ValueError('No detailed plate heat exchanger model defined: ' + self.name)

    def init_detailed_model(self):
        #Plate heat exchanger model constants
        self.phx_thick = 0.0005
        self.phx_thermalcond = 0.0249
        self.phx_flowgap = 0.00225
        self.phx_D_hyd = 0.003913
        self.phx_LtoW = 4
        self.phx_Length = 0.5
        self.phx_Factor = 1.15

    '''
    def area_solve(self):
        if self.duty > 0 and self.source_in.temperature > self.sink_out.temperature and self.source_out.temperature > self.sink_in.temperature:          
            Ch = self.duty / abs( self.source_in.temperature - self.source_out.temperature )
            Cc = self.duty / abs( self.sink_out.temperature - self.sink_in.temperature )
            self.c = min( Ch, Cc ) / max( Ch, Cc )
            
            if not ( self.c == 1 ):
                self.NTU = 1 / ( 1 - self.c )  * math.log(( 1 - self.hx_eff * self.c )  /  ( 1 - self.hx_eff ))
            else:
                self.NTU = self.hx_eff /  ( 1 - self.hx_eff )
            
            self.UA = self.NTU * min( Ch, Cc )
            
            self.sink_in.props_extend()
            self.sink_out.props_extend()
            self.source_in.props_extend()
            self.source_out.props_extend()

            Wp = self.phx_Length / self.phx_LtoW           
            self.U = 1.5 #Intitial estimate
            self.Area_per_plate = self.phx_Length * self.phx_Factor * Wp
            Nt = self.UA / ( self.U * self.Area_per_plate )
            Nt = 10 #initial guess for # plates
            Nt_0 = 0
            i = 0

            while Nt != Nt_0 and i < 100:
                Nt_0 = Nt
                G = self.source_in.mass_flow /  ( Wp * (Nt - 1) * self.phx_flowgap )

                #RE_h = ( G * self.phx_D_hyd )  /  ( Vis_h )
                RE_h = (G * self.phx_D_hyd )  /  ( self.source_in.viscosity)
                #Nu_h = 0.02 * ( Pr_h ** 0.67 ) * RE_h ** 0.4
                Nu_h = 0.02 * ( self.source_in.prandtl ** 0.67 ) * RE_h ** 0.4
                #HTC_h = ( k_h * Nu_h ) / self.phx_D_hyd
                HTC_h = ( self.source_in.thermal_conductivity * Nu_h ) / self.phx_D_hyd

                #RE_c = ( G * self.phx_D_hyd )  /  ( Vis_c )
                RE_c = ( G * self.phx_D_hyd )  /  ( self.sink_in.viscosity)
                #Nu_c = 0.02 * ( Pr_c ** 0.67 ) * RE_c ** 0.4
                Nu_c = 0.02 * ( self.sink_in.prandtl ** 0.67 ) * RE_c ** 0.4
                #HTC_c = ( k_c * Nu_c ) / self.phx_D_hyd
                HTC_c = ( self.sink_in.thermal_conductivity * Nu_c ) / self.phx_D_hyd              
               
                #Assuming no fouling for now (new HX)
                self.U = 1 / ( ( 1 / HTC_h ) + ( 1 / HTC_c ) + ( self.phx_thick / self.phx_thermalcond ) )
                Area_min_req = self.UA / self.U

                Nt = math.trunc(Area_min_req / self.Area_per_plate) + 1                                
                i += 1 

            self.Number_of_Plates = Nt
            self.Area = self.Area_per_plate * Nt
            
            #Calculates hotside pressure drop
            self.source_in.props_extend()
            den_avg_h = (self.source_in.density + self.source_out.density) /2
            f_h = 1.22 /  ( RE_h ** 0.252 )
            self.dP_h = ( ( 2 * f_h * self.phx_Length *  ( G ** 2 ) )  /  ( den_avg_h* self.phx_D_hyd ) )  / 1000
            
            #Calculates coldside pressure drop
            den_avg_c = (self.sink_in.density + self.sink_out.density) /2
            f_c = 1.22 /  ( RE_c ** 0.252 )
            self.dP_c = ( ( 2 * f_c * self.phx_Length *  ( G ** 2 ) )  /  ( den_avg_c * self.phx_D_hyd ) )  / 1000  
            test = 0           

    
    def area_solve(self):
        #Determines the design of the plate heat exchanger including thermal and hydralic performance
        if self.area_calc_on and self.duty > 0 and self.source_in.temperature > self.sink_out.temperature and self.source_out.temperature > self.sink_in.temperature:
            Ch = self.duty / abs( self.source_in.temperature - self.source_out.temperature )
            Cc = self.duty / abs( self.sink_out.temperature - self.sink_in.temperature )
            self.c = min( Ch, Cc ) / max( Ch, Cc )
            
            if not ( self.c == 1 ):
                self.NTU = 1 / ( 1 - self.c )  * math.log(( 1 - self.hx_eff * self.c )  /  ( 1 - self.hx_eff ))
            else:
                self.NTU = self.hx_eff /  ( 1 - self.hx_eff )
            
            self.UA = self.NTU * min( Ch, Cc )
            
            self.sink_in.props_extend()
            self.sink_out.props_extend()
            self.source_in.props_extend()
            self.source_out.props_extend()
            
            k_h = ( get_property('TCX', self.source_in.fluid, 'PT', self.source_in.pressure / 1000, ( self.source_in.temperature + self.source_out.temperature ) / 2, 0) ) / 1000
            Pr_h = get_property('prandtl', self.source_in.fluid, 'PT', self.source_in.pressure / 1000, ( ( self.source_in.temperature + self.source_out.temperature ) / 2 ), 0)
            Vis_h = ( get_property('VIS', self.source_in.fluid, 'PT', self.source_in.pressure / 1000, ( self.source_in.temperature + self.source_out.temperature ) / 2, 0) ) * ( 10 ** - 6 )
            
            k_c = ( get_property('TCX', self.sink_in.fluid, 'PT', self.sink_in.pressure / 1000, ( self.sink_in.temperature + self.sink_out.temperature ) / 2, 0) ) / 1000
            Pr_c = get_property('prandtl', self.sink_in.fluid, 'PT', self.sink_in.pressure / 1000, ( ( self.sink_in.temperature + self.sink_out.temperature ) / 2 ), 0)
            Vis_c = ( get_property('VIS', self.sink_in.fluid, 'PT', self.sink_in.pressure / 1000, ( self.sink_in.temperature + self.sink_out.temperature ) / 2, 0) ) * ( 10 ** - 6 )
            
            Wp = 0.5 * self.phx_Length
            
            self.U = 1.5 #Intitial estimate
            self.Area_per_plate = self.phx_Length * self.phx_Factor * Wp
            Nt = self.UA / ( self.U * self.Area_per_plate )
            Nt = 10 #initial guess for # plates
            Nt_0 = 0
            i = 0
            while Nt != Nt_0 and i < 100:
                Nt_0 = Nt
                G = self.source_in.mass_flow /  ( Wp * (Nt - 1) * self.phx_flowgap )

                #RE_h = ( G * self.phx_D_hyd )  /  ( Vis_h )
                RE_h = (G * self.phx_D_hyd )  /  ( self.source_in.viscosity)
                #Nu_h = 0.02 * ( Pr_h ** 0.67 ) * RE_h ** 0.4
                Nu_h = 0.02 * ( self.source_in.prandtl ** 0.67 ) * RE_h ** 0.4
                #HTC_h = ( k_h * Nu_h ) / self.phx_D_hyd
                HTC_h = ( self.source_in.thermal_conductivity * Nu_h ) / self.phx_D_hyd

                #RE_c = ( G * self.phx_D_hyd )  /  ( Vis_c )
                RE_c = ( G * self.phx_D_hyd )  /  ( self.sink_in.viscosity)
                #Nu_c = 0.02 * ( Pr_c ** 0.67 ) * RE_c ** 0.4
                Nu_c = 0.02 * ( self.sink_in.prandtl ** 0.67 ) * RE_c ** 0.4
                #HTC_c = ( k_c * Nu_c ) / self.phx_D_hyd
                HTC_c = ( self.sink_in.thermal_conductivity * Nu_c ) / self.phx_D_hyd              
               
                #Assuming no fouling for now (new HX)
                self.U = 1 / ( ( 1 / HTC_h ) + ( 1 / HTC_c ) + ( self.phx_thick / self.phx_thermalcond ) )
                Area_min_req = self.UA / self.U

                Nt = math.trunc(Area_min_req / self.Area_per_plate) + 1                                
                i += 1 

            self.Number_of_Plates = Nt
            self.Area = self.Area_per_plate * Nt
            
            #Calculates hotside pressure drop
            self.source_in.props_extend()
            #density_h = get_property('D', self.source_in.fluid, 'PT', self.source_in.pressure / 1000, ( ( self.source_in.temperature + self.source_out.temperature )  / 2 ), 0)
            den_avg_h = (self.source_in.density + self.source_out.density) /2
            #f_h = 1.22 /  ( RE_h ** 0.252 )
            f_h = 38 /  ( RE_h )
            #self.dP_h = ( ( 2 * f_h * self.phx_Length *  ( G ** 2 ) )  /  ( density_h * self.phx_D_hyd ) )  / 1000
            self.dP_h = ( ( 2 * f_h * self.phx_Length *  ( G ** 2 ) )  /  ( den_avg_h* self.phx_D_hyd ) )  / 1000
            
            #Calculates coldside pressure drop
            #density_c = get_property('D', self.sink_in.fluid, 'PT', self.sink_in.pressure / 1000, ( ( self.sink_in.temperature + self.sink_out.temperature )  / 2 ), 0)
            den_avg_c = (self.sink_in.density + self.sink_out.density) /2
            #f_c = 1.22 /  ( RE_c ** 0.252 )
            f_c = 38 /  ( RE_c)
            #self.dP_c = ( ( 2 * f_c * self.phx_Length *  ( G ** 2 ) )  /  ( density_c * self.phx_D_hyd ) )  / 1000
            self.dP_c = ( ( 2 * f_c * self.phx_Length *  ( G ** 2 ) )  /  ( den_avg_c * self.phx_D_hyd ) )  / 1000         
        '''
       