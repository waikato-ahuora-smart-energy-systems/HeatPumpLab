__author__ = 'Tim Walmsley'
__credits__ = ['']

from classes.streams.properties import REFPROP
from classes.streams.massStream import MassStream

class BaseHeatExc:

    def __init__(self) -> None:
        pass


    def init_stream_connections(self):
        #Set the stream connection points
        self.HotIn = MassStream()
        self.HotOut = self.HotIn
        self.ColdIn = MassStream()
        self.ColdOut = self.ColdIn    
    

    def init_performance_metrics(self):
        self.NTU = 0
        self.c = 0
        self.U = 0
        self.Area = 0
        self.Duty = 0


    def solve(self):
        if self.bypassed:
            self.HotOut = self.HotIn
            self.ColdOut = self.ColdIn
        else:
            if self.HotIn.Mass_Flow <= self.tol or self.ColdIn.Mass_Flow <= self.tol:
                self.HotOut.SetEqual(fld=self.HotIn)
                self.ColdOut.SetEqual(fld=self.ColdIn)          

            else:
                self.get_performance()


    def get_performance(self):
        if self.model == 'dTmin_terminal':
            self.get_dTmin_terminal_model()

        elif self.model == 'effectiveness':
            self.get_effectiveness_model()

        elif self.model == 'fixed cold duty':
            self.get_fixed_cold_duty_model()

        elif self.model == 'fixed hot duty':
            self.get_fixed_hot_duty_model()
        
        elif self.model == 'fixed cold temperature outlet':
            self.get_fixed_cold_T_outlet_model()
        
        elif self.model == 'fixed hot temperature outlet':
            self.get_fixed_hot_T_outlet_model()

        elif self.model == 'fixed cold delta temperature':
            self.get_fixed_cold_delta_T_model()
        
        elif self.model == 'fixed hot delta temperature':
            self.get_fixed_hot_delta_T_model()

        elif self.model == 'specialised':
            self.get_detailed_performance()

        else:
            raise ValueError('No heat exchanger model selected.')


    def get_dTmin_terminal_model(self):
        try:
            self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
            self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
            self.dH_h_max = self.HotIn.Mass_Flow * (self.HotIn.Enthalpy - REFPROP('H', self.HotOut.Fluid, 'PT', self.HotOut.Pressure / 1000, self.ColdIn.Temperature + self.dTmin, 0))
            self.dH_c_max = self.ColdIn.Mass_Flow * (REFPROP('H', self.ColdOut.Fluid, 'PT', self.ColdOut.Pressure / 1000, self.HotOut.Temperature - self.dTmin, 0) - self.ColdIn.Enthalpy )
            self.Duty = min(self.dH_h_max, self.dH_c_max)

            self.HotOut.Enthalpy = self.HotIn.Enthalpy - self.Duty / self.HotIn.Mass_Flow
            self.HotOut.PropsIn_mPH()
            self.ColdOut.Enthalpy = self.ColdIn.Enthalpy - self.Duty / self.ColdIn.Mass_Flow
            self.ColdOut.PropsIn_mPH()        
        except:
            raise ValueError('Warning: Phase change occuring in heat exchanger ' + self.Name)   

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False    


    def get_effectiveness_model(self):
        try:
            self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
            self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
            self.dT_max = abs( self.HotIn.Temperature - self.ColdIn.Temperature )
            self.dH_h_max = self.HotIn.Mass_Flow * (self.HotIn.Enthalpy - REFPROP('H', self.HotOut.Fluid, 'PT', self.HotOut.Pressure / 1000, self.HotIn.Temperature - self.dT_max * self.HX_Eff, 0))
            self.dH_c_max = self.ColdIn.Mass_Flow * (REFPROP('H', self.ColdOut.Fluid, 'PT', self.ColdOut.Pressure / 1000, self.ColdIn.Temperature + self.dT_max * self.HX_Eff, 0) - self.ColdIn.Enthalpy )
            self.Duty = min(self.dH_h_max, self.dH_c_max)

            self.HotOut.Enthalpy = self.HotIn.Enthalpy - self.Duty / self.HotIn.Mass_Flow
            self.HotOut.PropsIn_mPH()
            self.ColdOut.Enthalpy = self.ColdIn.Enthalpy - self.Duty / self.ColdIn.Mass_Flow
            self.ColdOut.PropsIn_mPH()
        
        except:
            raise ValueError('Warning: Phase change occuring in heat exchanger ' + self.Name)        

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False    

    def get_fixed_cold_duty_model(self):
        self.Duty = self.ColdOut.Heat_Flow - self.ColdIn.Heat_Flow
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.ColdOut.PropsIn_mPH()

        self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
        self.HotOut.Enthalpy = self.HotIn.Enthalpy - self.Duty / self.HotOut.Mass_Flow
        self.HotOut.PropsIn_mPH()

        if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature:
            print('Temperature cross occuring in the heat exchanger: ' + self.Name)

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False


    def get_fixed_hot_duty_model(self):
        self.Duty = self.HotIn.Heat_Flow - self.HotOut.Heat_Flow
        self.HotOut.Pressure = self.HotOut.Pressure - self.dP_h
        self.HotIn.PropsIn_mPH()

        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.ColdOut.Enthalpy = self.ColdIn.Enthalpy + self.Duty / self.ColdIn.Mass_Flow
        self.ColdOut.PropsIn_mPH()

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False


    def get_fixed_cold_T_outlet_model(self):
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.ColdOut.PropsIn_mPT()
        self.Duty = self.ColdOut.Heat_Flow - self.ColdIn.Heat_Flow
        
        self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
        self.HotOut.Enthalpy = self.HotIn.Enthalpy - self.Duty / self.HotOut.Mass_Flow
        self.HotOut.PropsIn_mPH()

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False


    def get_fixed_hot_T_outlet_model(self):
        self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
        self.HotOut.PropsIn_mPT()
        self.Duty = self.HotIn.Heat_Flow - self.HotOut.Heat_Flow
        
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.ColdOut.Enthalpy = self.ColdIn.Enthalpy + self.Duty / self.ColdIn.Mass_Flow
        self.ColdOut.PropsIn_mPH()

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False


    def get_fixed_cold_delta_T_model(self):
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.ColdOut.Temperature = self.ColdIn.Temperature + self.dT_c
        self.ColdOut.PropsIn_mPT()
        self.Duty = self.ColdOut.Heat_Flow - self.ColdIn.Heat_Flow
        
        self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
        self.HotOut.Enthalpy = self.HotIn.Enthalpy - self.Duty / self.HotOut.Mass_Flow
        self.HotOut.PropsIn_mPH()

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False 


    def get_fixed_hot_delta_T_model(self):
        self.HotOut.Pressure = self.HotIn.Pressure - self.dP_h
        self.HotOut.Temperature = self.HotIn.Temperature - self.dT_h
        self.HotOut.PropsIn_mPT()
        self.Duty = self.HotIn.Heat_Flow - self.HotOut.Heat_Flow
        
        self.ColdOut.Pressure = self.ColdIn.Pressure - self.dP_c
        self.ColdOut.Enthalpy = self.ColdIn.Enthalpy + self.Duty / self.ColdIn.Mass_Flow
        self.ColdOut.PropsIn_mPH()

        self.converged = True if self.HotOut.Temperature < self.ColdIn.Temperature or self.HotIn.Temperature < self.ColdOut.Temperature else False 


    def get_graph_data(self, divisions=50):
        H = []; T_hot = []; T_cold = []

        for i in range(divisions + 1):
            dh_hot = (self.HotIn.Enthalpy - self.HotOut.Enthalpy) * i / divisions
            dP_hot = (self.HotIn.Pressure - self.HotOut.Pressure) * i / divisions
            T_hot.append(REFPROP('T', self.HotOut.Fluid, 'PH', (self.HotOut.Pressure + dP_hot) / 1000, self.HotOut.Enthalpy + dh_hot))
                        
            dh_cold = (self.ColdOut.Enthalpy - self.ColdIn.Enthalpy) * i / divisions
            dP_cold = (self.ColdIn.Pressure - self.ColdOut.Pressure) * i / divisions
            T_cold.append(REFPROP('T', self.ColdIn.Fluid, 'PH', (self.ColdIn.Pressure - dP_cold) / 1000, self.ColdIn.Enthalpy + dh_cold))

            H.append(dh_hot * self.HotOut.Mass_Flow)

        return H, T_hot, T_cold
        