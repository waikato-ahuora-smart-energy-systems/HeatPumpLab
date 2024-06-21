__author__ = 'Tim Walmsley'
__credits__ = ['']

from streams import get_property_value, MassStream

class BaseHeatExchanger:
    def __init__(self) -> None:
        pass

    def init_stream_connections(self):
        #Set the stream connection points
        self.source_in = MassStream()
        self.source_out = self.source_in
        self.sink_in = MassStream()
        self.sink_out = self.sink_in    
    
    def init_performance_metrics(self):
        self.NTU = 0
        self.c = 0
        self.U = 0
        self.Area = 0
        self.duty = 0

    def solve(self):
        if self.bypassed:
            self.source_out = self.source_in
            self.sink_out = self.sink_in
        else:
            if self.source_in.mass_flow <= self.tol or self.sink_in.mass_flow <= self.tol:
                self.source_out.set_equal(fld=self.source_in)
                self.sink_out.set_equal(fld=self.sink_in)          

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

        elif self.model == 'saturated gas outlet evaporator':
            self.get_saturated_gas_sink_out_model()

        elif self.model == 'specialised':
            self.get_detailed_performance()

        else:
            raise ValueError('No heat exchanger model selected.')

    def get_dTmin_terminal_model(self):
        try:
            self.source_out.pressure = self.source_in.pressure - self.dP_h
            self.sink_out.pressure = self.sink_in.pressure - self.dP_c
            self.dH_h_max = self.source_in.mass_flow * (self.source_in.enthalpy - get_property_value('H', self.source_out.fluid, 'PT', self.source_out.pressure / 1000, self.sink_in.temperature + self.dTmin, self.RP))
            self.dH_c_max = self.sink_in.mass_flow * (get_property_value('H', self.sink_out.fluid, 'PT', self.sink_out.pressure / 1000, self.source_in.temperature - self.dTmin, self.RP) - self.sink_in.enthalpy )
            self.duty = min(self.dH_h_max, self.dH_c_max)

            self.source_out.enthalpy = self.source_in.enthalpy - self.duty / self.source_in.mass_flow
            self.source_out.props_in_mPH()
            self.sink_out.enthalpy = self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
            self.sink_out.props_in_mPH()        
        except:
            raise ValueError('Warning: Phase change occuring in heat exchanger ' + self.name)   

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False    

    def get_effectiveness_model(self):
        try:
            self.source_out.pressure = self.source_in.pressure - self.dP_h
            self.sink_out.pressure = self.sink_in.pressure - self.dP_c
            self.dT_max = abs( self.source_in.temperature - self.sink_in.temperature )
            self.dH_h_max = self.source_in.mass_flow * (self.source_in.enthalpy - get_property_value('H', self.source_out.fluid, 'PT', self.source_out.pressure / 1000, self.source_in.temperature - self.dT_max * self.HX_Eff, self.RP))
            self.dH_c_max = self.sink_in.mass_flow * (get_property_value('H', self.sink_out.fluid, 'PT', self.sink_out.pressure / 1000, self.sink_in.temperature + self.dT_max * self.HX_Eff, self.RP) - self.sink_in.enthalpy )
            self.duty = min(self.dH_h_max, self.dH_c_max)

            self.source_out.enthalpy = self.source_in.enthalpy - self.duty / self.source_in.mass_flow
            self.source_out.props_in_mPH()
            self.sink_out.enthalpy = self.sink_in.enthalpy - self.duty / self.sink_in.mass_flow
            self.sink_out.props_in_mPH()
        
        except:
            raise ValueError('Warning: Phase change occuring in heat exchanger ' + self.name)        

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False    

    def get_fixed_cold_duty_model(self):
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.props_in_mPH()

        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        self.source_out.props_in_mPH()

        # if self.source_out.temperature < self.sink_in.temperature or self.source_in.temperature < self.sink_out.temperature:
        #      print('temperature cross occuring in the heat exchanger: ' + self.name)

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False

    def get_fixed_hot_duty_model(self):
        self.duty = self.source_in.heat_flow - self.source_out.heat_flow
        self.source_out.pressure = self.source_out.pressure - self.dP_h
        self.source_in.props_in_mPH()

        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.enthalpy = self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
        self.sink_out.props_in_mPH()

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False

    def get_fixed_cold_T_outlet_model(self):
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.props_in_mPT()
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow
        
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        self.source_out.props_in_mPH()

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False

    def get_fixed_hot_T_outlet_model(self):
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        # Tsat = get_property_value('T', self.source_out.fluid, 'PQ', self.source_out.pressure / 1000, 0, RP=self.source_out.RP)
        # if abs(Tsat - self.source_out.temperature) < 1 or abs(Tsat - self.source_in.temperature) < 1:
        #     None
        self.source_out.props_in_mPT()
        self.duty = self.source_in.heat_flow - self.source_out.heat_flow
        
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.enthalpy = self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
        self.sink_out.props_in_mPH()

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False

    def get_fixed_cold_delta_T_model(self):
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.temperature = self.sink_in.temperature + self.dT_c
        self.sink_out.props_in_mPT()
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow
        
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        self.source_out.props_in_mPH()

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False 

    def get_fixed_hot_delta_T_model(self):
        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.temperature = self.source_in.temperature - self.dT_h
        self.source_out.props_in_mPT()
        self.duty = self.source_in.heat_flow - self.source_out.heat_flow
        
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.enthalpy = self.sink_in.enthalpy + self.duty / self.sink_in.mass_flow
        self.sink_out.props_in_mPH()

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False 

    def get_saturated_gas_sink_out_model(self):
        self.sink_in.props_in_mPH()
        self.sink_out.pressure = self.sink_in.pressure - self.dP_c
        self.sink_out.props_in_mPX()
        self.duty = self.sink_out.heat_flow - self.sink_in.heat_flow

        self.source_out.pressure = self.source_in.pressure - self.dP_h
        self.source_out.enthalpy = self.source_in.enthalpy - self.duty / self.source_out.mass_flow
        self.source_out.props_in_mPH()

        # if self.source_out.temperature < self.sink_in.temperature or self.source_in.temperature < self.sink_out.temperature:
        #      print('temperature cross occuring in the heat exchanger: ' + self.name)

        self.converged = True if self.source_out.temperature > self.sink_in.temperature and self.source_in.temperature > self.sink_out.temperature else False

    def get_graph_data(self, divisions=50):
        data = {'source': 
                    {'P': [], 'h': [], 'T': [], 's': [], 'v': [], 'H': []},
                'sink': 
                    {'P': [], 'h': [], 'T': [], 's': [], 'v': [], 'H': []}
                }

        for i in range(divisions + 1):
            dP_hot = (self.source_in.pressure - self.source_out.pressure) * i / divisions
            dh_hot = (self.source_in.enthalpy - self.source_out.enthalpy) * i / divisions
            data['source']['P'].append(self.source_in.pressure - dP_hot)
            data['source']['h'].append(self.source_in.enthalpy - dh_hot)            
            data['source']['T'].append(    get_property_value('T', self.source_out.fluid, 'PH', data['source']['P'][-1] / 1000, data['source']['h'][-1], self.RP))
            data['source']['s'].append(    get_property_value('S', self.source_out.fluid, 'PH', data['source']['P'][-1] / 1000, data['source']['h'][-1], self.RP))
            data['source']['v'].append(1 / get_property_value('D', self.source_out.fluid, 'PH', data['source']['P'][-1] / 1000, data['source']['h'][-1], self.RP))
            data['source']['H'].append(dh_hot * self.source_out.mass_flow)

            dP_cold = (self.sink_in.pressure - self.sink_out.pressure) * i / divisions
            dh_cold = (self.sink_in.enthalpy - self.sink_out.enthalpy) * i / divisions
            data['sink']['P'].append(self.sink_in.pressure - dP_cold)
            data['sink']['h'].append(self.sink_in.enthalpy - dh_cold)
            data['sink']['T'].append(    get_property_value('T', self.sink_in.fluid, 'PH', data['sink']['P'][-1] / 1000, data['sink']['h'][-1], self.RP))
            data['sink']['s'].append(    get_property_value('S', self.sink_in.fluid, 'PH', data['sink']['P'][-1] / 1000, data['sink']['h'][-1], self.RP))
            data['sink']['v'].append(1 / get_property_value('D', self.sink_in.fluid, 'PH', data['sink']['P'][-1] / 1000, data['sink']['h'][-1], self.RP))
            data['sink']['H'].append(dh_cold * self.sink_out.mass_flow)
        
        return data
        