__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams import get_property_value, MassStream

class Ejector:
    def __init__(self, name='Ejector', eff_is_ejc=1.0, tol=1e-6, RP=None):
        self.name = name
        self.tol = tol
        self.primary_flow_in = MassStream()
        self.sec_flow_in = MassStream()
        self.flow_out = MassStream()
        self.eff_is_ejc = eff_is_ejc
        self.dH_comp = 0
        self.power = 0
        self.pressure_ratio = 1
        self.m_ratio = 1
        self.converged = False
        self.RP = RP

    def solve(self):
        if self.eff_is_ejc > 0 and self.eff_is_ejc <= 1:
            self.get_work()
        else:
            self.flow_out.get_combined_stream(self.primary_flow_in, self.sec_flow_in)   
        self.ejector_mass_balance() 

    def ejector_mass_balance(self, i = 0):
        while self.converged == False and i < 100:
            #Solves the flow rates around the ejector and cycle
            m_ratio_1 = self.m_ratio
            self.flow_out.quality = get_property_value('QMASS', self.flow_out.fluid, 'PH', self.flow_out.pressure / 1000, self.flow_out.enthalpy, self.RP)
            self.m_ratio = ( 1 / self.flow_out.quality )  - 1
            h_out_1 = (self.sec_flow_in.enthalpy * self.m_ratio + self.primary_flow_in.enthalpy )  /  ( self.m_ratio + 1 )
            self.flow_out.enthalpy = ( h_out_1 + self.flow_out.enthalpy )  / 2

            if ( abs(self.flow_out.enthalpy / h_out_1 - 1) > self.tol or abs(self.m_ratio / m_ratio_1 - 1) > self.tol ):
                i += 1
            else:
                self.converged = True

        if i == 100:
            self.converged = False

        self.primary_flow_in.mass_flow = 1 / self.m_ratio
        self.flow_out.mass_flow = 1 / self.m_ratio + 1
        self.flow_out.props_in_mPH()

    def get_work(self):
        if abs(self.flow_out.pressure - self.sec_flow_in.pressure) / self.sec_flow_in.pressure < 0.0001:
            #Good estimate of the ejector outlet pressure
            #Use a large step for Newton's method
            y_min = ( get_property_value('H', self.sec_flow_in.fluid, 'PS', self.sec_flow_in.pressure / 1000, self.primary_flow_in.entropy, self.RP) - self.primary_flow_in.enthalpy )  * self.eff_is_ejc
            y_max = ( get_property_value('H', self.sec_flow_in.fluid, 'PS', self.primary_flow_in.pressure / 1000, self.sec_flow_in.entropy, self.RP) - self.sec_flow_in.enthalpy )  * self.m_ratio / self.eff_is_ejc
            self.flow_out.pressure = self.sec_flow_in.pressure + abs(y_min /  ( y_max - y_min )) *  ( self.primary_flow_in.pressure - self.sec_flow_in.pressure )
            #Set the step for use in Newton's method below
            dx = self.flow_out.pressure - self.sec_flow_in.pressure
        else:
            #Set the step for use in Newton's method below
            dx = 0.01
        
        i = 0
        while i < 100:
            #For a given self.primary_flow_in inlet pressure and fixed self.sec_flow_in pressure,
            #converge on point where the power of expansion equals the power of compression
            #within the ejector unit using Newton's method
            y_min = ( get_property_value('H', self.flow_out.fluid, 'PS', ( self.flow_out.pressure - dx )  / 1000, self.sec_flow_in.entropy, self.RP) - self.sec_flow_in.enthalpy )  * self.m_ratio / self.eff_is_ejc \
                    -  ( self.primary_flow_in.enthalpy - get_property_value('H', self.flow_out.fluid, 'PS', ( self.flow_out.pressure - dx )  / 1000, self.primary_flow_in.entropy, self.RP) )  * self.eff_is_ejc
            y_max = ( get_property_value('H', self.flow_out.fluid, 'PS', self.flow_out.pressure / 1000, self.sec_flow_in.entropy, self.RP) - self.sec_flow_in.enthalpy )  * self.m_ratio / self.eff_is_ejc \
                    -  ( self.primary_flow_in.enthalpy - get_property_value('H', self.flow_out.fluid, 'PS', self.flow_out.pressure / 1000, self.primary_flow_in.entropy, self.RP) )  * self.eff_is_ejc
            self.flow_out.pressure = self.flow_out.pressure -  ( y_max /  ( y_max - y_min ) )  * dx
            dh = abs(( get_property_value('H', self.flow_out.fluid, 'PS', self.flow_out.pressure / 1000, self.sec_flow_in.entropy, self.RP) - self.sec_flow_in.enthalpy )  * self.m_ratio / self.eff_is_ejc \
                    -  ( self.primary_flow_in.enthalpy - get_property_value('H', self.flow_out.fluid, 'PS', self.flow_out.pressure / 1000, self.primary_flow_in.entropy, self.RP) )  * self.eff_is_ejc) \
                    /  ( ( self.primary_flow_in.enthalpy - get_property_value('H', self.flow_out.fluid, 'PS', self.flow_out.pressure / 1000, self.primary_flow_in.entropy, self.RP) )  * self.eff_is_ejc )
            
            i += 1
            if dh < self.tol:
                break
        
        self.pressure_ratio = self.flow_out.pressure / self.sec_flow_in.pressure
        self.dH_comp = get_property_value('H', self.flow_out.fluid, 'PS', self.flow_out.pressure / 1000, self.sec_flow_in.entropy, self.RP) - self.sec_flow_in.enthalpy
        self.power = self.dH_comp * self.sec_flow_in.mass_flow / self.eff_is_ejc