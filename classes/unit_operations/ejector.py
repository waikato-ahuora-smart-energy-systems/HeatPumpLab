__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from classes.streams.properties import REFPROP
from classes.streams.massStream import MassStream

class Ejector:
    def __init__(self, Name='Ejector', Eff_is_ejc=1.0, tol=1e-8):
        self.Name = Name
        self.tol = tol
        self.PriFlowIn = MassStream()
        self.SecFlowIn = MassStream()
        self.FlowOut = MassStream()
        self.Eff_is_ejc = Eff_is_ejc
        self.dH_comp = 0
        self.Power = 0
        self.Pressure_Ratio = 1
        self.m_ratio = 1
        self.converged = False


    def solve(self):
        if self.Eff_is_ejc > 0 and self.Eff_is_ejc <= 1:
            self.get_work()
        else:
            self.FlowOut.Combine(self.PriFlowIn, self.SecFlowIn)   
        
        self.ejector_mass_balance() 


    def ejector_mass_balance(self, i = 0):
        while self.converged == False and i < 100:
            #Solves the flow rates around the ejector and cycle
            m_ratio_1 = self.m_ratio
            self.FlowOut.Quality = REFPROP('QMASS', self.FlowOut.Fluid, 'PH', self.FlowOut.Pressure / 1000, self.FlowOut.Enthalpy, 0)
            self.m_ratio = ( 1 / self.FlowOut.Quality )  - 1
            h_out_1 = (self.SecFlowIn.Enthalpy * self.m_ratio + self.PriFlowIn.Enthalpy )  /  ( self.m_ratio + 1 )
            self.FlowOut.Enthalpy = ( h_out_1 + self.FlowOut.Enthalpy )  / 2

            if ( abs(self.FlowOut.Enthalpy / h_out_1 - 1) > self.tol or abs(self.m_ratio / m_ratio_1 - 1) > self.tol ):
                i += 1
            else:
                self.converged = True

        if i == 100:
            self.converged = False

        self.PriFlowIn.Mass_Flow = 1 / self.m_ratio
        self.FlowOut.Mass_Flow = 1 / self.m_ratio + 1
        self.FlowOut.PropsIn_mPH()


    def get_work(self):
        if abs(self.FlowOut.Pressure - self.SecFlowIn.Pressure) / self.SecFlowIn.Pressure < 0.0001:
            #Good estimate of the ejector outlet pressure
            #Use a large step for Newton's method
            y_min = ( REFPROP('H', self.SecFlowIn.Fluid, 'PS', self.SecFlowIn.Pressure / 1000, self.PriFlowIn.Entropy, 0) - self.PriFlowIn.Enthalpy )  * self.Eff_is_ejc
            y_max = ( REFPROP('H', self.SecFlowIn.Fluid, 'PS', self.PriFlowIn.Pressure / 1000, self.SecFlowIn.Entropy, 0) - self.SecFlowIn.Enthalpy )  * self.m_ratio / self.Eff_is_ejc
            self.FlowOut.Pressure = self.SecFlowIn.Pressure + abs(y_min /  ( y_max - y_min )) *  ( self.PriFlowIn.Pressure - self.SecFlowIn.Pressure )
            #Set the step for use in Newton's method below
            dx = self.FlowOut.Pressure - self.SecFlowIn.Pressure
        else:
            #Set the step for use in Newton's method below
            dx = 0.01
        
        i = 0
        while i < 100:
            #For a given self.PriFlowIn inlet pressure and fixed self.SecFlowIn pressure,
            #converge on point where the Power of expansion equals the Power of compression
            #within the ejector unit using Newton's method
            y_min = ( REFPROP('H', self.FlowOut.Fluid, 'PS', ( self.FlowOut.Pressure - dx )  / 1000, self.SecFlowIn.Entropy, 0) - self.SecFlowIn.Enthalpy )  * self.m_ratio / self.Eff_is_ejc \
                    -  ( self.PriFlowIn.Enthalpy - REFPROP('H', self.FlowOut.Fluid, 'PS', ( self.FlowOut.Pressure - dx )  / 1000, self.PriFlowIn.Entropy, 0) )  * self.Eff_is_ejc
            y_max = ( REFPROP('H', self.FlowOut.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.SecFlowIn.Entropy, 0) - self.SecFlowIn.Enthalpy )  * self.m_ratio / self.Eff_is_ejc \
                    -  ( self.PriFlowIn.Enthalpy - REFPROP('H', self.FlowOut.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.PriFlowIn.Entropy, 0) )  * self.Eff_is_ejc
            self.FlowOut.Pressure = self.FlowOut.Pressure -  ( y_max /  ( y_max - y_min ) )  * dx
            dh = abs(( REFPROP('H', self.FlowOut.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.SecFlowIn.Entropy, 0) - self.SecFlowIn.Enthalpy )  * self.m_ratio / self.Eff_is_ejc \
                    -  ( self.PriFlowIn.Enthalpy - REFPROP('H', self.FlowOut.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.PriFlowIn.Entropy, 0) )  * self.Eff_is_ejc) \
                    /  ( ( self.PriFlowIn.Enthalpy - REFPROP('H', self.FlowOut.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.PriFlowIn.Entropy, 0) )  * self.Eff_is_ejc )
            
            i += 1
            if dh < self.tol:
                break
        
        self.Pressure_Ratio = self.FlowOut.Pressure / self.SecFlowIn.Pressure
        self.dH_comp = REFPROP('H', self.FlowOut.Fluid, 'PS', self.FlowOut.Pressure / 1000, self.SecFlowIn.Entropy, 0) - self.SecFlowIn.Enthalpy
        self.Power = self.dH_comp * self.SecFlowIn.Mass_Flow / self.Eff_is_ejc