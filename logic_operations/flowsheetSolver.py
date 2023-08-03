__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2023'

from scipy.optimize import brentq

class FlowSheetSolver:
    def __init__(self) -> None:
        pass
   
    def get_state_solution(self):
        self.converged = False
        i = 0; z = 1
        while z > 0.0 and i < 20: #Check for model convergance            
            z0 = [[s.mass_flow, s.pressure, s.enthalpy] for s in self.streams]

            for op in self.operations:
                op.solve()

            z1 = [[s.mass_flow, s.pressure, s.enthalpy] for s in self.streams]
            z = min(abs(z0 - z1))
            i += 1

        self.converged = True if i < 20 else False
        return self.converged
