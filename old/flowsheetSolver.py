__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2023'

import numpy as np

class FlowSheetSolver:
    def __init__(self) -> None:
        pass
   
    def get_state_solution(self, tol=1e-6):
        self.converged = False
        i = 0; z = 1
        while self.converged == False and i < 20: #Check for model convergance            
            z0 = np.array([[s.mass_flow, s.pressure, s.enthalpy] for s in self.streams])

            for op in self.operations:
                op.solve()

            z1 = np.array([[s.mass_flow, s.pressure, s.enthalpy] for s in self.streams])

            dz = abs(z0 - z1) - tol * z0
            if np.max(dz) < 0.0:
                self.converged = True

            i += 1

        self.converged = True if i < 20 else False
        return self.converged
