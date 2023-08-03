__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2023'

from scipy.optimize import brentq

class Solver:
    def __init__(self) -> None:
        pass

    def get_x_bounds(self, x):
        x_lb = x
        x_ub = x

        y_plus = self.get_state(x)
        y_neg = y_plus
        i = 0
        if y_neg > 0.0:
            while y_neg > 0.0 and i < 20:
                x_ub = x_lb
                x_lb *= 0.9
                y_neg = self.get_state(x_lb)
                i += 1
            
        else:
            while y_plus < 0.0 and i < 20:
                x_lb = x_ub
                x_ub *= 1.1
                y_plus = self.get_state(x_ub)
                i += 1
        
        feasible = int((-y_plus * y_neg) / abs(y_plus * y_neg))
        return x_lb, x_ub, feasible
   
    def get_state(self, x, x_obj, y_obj, y_target):
        self.converged = False
        x_obj = x
        i = 0; z = 1
        
        while z > 0.0 and i < 20: #Check for model convergance            
            z0 = [[s.mass_flow, s.pressure, s.enthalpy] for s in self.streams]

            for op in self.operations:
                op.solve()

            z1 = [[s.mass_flow, s.pressure, s.enthalpy] for s in self.streams]
            z = min(abs(z0 - z1))
            i += 1

        self.converged = True if i < 20 else False
        return y_obj - y_target

    def get_target(self, x, args):
        #Solve model
        r = [0]
        try:
            x_lb, x_ub, feasible = self.get_x_bounds(x=x)                     
            if feasible == 1:
                r = brentq(f=self.get_state, args=args, a=x_lb, b=x_ub, full_output=True, rtol=self.tol, xtol=self.tol)

        except:
            print('Model failed to solve.') 
            self.converged = False
                
        return r
