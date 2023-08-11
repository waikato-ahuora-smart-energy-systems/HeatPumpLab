__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2023'

from scipy.optimize import brentq

class Controller():
    def __init__(self) -> None:
        pass

    def get_target_values_difference(self, x, x_obj, y_obj, y_target):
        x_obj[0] = x
        self.converged = self.get_state_solution()
        return y_obj[0] - y_target[0]

    def get_x_bounds(self, x_obj, y_obj, y_target):
        x_lb = x_obj[0]
        x_ub = x_obj[0]

        self.get_state_solution()
        y_neg = y_obj[0]
        i = 0
        if y_neg > 0.0:
            while y_neg > 0.0 and i < 20:
                x_ub = x_lb
                x_lb *= 0.9
                x_obj[0] = x_lb
                y_neg = self.get_state_solution()
                i += 1
            
        else:
            while y_plus < 0.0 and i < 20:
                x_lb = x_ub
                x_ub *= 1.1
                x_obj[0] = x_ub
                y_plus = self.get_state_solution()
                i += 1
        
        feasible = int((-y_plus * y_neg) / abs(y_plus * y_neg))
        return x_lb, x_ub, feasible

    def get_target_value(self, x_init, x_obj, y_obj, y_target):
        args = {'x_obj'    : x_obj, 
                'y_obj'    : y_obj,
                'y_target' : y_target }
        
        r = [0]
        try:
            x_lb, x_ub, feasible = self.get_x_bounds(**args)                     
            if feasible == 1:
                r = brentq(f=self.get_target_values_difference, 
                           args=args, 
                           a=x_lb, 
                           b=x_ub, 
                           full_output=True, 
                           rtol=self.tol, 
                           xtol=self.tol)

        except:
            print('Controller block failed to converge.') 
            self.converged = False
                
        return r
