from scipy.optimize import shgo
import numpy as np

def f(x, p):  # (cattle-feed)
    return (24.55*x[0] + 26.75*x[1] + 39*x[2] + 40.50*x[3]) * p[0]

def g1(x, p):
    return (2.3*x[0] + 5.6*x[1] + 11.1*x[2] + 1.3*x[3] - 5) * p[0]  # >=0

def g2(x, p):
    return (12*x[0] + 11.9*x[1] +41.8*x[2] + 52.1*x[3] - 21
            - 1.645 * np.sqrt(0.28*x[0]**2 + 0.19*x[1]**2
                            + 20.5*x[2]**2 + 0.62*x[3]**2)
            ) * p[0] # >=0

def h1(x, p):
    return (x[0] + x[1] + x[2] + x[3] - 1) * p[0]  # == 0

p = [(2, 5)]

cons = ({'type': 'ineq', 'fun': g1, 'args': p},
        {'type': 'ineq', 'fun': g2, 'args': p},
        {'type': 'eq', 'fun': h1, 'args': p})
bounds = [(0, 1.0),]*4
res = shgo(f, args=p, bounds=bounds, n=150, constraints=cons)
print(res)