from scipy.optimize import shgo
res = shgo(
    lambda x, k: x[0]**k,
    bounds=[(0, 10)],
    args=(1,),
    constraints=({'type': 'ineq', 'fun': lambda x: x[0] - 1}),
    sampling_method='sobol',
)
print(res['x'])