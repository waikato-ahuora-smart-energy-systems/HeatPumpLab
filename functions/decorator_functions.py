from time import perf_counter

def time_func(func) -> any:
    def wrapper(*args, **kwargs) -> any:
        tic: float = perf_counter()
        r = func(*args, **kwargs)
        toc: float = perf_counter()
        return r, toc - tic
    return wrapper
