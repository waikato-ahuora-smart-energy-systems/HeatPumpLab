"""Utility decorators used across the codebase."""

from time import perf_counter


def timer_func(func) -> any:
    """Measure execution time of the decorated function."""
    def wrapper(*args, **kwargs) -> any:
        """Wrap the callable and return both result and elapsed time."""
        tic: float = perf_counter()
        r = func(*args, **kwargs)
        toc: float = perf_counter()
        return r, toc - tic
    return wrapper
