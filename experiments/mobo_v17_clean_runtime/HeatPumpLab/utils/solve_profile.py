"""Opt-in profiling helpers for heat-pump solve paths."""

from __future__ import annotations

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field
from functools import wraps
from time import perf_counter
from typing import Callable, Iterator, TypeVar

__all__ = [
    "SolveProfile",
    "active_solve_profile",
    "profile_section",
    "profile_solve",
    "profiled_section",
]

F = TypeVar("F", bound=Callable)

_ACTIVE_PROFILE: ContextVar["SolveProfile | None"] = ContextVar(
    "active_solve_profile", default=None
)


@dataclass
class SolveProfile:
    """Collect named timing and call-count data for one profiling scope."""

    timings: dict[str, float] = field(default_factory=dict)
    counts: dict[str, int] = field(default_factory=dict)

    @contextmanager
    def section(self, name: str) -> Iterator[None]:
        """Record elapsed time and call count for a named code section."""
        start = perf_counter()
        try:
            yield
        finally:
            self.record(name, perf_counter() - start)

    def record(self, name: str, elapsed: float) -> None:
        """Add one observation for a named section."""
        self.timings[name] = self.timings.get(name, 0.0) + elapsed
        self.counts[name] = self.counts.get(name, 0) + 1

    def as_dict(self) -> dict[str, dict[str, float] | dict[str, int]]:
        """Return a serialisable snapshot of recorded profile data."""
        return {"timings": dict(self.timings), "counts": dict(self.counts)}


@contextmanager
def profile_solve() -> Iterator[SolveProfile]:
    """Activate solve profiling inside the context."""
    profile = SolveProfile()
    token = _ACTIVE_PROFILE.set(profile)
    try:
        yield profile
    finally:
        _ACTIVE_PROFILE.reset(token)


def active_solve_profile() -> SolveProfile | None:
    """Return the active solve profile, if any."""
    return _ACTIVE_PROFILE.get()


@contextmanager
def profile_section(name: str) -> Iterator[None]:
    """Record a section only when solve profiling is active."""
    profile = active_solve_profile()
    if profile is None:
        yield
        return

    with profile.section(name):
        yield


def profiled_section(name: str) -> Callable[[F], F]:
    """Return a decorator that records function execution when profiling is active."""

    def decorator(func: F) -> F:
        @wraps(func)
        def wrapper(*args, **kwargs):
            with profile_section(name):
                return func(*args, **kwargs)

        return wrapper

    return decorator
