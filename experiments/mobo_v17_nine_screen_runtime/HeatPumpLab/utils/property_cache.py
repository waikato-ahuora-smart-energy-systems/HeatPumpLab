"""Scoped cache for thermodynamic property lookups."""

from __future__ import annotations

from collections.abc import Hashable, Iterator
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, field

__all__ = [
    "PropertyCacheKey",
    "PropertyCallCache",
    "active_property_call_cache",
    "property_call_cache",
]

PropertyCacheKey = tuple[Hashable, ...]

_ACTIVE_PROPERTY_CACHE: ContextVar["PropertyCallCache | None"] = ContextVar(
    "active_property_call_cache", default=None
)


@dataclass
class PropertyCallCache:
    """Store property results for a single solve scope."""

    values: dict[PropertyCacheKey, object] = field(default_factory=dict)
    hits: int = 0
    misses: int = 0

    def lookup(self, key: PropertyCacheKey) -> tuple[bool, object | None]:
        """Return a cached result and whether the lookup hit."""
        if key in self.values:
            self.hits += 1
            return True, self.values[key]

        self.misses += 1
        return False, None

    def store(self, key: PropertyCacheKey, value: object) -> None:
        """Store a result for later lookups in the same solve."""
        self.values[key] = value


@contextmanager
def property_call_cache() -> Iterator[PropertyCallCache]:
    """Cache property lookups inside the current context only."""
    cache = PropertyCallCache()
    token = _ACTIVE_PROPERTY_CACHE.set(cache)
    try:
        yield cache
    finally:
        _ACTIVE_PROPERTY_CACHE.reset(token)


def active_property_call_cache() -> PropertyCallCache | None:
    """Return the active property cache, if one has been enabled."""
    return _ACTIVE_PROPERTY_CACHE.get()
