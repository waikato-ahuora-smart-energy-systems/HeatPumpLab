"""Observation validation must tolerate unit-conversion round-off only."""

import numpy as np
import pytest

from rhp_mobo_search import ParetoBOSearch


def test_observation_accepts_ulp_overshoot_without_changing_coordinates():
    search = ParetoBOSearch(123.33333333333333)
    x = np.asarray([76.66171575, 123.33333333333333-76.66171575, .39451915, .82349311, .17650689])
    x[1] = np.nextafter(x[1], np.inf)
    encoded = search.encode(x)
    search.tell(x, 2.4, 3., True)
    assert np.array_equal(search.x[0], encoded)
    assert len(search._admissible(np.asarray([[.8, .43333333334, .5, .5]]))) == 0


def test_material_domain_violation_is_still_rejected():
    search = ParetoBOSearch(123.33333333333333)
    with pytest.raises(ValueError, match="outside"):
        search.tell([80., 43.333334, .5, .5, .5], 2.4, 3., True)
