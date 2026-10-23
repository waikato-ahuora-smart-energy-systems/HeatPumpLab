import numpy as np
import pytest


@pytest.mark.parametrize("components", [2, 3])
def test_min_fraction_maps_whole_search_domain(components):
    search_mod = pytest.importorskip("rhp_mobo_search")
    s = search_mod.ParetoBOSearch(100.0, 1, components=components, min_fraction=0.01)
    # _candidates() measures distance to observed designs, so tell one first (as in a real run)
    s.tell(s.decode(s.initial[0]), 2.5, 2.0, True)
    pool, _ = s._candidates()
    for u in np.vstack([s.initial, pool]):
        x = s.decode(u)
        w = x[3:]
        assert len(w) == components
        assert w.min() >= 0.01 - 1e-12
        assert w.sum() == pytest.approx(1.0, abs=1e-12)
        assert np.allclose(s.encode(x), u, atol=1e-12)
    # pure-fluid anchor (encoded vertex) becomes 99/1 or 98/1/1
    vertex = s.decode(s.initial[0])[3:]
    assert vertex.max() == pytest.approx(1 - 0.01 * (components - 1))


def test_zero_min_fraction_unchanged():
    search_mod = pytest.importorskip("rhp_mobo_search")
    a = search_mod.ParetoBOSearch(100.0, 1, components=3)
    for u in a.initial:
        x = a.decode(u)
        assert np.allclose(x[3:5], u[-2:]) and x[5] == pytest.approx(max(0.0, 1 - u[-2:].sum()))


def test_min_fraction_too_large():
    search_mod = pytest.importorskip("rhp_mobo_search")
    with pytest.raises(ValueError):
        search_mod.ParetoBOSearch(100.0, 1, components=3, min_fraction=0.34)
