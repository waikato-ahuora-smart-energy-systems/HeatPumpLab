"""Local anchors survive a front that pure-fluid designs almost entirely dominate."""
import json
from pathlib import Path

import pytest

from rhp_local_pareto_search import choose_anchors, front_indices


def record(cop, vhc, comp):
    return dict(x=[1., 50., .5, *comp], COP=cop, VHC_MJ_per_m3=vhc, success=True)


def test_mixed_layers_fill_anchors_when_pure_designs_dominate():
    pure = [record(2.6-.05*k, 1.2+.1*k, [1., 0., 0.]) for k in range(8)]
    mixed_front = [record(2.62-.04*k, 1.0+.05*k, [.5, .5, 0.]) for k in range(3)]
    dominated = [record(2.3-.03*k, 1.0+.08*k, [.4, .3, .3]) for k in range(6)]
    records = pure+mixed_front+dominated
    anchors = choose_anchors(records)
    assert len(anchors) == 6 and len(set(anchors)) == 6
    assert all(max(records[i]['x'][3:]) <= .99 for i in anchors)
    # Only the first mixed design is on the overall front; the rest come from mixed layers.
    on_front = set(front_indices([[r['COP'], r['VHC_MJ_per_m3']] for r in records]))
    assert on_front & set(range(8, len(records))) == {8} and 8 in anchors


def test_too_few_mixed_designs_still_raises():
    records = [record(2.6-.05*k, 1.2+.1*k, [1., 0., 0.]) for k in range(8)]
    records += [record(2.3, 1.0+.1*k, [.5, .5, 0.]) for k in range(5)]
    with pytest.raises(ValueError, match='Too few mixed-fluid'):
        choose_anchors(records)


RUN = Path('results/RHPfixed200_v24/Butane_Isopentane_Neopentane/result.json')


@pytest.mark.skipif(not RUN.exists(), reason='sweep result not present')
def test_failed_neopentane_refresh_now_has_six_anchors():
    records = json.loads(RUN.read_text())['records'][:406]
    assert len(choose_anchors(records)) == 6
