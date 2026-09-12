"""Observed hypervolume and replay-rejection bookkeeping."""

from main_RHPmobo import frontier_indices
from main_RHPmobo_report import hypervolume
import numpy as np
import pytest


def record(cop, vhc, **extra):
    return dict(COP=cop, VHC_MJ_per_m3=vhc, success=True, **extra)


def test_hypervolume_uses_both_maximisation_objectives_and_ignores_invalid_points():
    records = [record(2., 4.), record(3., 3.), record(1., 1.),
               dict(COP=None, VHC_MJ_per_m3=None, success=False)]
    assert hypervolume(records) == 11.
    assert hypervolume(records+[record(2., 4.)]) == 11.
    assert hypervolume(records+[record(4., 4.)]) == 16.


def test_rejected_replay_exposes_previously_dominated_candidate():
    records = [record(3., 4.), record(2., 3.), record(4., 2.)]
    assert frontier_indices(records) == [0, 2]
    records[0]["validation_rejected"] = True
    assert frontier_indices(records) == [1, 2]


def test_fast_hypervolume_matches_frontier_integration_with_ties_and_rejections():
    rng = np.random.default_rng(14)
    records = [record(*values) for values in rng.uniform(.1, 8., (100, 2))]
    records += [records[0].copy(), record(20., 20., validation_rejected=True),
                record(float('nan'), 1.), record(-1., 5.)]
    front = sorted((records[i] for i in frontier_indices(records)), key=lambda r: r['VHC_MJ_per_m3'])
    expected, previous = 0., 0.
    for r in front:
        expected += (r['VHC_MJ_per_m3']-previous)*r['COP']
        previous = r['VHC_MJ_per_m3']
    assert hypervolume(records) == pytest.approx(expected)
    assert hypervolume(list(reversed(records))) == pytest.approx(expected)
