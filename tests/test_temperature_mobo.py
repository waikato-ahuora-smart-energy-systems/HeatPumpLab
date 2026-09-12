"""Joint discrete-temperature objective geometry, transfer and milestone tests."""
import copy

import numpy as np
import pytest
import torch
from botorch.utils.multi_objective.hypervolume import Hypervolume

from rhp_temperature_mobo import (
    TEMPERATURES, TemperatureGridPartition, TemperatureParetoSearch,
    add_known_temperature, admissible, decode, encode, front_indices,
    hypervolume_3d, objectives,
)
from main_RHPtemperature_mobo import milestone_snapshot, read_journal, append, write_milestones


def record(temperature, cop=3., vhc=2., x=None, success=True):
    return dict(T_air_tar=float(temperature), COP=cop, VHC_MJ_per_m3=vhc,
                x=x or [0., 0., .1, .5, .5], success=success)


def test_temperature_is_exact_in_each_posterior_draw():
    X = torch.tensor([[[0., 0., 0., .5, 0.]], [[0., 0., 0., .5, 1.]]], dtype=torch.double)
    samples = torch.randn(64, 2, 1, 2, dtype=torch.double)
    output = add_known_temperature(samples, X)
    assert torch.equal(output[..., :2], samples)
    assert torch.all(output[:, 0, :, 2] == .04)
    assert torch.all(output[:, 1, :, 2] == 1.)
    assert torch.all(output[..., 2].var(dim=0) == 0.)


def test_temperature_grid_and_temperature_dependent_lift():
    for temperature in TEMPERATURES:
        u = encode([20., 30., .4, .3, .7], temperature)
        assert decode(u)[1] == temperature
        assert admissible(u)[0]
    assert not admissible(encode([40., 40., .4, .3, .7], 100))[0]
    assert admissible(encode([40., 40., .4, .3, .7], 220))[0]
    with pytest.raises(ValueError):
        encode([0., 0., .1, .5, .5], 102)
    with pytest.raises(ValueError):
        decode([0., 0., 0., .5, .01])


def test_joint_dominance_compares_temperatures_and_keeps_tradeoffs():
    records = [record(100, 3, 2), record(105, 3, 2), record(220, 2, 1),
               record(220, 9, 9, success=False), record(105, 3, 2)]
    assert front_indices(records) == [1, 2, 4]


def test_grid_hypervolume_and_cells_match_independent_3d_geometry():
    rng = np.random.default_rng(124)
    points = np.column_stack((rng.uniform(.2, 2, (80, 2)), rng.integers(1, 26, 80)/25))
    expected = Hypervolume(torch.zeros(3, dtype=torch.double)).compute(torch.tensor(points))
    assert hypervolume_3d(points) == pytest.approx(expected, rel=1e-12)
    partition = TemperatureGridPartition(points)
    lo, hi = partition.bounds.numpy()
    for candidate in np.column_stack((rng.uniform(.2, 3, (30, 2)), rng.integers(1, 26, 30)/25)):
        cell_gain = np.maximum(np.minimum(candidate, hi)-lo, 0).prod(1).sum()
        expected_gain = hypervolume_3d(np.vstack((points, candidate)))-expected
        assert cell_gain == pytest.approx(expected_gain, abs=1e-12)


def test_shared_training_subset_uses_observations_across_temperatures():
    records = [record(t, x=[0., 0., .1, .2+i*.02, .8-i*.02]) for i, t in enumerate(TEMPERATURES)]
    search = TemperatureParetoSearch(records, subset_size=8)
    indices = search.training_indices(encode([0., 0., .1, .44, .56], 160))
    assert len(indices) == 8
    assert len({records[i]['T_air_tar'] for i in indices}) > 1
    # Same equipment design at a different temperature is a distinct evaluation.
    search.tell(record(100, x=[0., 0., .1, .44, .56]))
    with pytest.raises(ValueError, match='Duplicate'):
        search.tell(record(100, x=[0., 0., .1, .44, .56]))


def test_milestones_use_only_prefix_and_fixed_reference(tmp_path):
    records = [record(100), record(105), record(220, 4, 3)]
    write_milestones(tmp_path, records, [1, 2, 3])
    import json
    summaries = json.loads((tmp_path/'milestone_comparison.json').read_text())
    assert [r['evaluations'] for r in summaries] == [1, 2, 3]
    assert summaries[0]['hypervolume_normalized'] == pytest.approx(.75*.2*.04)
    assert summaries[1]['hypervolume_normalized'] == pytest.approx(2*summaries[0]['hypervolume_normalized'])
    assert summaries[2]['hypervolume_gain_percent_since_first_mark'] > 0
    assert len(milestone_snapshot(records)['per_temperature']) == 25


def test_torn_journal_recovers_only_incomplete_final_line(tmp_path):
    path = tmp_path/'records.jsonl'
    append(path, record(100))
    with path.open('ab') as stream:
        stream.write(b'{"T_air_tar":')
    assert read_journal(path) == [record(100)]
    assert path.read_bytes().endswith(b'\n')
    with path.open('ab') as stream:
        stream.write(b'broken\n')
    with pytest.raises(ValueError):
        read_journal(path)


def test_actual_joint_bo_proposal_and_resume():
    search = TemperatureParetoSearch(seed=812)
    for i, u in enumerate(search.initial):
        x, temperature = decode(u)
        search.tell(record(temperature, cop=4-.008*(temperature-100)+.3*x[3],
                           vhc=1+x[3]+.005*(temperature-100), x=x.tolist()))
    restored = TemperatureParetoSearch(search.records, seed=812, state=copy.deepcopy(search.state))
    x, temperature = search.ask()
    assert admissible(encode(x, temperature))[0]
    assert search.diagnostic['method'] == 'joint_qLogEHVI'
    assert search.diagnostic['training_points'] == 128
    assert len(search.diagnostic['training_temperature_levels']) > 1
    other_x, other_temperature = restored.ask()
    assert temperature == other_temperature
    np.testing.assert_allclose(x, other_x, rtol=0, atol=1e-10)
