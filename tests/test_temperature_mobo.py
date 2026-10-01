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
                x=x or [0., 0., .5, .5, .5], success=success)


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
        encode([0., 0., .5, .5, .5], 102)
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
    records = [record(t, x=[0., 0., .5, .2+i*.02, .8-i*.02]) for i, t in enumerate(TEMPERATURES)]
    search = TemperatureParetoSearch(records, subset_size=8)
    indices = search.training_indices(encode([0., 0., .5, .44, .56], 160))
    assert len(indices) == 8
    assert len({records[i]['T_air_tar'] for i in indices}) > 1
    # Same equipment design at a different temperature is a distinct evaluation.
    search.tell(record(100, x=[0., 0., .5, .44, .56]))
    with pytest.raises(ValueError, match='Duplicate'):
        search.tell(record(100, x=[0., 0., .5, .44, .56]))


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
    assert search.diagnostic['method'] == 'temperature_qLogEHVI'
    assert search.diagnostic['training_points'] == 128
    assert len(search.diagnostic['training_temperature_levels']) > 1
    other_x, other_temperature = restored.ask()
    assert temperature == other_temperature
    np.testing.assert_allclose(x, other_x, rtol=0, atol=1e-10)


def test_schedule_visits_every_temperature_once_per_25_steps():
    from rhp_temperature_mobo import LEVELS, SCHEDULE_STRIDE
    assert sorted((SCHEDULE_STRIDE*k) % LEVELS for k in range(LEVELS)) == list(range(LEVELS))
    # Consecutive proposals in a batch of four land on well-separated temperatures.
    assert [(SCHEDULE_STRIDE*k) % LEVELS for k in range(4)] == [0, 7, 14, 21]


def test_fronts_and_hypervolume_are_per_temperature():
    from rhp_temperature_mobo import (mean_temperature_hypervolume, per_temperature_front_indices,
                                      per_temperature_hypervolumes)
    # A 220 C design no longer dominates an identical 100 C design, and within one temperature
    # dominance is 2D.
    records = [record(100, 3, 2), record(220, 3, 2), record(100, 2, 1, x=[1., 0., .5, .5, .5]),
               record(100, 2, 3, x=[2., 0., .5, .5, .5])]
    assert per_temperature_front_indices(records) == [0, 1, 3]
    hv = per_temperature_hypervolumes(records)
    assert hv[220.] == pytest.approx(.75*.2)
    assert hv[100.] == pytest.approx(.5*.3+.25*.2)
    assert mean_temperature_hypervolume(records) == pytest.approx((hv[100.]+hv[220.])/25)


def test_sparse_temperature_borrows_centres_from_neighbours():
    search = TemperatureParetoSearch()
    for i, t in enumerate((100, 105)):
        for j in range(3):
            search.tell(record(t, cop=3-.1*j, vhc=1+.5*j, x=[10.*j, 5., .45, .2+.1*j+.01*i, .8-.1*j-.01*i]))
    centres = search.level_centres(24)          # nothing observed at 220 C yet
    assert len(centres) == 4 and np.all(centres[:, 4] == 1.)
    assert admissible(centres).all()


def test_each_zero_lift_face_gets_a_region_and_face_steps_stay_on_the_face():
    from rhp_temperature_mobo import face_of
    search = TemperatureParetoSearch(seed=3)
    # Interior designs at 200 C, and dominated designs on each zero-lift face.
    for j in range(6):
        search.tell(record(200, cop=3-.1*j, vhc=1+.5*j, x=[10.+j, 20.+j, .5, .2+.1*j, .8-.1*j]))
    search.tell(record(200, cop=2.5, vhc=2.0, x=[70., 0., .5, .8, .2]))
    search.tell(record(200, cop=2.6, vhc=1.9, x=[0., 40., .5, .6, .4]))
    faces = {face for _, face in search.face_centres(20)}
    assert faces == {'rhx1_zero', 'rhx2_zero'}
    # A temperature with no face designs borrows them from nearby levels, scaled into its lift limit.
    assert {face for _, face in search.face_centres(17)} == {'rhx1_zero', 'rhx2_zero'}
    assert search.face_centres(0) == []        # 100 C is more than four levels away
    assert face_of([0., .3, .5, .5, 1.]) == 'rhx1_zero' and face_of([.3, 0., .5, .5, 1.]) == 'rhx2_zero'


def test_boundary_candidates_cover_every_bound_and_the_origin():
    from rhp_temperature_mobo import boundary_candidates
    pts = np.random.default_rng(1).uniform(.2, .8, (20, 5))
    b = boundary_candidates(pts)
    assert len(b) == 9*16
    for axis in range(4):
        assert (b[:, axis] == 0.).any() and (b[:, axis] == 1.).any()
    assert ((b[:, 0] == 0.) & (b[:, 1] == 0.)).any()


def test_split_bounds_restrict_proposals_and_observations():
    search = TemperatureParetoSearch(seed=4, split_bounds=(.4, .6))
    for _ in range(12):
        x, temperature = search.ask()
        assert .4-1e-12 <= x[2] <= .6+1e-12
        search.tell(record(temperature, x=x.tolist()))
    with pytest.raises(ValueError, match='outside'):
        search.tell(record(150, x=[0., 0., .3, .5, .5]))
    u = encode([10., 5., .55, .3, .7], 150, (.4, .6))
    assert u[2] == pytest.approx(.75) and decode(u, (.4, .6))[0][2] == pytest.approx(.55)
    with pytest.raises(ValueError, match='Split bounds'):
        TemperatureParetoSearch(split_bounds=(.05, .6))
