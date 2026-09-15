"""A converged cycle is not an admissible design if a compressor PR exceeds its limit."""
import numpy as np
import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import OuterDesignOptimisationConfig
from tests.test_outer_design_optimisation import FakeHeatPump


@pytest.mark.parametrize('n_comp,ratios,success,violating_stage', [
    (1, [10., 99., 99.], True, None),
    (2, [10., 10., 99.], True, None),
    (2, [10.000001, 1.1, 1.], False, 1),
    (2, [1.1, 10.000001, 1.], False, 2),
    (3, [2., 3., 10.000001], False, 3),
    (3, [10., 10., 10.], True, None),
    (2, [5., None, 1.], False, 2),
    (2, [np.nan, 5., 1.], False, 1),
    (2, [np.inf, 5., 1.], False, 1),
    (2, [.9, 5., 1.], False, 1),
])
def test_each_active_stage_is_checked_and_failures_remain_in_history(n_comp, ratios, success, violating_stage):
    def analyse(hp_cls, case, **kwargs):
        hp = FakeHeatPump(.5, case)
        for i, ratio in enumerate(ratios, 1):
            if ratio is None:
                hp.output.pop(f'PR_{i}')
            else:
                hp.output[f'PR_{i}'] = ratio
        return hp, .01

    objective = _OuterDesignObjective(max_evals=1, base_case={'n_comp': n_comp},
                                     analyser=analyse, molar_masses=(44., 58.))
    x = [0., 0., .5, .5, .5]
    result = objective.evaluate(x)
    assert result.success is success
    assert objective.evaluate(x) is result
    assert objective.exact_evaluations == 1
    if not success:
        assert result.cop is None and result.sec is None
        assert result.objective >= 1e6
        assert any(k.startswith(f'stage_{violating_stage}_') for k in result.constraint_errors)
        assert result.report['solve_failure_kind'] == 'constraint_violation'
        assert result.report['solve_error_stage'] == 'compressor_pressure_ratio'
    else:
        assert result.objective == .5


@pytest.mark.parametrize('limit', [0., 1., -1., np.nan, np.inf])
def test_invalid_limits_are_rejected(limit):
    with pytest.raises(ValueError, match='max_stage_pressure_ratio'):
        OuterDesignOptimisationConfig(max_stage_pressure_ratio=limit)
    with pytest.raises(ValueError, match='max_stage_pressure_ratio'):
        _OuterDesignObjective(max_evals=1, max_stage_pressure_ratio=limit)


@pytest.mark.parametrize('limit,success', [(10., True), (8., False)])
def test_native_isolated_evaluation_propagates_limit(tmp_path, limit, success):
    # Previously solved best-COP Propane/Neopentane point: stage PRs 1.2693 and 8.5533.
    objective = _OuterDesignObjective(max_evals=1, refrigerant='Propane;Neopentane|0.5;0.5|mass',
        results_dir=tmp_path, max_stage_pressure_ratio=limit,
        exact_evaluation_timeout_seconds=90, exact_evaluation_timeout_retries=0,
        exact_evaluation_worker_max_tasks=1)
    try:
        result = objective.evaluate([26.963807397739426, 59.560832120665, .9,
                                     .42831779308637935, .5716822069136207])
        assert result.success is success, result.message
        assert result.report['max_stage_pressure_ratio_limit'] == limit
        assert result.report['PR_2'] == pytest.approx(8.553322236608999, rel=1e-6)
        for stage in (1, 2):
            for end in ('in', 'out'):
                assert result.report[f'Comp_{stage}_phase_{end}'] in {
                    'saturated_vapour', 'superheated_vapour', 'supercritical'}
                assert result.report[f'Comp_{stage}_q_refprop_{end}'] >= 1
        if not success:
            assert 'stage_2_pressure_ratio_excess' in result.constraint_errors
    finally:
        objective.close()


def test_previous_best_cop_pentane_design_is_now_rejected(tmp_path):
    objective = _OuterDesignObjective(max_evals=1, refrigerant='Propane;Pentane|0.5;0.5|mass',
        results_dir=tmp_path, exact_evaluation_timeout_seconds=90,
        exact_evaluation_timeout_retries=0, exact_evaluation_worker_max_tasks=1)
    try:
        result = objective.evaluate([1.2546887442664015e-05, 14.183255155228972,
            .1400290210468387, .18354082715025566, .8164591728497443])
        assert not result.success
        assert result.report['PR_1'] == pytest.approx(26.46035392941161, rel=1e-6)
        assert result.report['COP'] == pytest.approx(2.8293112878483484, rel=1e-6)
        assert 'stage_1_pressure_ratio_excess' in result.constraint_errors
    finally:
        objective.close()
