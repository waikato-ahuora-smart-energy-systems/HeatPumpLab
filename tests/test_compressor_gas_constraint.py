"""Gas-only constraints use native phase evidence at all active compressor ends."""

import math

import pytest

from HeatPumpLab.components.stream import Stream
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
from HeatPumpLab.utils.properties_package import PropertyPackage, phase_from_refprop_quality
from tests.test_outer_design_optimisation import FakeHeatPump


@pytest.mark.parametrize('quality,phase', [
    (None, 'unknown'), (math.nan, 'unknown'), (math.inf, 'unknown'),
    (-9999990., 'unknown'), (1000., 'unknown'),
    (-998., 'liquid'), (-.1, 'liquid'), (0., 'saturated_liquid'),
    (.99, 'two_phase'), (1., 'saturated_vapour'),
    (1.001, 'superheated_vapour'), (998., 'superheated_vapour'),
    (999., 'supercritical'),
])
def test_native_phase_flags(quality, phase):
    assert phase_from_refprop_quality(quality) == phase


@pytest.mark.parametrize('n_comp,stage', [(1, 1), (2, 1), (2, 2), (3, 3)])
@pytest.mark.parametrize('end', ['in', 'out'])
@pytest.mark.parametrize('quality,success', [
    (1., True), (998., True), (999., True),
    (.999999999, False), (0., False), (-998., False),
    (None, False), (math.nan, False), (-9999990., False),
])
def test_each_active_compressor_end(n_comp, stage, end, quality, success):
    def analyse(hp_cls, case, **kwargs):
        hp = FakeHeatPump(.5, case)
        field = f'Comp_{stage}_q_refprop_{end}'
        if quality is None:
            hp.output.pop(field)
        else:
            hp.output[field] = quality
        # A stale/inconsistent display label must never override raw evidence.
        hp.output[f'Comp_{stage}_phase_{end}'] = 'superheated_vapour'
        for inactive in range(n_comp+1, 4):
            hp.output[f'Comp_{inactive}_q_refprop_in'] = .5
            hp.output[f'Comp_{inactive}_q_refprop_out'] = .5
        return hp, .01

    objective = _OuterDesignObjective(max_evals=1, base_case={'n_comp': n_comp},
                                     analyser=analyse, molar_masses=(44., 58.))
    result = objective.evaluate([0., 0., .5, .5, .5])
    assert result.success is success
    assert objective.exact_evaluations == 1
    assert objective.evaluate(result.x) is result
    if not success:
        assert result.cop is None and result.sec is None
        assert result.objective >= 1e6
        assert result.constraint_errors == {f'stage_{stage}_{end}_not_gas': 1.}
        assert result.report['solve_error_stage'] == 'compressor_phase'
        assert result.report['solve_failure_kind'] == 'constraint_violation'
        assert result.report['COP'] == 2.5


@pytest.mark.parametrize('pair,a,b,phase', [
    ('PQ', .5, .99, 'two_phase'),
    ('PQ', .5, 1., 'saturated_vapour'),
    ('PT', .1, 20., 'superheated_vapour'),
    ('PT', 1., 0., 'liquid'),
    ('PT', 5., 150., 'supercritical'),
])
def test_native_phase_survives_property_conversion_stream_update_and_copy(pair, a, b, phase):
    pp = PropertyPackage('Propane')
    props = pp.get_std_properties(pair, a, b, pp.fld_composition)
    assert props['phase'] == phase
    assert phase_from_refprop_quality(props['q_refprop']) == phase
    stream = Stream()
    stream._update_props(props)
    copied = Stream()
    copied.copy_fld_props(stream)
    assert copied.phase == phase
    assert copied.refprop_quality == props['q_refprop']
    # Preserve callers' legacy convention for undefined single-phase quality.
    assert copied.quality == (props['q_refprop'] if 0 <= props['q_refprop'] <= 1 else pp.DEFAULT_VAL)
