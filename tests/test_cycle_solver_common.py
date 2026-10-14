from types import SimpleNamespace

import numpy as np
import pytest

from HeatPumpLab.cycles.common import (
    CycleFailure,
    calculate_cop,
    calculate_sec,
    cycle_convergence_failure,
    cycle_has_converged,
    find_lowest_log_bracketed_target,
    find_sampled_log_target,
)


def test_cycle_failure_keeps_only_finite_nonnegative_constraint_errors():
    failure = CycleFailure(
        "failed",
        {"negative": -2.0, "valid": 3.0, "invalid": np.inf},
    )

    assert failure.message == "failed"
    assert failure.constraint_errors == {"negative": 0.0, "valid": 3.0}


def test_cycle_convergence_requires_feasibility_operations_and_positive_work():
    hp = SimpleNamespace(
        feasible=True,
        iteration_converged=True,
        operations=[SimpleNamespace(converged=True)],
        compressors=[SimpleNamespace(power=2.0)],
    )

    assert cycle_has_converged(hp) is True
    assert cycle_convergence_failure(hp) is None


def test_cycle_convergence_reports_heat_exchanger_reversal():
    exchanger = SimpleNamespace(
        converged=False,
        source_in=SimpleNamespace(temperature=40.0),
        source_out=SimpleNamespace(temperature=42.0),
        sink_in=SimpleNamespace(temperature=20.0),
        sink_out=SimpleNamespace(temperature=30.0),
    )
    hp = SimpleNamespace(
        feasible=True,
        iteration_converged=True,
        operations=[exchanger],
        compressors=[SimpleNamespace(power=2.0)],
    )

    failure = cycle_convergence_failure(hp)

    assert failure is not None
    assert failure.constraint_errors == {
        "heat_exchanger_temperature_reversal": pytest.approx(2.0)
    }


def test_common_metrics_preserve_converged_and_failed_cycle_behaviour():
    assert calculate_cop(10.0, 2.0, converged=True) == pytest.approx(5.0)
    assert calculate_sec(10.0, 2.0, converged=True) == pytest.approx(0.25)
    assert calculate_cop(10.0, 2.0, converged=False) == 0.0
    assert calculate_sec(10.0, 2.0, converged=False) == np.inf


def test_lowest_log_bracketed_target_finds_entry_into_tolerance_band():
    result = find_lowest_log_bracketed_target(
        lambda value: -0.2 if value < 1.2 else 0.01,
        (1.0, 2.0),
        (-0.2, 0.01),
        target_tolerance=0.05,
        relative_tolerance=1e-6,
    )

    assert result.success is True
    assert result.value == pytest.approx(1.2, rel=1e-5)
    assert result.residual <= 0.05


def test_sampled_log_target_finds_hidden_non_monotonic_target():
    def residual(value):
        log_distance = np.log(value / 3.0)
        return 0.02 - (log_distance / 0.04) ** 2

    result = find_sampled_log_target(
        residual,
        (1.0, 10.0),
        4.0,
        target_tolerance=0.05,
        x_tolerance=1e-6,
    )

    assert result.success is True
    assert 2.8 < result.value < 3.2
    assert result.residual <= 0.05


def test_sampled_log_target_skips_failed_evaluations_and_reports_shortfall():
    def residual(value):
        if value < 2.0:
            raise ValueError("property failure")
        return -0.75

    result = find_sampled_log_target(
        residual,
        (1.0, 10.0),
        4.0,
        target_tolerance=0.05,
        x_tolerance=1e-6,
    )

    assert result.success is False
    assert result.value == pytest.approx(4.0)
    assert result.residual == pytest.approx(0.75)


def test_sampled_log_target_rejects_a_root_on_a_discontinuity():
    # the residual jumps from -5 to +5 at 2.0 (no root there) and enters the band properly near 6.0
    def residual(value):
        if value < 2.0:
            return -5.0
        if value < 5.0:
            return 5.0
        return 0.02 * (value - 6.0)

    result = find_sampled_log_target(residual, (1.0, 10.0), 4.0, target_tolerance=0.05, x_tolerance=1e-6)

    assert result.success is True
    assert result.value > 4.9 and result.residual <= 0.05


def test_sampled_log_target_propagates_programming_errors():
    def residual(value):
        raise TypeError("bad call")

    with pytest.raises(TypeError):
        find_sampled_log_target(residual, (1.0, 10.0), 4.0, target_tolerance=0.05, x_tolerance=1e-6)
