# Recuperative heat pump: development notes

Notes written while developing the recuperative heat pump (RHP) solver and its design optimiser, kept as a record
of what was tried and why. The user documentation is in
[`HeatPumpLab/cycles/RecuperativeHeatPump/`](../../../HeatPumpLab/cycles/RecuperativeHeatPump/)
([ABOUT.md](../../../HeatPumpLab/cycles/RecuperativeHeatPump/ABOUT.md) and
[BAYESIAN_OPTIMISATION.md](../../../HeatPumpLab/cycles/RecuperativeHeatPump/BAYESIAN_OPTIMISATION.md)).

Some notes link to result files under `results/` and to runner scripts that are not part of the public
repository; those links do not resolve here.

| Topic | Notes |
|---|---|
| Optimiser benchmarks | [OPTIMISER_COMPARISON.md](OPTIMISER_COMPARISON.md), [ORIGIN_OPTIMISER_COMPARISON.md](ORIGIN_OPTIMISER_COMPARISON.md), [OPTIMISER_RETEST_V13_400.md](OPTIMISER_RETEST_V13_400.md) |
| BoTorch versus the custom BO | [BOTORCH_PROTOCOL.md](BOTORCH_PROTOCOL.md), [BOTORCH_COMPARISON.md](BOTORCH_COMPARISON.md), [BOTORCH_FAILURE_INVESTIGATION.md](BOTORCH_FAILURE_INVESTIGATION.md) |
| Trust-region methods | [TRUST_REGION_PROTOCOL.md](TRUST_REGION_PROTOCOL.md), [TRUST_REGION_COMPARISON.md](TRUST_REGION_COMPARISON.md) |
| Solver robustness | [ROBUSTNESS_INVESTIGATION.md](ROBUSTNESS_INVESTIGATION.md), [REMAINING_ISSUES_INVESTIGATION.md](REMAINING_ISSUES_INVESTIGATION.md), [NEAR_ZERO_INVESTIGATION.md](NEAR_ZERO_INVESTIGATION.md), [PH_INTERPOLATION_INVESTIGATION.md](PH_INTERPOLATION_INVESTIGATION.md), [SMOOTHNESS_FIXES.md](SMOOTHNESS_FIXES.md) |
| Solver fixes by evaluation version | [V14_SOLVER_FIXES.md](V14_SOLVER_FIXES.md), [V15_SOLVER_FIXES.md](V15_SOLVER_FIXES.md), [V17_INTERNAL_PH_FALLBACK.md](V17_INTERNAL_PH_FALLBACK.md) |
