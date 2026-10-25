# v15 numerical-failure fixes

The implemented fixes recover **19 of the 44 previously failed candidates**. All four timeout candidates now solve successfully, and none of the 82 final validation evaluations timed out. The remaining 25 candidates return explicit REFPROP PH error 226; they remain failed evaluations, without a claim of physical infeasibility.

| Original failure | Candidates | Recovered |
|---|---:|---:|
| Native PH error 226 | 36 | 11 |
| Native timeout | 4 | 4 |
| Negligible-duty temperature reversal | 2 | 2 |
| Small final pinch shortfall | 2 | 2 |
| Total | 44 | 19 |

## Implemented behaviour

**Validate sink reconstruction before native calls.** Nonfinite/nonpositive total duty, materially negative individual duty, nonpositive enthalpy rise, and invalid mass flow raise a structured `InvalidSinkStateError`. Each reconstructed sink enthalpy is also checked against the inlet/outlet enthalpy interval before PH evaluation. This additional check catches tiny tolerated negative duties amplified by near-zero net flow. The sampled pressure fallback skips invalid probes; the cycle boundary retains the actual constraint violation.

**Use bounded numerical retries.** Final mixture pinch shortfalls up to 0.05 K beyond the existing allowance receive one fresh retry at cycle tolerance `1e-8`, replacing the previous 0.001 K retry window. Final physical acceptance still requires the original target minus its existing 0.05 K allowance. Whole-case numerical failures can then try fresh pressure seeds of 3 and 7 MPa, skipping the initial seed if it is already one of those values. Successful and clearly physical-failure cases are not automatically retried. Retry history is retained, and an exhausted retry sequence preserves the original failure cause.

**Preserve negligible-duty states.** At zero pressure drop, nonnegative heat duty below the component's existing absolute duty tolerance can preserve both inlet states instead of reflashing them. Fixed temperature-rise/drop models additionally require the requested temperature change to lie within the component's temperature tolerance. Both streams conserve energy exactly. Negative duty, pressure drops and finite heat transfer continue through their normal validation paths.

**Refine successful binary two-phase PH results when necessary.** The boundary checks exposed native PH calls reporting success while missing the requested enthalpy by roughly 1.6–4.1 × 10⁻⁵ kJ/kg. A guarded local refinement solves both phase pressures, two fugacity equalities, material balance and enthalpy balance. It starts only from a successful native two-phase state; it never seeds from a failed flash. Densities remain positive and compositions/quality bounded through variable transformations. The refinement has a 60-function-evaluation limit, must remain near the original state, and must pass the same independent phase checks used by the existing recovery path. If it cannot, the original successful native state is retained. No phase-equilibrium validation tolerance was relaxed.

**Retain diagnostics.** Failed objective records preserve the underlying message, failure kind, solver stage, property fluid/input pair and inputs, and numerical retry history. The benchmark serializer now includes these fields too.

The objective version is now `isolated-refprop-pattern-search-v15`, invalidating completed signatures from the previous model. Historical benchmark files have not been rewritten.

## Final validation

The final source was exercised in **82 fresh-process evaluations** using the original design variables and compositions, with a 200-second diagnostic deadline per complete evaluation:

- 44 original failed candidates: 19 successes, 25 explicit PH failures, zero timeouts.
- 18 original benchmark winners: all successful.
- 12 sampled successful controls, four per composition: all successful.
- Eight additional near-zero probes around the two boundary candidates: all successful.

The four previous timeout cases now solve in approximately **2.2, 22.6, 52.4 and 1.35 seconds**, respectively. These are instrumented concurrent replay times; the old benchmark spent approximately 362 seconds per timeout including its retry. The slowest final replay, including bounded numerical retries, took about 114.3 seconds.

| Metric | Largest relative change among 18 winners |
|---|---:|
| SEC | 5.74 × 10⁻⁸ |
| COP | 3.48 × 10⁻⁸ |
| Discharge pressure | 7.12 × 10⁻⁸ |
| VHC | 1.26 × 10⁻⁷ |

All winner changes remain below `1e-6` relative tolerance. The largest relative COP change is **0.00000348%**. The twelve sampled controls also remain within `1e-6` relative tolerance for all four metrics.

![Near-zero validation](../../../results/RHPcomposition/v15_failure_fixes/verified/near_zero.png)

The near-zero sweep tests RHX2 lifts of 0, 10⁻⁸, 10⁻⁶ and 10⁻⁵ K at each boundary candidate. COP variation across each sweep is less than 7 × 10⁻⁸ in absolute COP. These checks demonstrate recovery of the sampled failure regions, not a proof of global differentiability.

The full suite passes: **544 tests**. After extending the benchmark serializer, **30 targeted diagnostic and resume tests** also pass; these overlap the full suite. Ruff passes on all files changed for this work. The sole pytest warning is an existing dependency warning about PyTorch's `torch.jit.script` deprecation.

## Remaining limitations

Twenty-five candidates still fail validated native PH calculations after the bounded retries. The new successful-PH refinement does not accept or polish failed native output, and existing failed-PH recovery checks remain in force. Further recovery work would need another validated flash or pressure-search path. Failed samples should continue to be treated as unavailable numerical results rather than automatic proof of physical infeasibility.

The bounded retries increase effort on difficult evaluations. The normal isolated objective timeout still applies to the complete candidate evaluation, including retries. No new 400-evaluation optimiser comparison was run in this change; method rankings should be measured again against v15 before interpreting performance changes.

## Evidence and implementation

- [Final numerical summary](../../../results/RHPcomposition/v15_failure_fixes/verified/summary.json)
- [All 44 candidate outcomes](../../../results/RHPcomposition/v15_failure_fixes/verified/candidate_results.md)
- [Full tests](../../../results/RHPcomposition/v15_failure_fixes/verified/full_tests.log)
- [Serialization and resume checks](../../../results/RHPcomposition/v15_failure_fixes/verified/serialization_tests.log)
- [Frozen thermodynamic source](../../../results/RHPcomposition/v15_failure_fixes/verified/source_manifest.json)
- [Evidence manifest](../../../results/RHPcomposition/v15_failure_fixes/verified/evidence_manifest.json)
- [Sink guards](../../../HeatPumpLab/cycles/RecuperativeHeatPump/solver/_cycle_analysis.py), [cycle retries](../../../HeatPumpLab/cycles/RecuperativeHeatPump/solver/cycle_solver.py), [pressure-seed retries](../../../HeatPumpLab/cycles/RecuperativeHeatPump/utils/analysis.py)
- [Exchanger state handling](../../../HeatPumpLab/components/common/base_heat_exchanger.py), [validated PH refinement](../../../HeatPumpLab/utils/properties_package.py)
- [Objective diagnostics](../../../HeatPumpLab/cycles/RecuperativeHeatPump/optimisation/_objective.py), [benchmark serializer](../../../main_RHPoptimiser_benchmark.py)

The `verified` directory contains the final measurements and frozen source. Earlier files under `v15_failure_fixes` and its `final` subdirectory are intermediate diagnostic runs, superseded by this report. Raw requests, results, trace logs and statuses are retained. Rebuild the final summary with `python main_RHPv15_report.py` from the repository root.
