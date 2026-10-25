# Nearby-state PH interpolation experiment

**Yes: nearby successful flashes can recover these errors.** Of 25 representative failed PH states, 14 admit a directly validated interpolation and 20 admit recovery when interpolation seeds the existing local equilibrium refinement. Four of five selected complete heat-pump replays also succeed with this experimental fallback.

This is a diagnostic experiment. Production v15 recovery and the historical benchmark results are unchanged.

## Property-state results

One stored failed state was selected from each of the 25 remaining failed candidates: the first property context in its v15 retry history. All 25 reproduce native PH error 226 in a fresh process.

| Check | States out of 25 |
|---|---:|
| Independently valid nearby states bracket the target | 21 |
| Interpolation alone passes all existing phase checks | 14 |
| Interpolation with local equilibrium refinement recovers the target | 20 |
| Still unrecovered by the tested procedure | 5 |

The interpolation-only cases are included in the 20 recoveries. These counts refer to individual PH states, not completed cycle evaluations.

At constant pressure, the experiment samples enthalpy offsets of ±10⁻⁵, ±10⁻⁴, ±0.001, ±0.01, ±0.1, ±1 and ±5 kJ/kg. At constant enthalpy it samples relative pressure offsets of ±10⁻⁶, ±10⁻⁵, ±10⁻⁴, ±0.001 and ±0.01. Both symmetric and unequal-distance bracketing pairs are tested. Symmetric pairs alone recover 9 states by interpolation and 17 with refinement; unequal-distance pairs increase those counts to 14 and 20.

Only successful native two-phase endpoints that pass independent validation are used. Successful endpoints may first receive the existing local precision refinement. The experiment interpolates temperature, phase densities, phase compositions, quality, enthalpy, entropy and internal energy. Bulk density is reconstructed from the phase specific volumes. Weights place the result at the original requested pressure and enthalpy; there is no extrapolation.

Every candidate must satisfy the existing pressure, fugacity, material, density and enthalpy checks. If interpolation does not pass, it seeds the existing binary equilibrium solver at the exact requested pressure and enthalpy. Its iteration limit, close-state guard and validation tolerances remain unchanged.

Four states have no independently valid bracketing pair in this search. One has a bracket but yields no accepted interpolation/refinement. Native success alone is insufficient: some nearby native calls return success but fail the phase checks.

## Why raw interpolation needs checks

Interpolation does not generally preserve nonlinear phase equilibrium. For example, a ±1% pressure pair around the failed state at 3.927491 MPa and 645.631920 kJ/kg produces an interpolated state with a **0.30677 kJ/kg phase enthalpy-balance error**, a **0.03846% phase-pressure error**, and a material-balance error of about **1.77 × 10⁻⁴**. This result is rejected. Close pairs can pass, but a numerical value alone is not sufficient evidence of recovery.

## Complete-cycle checks

Five targeted candidates were replayed in fresh processes with a 200-second deadline. They include both propane/pentane 60% mass-fraction candidates, a propane/isobutane candidate requiring unequal bracketing, one requiring pressure bracketing, and one whose original failed state remains unrecovered. This selection is diagnostic rather than a random estimate of recovery rate.

| Candidate | Outcome | COP | Solve seconds | Recovered PH calls / attempts |
|---|---|---:|---:|---:|
| `6da99772a7_botorch_log_ei_20260910_295` | Success | 2.46483974 | 28.6 | 14 / 14 |
| `6da99772a7_gp_ei_20260909_319` | Success | 2.42164483 | 61.7 | 32 / 33 |
| `b640ecc9cf_botorch_log_ei_20260909_128` | Success | 2.50207350 | 44.8 | 7 / 8 |
| `b640ecc9cf_botorch_log_ei_20260910_331` | Success | 2.51691505 | 27.1 | 7 / 7 |
| `b640ecc9cf_botorch_log_ei_20260911_275` | PH error 226 | — | 64.6 | 4 / 10 |

All four successful cycles pass the existing feasibility checks, including the gas-cooler minimum approach. None of the five replays times out. Timings are instrumented concurrent runs, not standalone performance benchmarks. The fifth case advances through several recovered property calls but ultimately still fails.

## Recommended implementation direction

Use nearby interpolation as a bounded fallback and as a starting point for an equilibrium correction. Retain strict acceptance checks, prefer close valid brackets, record the recovery path, and cap the number of native calls. The experiment supports this mechanism; it does not establish recovery of all 25 complete cycles or global smoothness. Wider validation and regression checks are needed before enabling it in production.

## Evidence

- [Diagnostic implementation](../../../main_RHPph_interpolation_probe.py)
- [Numerical summary and cycle outcomes](../../../results/RHPcomposition/ph_interpolation_probe/summary.json)
- [All 25 state outcomes](../../../results/RHPcomposition/ph_interpolation_probe/state_outcomes.md)
- [Source manifest](../../../results/RHPcomposition/ph_interpolation_probe/source_manifest.json)

Raw endpoint results and balance residuals are retained beside the summary. Cycle requests, results, traces, logs, deadlines and fallback counts are in its `cycles` subdirectory. Run `python main_RHPph_interpolation_probe.py` to repeat the 25 state scans, or `python main_RHPph_interpolation_probe.py --cycle REQUEST_PATH` to repeat a saved cycle request.
