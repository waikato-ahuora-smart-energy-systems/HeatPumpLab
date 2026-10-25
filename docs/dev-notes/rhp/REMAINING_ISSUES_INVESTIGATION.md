# Investigation following the v13 optimiser retest

The pressure-initialisation sensitivity is an inner-solver bug. The BO failure is a native two-phase flash failure encountered along an intermediate solve path. Replays also uncovered an independent problem: optional LMTD reporting can invalidate an otherwise converged, pinch-checked cycle. Experimental remedies are promising, but an unconditional evaporator seed introduces regressions and should not be shipped.

Production code was held at the benchmark model hash `75aff32de197b6a1daa23d8453685c0e15deae167a9ae68d74649dd4bee42b11`. All changes for this investigation are isolated diagnostic scripts; the earlier optimiser results have not been overwritten or rescored.

## Pressure minimisation stops before finding the minimum

For the 60/40 propane/pentane CMA-ES winner, the existing solver returns:

| Initial pressure, MPa (0 = automatic) | Selected pressure, MPa | COP |
|---:|---:|---:|
| 0 | 3.883477 | 2.538172 |
| 3 | 3.883477 | 2.538172 |
| 5 | 4.169764 | 2.495057 |
| 7 | 5.822657 | 2.321442 |
| 10 | 8.301996 | 2.175350 |

`expand_discharge_pressure_bounds` stops once the pinch residual enters the ±0.05 K acceptance band. `_select_accepted_endpoint` can then accept a higher pressure without a lower infeasible point having been bracketed. With a 7 MPa guess, residuals are positive but only about 0.003–0.006 K. These are feasible points, but feasibility does not establish minimum pressure. Cached bounds and the zero-residual bracketing test also differ from the interpolation target at the lower tolerance boundary.

The experimental `band` mode searches for a bracket around the actual target, −0.04999 K, rebuilds it after coupled-state changes, and uses Brent interpolation in log pressure. This retains the existing 0.05 K physical acceptance allowance. It gives COP 2.538172 for all five guesses in the affected 60/40 design.

Evidence: `screen_v2_results.json`; previous `optimiser_v13_400/pressure_trace_0.json`, `pressure_trace_7.json`, and `pressure_followup.json`.

## BO encounters an invalid intermediate PH flash

The propane/isobutane BO winner fails from a 3 MPa starting guess inside the gas-cooler internal pinch minimisation. The captured native inputs are:

- Pressure: 3869.8048796042 kPa.
- Enthalpy: 30438.25761068768 J/mol.
- Propane/isobutane mole fractions: 0.36098059287265455 / 0.6390194071273454.
- REFPROP error 226: two-phase iteration did not converge.

The failed flash returns enthalpy 30397.53377282915 J/mol, about 40.724 J/mol below the requested value. Reusing that output would violate the requested state. The native documentation identifies 226 as a two-phase convergence error; it does not establish physical infeasibility. [REFPROP DLL documentation](https://refprop-docs.readthedocs.io/en/latest/DLL/legacy.html).

This pressure is encountered before evaporator alignment. The final successful operating point is around 3.94324 MPa. Tightening pressure bracketing alone makes all five tested BO starting guesses encounter a failing intermediate flash. Starting at the existing superheater temperature limit minus the alignment tolerance avoids this path for the winner. It also makes the affected 60/40 case independent of the pressure guess.

A direct, fresh-process PH reproduction fails with the same enthalpy discrepancy. A naive quality-based inversion fails after three PQ calls with error 226 as well. An unchecked PQ inversion or forced single-phase flash is therefore not an established recovery method.

Evidence: `screen_v2_results.json`, `seeded_results.json`, `native_flash_probe.json`; scripts `investigate_v13_remaining_issues.py` and `probe_v13_native_flash.py` in `paper/analysis`.

## Full replay separates cycle failures from reporting failures

All 75 failed benchmark candidates fail again under fresh-process baseline replays. Their terminal failure locations are:

| Failure location | Baseline | Pressure bracket + upper evaporator seed |
|---|---:|---:|
| Internal pinch calculation | 38 | 17 |
| Cycle calculation | 13 | 3 |
| LMTD reporting | 20 | 9 |
| Final pinch tolerance validation | 4 | 0 |
| **Total failures** | **75** | **29** |

The baseline native failures comprise 70 PH errors and one PT error, all code 226. The four other failures miss the final pinch acceptance boundary by approximately 0.000001–0.000098 K.

The combined prototype recovers 46 candidates, including the four final-pinch failures. Nine remaining candidates have converged operations and a valid final pinch, then fail in `populate_full_metrics → get_LMTD_ave → get_graph_data`, which samples additional PH states for reporting. A separate diagnostic catches only REFPROP errors from this LMTD calculation and marks LMTD unavailable. All nine then produce successful objective evaluations, with minimum approaches at least 19.95000999975 K. No failed state is reused and no replacement LMTD value is invented.

Thus 55 of the 75 originally failed candidates can be evaluated by these experimental changes. The remaining 20 fail during internal pinch (17) or cycle calculations (3); those are still unresolved native numerical failures. These counts describe a diagnostic recovery experiment, not a validated production fix or a rerun of the optimisation algorithms.

Evidence: `failed_candidates_results.json`, `optional_metrics_results.json`; terminal stacks and failed native inputs are included in each record.

## Regression checks reject unconditional upper-temperature seeding

The combined prototype passes **48/48 winner checks**: all 12 method/composition winners, each from pressure guesses 0, 3, 7 and 10 MPa. The largest within-design COP spread is 0.0000000502%. This demonstrates that the pressure sensitivity can be removed at the tested winners.

However, replaying the saved **314-point smoothness set** gives **307 successful solves**, compared with 312 under production v13:

- One of the two existing native failure points now succeeds; one still fails.
- Six previously successful propane/isobutane points now fail, at centres `2_15` and `2_16` (three nearby samples each).
- At those six points, the upper evaporator seed produces an RHX2 temperature crossing before evaporator alignment. The outer solve rejects the intermediate state before alignment can reach the valid state.
- Repeating all seven failures with the pressure-bracketing change and the normal evaporator seed recovers all six regressions. The remaining original native flash failure persists.

All 12 saved discontinuity-bracket endpoints succeed; the largest COP difference across those narrow brackets is 0.0000024061%. All 14 saved feasibility-gap points succeed. These finite probes do not establish smoothness everywhere, and the six regressions prohibit adopting the combined prototype unchanged.

Evidence: `winners_results.json`, `smoothness_results.json`, `regression_results.json`, `summary.json`.

## Optimiser behaviour is a separate issue

The deterministic line-search stops after 66/89/89 actual calls because its stages finish. `max_exact_evaluations=400` is a ceiling: line and split budgets are increased only to 12, while polish and lift-feedback stages retain their own caps. Spare evaluations are not repeatedly allocated to another joint search. Its result is therefore not an equal-400-call comparison. A stronger 400-call variant requires an explicit budget allocation and continuation policy.

The weak 30/70 propane/pentane CMA-ES run (seed 20260910) makes no failed evaluations and records no native stopping/restart event. Its incumbent approaches RHX1 = RHX2 = 0 K, compressor fraction = 0.9, COP 2.588390. At 200 calls its best COP is already 2.588375; the last 200 calls gain only about 0.000015 COP. Other seeds reach 2.6723. The history supports stagnation in a weaker region; it does not prove a local or global optimum. Native-stop-only restarts do not address this within a 400-call budget.

## Recommended implementation order

1. Make the mixture pressure bracketing and interpolation target consistent, and validate the refreshed final state against the same tolerance. Preserve the pure-fluid path and error diagnostics.
2. Keep the normal evaporator initialisation. Add a bounded alternate-seed retry for appropriate numerical failures using a fresh model state, rather than universally imposing the upper-temperature seed. Test both winner cases and the six recovered regressions. The retry policy itself still needs implementation and broader validation.
3. Separate optional LMTD reporting failure from cycle/objective validity. Surface the unavailable metric and its native error; retain failures in required cycle and pinch calculations.
4. Investigate the remaining two-phase states with independently checked phase and material/energy balances before introducing a general flash fallback. Current evidence does not validate one.
5. Once the objective is stable, add an explicit 400-call continuation policy for line-search and budget-aware stagnation/restarts for CMA-ES, then repeat the benchmark with new versioned outputs.

Raw artifacts live in `results/RHPcomposition/v13_remaining_issues`. `manifest.json` records source and result hashes, and `experiment_sources` preserves the experimental drivers. The initial `screen_results.json` is excluded: it contains harness serialization failures from the first attempt; use `screen_v2_results.json` instead.
