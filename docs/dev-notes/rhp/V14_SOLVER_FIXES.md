# v14 solver fixes and failed-candidate investigation

The minimum-pressure search now brackets the lower pinch-tolerance boundary instead of accepting an arbitrary feasible pressure. Mixture searches rebuild bounds after coupled-state changes and use the same terminal and internal pinch candidates as final validation. Pure-fluid pressure targeting retains its existing path.

Native property errors trigger at most one fresh-model retry at the superheater-limited evaporator seed. The normal seed remains the first choice, preserving the six previously identified regression cases. A final pinch shortfall no larger than 0.001 K can trigger one fresh-model retry at a tighter `1e-8` convergence tolerance. This does not relax the physical pinch limit.

An LMTD reporting flash failure now leaves `LMTD_gcs` unavailable (`NaN`) and records its explanation in `metric_errors`. It does not invalidate the independently checked cycle. Required cycle and pinch failures still reject the candidate. Reports also expose numerical retry reasons.

## Guarded two-phase PH recovery

For a mixture PH error 226, the property wrapper may invert PT between checked bubble/dew temperatures. It accepts a result only when all of these checks pass:

- The requested enthalpy lies between the saturation endpoint enthalpies, and the result has a two-phase quality and distinct positive phase densities.
- Requested and phase-weighted enthalpy agree within `1e-6 J/mol`.
- Material balance, phase composition sums and density lever rule agree within `1e-8`.
- Both independently evaluated phase pressures agree with the requested pressure within `1e-7` relative error.
- Phase fugacities agree within `1e-7` in log ratio.

The failed PH output is never reused. Native errors remain recorded, and rejected recovery states are not cached. Tests include a captured dense-profile failure and deliberate corruption of phase pressure to verify rejection.

Direct tests of 21 captured PH states found only one PT inversion that passed all independent checks: ten failed balance checks, five encountered another native convergence error, and five did not return a two-phase state. A separate captured PT error has no PH inversion to apply. This is why the recovery is guarded rather than treating native success or an enthalpy match alone as sufficient. The additional isopentane/neopentane profile state passed every check and is covered by a regression test. Error 226 denotes a two-phase iteration failure. [REFPROP documentation](https://refprop-docs.readthedocs.io/en/latest/DLL/legacy.html).

## Optimiser policies

The line-search benchmark now defaults to `--budget-policy continue`: it retains the staged production search, then spends the remaining allowance on joint coordinate polls with deterministic restarts. Previously evaluated designs are excluded from continuation proposals. `--budget-policy staged` preserves the old ceiling-only experiment. The general production optimiser's default stage allocation is unchanged; the benchmark's continuation policy is explicit.

CMA-ES now also restarts after ten generations without useful objective improvement, using a `1e-4` relative threshold. Restarts and their reasons remain reproducible from recorded observations.

For propane/pentane 30/70 at 200°C, the targeted validation runs used exactly 400 calls each:

| Method | Previous COP | Updated COP | Restarts |
|---|---:|---:|---:|
| Line-search with continuation | 2.646941 | 2.671131 | 4 continuation restarts |
| CMA-ES, seed 20260910 | 2.588390 | 2.671650 | 2 stagnation restarts |

Both winners were replayed in fresh processes on the final model. These are targeted policy checks, not a replacement for the complete three-mixture, three-seed comparison. The 400-call runs preceded the final small-pinch retry; their exact source is preserved in `pre_pinch_retry_source.zip`. Neither run had a failed evaluation. The line-search artifact's reported total was corrected to include its continuation calls; its 400 recorded observations and replay are unchanged, and the executed driver is preserved.

## Remaining candidate behaviour

Final validation on model hash `812b3bebdaa94bd653a110f3cb2120c424bbdce01c886b4768dfb38ee1bea087`:

| Check | Result |
|---|---:|
| Automated tests | 500 passed |
| Saved smoothness points, brackets and gaps | 314/314 solved |
| Twelve benchmark winners × four pressure guesses | 48/48 solved |
| Originally failed benchmark candidates | 57/75 recovered |

The remaining 18 candidates comprise 16 PH errors, one PT error (all native code 226), and one 180-second timeout. They remain explicit failures. The six regressions from unconditional evaporator seeding are absent. Every successful replay passes the final 19.95 K minimum-approach requirement.

The small final-pinch rejection had residual changes of roughly `0.00011 K` over a very narrow pressure interval. It solved from a 7 MPa starting guess and with convergence tolerance `1e-8`, while `1e-7` still failed. This motivated the bounded final-validation retry.

The slow propane/isobutane BO candidate (`b640ecc9cf_gp_ei_20260911:344`) reaches the 180-second evaluation limit during sampled pressure fallback inside evaporator alignment. A separate 40-second trace completed about 22,000 property calls. It is excessive fallback work, not evidence of one permanently frozen native call. The evaluation remains a timeout rather than a fabricated infeasibility penalty or a successful solve.

Other unresolved native errors remain failed evaluations; the investigation does not establish physical infeasibility. The guarded recovery does not cover every two-phase state.

## Reproducibility

Production cohort replays, native traces and summaries are in `results/RHPcomposition/v14_fixes`. The two 400-call searches are in `results/RHPcomposition/v14_search_validation`. `validate_v14_fixes.py` launches each candidate in a fresh process, records its full model source hash, and retains native errors and timeout stacks. Earlier investigation and benchmark artifacts remain intact.

Final validation counts and hashes are recorded in `validation_summary.json` and `manifest.json` beside the raw results.
