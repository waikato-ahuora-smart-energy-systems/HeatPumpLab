# v13 optimiser retest: 400-call allowance

Completed 9 September 2026 at a 200°C sink target, with fixed mass composition and origin RHX1 = RHX2 = 0 K, compressor 2 power fraction = 0.5. RBF, BO and CMA-ES each ran three seeds (20260909–20260911), with exactly 400 model calls per seed. The unchanged deterministic production line-search had a 400-call ceiling and finished after 66, 89 and 89 calls. Its internal stage caps were retained.

## Median COP under default numerical initialisation

| Mixture by mass | Line-search | RBF | BO (GP-EI) | CMA-ES |
|---|---:|---:|---:|---:|
| Propane/pentane 30/70 | 2.646941 | 2.672332 | 2.671888 | 2.672328 |
| Propane/pentane 60/40 | 2.522823 | 2.533022 | 2.537141 | 2.538127 |
| Propane/isobutane 30/70 | 2.511256 | 2.538312 | 2.540916 | 2.541161 |

RBF and CMA-ES have essentially the same median on propane/pentane 30/70, but CMA-ES has a much wider seed range (2.588390–2.672346 versus RBF 2.672285–2.672337). CMA-ES has the highest median on the other two mixtures. BO has a higher worst-seed COP than CMA-ES on both 30/70 mixtures. These observations describe the tested configurations and three seeds.

## Remaining numerical issues

- The best CMA-ES propane/pentane 60/40 design changes from COP 2.538172 at 3.883477 MPa under automatic initialisation to COP 2.321442 at 5.822657 MPa when starting the pressure search at 7 MPa. All three CMA-ES winners for this mixture show the same sensitivity.
- Instrumented direct solves show pressure-bound expansion and endpoint acceptance stopping inside the ±0.05 K pinch tolerance band, before a lower infeasible pressure is bracketed. Both final states pass the pinch check. A lower feasible pressure has already been demonstrated, so the accepted higher-pressure state does not fulfil the intended minimum-pressure search.
- Disabling the fast inner objective at the 7 MPa guess recovers COP 2.537801, pressure 3.882546 MPa and a 0.045038 K gas-cooler offset. The solver should consistently locate the lower edge of the feasible pressure band before treating its output as an initialisation-independent optimum.
- The best BO propane/isobutane design passes its default and 7 MPa replays but encounters native REFPROP PH error 226 from a 3 MPa pressure guess.

## Validation and artifacts

Validated 11,044 optimisation calls, including 75 failed calls counted against budgets, and all 30 fresh winner replays. All stochastic budgets, origins, compositions and bounds passed validation. Of 24 alternate-pressure checks, 23 solved. Twelve further objective calls and two direct cycle traces investigate the pressure sensitivity; none trained the optimisers. The optimiser tests also passed (16 tests), plus 400-proposal feasibility/uniqueness checks for each stochastic method.

[Detailed report, all ranges, selected designs and diagnostics](../../../results/RHPcomposition/optimiser_v13_400/retest_report.md)
[Convergence plot](../../../results/RHPcomposition/optimiser_v13_400/convergence.png)
[All trial results](../../../results/RHPcomposition/optimiser_v13_400/trial_results.csv)
[Validation](../../../results/RHPcomposition/optimiser_v13_400/validation.json)
[Reproduction manifest](../../../results/RHPcomposition/optimiser_v13_400/benchmark_manifest.json)

All runs used the same v13 model source hash: `75aff32de197b6a1daa23d8453685c0e15deae167a9ae68d74649dd4bee42b11`. The exact model and benchmark drivers are archived in `source_snapshot.zip` alongside the results. Historical input rows are preserved in `reference_cases.csv`; their objective values were not reused.
