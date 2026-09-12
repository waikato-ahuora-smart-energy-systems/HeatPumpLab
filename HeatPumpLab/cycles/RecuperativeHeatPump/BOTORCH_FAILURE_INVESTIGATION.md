# Failed-candidate investigation: BoTorch and custom BO

All **44 failed candidates reproduce** from fresh processes. Most are numerical failures, and **14 solve successfully when only the initial pressure changes**. The original benchmark results remain unchanged. Investigation date: 10 September 2026.

| Cause | BoTorch | Custom BO | Total |
|---|---:|---:|---:|
| REFPROP PH flash error 226 | 19 | 17 | 36 |
| Native timeout | 2 | 2 | 4 |
| Tiny recuperator temperature reversal | 2 | 0 | 2 |
| Final gas-cooler pinch shortfall | 0 | 2 | 2 |
| Total | 23 | 21 | 44 |

Thirty-four failures are in Propane/Isobutane 30/70 by mass; 27 of those occur at compressor-2 fraction 0.9. The late BoTorch cluster is predominantly around RHX1 55–63 K and RHX2 19–24 K. These are numerical failure regions; unsuccessful evaluations alone do not establish physical infeasibility.

![Failed candidate locations](../../../results/RHPcomposition/botorch_v14_failures/failure_locations.png)

This is a two-dimensional projection across all compressor fractions, not a feasibility map at one fixed fraction. Green circles identify candidates recovered by pressure-seed or tolerance sensitivity checks.

## 1. All four timeouts reach invalid air property inputs

The unbracketed pressure fallback probes very low discharge pressures, near suction pressure. In these cases, the sequential component solve produces a negative gas-cooler duty that nearly cancels the other duties. Despite the unconverged/invalid state, `set_gas_cooler_sink_temperatures` divides by the resulting tiny sink mass flow and calls the air PH flash.

Captured examples have gas-cooler duty around −10 to −11 kW, sink flow around 4–5 × 10⁻⁷ kg/s, and air enthalpy around 2.2–2.65 × 10⁷ kJ/kg. All four eventually enter `ABFLSHdll("PH")` for Air at 101 kPa with native molar enthalpy approximately 6.26–7.67 × 10⁸ J/mol. The call does not return before the diagnostic deadline. The low-pressure probes are approximately 0.063 MPa for the pentane case and 0.345 MPa for the isobutane cases.

This is a missing validation step before a native call. Increasing the timeout or repeating the same initialisation does not address it. The original benchmark spent about 362 seconds on each candidate, including its retry.

A **diagnostic-only guard** rejects nonpositive total duty or gas-cooler duty below its existing negative-duty tolerance before sink reconstruction. The sampled pressure search can then skip that invalid probe. All four original stalls disappear:

| Candidate | Guarded solve time | Remaining outcome |
|---|---:|---|
| Isobutane, BoTorch seed 20260909, evaluation 308 | 6.6 s | Final pinch shortfall 0.004141 K |
| Isobutane, BoTorch seed 20260911, evaluation 262 | 18.5 s | Final pinch shortfall 0.033616 K |
| Isobutane, custom BO seed 20260911, evaluation 37 | 88.6 s | Final pinch shortfall 0.00000329 K |
| Pentane 30/70, custom BO seed 20260910, evaluation 162 | 7.8 s | Final pinch shortfall 0.002411 K |

Shortfalls in this table are beyond the existing 0.05 K pinch allowance. The guard prevents the unsafe intermediate calculation; it does not by itself make these designs converge. Three known BoTorch winners, one per composition, retain **exactly identical SEC, COP, pressure and VHC** with the guard. This is a promising prototype, not full regression validation or a production change.

Relevant code: [sink reconstruction](solver/_cycle_analysis.py), [pressure fallback](solver/_discharge_pressure.py), and [sampled scalar search](../common/scalar_search.py).

## 2. Native PH errors create numerical failure regions

All 36 ordinary native failures repeat error 226 in a fresh process. The traces show failed refrigerant PH calculations during internal gas-cooler pinch evaluation, including pressure interpolation and coupled evaporator realignment. The existing evaporator-seed retry does not rescue these default-seed cases.

Detailed property tracing of two representatives explains why the existing PH recovery rejects the result:

- Propane/Pentane 60/40, BoTorch seed 20260910, evaluation 295: one saturation-bounded PT recovery fails the independent equilibrium checks. Normalised vapor-pressure and fugacity errors reach approximately 5.27 and 13.94 times their acceptance limits. Another attempt fails with a native PT error 226.
- Propane/Isobutane 30/70, BoTorch seed 20260911, evaluation 376: the temperature root finder reports convergence, but the returned PT state is single phase (`q=998`) rather than the required two-phase state. It is correctly rejected.

**Keep those thermodynamic validation checks.** Numerical root convergence or a native success code alone is insufficient. Candidate fixes should change the route to a validated solution—for example, better pressure bracketing, fresh coupled state at retries, or another independently validated flash formulation.

With unchanged composition, design variables and physical tolerances, pressure seeds of 3 and 7 MPa recover 7 and 13 candidates respectively, **14 distinct candidates** in total. These include nine PH-error candidates, three timeouts and both final-pinch failures. Eleven belong to custom BO and three to BoTorch. Other cases continue to fail or hit the sensitivity deadline; this limited seed sweep does not prove them infeasible.

Every recovered result passes the model's final convergence checks and has minimum gas-cooler approach at least 19.95 K. None improves the best observed COP for its composition. Correcting their evaluations could nevertheless alter the subsequent optimisation trajectories, so the archived method rankings should not be silently recomputed using these replays.

Relevant code: [PH recovery and validation](../../utils/properties_package.py), [cycle retry](solver/cycle_solver.py), and [pinch property calculations](solver/_gas_cooler_pinch.py).

## 3. Two pinch failures recover with tighter coupled convergence

The custom-BO Propane/Pentane 30/70 candidates at evaluations 77 (seed 20260910) and 41 (seed 20260911) miss the permitted pinch boundary by 0.001804 K and 0.007703 K. Repeating with cycle tolerance `1e-8` recovers both, giving COP 2.54853367 and 2.51780961, with final approach about 19.95001 K.

The current automatic tighter-tolerance retry only triggers for shortfall at most 0.001 K, so both cases fall outside that recovery window. A wider, bounded numerical retry is worth testing. The final physical acceptance limit should remain unchanged.

## 4. Two almost-zero-duty recuperators fail direction checks

The BoTorch Propane/Pentane 30/70 candidates at evaluations 106 (seed 20260909) and 137 (seed 20260911) have recuperator duties of approximately 3.44 × 10⁻¹¹ and 1.73 × 10⁻⁷ kW. Their source outlet temperature rises by 1.2347 × 10⁻⁶ and 2.4590 × 10⁻⁶ K, exceeding the component's `1e-6` temperature tolerance.

The recorded outlet enthalpy also rises despite negligible positive duty, indicating that the PH property round trip overwhelms the tiny intended enthalpy change. Tightening the outer cycle tolerance does not cure either failure. Investigate a consistent near-zero-duty state treatment or improved PH consistency, with energy-balance checks. A blanket relaxation of exchanger direction checks would conceal larger physical errors.

Relevant code: [heat-transfer direction validation](../../components/common/base_heat_exchanger.py).

## 5. Failure reporting loses the useful cause

`_OuterDesignObjective._record_from_hp` replaces the model's `solve_error` with a generic “infeasible, unconverged, or missing SEC” message. Consequently, the benchmark originally grouped native errors with physical failures. Preserve the original exception/message, solver stage, property fluid and input pair, retry history, and actual constraint violations in structured records.

The two BoTorch GP-fit exceptions recorded in the benchmark are separate events from these 44 thermodynamic candidate failures.

## Recommended implementation order

1. Add validated duty/flow checks before sink PH reconstruction, and treat invalid pressure probes as unavailable samples. Preserve detailed failure causes.
2. Extend the bounded tighter-convergence retry for small final pinch shortfalls, retaining the same final physical limit.
3. Repair the negligible-duty property/state inconsistency with focused energy and heat-direction tests.
4. Improve deterministic pressure/evaporator initialisation and validated flash recovery for the remaining PH failures. Do not train a physical-feasibility classifier on these numerical-failure labels indiscriminately.

The guard experiment above is confined to the diagnostic harness. **No production thermodynamic or optimiser source was changed.**

## Evidence and reproducibility

The investigation ran **155 fresh diagnostic processes**: 44 baseline replays; 88 pressure-seed checks; four tighter-tolerance checks; six detailed flash traces; four extra timeout-state traces; four guarded timeout cases; three guarded winner controls; and two detailed boundary cases. All processes either returned normally or reached their explicit diagnostic deadline. There were 17 deadline terminations across these overlapping experiments, not 17 additional distinct benchmark failures.

Baseline replay deadlines were 200 seconds for one attempt, compared with the original 180-second limit plus one retry. Seed and guard checks used 120 seconds; detailed flash traces 90 seconds; extra timeout-state traces 35 seconds; boundary traces 60 seconds. The shorter extra trace of the fourth timeout stops before its eventual air stall; the full 90-second trace captures it. Diagnostic timings include tracing/profiling overhead and are not optimiser speed comparisons.

All 77 frozen benchmark source files still match their archived hashes. The original 7,200 search records and winner replays remain untouched. Ruff passes for both diagnostic scripts. The three guarded winner checks are targeted controls, not a substitute for the full regression suite before implementing a production fix.

- [All 44 candidates and sensitivity outcomes](../../../results/RHPcomposition/botorch_v14_failures/candidate_table.md)
- [Complete numerical summary, states and profiles](../../../results/RHPcomposition/botorch_v14_failures/summary.json)
- [Diagnostic harness](../../../main_RHPfailure_investigation.py)
- [Report generator](../../../main_RHPfailure_report.py)
- [Evidence manifest](../../../results/RHPcomposition/botorch_v14_failures/manifest.json)

Raw `.trace.jsonl`, `.log`, `.request.json`, `.result.json`, and `.status.json` files accompany the summary. Run the harness from the repository root to reuse the completed baseline checkpoints; each sensitivity batch has its own saved requests file.
