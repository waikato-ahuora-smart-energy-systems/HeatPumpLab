# RHP optimisation convergence investigation

Date: 2026-07-30

## Cases investigated

The current results file contains 1,702 rows and 1,702 unique case
signatures. It therefore has no duplicate signatures. It has 1,646 unique
physical `(input_refrigerant, T_air_tar)` cases because 56 cases have both a
v4 and v5 result. Those are deliberate cross-version reruns, not duplicate
work within one optimiser version.

Four of the 1,650 default physical cases were absent:

| Refrigerant | Target |
|---|---:|
| Pentane | 105 °C |
| Butane/CO2, 50/50 mass | 170 °C |
| Propane/Butane/CO2, equal mass | 170 °C |
| Isobutane/Butane/CO2, equal mass | 170 °C |

The exact v5 case signatures are `77be17609e714345`,
`3704a3fdf3cb2e3b`, `303be9c1902415bd`, and `432747d90a18ff7a`,
respectively.

## Root causes

There is no unbounded Python loop in these four cases. Stack sampling during
the apparent hang consistently located the active worker in REFPROP
`ABFLSHdll`, reached through a P-h flash. The excessive runtime came from two
separate bounded-search defects.

### Pentane: infeasible boundary prevents offset optimisation

Pentane at these targets is infeasible at zero gas-cooler outlet-temperature
offset but feasible at a small positive offset. The failed boundary solve
returns an infinite SEC. Finite differencing then forms a non-finite gradient
and exits without entering the feasible interval.

The retained geometric bracket and eight bisections find the first feasible
offset before estimating the boundary gradient. Non-finite local objectives
are also converted to a finite penalty.

### Mixtures: seed-dependent pressure branch and nested flash amplification

A discharge pressure copied from the preceding outer-design point is not a
safe continuation variable for mixtures. The gas-cooler pinch residual can
have more than one accepted pressure root. The existing local bracket stops
at the root around its seed, so identical thermodynamic designs can return
different pressures and SEC values:

| Initialisation for the same Butane/CO2 split point | Time | SEC | Pressure |
|---|---:|---:|---:|
| Cold | 7.18 s | 0.759701 | 7.927 MPa |
| Previous corner pressure | 4.02 s | 0.830261 | 9.593 MPa |
| Previous split-line pressure | over 90 s in one trial | — | — |

The higher root is both thermodynamically worse and capable of sending the
nested property calculations into a very expensive region.

The second multiplier is the internal gas-cooler pinch minimisation. A
profile of one cold Butane/CO2 edge took 71.7 seconds. Its 56 pressure
residual evaluations generated 8,003 uncached REFPROP calls; 71.4 seconds
were inside REFPROP. The internal pinch location used the global `1e-6`
cycle tolerance even though it only supports the explicitly checked 0.05 K
pinch-acceptance band.

## Prototypes trialled

Times are local wall-clock measurements and are useful as relative, not
portable, benchmarks.

### Pentane offset recovery

| Prototype | Result | Decision |
|---|---|---|
| Ten-point fixed offset grid | Fast but only resolved the offset to the grid spacing and always spent ten solves | Rejected |
| Geometric bracket plus eight bisections | Pentane 100–165 °C all solved; deterministic finite work | Kept |
| Bounded offset minimisation with finite penalty | About 19 calls for no material SEC improvement | Rejected |

### Mixture pressure and pinch search

| Prototype | Result | Decision |
|---|---|---|
| Reuse pressure and evaporator temperature | Fast on a selected edge, but seed-dependent and sometimes chose the wrong pressure root | Rejected |
| Global 8-point pressure scan to 40 MPa | Exceeded 120 s on the pathological point because the scan visited expensive flash regions | Rejected |
| Disable mixture pressure seed only | Restored the lower root, but the 21-point trajectory exceeded 180 s before completion | Necessary but insufficient |
| Internal pinch duty tolerance 1% | Numerically unstable trajectory; one isolated run exceeded 90 s | Rejected |
| Internal pinch duty tolerance 0.1% | Correct result, 40.5 s on the profiled edge | Rejected in favour of faster matched tolerance |
| Internal pinch duty tolerance 0.3% | Same SEC and pressure as the strict solve, 33.2 s versus 71.7 s | Kept |
| No mixture pressure seed plus 0.3% pinch tolerance | All 21 formerly problematic 50/50 search points completed in 234.5 s; zero failures | Kept |
| Per-evaluation subprocess watchdog | Would bound a single REFPROP call but adds process startup/state-transfer cost to every exact evaluation | Not needed after complete-case validation |

The kept pinch tolerance is paired with a 40-iteration scalar cap. Pressure
bound expansion is capped at 100 steps instead of 1,000. Existing pressure
root, evaporator alignment, offset, line-search, polish, and total exact
evaluation budgets remain finite.

## Retained implementation

- Find and refine a feasible gas-cooler outlet offset before differentiating
  when the zero-offset solve fails.
- Convert non-finite local objective values to a finite penalty and cap the
  local offset optimiser.
- Evaluate every outer-design vector from its own case data. Do not carry
  discharge pressure or evaporator-outlet temperature from another vector,
  because either continuation variable can select a different successful
  inner-solver branch and make the objective evaluation-order-dependent.
- Use `1e-3` relative pressure-location tolerance for mixtures, still subject
  to the existing 0.05 K residual acceptance check.
- Resolve internal gas-cooler pinch location to 0.3% of duty with a
  40-iteration cap.
- Cap pressure-bound expansion at 100 iterations.

## Full end-to-end validation

All four cases were run through the public outer optimiser with the default
budgets. Mixture composition starts and polishes used five workers.

| Case | Wall time | Exact evaluations | Failed | Best SEC |
|---|---:|---:|---:|---:|
| Pentane, 105 °C | 54.3 s | 56 | 0 | 0.466352 |
| Butane/CO2, 170 °C | 345.3 s | 137 | 0 | 0.523900 |
| Propane/Butane/CO2, 170 °C | 424.6 s | 181 | 0 | 0.525044 |
| Isobutane/Butane/CO2, 170 °C | 370.1 s | 169 | 0 | 0.546811 |

The slowest individual composition starts were 235, 244, and 240 seconds for
the three mixture cases. These remain expensive high-temperature CO2 mixture
calculations, but they now complete with no manually terminated workers, no
silently discarded composition starts, and no failed exact evaluations.

Focused solver, pressure, profiling, and outer-optimiser tests pass:
153 passed. The full suite has 392 passing tests and one unrelated,
pre-existing failure:
`test_main_rhpopt_run_many_delegates_to_public_batch_optimiser` expects
unreversed default lists even though `run_many(reverse=True)` deliberately
reverses them.

## v4 versus v6 random sample

A reproducible five-case sample was drawn from the v4 result rows with random
seed `20260730`. Each v6 case used the default optimisation budget and one
inner worker, matching the way a case is executed inside `run_many`.

| Case | v4 SEC | v6 SEC | Relative change | Runtime speed-up |
|---|---:|---:|---:|---:|
| Propane, 205 °C | 0.702182 | 0.711817 | +1.372% | 1.73x |
| Butane/Neopentane/CO2, 170 °C | 0.551651 | 0.543385 | -1.498% | 12.36x |
| Isobutane/CO2, 210 °C | 0.676096 | 0.699224 | +3.421% | 8.65x |
| Neopentane/Pentane, 210 °C | 0.680042 | 0.680521 | +0.071% | 7.48x |
| Butane/Isopentane/Pentane, 205 °C | 0.650236 | 0.650073 | -0.025% | 5.39x |

The versions therefore do not return the same final optimum for this sample.
Two of five agree within 0.1% SEC; v6 finds a 1.50% better point in one case
and a 1.37% and 3.42% worse point in two cases.

Cold v6 replay of each stored v4 optimum changes its SEC by no more than
0.012%. The thermodynamic objective at those designs is effectively
unchanged. The larger differences arise because the edge-first outer search
follows a different path and selects different design vectors. In particular,
the stored v4 design remains superior when replayed for the Propane and
Isobutane/CO2 cases. The v6 runtime is substantially lower and all five cases
finish with zero failed exact evaluations; the corresponding v4 rows contain
16 failed evaluations in total.

## Why v6 missed the Propane and Isobutane/CO2 v4 basins

Controlled factorial runs isolate the regression to state carried between
otherwise independent outer-design evaluations. The outer objective is
therefore history-dependent: evaluating a design vector after a neighbouring
point can select a different successful inner-solver branch than evaluating
the same vector cold. The existing fallback does not protect against this
because it retries without the seed only when the seeded solve fails; the
wrong branches in these two cases report success.

### Propane, 205 °C

| Variant | Time | Exact evaluations | Best SEC |
|---|---:|---:|---:|
| v6 pressure and temperature continuation | 21.4 s | 56 | 0.711817 |
| Temperature continuation only | 23.3 s | 58 | 0.702182 |
| No state continuation | 23.3 s | 58 | 0.702182 |
| Strict v4 pinch tolerance, with both seeds | 29.3 s | 56 | 0.711817 |
| Strict v4 pinch tolerance, no seeds | 32.5 s | 58 | 0.702182 |

Discharge-pressure continuation is the sole causal factor for this case.
Removing it exactly recovers the stored v4 design and SEC; changing the
internal pinch tolerance does not.

At the same edge design, `(dT1, dT2, split) = (66.667, 0, 0.5)`, a cold
evaluation returns SEC 0.710132 and 6.314 MPa. Following the actual seeded
edge trajectory carries a pressure of 12.472 MPa into that point and returns
SEC 0.912119 at 11.396 MPa. Both pressure branches are accepted as successful.
The falsely poor edge value changes the split source and subsequent polish,
so the v4 basin is never explored.

### Isobutane/CO2, 210 °C

Pressure continuation is already disabled for mixtures, but evaporator-outlet
temperature continuation has the same class of effect:

| Variant | Time | Exact evaluations | Best SEC |
|---|---:|---:|---:|
| v6 temperature continuation | 412.1 s | 139 | 0.699224 |
| No state continuation | 701.4 s | 138 | 0.676126 |
| Stored v4 result | 3660.7 s | 138 | 0.676096 |

The no-continuation run uses all current v6 pressure and pinch tolerances yet
recovers the v4 basin to 0.0044% SEC while remaining 5.2 times faster than the
stored v4 run.

The decisive divergence occurs in the 75/25 Isobutane/CO2 edge search. At
`dT2 = 23.607 K`, temperature continuation returns SEC 0.808255 and
18.244 MPa, while a cold evaluation returns SEC 0.790421 and 17.290 MPa.
That difference changes the bounded scalar search:

- With continuation, it turns back towards 33 K and selects an edge point at
  33.333 K.
- Without continuation, it probes 14.590 K and selects that lower edge point.

The edge choice changes the split-line length. Only the cold trajectory
constructs the split point `(47.214, 29.180, 0.5, 0.75, 0.25)`, which enters
the three-source polish and reaches the stored v4 basin at approximately
`(47.214, 29.180, 0.6, 0.85, 0.15)`.

### Conclusion and proposed correction

The relaxed internal pinch-location tolerance and mixture pressure tolerance
are not responsible for these two misses. No-seed searches using those current
tolerances recover the v4 basins, while strict pinch searches with the
problematic seed retain the wrong Propane basin.

For a deterministic and robust objective, do not carry either
`P_H_est` or `Tref_evap_outlet` between outer-design vectors. A successful
warm solve cannot be trusted solely on its success flag because the inner
problem has multiple accepted local branches. Cold verification of every
warm result would also restore determinism, but would retain the complexity
while losing most of the speed benefit.

The measured cost of removing both continuation variables is modest for
Propane (21.4 to 23.3 seconds) and larger for Isobutane/CO2 (412 to
701 seconds). The latter is still substantially faster than v4 and buys the
correct basin. This correction was subsequently implemented in optimiser v7.

## Additional robustness investigation and optimiser v8

The v7 outer line search evaluated four equally spaced points with the default
eight-point budget, then started a separate whole-interval bounded scalar
minimisation. The scalar minimiser did not use the values of the initial
samples and assumed that the full interval was unimodal. Consequently, half
of the budget provided coverage information that had no influence on the
adaptive probes.

### Equal-budget scalar-search prototypes

Four policies were compared on 4,000 smooth quadratic surfaces and 4,000
random multi-basin surfaces with occasional infeasible intervals. Each policy
received eight evaluations.

| Policy | Multi-basin median regret | 90th-percentile regret | Within 0.05 of reference | Decision |
|---|---:|---:|---:|---|
| v7 independent global bounded search | 0.01551 | 0.26249 | 66.9% | Replaced |
| Eight-point uniform grid | 0.02681 | 0.22982 | 63.0% | Rejected |
| Adaptive interval subdivision | 0.00873 | 0.15335 | 78.5% | Rejected |
| Five-point screen plus best-bracket refinement | 0.00670 | 0.13537 | 81.3% | Kept |

On the ideal quadratic set, the kept policy remained within 0.01 of the
reference in 100% of trials. Its median regret was 0.00026 versus zero for
the v7 policy, an acceptable local-precision trade-off for substantially
better multi-basin coverage.

The kept policy evaluates five equally spaced points, selects the best
successful observation, and spends the remaining three calls only inside its
neighbouring sample bracket. Failed samples remain finite penalties and cannot
attract refinement over a successful point.

### Real-case validation

| Case | v7 SEC | Bracketed SEC | Relative change | v7 time | Bracketed time | Failed |
|---|---:|---:|---:|---:|---:|---:|
| Propane, 205 °C | 0.702182 | 0.698440 | -0.533% | 22.6 s | 22.8 s | 0 |
| Pentane, 105 °C | 0.466352 | 0.466035 | -0.068% | 84.8 s | 64.1 s | 0 |
| Isobutane/CO2, 210 °C | 0.676126 | 0.676044 | -0.012% | 701.4 s | 677.7 s | 0 |

The difficult Isobutane/CO2 result is also slightly better than the stored v4
SEC of 0.676096. The bracketed policy improved all three cases without
increasing the exact-evaluation budget and was retained as optimiser v8.

### Joint recuperator-lift polish

Allowing the local pattern search to reopen both recuperator lifts improved
the bracketed Propane result from SEC 0.698440 to 0.695263 using 23 additional
evaluations and 8.7 seconds. It was not safe as an unconditional default,
however. The corresponding Isobutane/CO2 trial was terminated after more than
nine minutes before its 24-evaluation allowance completed. The interrupted
stack was inside REFPROP `ABFLSHdll`, reached through the internal gas-cooler
pinch calculation.

This prototype is deferred as a possible opt-in verification stage after a
hard per-evaluation wall-time limit exists. An evaluation-count limit alone
does not control its runtime.

### Composition coverage

A ternary case currently produces seven distinct starts: the centroid, three
pure vertices, and three binary-edge midpoints. Therefore
`max_composition_starts=8` does not by itself add an asymmetric interior
start. On a dense simplex grid, the worst nearest-start distance is 0.4042.
Adding only one asymmetric 50/25/25 point leaves that worst distance
unchanged because two symmetric holes remain. Adding all three permutations
reduces it to 0.3536 but raises the full-search count from seven to ten.

A single biased eighth point was rejected. Symmetric interior coverage should
instead be offered as a higher-budget robust mode, or via a two-stage
composition screen that promotes only promising starts.

### Remaining operational safeguards

Iteration caps now bound Python and SciPy loops, but cannot interrupt one slow
or hung native property call. `ProcessPoolExecutor.future.result()` is used
without a timeout at both the batch and composition levels, and normal
executor shutdown waits for running workers. A reliable watchdog therefore
requires disposable worker processes that the parent can terminate, not a
thread timeout around REFPROP.

The recommended next operational increment is:

1. Persist each completed composition search and polish stage so a killed
   case can resume without repeating completed work.
2. Run a whole case or composition search in a disposable process with a
   configurable wall-time limit.
3. On timeout, terminate that process, record a structured timeout result,
   and continue the remaining batch.
4. Record the slow design vector and phase for later isolated replay.

Pentane recovery originally emitted many expected inner trial-failure messages
while the final outer evaluation succeeded. Non-recording optimisation solves
now set a separate `report_failed_solve=False` flag, suppressing those expected
messages without changing failed-case persistence. Unexpected exceptions and
normal standalone failures remain visible. A live Pentane 105 °C replay was
silent and reproduced SEC 0.466035 with zero failed outer evaluations.

## Second random v4 versus v8 sample

A second reproducible sample of five successful v4 rows was drawn with random
seed `20260731`. The five cases used in the earlier random sample were excluded,
leaving 1,628 eligible rows. Each stored v4 optimum was first evaluated cold
with the current solver, then the complete v8 optimisation was run with default
budgets and one inner worker. The five cases ran concurrently.

SEC changes below compare v8 with the cold current-solver replay of the stored
v4 design, isolating search effectiveness from solver-model changes.

| Case | Replayed v4 SEC | v8 SEC | SEC change | Runtime speed-up | Failed evaluations v4 → v8 |
|---|---:|---:|---:|---:|---:|
| Isobutane/Neopentane, 105 °C | 0.443546 | 0.441941 | -0.362% | 3.88x | 0 → 0 |
| Pentane/CO2, 130 °C | 0.546652 | 0.448501 | -17.955% | 4.31x | 25 → 0 |
| Water/Ammonia, 220 °C | 0.751928 | 0.751928 | 0.000% | 1.38x | 21 → 9 |
| Propane/Isobutane/Butane, 150 °C | 0.507937 | 0.507606 | -0.065% | 2.58x | 1 → 0 |
| Isobutane/Pentane/CO2, 175 °C | 0.536158 | 0.536401 | +0.045% | 3.78x | 14 → 0 |

The sample contains three improved results, one identical result, and one
0.045% regression within the same design basin. Median SEC change is -0.065%
and median runtime speed-up is 3.78x. Total exact evaluations are nearly
unchanged (705 for v4 versus 703 for v8), while failed evaluations fall from
61 to 9, an 85% reduction. Summed CPU time falls from 10,992 to 5,099 seconds,
an aggregate 2.16x speed-up.

Pentane/CO2 is the strongest basin correction. The v4 design is pure CO2 at
16.3 MPa; v8 selects 90/10 Pentane/CO2 at 1.42 MPa and reduces SEC by 18%.
Isobutane/Pentane/CO2 is the only regression: v8 changes the second recuperator
lift from 51.3 to 55.9 K and compressor split from 0.45 to 0.50 while retaining
the same 45/50/5 composition, so this is a small line/polish precision
difference rather than a missed remote basin.

Water/Ammonia remains the operational outlier. It returns the identical pure
Ammonia design but takes 3,049 seconds and still has nine failed outer
evaluations. Its invalid trial points raise an LMTD temperature-cross
`ValueError`, currently classified as an unexpected exception, so their
messages bypass expected-failure console suppression. This should be
reclassified as a physical `CycleFailure`; the case also reinforces the need
for process-level timeouts and checkpoints.

## No-output v8 batch investigation and optimiser v9

A new default `run_many()` attempt produced many apparent failure messages but
no result row for more than an hour. The existing
`RHPopt_ncomponent_results_v6.csv` was not evidence from that run: all 1,702
rows in it are v4 or v5. No v8 case had completed and reached its temp row.

### Why the batch appeared stuck

Three effects aligned:

1. `run_many(reverse=True)` reversed both dimensions of the case grid. The
   first eight workers therefore received Water/Ammonia targets from 220 down
   to 185 °C. Historical v4/v5 runtimes for these cases are 2,541–7,042
   seconds per case. The sampled v8 220 °C case took 3,049 seconds.
2. Completed worker rows were written to per-case temp files, but the parent
   merged those files into the master CSV only after every submitted future
   returned. The visible master file could therefore remain unchanged for the
   complete 1,702-case run.
3. A Water/Ammonia LMTD temperature cross was caught by the generic exception
   handler. It was labelled unexpected, so `report_failed_solve=False` did not
   suppress it. The repeated messages looked like whole-case failures even
   though they were infeasible internal trial points.

There is still no unbounded Python loop in this path. Both interrupted
Water/Ammonia trials were inside REFPROP `ABFLSHdll`; one was reached through
the internal gas-cooler pinch search.

### Offset-recovery amplification

The Pentane recovery was being applied to every failed zero-offset seed,
regardless of its physical failure class. A live offset sweep separated the
two modes:

| Candidate | Zero-offset result | First feasible offset | Feasible SEC |
|---|---|---:|---:|
| Pentane, 105 °C | cycle iteration failure in 0.05 s | about 0.8 K | 0.4800 |
| Water/Ammonia 50/50, 220 °C | LMTD cross in 3.75 s | 40 K | 1.0870 |

For Water/Ammonia, offsets 0, 0.1, 0.2, 0.4, 0.8, 1.6, and 3.2 K all produced
the same temperature-cross class at roughly 3.7–4.1 seconds each. The
geometric bracket plus eight refinements could therefore turn one failed outer
point into roughly 18 full solves, only to recover a very poor extreme-bound
design.

Three alternatives were considered:

| Prototype | Decision |
|---|---|
| Suppress all unexpected exceptions during optimisation | Rejected because it would hide property and programming faults |
| Keep recovering every physical failure up to 40 K | Rejected because temperature-cross designs caused large, low-value retry amplification |
| Recover only the observed `cycle_iteration` boundary-failure class | Kept; it preserves Pentane recovery and returns structural temperature crosses immediately |

A dedicated `TemperatureCrossError(ValueError)` now carries the approach
violation. The cycle solver converts it to an expected `CycleFailure` with a
`heat_exchanger_temperature_cross` diagnostic. Non-recording optimiser trials
are consequently quiet, while genuine unexpected exceptions remain visible.

Live post-change checks:

- Pentane 105 °C recovered at 0.7953 K in 0.70 seconds.
- Water/Ammonia 50/50 at 220 °C returned once in 3.76 seconds with a
  132.07 K temperature-cross diagnostic and no console output.

### Reporting-only solve discovered in every candidate

Stack sampling of a bounded Water/Ammonia search found another multiplier:
after every successful candidate, `analyse_RHP_case()` called `_get_dQdW()`.
That reporting metric launches a second full heat-pump solve at a perturbed
air target. Outer optimisation does not use `dQdW`, and its calls explicitly
set `record_result=False`.

The follow-up solve is now skipped for non-recording optimiser candidates and
retained for normal standalone/reporting analyses. A repeat 25-evaluation
Water/Ammonia trial advanced from the second edge search into the split search
within the same three-minute observation window. It still did not finish:
successful Water/Ammonia REFPROP mixture calculations remain intrinsically
expensive, so a process-level watchdog is still warranted.

### Batch progress and persistence

The retained batch changes are:

- `run_many()` no longer reverses the case grid by default. It starts with
  ordinary low-temperature Propane cases instead of the slowest
  Water/Ammonia block. Explicit `reverse=True` remains available.
- Parallel futures are consumed in completion order. Returned rows are kept
  in original case order for API compatibility.
- The parent merges only the completed worker's temp file into the master CSV
  immediately, avoiding races with other workers that may still be writing.
- One concise completion line reports progress and success status.

A one-evaluation live batch confirmed that the master CSV existed with one row
as soon as the case returned and that its temp file was removed. These
behavioural changes are recorded as optimiser v9 so they cannot be confused
with prior v8 signatures.

Focused batch, outer-optimiser, analysis, component, and solver tests pass:
125 passed. The complete project suite passes: 396 passed.

## Multibasin screened optimiser v10

Incumbent replay was prototyped first, but it only prevented a regression; it
did not make the search find a better basin. It is retained as a final
guardrail, not as the v10 search improvement.

### Scale and root cause

Of the first 746 comparable successful v9 rows, 227 had a stored SEC more than
0.01% above v4, 43 were more than 0.5% above v4, and 14 were more than 1%
above v4. Some differences include model drift, but live cold checks confirmed
the most important missed basins.

The systematic search weakness was:

1. Recuperator lines were searched only at compressor split 0.5 and a small
   set of fixed compositions.
2. Each eight-evaluation line screened five points and refined only the
   lowest observed basin.
3. The final pattern search could vary compressor split and composition but
   froze both recuperator lifts.
4. The nominal 250-evaluation budget was often badly underused: a pure-fluid
   case used roughly 40–60 evaluations and a binary case roughly 120–160.

The Propane/Isopentane 180 °C case illustrates the interaction. At 25/75 mass
and compressor split 0.5, the five-point recuperator screen returned SEC
0.549419 at zero lift and 0.549790 at 50 K. Refining only the nominally better
zero-lift point discarded the second basin. Refining the separate 50 K basin
with the same thermodynamic model sampled 55.90 K and SEC 0.541127. A later
compressor/composition screen then improved it further.

### Prototypes

| Prototype | Finding | Decision |
|---|---|---|
| Historical incumbent replay only | Prevents regressions but does not improve independent basin discovery | Retained only as a guardrail |
| Spend all three adaptive line calls on a near-best interior basin | Recovered Propane/Isopentane but lost Neopentane's shallow endpoint basin | Rejected |
| Dense joint grid over lift, compressor split, and composition | Revealed the basin interaction but 32 cold points took about 40 s for one slice | Rejected as the default |
| Budget-aware two-basin line refinement | Preserves endpoint and interior basins using otherwise idle budget | Retained |
| Full-range compressor-split screen before local polling | Finds remote 0.1 or 0.9 split basins that 0.05 local steps cannot reach | Retained |
| Select polish sources only after collapsing each direction globally | Discarded promising composition-specific basins in Butane/CO2 | Rejected |
| Rank polish sources across every composition and direction | Preserves strong alternative composition/lift combinations | Retained |
| Joint recuperator polish for every mixture | Previously exceeded nine minutes for one Isobutane/CO2 trial | Rejected pending a hard wall-time watchdog |
| Joint recuperator polish for pure fluids only | Recovers shallow coupled lift/split basins cheaply | Retained |
| Feed polished mixture coordinates back through a seven-call lift line | Recovers coupled lift/composition basins using only unused case budget | Retained |

### Retained search

- Composition coverage is selected first.
- The remaining case budget is redistributed across those starts. With the
  default 250 evaluations this expands a binary start from 25 to 35
  evaluations and a pure start to 37, while seven-start ternary coverage
  remains at 25 evaluations per start.
- Expanded lines screen seven values and refine the two strongest observed
  local basins rather than forcing one unimodal bracket.
- Each selected polish source screens compressor split at 0.1, 0.3, 0.5, 0.7,
  and 0.9 before local 0.05/0.025 refinement.
- Polish sources are ranked across all composition starts and line directions;
  a good secondary composition is no longer discarded by the global
  per-direction reduction.
- The polish allowance increases from 12 to 24 exact evaluations.
- Pure-fluid polishing can also reopen both recuperator lifts with a 5 K
  pattern step. Mixtures retain fixed-lift polishing to bound REFPROP runtime.
- Successful mixture polishes on a nonzero lift are processed worst-first
  through a seven-call one-dimensional lift feedback search, using only the
  evaluations still available inside the 250-evaluation case budget.
- Partial `base_case` inputs now merge with the standard defaults, making
  direct optimisation consistent with batch case generation.

Historical incumbent discovery and current-solver replay remain enabled after
the new search. They are outside the 250-evaluation search allowance and are
reported separately. V10 writes to
`RHPopt_ncomponent_results_v10.csv`.

### No-incumbent live benchmarks

All results below were produced with historical incumbent preservation
disabled. V10 therefore had to discover the result through its own search.

| Case | v4 SEC | v9 SEC | v10 SEC | Change vs best prior |
|---|---:|---:|---:|---:|
| CO2, 135 °C | 0.558035346 | 0.564852778 | 0.557546108 | -0.088% |
| Propane, 185 °C | 0.630457683 | 0.635234536 | 0.630457668 | approximately 0% |
| Neopentane, 180 °C | 0.622825590 | 0.589514920 | 0.589127563 | -0.066% |
| Propane/Isopentane, 180 °C | 0.529159316 | 0.542104263 | 0.524798759 | -0.824% |
| Propane/Pentane, 195 °C | 0.566247892 | 0.577140420 | 0.557678017 | -1.514% |
| Isobutane/Pentane, 195 °C | 0.604984647 | 0.612774014 | 0.593034065 | -1.975% |
| Butane/CO2, 210 °C | 0.656536845 | 0.670820154 | 0.650076687 | -0.984% |

This set deliberately includes cases where v9 missed v4 and a case where v9
was already 5.35% better than v4. The new search matched or improved the best
prior result in all seven without using historical designs. The Butane/CO2
case was validated in stages to avoid repeating its complete nine-minute
screen: the retained composition-diverse polish produced
`(dT1, split, composition) = (68.675 K, 0.896875, 84.6875/15.3125)` and the
retained seven-call feedback line independently reached 84.262 K and SEC
0.650076687.
