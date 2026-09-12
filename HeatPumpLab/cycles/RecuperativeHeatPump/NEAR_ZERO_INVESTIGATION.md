# Near-zero RHX1 discontinuity investigation

The large near-zero COP peaks in the R290/R601 30/70 mass benchmark are caused
by premature cycle convergence after an evaporator-pressure adjustment. The
selected high-COP states violate the required gas-cooler approach. A diagnostic
guard that waits for compressor-1 suction pressure and enthalpy to settle removes
the peaks in the tested cases. This establishes a solver defect behind these
results; it is not evidence of a thermodynamic benefit from millikelvin RHX1 lifts.

Investigation date: 9 September 2026. The initial investigation used a diagnostic
guard only. The subsequent production correction is described below; the tables
in the investigation retain the original and diagnostic results for provenance.

## Production correction

The cycle now requires pressure, enthalpy, and mass flow to converge at **every
compressor inlet**, in addition to the existing heat/power closure checks. The
existing ten-sweep limit remains bounded. After pressure targeting, the solver
recomputes the final gas-cooler approach, including all terminal points, and
rejects nonfinite values or a shortfall beyond the 0.05 K tolerance. Fixed-pressure
offset-gradient probes retain their approximate objective; their selected design
still passes through the exact validated solve.

Pinch-property refresh also now multiplies the internal pinch fraction by the
exchanger duty before evaluating its stream properties. Previously the fraction
was passed to a helper expecting heat duty.

All **29 production replays without a diagnostic guard pass** nominal convergence
and the independent final-profile check. The captured near-zero case now has
COP 2.622839, with a 19.951298 K minimum approach. The earlier 17.5/82.5 mixture
and R600 VHC results retain their quoted values. Production replay inputs and
outputs are in `results/RHPcomposition/near_zero_fix/results.json`.

The numerical objective is versioned as `isolated-refprop-pattern-search-v12`.
Completed v11 signatures are no longer reused, and default sweep filenames use
v12 to keep new results separate. Both benchmark runners include the evaluation
version in checkpoint metadata and reject old checkpoints before evaluating.
Existing optimisation records are historical outputs; rerun a benchmark into a
new output directory to compare methods on the corrected objective.

Regression tests cover the late inlet pressure/enthalpy/flow change, final pinch
rejection and tolerance acceptance, pinch-location units, the near-zero lift
sequence, pressure seeds for the final candidate, and old-result invalidation.
To reproduce the fresh-process production replays:

```sh
.venv/bin/python results/RHPcomposition/near_zero_fix/verify_fix.py
```

## Reproduction and numerical evidence

The main case uses a 200 °C air target, 15 °C air inlet, 40 °C source inlet,
10 °C maximum refrigerant evaporator inlet, and a 20 K gas-cooler approach target.
RHX2 lift is 31.029212817770762 K and compressor-2 power fraction is
0.2947185847601135. The gas-cooler outlet offset is fixed at zero, as it was in
the archived high/low results. The direct cycle solves reproduce those results
without involving an outer optimisation method or offset approximation.

| Case | RHX1 lift, K | COP | Discharge pressure, MPa | Minimum GC approach, K |
|---|---:|---:|---:|---:|
| Original, zero | 0 | 2.622854 | 5.784563 | 19.954337 |
| Original, captured peak | 0.005809540466 | 2.742412 | 4.339044 | **2.618325** |
| Suction-state guard, same design | 0.005809540466 | 2.622839 | 5.784197 | 19.951298 |
| Original, final CMA-ES winner | 0.016142113040 | 2.769876 | 5.120857 | **1.528072** |
| Suction-state guard, same final design | 0.016142113040 | 2.664071 | 6.906473 | 19.969185 |

The final CMA-ES design has RHX2 lift 16.268979769802385 K and compressor-2 power
fraction 0.16987922757884472. The two peak states fail by approximately 17.38 K
and 18.47 K. The model's permitted pressure-target residual is 0.05 K, so these
are substantial violations, not borderline tolerance decisions.

Each final gas-cooler profile was checked independently using all terminal
points and 101 evenly spaced duty fractions, as well as the existing bounded
internal-minimum calculation. The initial captured-case check used 501 fractions.
The decisive minimum is the hot terminal, so its detection does not depend on
the internal sampling resolution or local minimiser.

![COP, pressure, and final pinch across the near-zero lift sweep](../../../results/RHPcomposition/near_zero_investigation/near_zero_diagnostic.png)

## Exact cause

1. Evaporator alignment raises its outlet from 10 °C to the superheater-limited
   34.9999 °C. Suction pressure must rise from 0.062993107 to 0.161897977 MPa.
2. In `rhp_helpers.py::_get_RecuperativeHP_solver_pathway`, the compressor and
   intermediate-pressure relationship execute before the hot-side RHX1 solve.
   `set_rhx1_dT_cold` updates temperature at the existing outlet pressure.
   RHX1's later `get_fixed_cold_delta_T_model` propagates the new pressure to
   its cold outlet, which is compressor-1 suction.
3. In the first sweep after alignment, compressor 1 consequently still sees
   0.062993107 MPa. RHX1 later changes that inlet to 0.161897977 MPa, after the
   compressor has run. Its inlet enthalpy also changes from 445.398362 to
   441.292907 kJ/kg.
4. `_cycle_analysis.py::run_sequential_cycle_analysis` checks expansion-valve
   inlet enthalpy, high-side pressure drop, total compressor power, and total
   gas-cooler duty. It does not check compressor suction-state changes. In the
   captured transition, power remains 282.283350 kW and heat changes by only
   −0.000203467 kW, less than its 0.000624437 kW tolerance. The sweep is marked
   converged even though suction pressure changed by a factor of 2.57009.
5. Pressure targeting receives a stale residual of +0.00406883 K at
   4.339043760 MPa. Subsequent pressure trials propagate the changed state and
   return very negative residuals. The scalar search can retain the initially
   accepted upper endpoint and its cached residual. It returns that pressure.
6. The final cycle reruns settle the compressor state, giving the inflated COP
   at an insufficient discharge pressure. `solve_at_discharge_pressure` updates
   pinch properties but returns `cycle_convergence_failure`; that validation
   checks operation convergence and positive work, not the final GC pinch target.
   The result is therefore still reported as converged.

The small positive lift adds RHX1 to the sequential dependency path. Whether
another sweep happens is governed by numerical convergence residuals. This
explains why the false branch appears over a narrow positive-lift interval
rather than at every positive lift or only at exactly zero.

The detailed operation trace is in
`results/RHPcomposition/near_zero_investigation/small_ops.json`; the first sweep
after evaporator adjustment is identifiable by `calls == 9` and
`Tevap_out > 34`. `small.json` and `zero.json` contain pressure-trial traces and
501-point gas-cooler profiles.

## Controlled diagnostic and limits

The diagnostic repeats the existing sequential cycle routine until compressor-1
suction pressure and enthalpy remain within the existing relative cycle
tolerance, with a finite ten-call limit. It is installed only inside the fresh
diagnostic worker. No property correlations, physical design inputs, pinch
targets, or optimiser settings are changed.

The experiment comprises **58 fresh-process solves**, 29 original and 29 guarded:

- Seventeen lifts from 0 to 0.1 K, including the captured design.
- Both captured and final designs with pressure guesses automatic, 3, 4.3,
  5.5, and 7 MPa.
- The earlier 17.5/82.5 propane/pentane and R600 VHC comparison designs.

All 58 return nominal convergence. **Eleven of the 29 original results fail
the independent pinch check; all 29 guarded results pass.** In the lift sweep,
seven original peaks between the sampled lifts 0.003 and 0.01 K disappear.
These sampled bounds are not a precise location of the branch transitions.

The final design's COP range across pressure guesses shrinks from
2.664030–2.769876 to 2.664030–2.664175. Some residual seed dependence remains:
the earlier design gives guarded COP 2.621483–2.622839, and small steps remain
in the lift curve. The present pressure tolerance and endpoint-selection policy
have not been redesigned. This investigation therefore does not claim exact
initialisation independence or certify a global optimum.

The previous VHC comparison remains valid at the quoted precision: guarded VHC
is 1.83609909 MJ/m³ for the 17.5/82.5 mixture and 2.02076583 MJ/m³ for R600,
still 9.14% lower for the mixture.

## Recommended production correction

Include the changing compressor inlet states in cycle convergence, ideally
checking pressure, enthalpy, and flow across the coupled stream states rather
than relying only on aggregate heat and power. Ensure the pressure residual
is evaluated on a settled cycle and recheck the selected pressure's final
pinch after all evaporator adjustments. An accepted endpoint with a stale
residual must not be sufficient for success. Pinch-stream property refresh
alone does not constitute that validation.

Regression coverage should include this delayed-pressure transition, the
captured and final designs across pressure guesses, and the zero/near-zero
lift sequence. Rerun optimiser comparisons after correcting the production
objective; the archived high-COP CMA-ES winners above are infeasible under the
specified 20 K pinch constraint. The diagnostic guard is a causal experiment,
not a fully qualified replacement cycle solver.

## Reproduce

From the repository root:

```sh
MPLCONFIGDIR=/private/tmp/hpl-mpl MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python paper/analysis/investigate_near_zero.py
```

The script creates fresh workers with a 180-second timeout and at most two
concurrent cases. Results, inputs, profiles, and plots are written to
`results/RHPcomposition/near_zero_investigation/`. The full run is in
`results.json`, with individual `case_*.json` files and `run.log`. To repeat the
operation trace separately:

```sh
MPLCONFIGDIR=/private/tmp/hpl-mpl MPLBACKEND=Agg .venv/bin/python results/RHPcomposition/near_zero_investigation/trace_probe.py --lift 0.005809540466003726 --ops --out /private/tmp/near_zero_trace.json
```
