# Fixed-composition optimisation investigation — 8 September 2026

The subsequent [200-evaluation origin-start study](ORIGIN_OPTIMISER_COMPARISON.md)
compares BO, RBF and CMA-ES across repeated seeds. Its pressure-initialisation
diagnostics show that identical physical designs can return different COPs.
Read the nominal improvements below with that inner-solver limitation in mind.

The v11 search misses better designs at fixed composition. This is supported
by exact thermodynamic evaluations, not just the irregular composition curve.
Joint RBF and Gaussian-process Bayesian optimisation are credible alternatives;
the experiment below is a pilot, not a global-optimality certification.

## What is going wrong

The initial searches in `optimisation/_multistart.py` examine recuperator edges
and a fixed-total-lift split line at a compressor power fraction of 0.5.
In `optimisation/_polish.py`, `_run_polish` enables joint recuperator polishing
only when the refrigerant has one component. Mixtures keep both lifts fixed
while compressor split is polished. A seven-call lift feedback line provides
limited correction afterwards. These stages can miss coupled lift/split basins.

The restriction was originally a runtime safeguard: the previous robustness
investigation recorded a slow joint mixture trial. V11 now has subprocess
timeouts and worker recycling, allowing a bounded investigation of joint search.
The 250-call setting is a ceiling, not a convergence guarantee: 32 of the 41
fixed-composition cases used 89 evaluations. Staged allowances and local stopping
leave much of the nominal budget unused.

There is a direct counterexample for R290/R600A at 30/70 by mass. The saved
design has `(RHX1, RHX2, split) = (79.398867 K, 0 K, 0.9)` and COP 2.511081.
Changing only RHX2 to 5 K gives COP 2.534672, a 0.940% increase. Reopening both
lifts with the existing pattern search eventually gives COP 2.538510 after
40 additional evaluations. This improvement reproduced in a fresh process.

Local repair alone is insufficient. The same warm joint polish did not improve
the saved R290/R601 30/70 or 60/40 designs, whereas both cold surrogate searches
did. For the 30/70 case, the old and GP-EI designs are:

| Quantity | v11 | GP-EI |
|---|---:|---:|
| RHX1 cold lift, K | 0 | 0 |
| RHX2 cold lift, K | 62.732200 | 18.103164 |
| Compressor 2 power fraction | 0.493750 | 0.163545 |
| COP | 2.646844 | 2.671504 |
| VHC, MJ/m³ | 2.205391 | 2.179909 |
| Maximum pressure, MPa | 3.602665 | 6.709670 |

An 11-point exact evaluation along the straight line between these designs
falls to COP 2.623279 before rising into the better region. This tested path
contains a performance valley; it is consistent with the observed failure of
local polishing to escape. It does not prove that every possible connecting
path has a valley. Disabling the fast inner objective at both endpoints gives
exactly the same COP values, so this particular improvement is not caused by
that approximation. The pressure increase and small VHC reduction are material
trade-offs despite the higher COP.

## Pilot design

Three fixed mixtures at a 200°C sink target were selected: R290/R601 30/70,
R290/R601 60/40 and R290/R600A 30/70, all by mass. Source, sink, compressor and
pinch settings were loaded directly from their original sweep rows.

The cold-start methods each receive 64 evaluations, including the identical
24-point initial design. The seed is 20260908. Historical designs are used
only for comparison. The three optimised coordinates are both recuperator
lifts and compressor power fraction; composition is never a search coordinate.
The normalised candidates satisfy the individual bounds and combined lift
limit before calling the model, avoiding projection-induced duplicate designs.

- **RBF:** regularised cubic interpolation, with candidate selection balancing
  predicted SEC against distance from previous evaluations. This adopts the
  response-value/distance idea of
  [Regis and Shoemaker](https://pubsonline.informs.org/doi/10.1287/ijoc.1060.0182).
- **GP-EI Bayesian optimisation:** a Matérn 5/2 Gaussian process with three fitted
  length scales, standardised SEC and expected improvement over a candidate pool.
- **Space-filling control:** maximin selection over the same candidate-generation
  scheme, including local candidates around the best observations.

The surrogate methods include explicit boundary candidates, local perturbations
and periodic maximin exploration. Failed solves are excluded from the SEC
response model; a separate distance heuristic downweights nearby candidates.
This is not a calibrated feasibility model. These lightweight NumPy/SciPy
prototypes are neither BoTorch nor full implementations of published algorithms.

Every winning design is replayed in a fresh process outside the search budget.
The warm pattern-search diagnostic uses the original v11 optimum and an
additional 40-call ceiling, so its costs are reported separately. Candidate
counts include unsuccessful evaluations; the existing timeout and retry
policy remains active.

## Results and evidence

| Fixed composition by mass | v11 COP | RBF COP | GP-EI COP | Space-filling COP |
|---|---:|---:|---:|---:|
| R290/R601 30/70 | 2.646844 | 2.670062 | 2.671504 | 2.619801 |
| R290/R601 60/40 | 2.522545 | 2.536492 | 2.532893 | 2.523878 |
| R290/R600A 30/70 | 2.511081 | 2.541531 | 2.540457 | 2.536549 |

All nine cold trials completed their 64-call budgets. The original v11 runs
used 66, 89 and 89 calls respectively. All 576 cold candidate evaluations
were successful. The three warm diagnostics used 101 additional calls in
total. All 12 selected designs reproduced exactly in fresh processes.
RBF improved COP by 0.553–1.213%; GP-EI improved it by 0.410–1.170%.
RBF finished ahead on two cases and GP-EI on one. This single-seed result
is insufficient to infer a reliable ranking.

The generated `results/RHPcomposition/optimiser_investigation/benchmark_report.md`
contains the final comparison table, evaluation counts, replayed designs,
pressure, VHC and failure counts. `optimiser_convergence.png` and its SVG show
best COP against evaluations. `pentane30_landscape.png` shows the interpolated
design path. Individual JSON checkpoints retain every sampled design and result;
the original composition sweep is preserved as the comparator.

## Recommended approach

Use a **joint global surrogate search followed by a bounded joint local polish**.
The evidence supports changing which variables can move together and providing
global exploration, rather than merely increasing the present line-search
budgets. RBF is a practical first candidate using the current dependencies;
GP-EI is a competitive alternative. A hybrid with final polish is a proposed
next experiment: this pilot tested cold surrogate search and warm v11 polishing
separately, not the complete hybrid.

Before choosing a production default, compare several seeds across all four
binary families and pure endpoints, at matched evaluation budgets. Preserve
and replay the best known design, but distinguish improvements obtained from
those warm starts from independent discovery. Extend the benchmark to pressure
or discharge-temperature constraints if those are part of the intended design
requirements; the current experiments retain the existing model constraints.

For a complete COP–VHC Pareto front, use a multiobjective search at each fixed
composition. Selecting the highest-COP design and plotting its VHC only exposes
one design per composition. BoTorch's
[multiobjective acquisition functions](https://botorch.org/docs/multi_objective),
including qLogNEHVI, provide a suitable route to sampling the trade-off directly.
For solver feasibility, its
[outcome-constraint models](https://botorch.org/docs/constraints) are more
principled than the prototype's distance heuristic. These are recommendations
for a subsequent implementation, not methods benchmarked here.

This pilot does not establish global optima, a robust RBF-versus-BO ranking, or
the accuracy of the thermodynamic model. The nested gas-cooler offset search is
still local. Repeatable selected designs and the two exact-inner endpoint
checks narrow the diagnosis but do not validate every solver branch. The
existing composition curves should therefore be described as **best-found
designs**, not certified optima.

## Reproduction and validation

From the repository root, with the original fixed-composition sweep available:

```sh
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python main_RHPoptimiser_benchmark.py --n-cpu 6
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python main_RHPoptimiser_benchmark.py --method joint_polish --budget 40 --n-cpu 3
.venv/bin/python main_RHPoptimiser_benchmark.py --report-only
```

Cold searches resume matching checkpoints. A changed budget needs a new output
directory. `--seed` and repeated `--refrigerant` / `--method` options support
further trials. The straight-line diagnostic is reproducible from the two
designs tabulated above using 11 equally spaced interpolation fractions and
fresh `_OuterDesignObjective` evaluations; its raw data also contains the two
`use_fast_objective=False` endpoint checks.

The full suite passed with 430 tests. Nine new behavioural tests cover coupled
optimisation, the feasible domain, duplicate avoidance, failed-observation
handling and deterministic checkpoint restoration. No production default was
changed by this investigation.
