# Frozen v14 comparison: multistart COBYQA, TuRBO-1 and BO

Protocol fixed before the main runs, 9 September 2026 (New Zealand time).

## Common experiment

- Fixed mass fractions: Propane/Pentane 0.30/0.70 and 0.60/0.40;
  Propane/Isobutane 0.30/0.70. Sink target 200 °C. All other physical
  settings loaded by `load_case` from the existing v11 composition table;
  historical scores and optimised designs are not used to guide these searches.
- Minimise SEC in three coordinates: RHX1 and RHX2 lifts in [0,100] K;
  compressor-2 fraction in [0.1,0.9]. Combined lift <= the common objective's
  physical limit (123.333333333 K). All methods operate in normalised coordinates.
- Each run starts at [0 K,0 K,0.5]. Different subsequent initialisation policies
  are part of the methods being compared, and their evaluations count.
- 400 unique expensive candidate evaluations per run, including failed solves;
  seeds 20260909, 20260910, 20260911: 27 runs, 10,800 search evaluations.
  Cheap known-domain rejections and cache hits do not consume expensive calls.
- Same isolated native REFPROP oracle: timeout 180 seconds, one timeout retry,
  worker recycled after 20 candidate calls. Retries are internal to one budgeted
  candidate, so 400 candidates need not mean 400 underlying cycle solve attempts.
- Frozen production model SHA256:
  `812b3bebdaa94bd653a110f3cb2120c424bbdce01c886b4768dfb38ee1bea087`.
  Experimental code lives outside `HeatPumpLab` to preserve this model hash.
- SciPy 1.18.1, verified against PyPI and imported from the project environment;
  other versions and complete model/search hashes recorded in every checkpoint.
- All 27 runs use fresh v14 evaluations. Earlier v13/v14 search results are not
  substituted for the BO comparator.
- Fresh worker replay of every SEC-selected run winner, outside search budget.
  Match criterion: relative SEC tolerance 1e-6, absolute 1e-9.
- Checkpoint every expensive evaluation. Resumption reconstructs the complete
  proposal trajectory and rejects incompatible source/settings or proposals.

## Methods and predetermined settings

### Multistart COBYQA

SciPy `minimize(method="COBYQA")`, explicit box and linear combined-lift
constraints. First start is the origin. Subsequent starts maximise distance
to all evaluated points among 1,024 seeded scrambled Sobol proposals filtered
to the allowed domain. Each local start gets up to 80 new expensive evaluations;
early convergence releases unused budget to additional starts. Every run spends
exactly 400 evaluations. Initial normalised trust radius 0.15, final radius
1e-5, feasibility tolerance 1e-10. Failed solves return a fixed 1e6 penalty;
their SEC remains missing in the stored data. No post-hoc polish or tuning.

### TuRBO-1 (NumPy/SciPy adaptation)

Sequential single-region TuRBO with restart-local GP data and joint Thompson
sampling. First bootstrap is the origin plus five accepted Latin-hypercube
proposals; restarts use six accepted global proposals. Rejection means the
accepted constrained subset is no longer a strict Latin hypercube.

Exact Matérn-5/2 GP with ARD length scales, fitted constant mean, output variance
and observation noise. Responses standardised by median and standard deviation.
Length-scale bounds [0.005,2], output variance [0.05,20], noise variance
[0.0005,0.2]. SciPy L-BFGS-B maximum likelihood fit every iteration, up to 50
iterations, warm-started within each region. This differs from the authors'
GPyTorch/Adam training implementation; no upstream code is vendored.

Trust length initially 0.8, capped at 1.6. Three consecutive meaningful
improvements double it; four unsuccessful steps halve it; below 0.5^7 the
region restarts. Meaningful SEC decrease is >0.001*abs(previous best SEC).
ARD weights have geometric mean one. Select up to 300 admissible, previously
unevaluated candidates from a scrambled Sobol pool inside the trust region.
Sample the full correlated noisy posterior, rather than independent marginal
normals. In three dimensions all axes are perturbed. No explicit boundary
snapping or extra global steps beyond restarts.

Domain constraints and failures require adaptations to unconstrained TuRBO:
failed solves never enter its SEC GP, count as unsuccessful trust-region steps,
and contribute an additive standardised-score penalty (1-feasibility) using
the existing BO distance heuristic. If fewer than three local successful
observations exist after bootstrap, use global maximin proposals until enough
data or a restart. The heuristic is not a calibrated feasibility classifier.

### Existing BO (GP-EI)

Unchanged `OriginSurrogateSearch(method="gp_ei")`: origin plus four local
bootstrap probes; global and incumbent-centred candidate pool with explicit
box-face samples; exact Matérn-5/2 ARD GP; expected improvement weighted by the
existing failed-region distance heuristic; every eighth adaptive step maximin.
Its existing fitting schedule/noise settings are preserved. Failed observations
are excluded from the SEC GP. This compares complete search implementations,
including their initialisation, fitting and failure adapters, not just their
acquisition functions.

## Reporting

Primary comparison: distribution of final incumbent COP (incumbent selected by
minimum SEC), median and range over three seeds. Also report SEC, full winning
designs, failures, convergence at 100/200/400 evaluations, proposal/oracle time,
restarts and replay results. Three seeds describe variability and do not support
strong statistical claims. Wall times are measured under concurrent execution.

## Primary references

- [SciPy COBYQA interface](https://docs.scipy.org/doc/scipy/reference/optimize.minimize-cobyqa.html)
- [TuRBO paper](https://arxiv.org/abs/1910.01739)
- [TuRBO reference implementation](https://github.com/uber-research/TuRBO)
- [SciPy release on PyPI](https://pypi.org/project/scipy/1.18.1/)
