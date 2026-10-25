# Origin-start BO, RBF and CMA-ES comparison

This follow-up uses the user-confirmed origin `(RHX1, RHX2, compressor 2 power
fraction) = (0 K, 0 K, 0.5)`. Composition remains fixed by mass. The three cases
are R290/R601 30/70, R290/R601 60/40 and R290/R600A 30/70, at a 200°C sink target.
Each method receives 200 candidate evaluations for each of three random seeds
(20260909, 20260910, 20260911). Previous designs are comparison values only.

## Protocol

The origin is evaluation 1 for every run. BO and RBF then use four local
observations: each recuperator raised separately by 5 K and compressor power
fractions 0.54 and 0.46. This replaces the previous global 24-point prescreen.
Their subsequent proposals can explore the full feasible domain. Both minimise
SEC, which is equivalent to maximising COP over the feasible results here.

CMA-ES uses official pycma 4.4.4, starting its distribution at the same origin,
with a normalised standard deviation of 0.25 and population size seven. It uses
the library's box-bound transform and resamples points outside the known
combined-lift limit before model evaluation. Native stopping permits a fresh
origin-centred restart within the total budget. The final partial population is
evaluated and eligible to win, but does not update the sampling distribution.
These choices follow the
[official ask/tell interface](https://cma-es.github.io/apidocs-pycma/cma.evolution_strategy.CMAEvolutionStrategy.html).

Initialisation, failed model evaluations and restarts count against the same
200-candidate budget. A native timeout retry does not create a new candidate
entry. Every selected winner receives an additional fresh-process replay.
The inner model and production optimiser are unchanged.

## Results

<!-- ORIGIN_SUMMARY -->

The complete results, each seed's design, COP, VHC and pressure, and the
convergence plot are generated in
`results/RHPcomposition/origin_200/origin_report.md`. The shaded convergence
bands show the observed minimum–maximum across seeds, not confidence intervals.
The validation checks budgets, the first evaluation, constant requested
composition, bounds, duplicate candidates and fresh winner replays.

## A separate problem with the inner objective

The high-scoring first CMA-ES design for R290/R601 30/70 is
`(0.016142113 K, 16.268979770 K, 0.169879228)`. Its default COP of 2.769876
reproduced in a fresh process. However, changing only the numerical pressure
starting guess changes its result substantially:

| Initial pressure guess, MPa | CMA-ES-selected design COP | BO-selected design COP |
|---|---:|---:|
| Automatic | 2.769876 | 2.671786 |
| 3.0 | 2.664067 | 2.671792 |
| 4.3 | 2.664030 | 2.671755 |
| 5.5 | 2.664053 | 2.671775 |
| 7.0 | 2.691406 | 2.671725 |

Physical design inputs remain identical within each column. The CMA-ES design
selects pressures from 5.121 to 6.907 MPa; the BO design remains near 6.638 MPa.
This shows that default replay reproducibility is insufficient to establish
independence from numerical initialisation. It does not establish that the
higher-COP solution is physically invalid.

An earlier CMA-ES candidate provides another sensitivity check. With RHX2 and
compressor split fixed, RHX1 lifts around 0.0058–0.01 K give COP 2.7424, while
0, 0.001 and 0.1 K give about 2.623. Tiny positive lifts also return the lower
result, so the effect is not simply switching the recuperator on at zero.
Disabling the fast inner objective reproduces the captured high result and
the zero-lift result. The pressure/evaporator solution path needs investigation
beyond the fast gas-cooler offset approximation.

The current GP length-scale lower bound and the smallest local surrogate
candidate standard deviation correspond to 2 K in the recuperator coordinates.
These settings are poorly suited to features at hundredths or thousandths of a
kelvin. RBF can encounter them by chance; CMA-ES adapts its sampling scale.
Results therefore compare these implementations and settings, not all possible
BO, RBF and CMA-ES configurations.

## Interpretation

A robust production solution needs a consistent inner objective as well as a
good outer search. Validate pressure/evaporator branch selection against
different initial guesses and small design perturbations, then repeat the
comparison on that controlled objective. Treat the current scores as nominal
model results, not certified global thermodynamic optima. This study does not
establish a full COP–VHC Pareto front.

## Reproduction

From the repository root, with the previous fixed-composition sweep available:

```sh
uv sync --group benchmark
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python main_RHPorigin_benchmark.py --n-cpu 6
.venv/bin/python main_RHPorigin_benchmark.py --report-only
```

The runner checkpoints every evaluation and resumes matching settings and
search-code hashes. `--seed`, `--method`, `--refrigerant`, `--budget` and
`--output` support further trials. Keep different budgets in separate output
directories. The benchmark dependency group adds pycma without changing the
core model dependencies.

The pressure-guess diagnostic is reproducible using the two selected designs
in `df3fb5a45e_cma_es_20260909.json` and `df3fb5a45e_gp_ei_20260909.json`:
evaluate each with base-case `P_H_est` values 0, 3, 4.3, 5.5 and 7 MPa in fresh
workers. `pressure_seed_diagnostic.json` contains the raw results.
`boundary_diagnostic.json` preserves its captured design and every perturbation.
These additional diagnostic evaluations were never fed back into the searches.

The full test suite passes with 437 tests. Seven new tests exercise origin
initialisation, a 200-evaluation coupled problem, physical bounds, uniqueness,
local bootstrap rank and deterministic restoration of each method's history.
