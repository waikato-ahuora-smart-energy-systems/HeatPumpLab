# The Heat Pump Lab

The Heat Pump Lab by the Ahuora Centre is a digital experimentation package for constructing and analysing high temperature heat pump cycles. At present, the lab contains one new heat pump cycle, the Recuperative Heat Pump, which is described below.

The [v14 solver notes](HeatPumpLab/cycles/RecuperativeHeatPump/V14_SOLVER_FIXES.md)
describe the pressure-targeting fixes, guarded flash recovery, remaining failed
candidates, and the historical 400-call optimiser comparisons.

## Recuperative heat pump model for high temperature glide heating

The idea of heat recuperation within vapour-compression and transcritical heat pump cycles is not new. Internal heat exchangers (IHXs) are widely used to reduce flash gas generation after expansion, preheat suction gas, provide greater superheat after compression, and enable slightly lower high-side pressures (e.g. when producing hot water). In conventional cycles, the IHX transfers heat from the warm refrigerant immediately before expansion to the evaporator pressure. Pinch analysis states that heat is in excess below the pinch temperature and in deficit above it. A heat pump should ideally absorb heat below the pinch and reject upgraded heat above it. Applying this logic to the heat pump cycle itself suggests that repositioned IHXs could deliver thermodynamic benefits. This model calculates heat pump performance for alternate IHX placements. To distinguish a conventional IHX from a repositioned one, this study uses the term recuperative heat exchanger (RHX). If you wish to learn more, we refer you to the paper below.

Cite: T.G. Walmsley, K.-M. Adamson, C. Arpagaus, F. Schlosser, Recuperative heat pump cycle concept for high-temperature glide process heating applications, Applied Thermal Engineering 301 (2026) 131660. https://doi.org/10.1016/j.applthermaleng.2026.131660.

Contact: Tim Walmsley, tim.walmsley@waikato.ac.nz

## Production Bayesian optimisation

`main_RHPopt.py` and the public `optimise_recuperative_heat_pump` API now use
joint BoTorch LogEI (v18). Run `uv sync` to install the core dependencies.

```python
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    OuterDesignOptimisationConfig, optimise_recuperative_heat_pump,
)

result = optimise_recuperative_heat_pump(
    refrigerant="Propane;Butane|0.3;0.7|mass",
    config=OuterDesignOptimisationConfig(
        max_exact_evaluations=200, seed=20260910,
        vary_composition=False, preserve_incumbent=False,
    ),
)
```

The origin is `(0 K, 0 K, 0.5)` for free cycle coordinates; fixed settings stay
fixed. The legacy field `comp2_power_frac` represents logarithmic pressure-ratio
allocation in the current model. Both lifts, stage allocation and all enabled
independent mass fractions are searched jointly. Fixed variables are removed
from the GP; variable mixtures use an explicit simplex for any component count.

One shared 200-candidate budget includes initial observations, failures and any
current-solver historical replays. The default bootstrap has 16 observations,
with boundary-aware sampling and periodic maximin exploration. Failed penalties
never enter the response GP. Batch cases run in parallel; selection within a
case is sequential. The edge/split/feedback searches and pattern-search polish,
their configuration fields and their executable benchmark have been removed.

New signatures and `*_bo_v18` output filenames prevent old completed searches
from being mistaken for BO results. The manuscript methods and diagrams are in
`paper/article_outline.tex`; its existing v11 result archive still needs a BO
rerun before the rankings can be attributed to the new method. See
[the production BO protocol](HeatPumpLab/cycles/RecuperativeHeatPump/BAYESIAN_OPTIMISATION.md)
for settings and numerical limitations.

## Fixed-composition refrigerant analysis

Run `.venv/bin/python main_RHPcomposition.py` to compare R600 (butane) with
R290/R600, R290/R601A (isopentane), R290/R601 (pentane), and R290/R600A
(isobutane) at a 200°C sink
target. Each binary uses 0–100% R290 by mass in 10 percentage-point increments.
Composition stays fixed throughout optimisation; recuperator temperature lifts
and compressor-stage allocation remain free. Other cycle inputs retain the standard
defaults (40°C source inlet, 10°C evaporator inlet, 15°C sink inlet).

The sweep solves 41 unique cases, reusing pure-fluid endpoints, and saves
resumable results under `results/RHPcomposition/`. Historical incumbent replay
is disabled for this independent sweep. Use `--dry-run` to inspect the grid,
`--n-cpu 6` to select the worker count, or `--target 180` to change the sink target.
On completion, a comparison CSV and plots are written to `analysis_200C/`
inside the results directory. Use `--report-only` to regenerate them.
Power savings compare electricity at equal delivered heat duty; the separate
SEC column uses the model definition `W / (Q_hot - W)`.

The COP–VHC plot labels R290 mass percentages and highlights designs that are
non-dominated when maximising both metrics. It compares the highest-COP design
found at each sampled composition, rather than performing multi-objective design
optimisation. Fresh isolated replays supply paired COP and VHC values, cached
by case signature and model version. Total VHC is heating duty divided by the
sum of active compressor suction volume flows: `Q_H_total / (V_suc_1 + V_suc_2)`
for two compressors. It is converted from the model's kJ/m³ to MJ/m³ for plotting.

## Experimental optimiser comparison

`main_RHPoptimiser_benchmark.py` compares joint RBF, Gaussian-process expected
improvement (Bayesian optimisation), and space-filling searches on three fixed
mixtures at 200°C. All use the same 24 initial designs and a 64-evaluation
budget, without seeding from the previous sweep. It requires the existing
fixed-composition results as comparators. Run with the repository environment:

```sh
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python main_RHPoptimiser_benchmark.py --n-cpu 6
.venv/bin/python main_RHPoptimiser_benchmark.py --report-only
```

Historical local-polish results remain readable, but that diagnostic is no longer executable.
Results, evaluation checkpoints, fresh-process winner replays and
convergence plots are saved under `results/RHPcomposition/optimiser_investigation/`.
Use `--seed`, `--refrigerant`, `--method` and `--output` for further trials. Cold
searches resume matching checkpoints; a changed budget requires a new output
directory. These experiments do not change the production optimiser or replace
the original composition sweep. See the
[investigation](HeatPumpLab/cycles/RecuperativeHeatPump/OPTIMISER_COMPARISON.md)
for findings and limitations.

The repeated origin-start comparison is available in `main_RHPorigin_benchmark.py`:

```sh
uv sync --group benchmark
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python main_RHPorigin_benchmark.py --n-cpu 6
.venv/bin/python main_RHPorigin_benchmark.py --report-only
```

This runs BO, RBF and official pycma CMA-ES from `(0 K, 0 K, 0.5)` on the three
pilot compositions, using three seeds and 200 candidate evaluations per run.
Initialisation counts against the budget; no historical designs seed the search.
Composition remains fixed. Results and replays are checkpointed under
`results/RHPcomposition/origin_200/`. Use `--seed`, `--method`, `--refrigerant`
and `--output` to select further trials. Resume requires matching settings and
search-code hashes. The report validates complete budgets, origin, composition,
unique candidates and winner replays before aggregating the repeats.
See the [origin-start study](HeatPumpLab/cycles/RecuperativeHeatPump/ORIGIN_OPTIMISER_COMPARISON.md)
for the repeated-trial results and the inner solver's sensitivity to numerical
pressure initialisation, which limits interpretation of the nominal COP peaks.

## Tests and archived experiments

Run the maintained test suite from the repository root with
`.venv/bin/python -m pytest`. The `experiments/` directory is optional: maintained
search modules, reporting entry points and tests import the root source and
`HeatPumpLab/` package directly. Local Pareto continuation checkpoints are
restored by `restore_search` in `rhp_local_pareto_search.py`.

Archived experiment runners retain their frozen implementations for historical
reproduction. Run any archive-specific tests in a separate process from that
archive's directory; they are not dependencies of the maintained test suite.

The two historical paper-figure audits run only when the local figure generator
under `paper/figures/` and the saved v11 dataset under `results/RHPopt/` are
available. They explicitly skip otherwise, so a source checkout does not require
those ignored research artifacts to collect or run the suite.
