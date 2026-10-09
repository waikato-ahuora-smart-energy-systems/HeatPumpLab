# The Heat Pump Lab

The Heat Pump Lab by the Ahuora Centre is a digital experimentation package for constructing and analysing high temperature heat pump cycles. At present, the lab contains one new heat pump cycle, the Recuperative Heat Pump, which is described below.

The [v14 solver notes](HeatPumpLab/cycles/RecuperativeHeatPump/V14_SOLVER_FIXES.md)
describe the pressure-targeting fixes, guarded flash recovery, remaining failed
candidates, and the historical 400-call optimiser comparisons.

The `main_*.py` runner scripts and their tests live in `local_scripts/` (gitignored, not part of
the package). Run them from the repository root, e.g. `.venv/bin/python local_scripts/main_RHPopt.py`;
their tests run with `.venv/bin/python -m pytest local_scripts/tests`.

## Recuperative heat pump model for high temperature glide heating

The idea of heat recuperation within vapour-compression and transcritical heat pump cycles is not new. Internal heat exchangers (IHXs) are widely used to reduce flash gas generation after expansion, preheat suction gas, provide greater superheat after compression, and enable slightly lower high-side pressures (e.g. when producing hot water). In conventional cycles, the IHX transfers heat from the warm refrigerant immediately before expansion to the evaporator pressure. Pinch analysis states that heat is in excess below the pinch temperature and in deficit above it. A heat pump should ideally absorb heat below the pinch and reject upgraded heat above it. Applying this logic to the heat pump cycle itself suggests that repositioned IHXs could deliver thermodynamic benefits. This model calculates heat pump performance for alternate IHX placements. To distinguish a conventional IHX from a repositioned one, this study uses the term recuperative heat exchanger (RHX). If you wish to learn more, we refer you to the paper below.

Cite: T.G. Walmsley, K.-M. Adamson, C. Arpagaus, F. Schlosser, Recuperative heat pump cycle concept for high-temperature glide process heating applications, Applied Thermal Engineering 301 (2026) 131660. https://doi.org/10.1016/j.applthermaleng.2026.131660.

Contact: Tim Walmsley, tim.walmsley@waikato.ac.nz

## Production Bayesian optimisation

`local_scripts/main_RHPopt.py` and the public `optimise_recuperative_heat_pump` API now use
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

Run `.venv/bin/python local_scripts/main_RHPcomposition.py` to compare R600 (butane) with
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

`local_scripts/main_RHPoptimiser_benchmark.py` compares joint RBF, Gaussian-process expected
improvement (Bayesian optimisation), and space-filling searches on three fixed
mixtures at 200°C. All use the same 24 initial designs and a 64-evaluation
budget, without seeding from the previous sweep. It requires the existing
fixed-composition results as comparators. Run with the repository environment:

```sh
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python local_scripts/main_RHPoptimiser_benchmark.py --n-cpu 6
.venv/bin/python local_scripts/main_RHPoptimiser_benchmark.py --report-only
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

The repeated origin-start comparison is available in `local_scripts/main_RHPorigin_benchmark.py`:

```sh
uv sync --group benchmark
MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 .venv/bin/python local_scripts/main_RHPorigin_benchmark.py --n-cpu 6
.venv/bin/python local_scripts/main_RHPorigin_benchmark.py --report-only
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

**Ternary mixtures.** `local_scripts/main_RHPfixed_mobo.py --third <fluid>` searches the
three-component mixture `--first`/`--partner`/`--third` (e.g. Propane/Isobutane/Isopentane):
the composition is encoded as the first two mass fractions with their sum at most 1,
initial and boundary designs cover the triangle's vertices, edges and centroid, and ridge
seeding solves lines at the vertices, edge midpoints and centroid (about 210 evaluations,
no refinement lines). Binary runs are unchanged.

**Compressor efficiency.** From optimiser version v25 the production base case uses
`comp_model_index = 4`: each stage's isentropic efficiency follows the Pareto front of
published efficiency against pressure ratio, eta_is(PR) = 0.87096 - 0.01991 (ln PR)^2.63312
(`pareto_front_efficiency` in `components/compressor.py`; about 0.86 at PR 2, 0.80 at PR 5,
0.69 at PR 10, projected to 0.51 at PR 20 beyond the fitted data), with a 96% motor efficiency. Up to v24 the constant model 1 was used
(0.729 isentropic, 0.70 overall). Results from the two are not comparable.

**Compressor power split constraint.** Every optimisation search limits the
compressor-2 share of total compressor power to 0.4-0.6
(`COMP2_POWER_FRAC_BOUNDS` in `optimisation/_types.py`, which also sets
`DEFAULT_BOUNDS`): the production BO, the fixed-split and temperature Pareto
runners, the surrogate and trust-region searches, and the objective's design
repair. The cycle model itself still solves 0.1-0.9. Earlier results (up to
optimiser version v23) used 0.1-0.9; designs outside the new range should be
excluded when they are compared with new runs.

The fixed-200°C Pareto runner `local_scripts/main_RHPfixed_mobo.py` accepts `--comp-frac 0.5`
to fix the two-stage pressure allocation. Both global BO and local-128 BO then
search only the two recuperator lifts and composition, with three-dimensional
GP inputs. The fixed setting is recorded in checkpoints and retained during
resume and final replay. At zero pressure loss, 0.5 gives equal stage pressure
ratios. Use a fresh output directory when changing this setting.

Optimisation enforces a maximum pressure ratio of 20 per active compressor stage (10 up to
optimiser v24; raised with compressor model 4, whose efficiency falls with PR)
by default (`OuterDesignOptimisationConfig.max_stage_pressure_ratio`). This check
uses the solved stage ratios, including during isolated-worker evaluations and
historical-candidate replays. Missing or invalid active-stage ratios fail the
constraint. Rejected candidates retain their report and consume their evaluation
budget, but cannot become successful incumbents or Pareto points. Historical
results produced without this constraint require re-evaluation; the optimiser
version and run source signatures prevent silent reuse under the new settings.

`local_scripts/main_RHPfixed_mobo.py` seeds each run with pinch-balance ridge solves by default
(`rhp_ridge.py`; disable with `--no-ridge-seed`). At a fixed lift share
theta = RHX1/(RHX1+RHX2), COP peaks where two adjacent gas coolers are pinched at
once; each solve reports the per-gas-cooler minimum approaches (`gc_approaches`),
and the seeding solves A[k+1] = A[k] along 12 composition/theta lines (about 130
evaluations) before BO starts. A second pass adds one line per theta at the
composition midway across the widest COP-VHC gap between neighbouring solved
ridge points (for Propane/Pentane, 5/6 propane, near the 80% ridge-B design). In a
3-seed Propane/Pentane benchmark the first pass alone found the RHX2 = 0 branch of
the front in every seed (2 of 3 without it). `run_family()` defaults to seeding and
needs a budget of at least 200 with it.

The temperature runner `local_scripts/main_RHPtemperature_mobo.py` finds the COP-VHC front at
each delivery temperature on the 100:5:220 °C grid. Temperature is a condition,
not an objective: quality is the mean of the 25 per-temperature 2D hypervolumes
(COP/4, VHC/10 against the origin). One GP (five inputs including temperature, two
outputs) is shared across temperatures; each proposal targets one temperature on a
fixed rotation (adaptive step k takes level 7k mod 25) and uses qLogEHVI against
that temperature's front. Each temperature has four local centres along its front
(borrowed from neighbouring temperatures while a front is sparse) plus, as in the
fixed-split runner, one face-restricted region per zero-lift face (RHX1 = 0,
RHX2 = 0); every candidate pool also includes points on each design bound. `--split-bounds` can narrow
the compressor-2 power fraction further within the global limit below. Without
the face regions (v25), the 200 °C front of a 23,000-evaluation Propane/Pentane run
missed the RHX2 = 0 branch and trailed the fixed-200 °C searches by 1-4% in COP. It uses the same machinery as
the fixed-split runner: batches of `--batch-size` proposals (default 4) with a
pending-point penalty, the constraint-margin GP with the distance penalty for
numerical failures only, `--min-discharge-superheat`, ridge seeding at 100, 120,
..., 220 °C (`--no-ridge-seed`; budget at least 2,000), and a stop rule on the
mean hypervolume (relative change below 2e-4 over three blocks of 1,000
evaluations, not before 5,000; `--no-stop-rule`). Up to v24 temperature was a
third objective, which concentrated BO at 220 °C (454 of 645 steps in a
2,000-evaluation trial); each run version needs a fresh output directory.

`local_scripts/main_RHPehvi_speed_check.py` reports whether BoTorch's fused qLogEHVI C++ kernel
loaded (it needs the `ninja` package and a C++ compiler) and times it against the
pure-Python path.

Each solved design also reports continuous constraint margins
(`report["constraint_margins"]`): per-stage discharge and suction superheat,
`T − T_dew(P)` in K (negative when wet), and pressure-ratio headroom `10 − PR`.
The searches fit a GP to these margins and multiply the acquisition by the
probability that all are non-negative; the old distance penalty now applies
only to numerical failures, which report no margins. Set
`OuterDesignOptimisationConfig.min_discharge_superheat_K` (or
`local_scripts/main_RHPfixed_mobo.py --min-discharge-superheat`) to require a discharge
superheat margin; supercritical discharge is exempt. The default, 0 K, keeps
the previous acceptance rule.

Every active compressor's suction and actual discharge must also be dry vapour
or single-phase supercritical fluid. Saturated vapour (`q=1`) is accepted without
an added superheat margin; any liquid fraction, liquid state, or missing phase
evidence rejects the candidate. The property wrapper preserves REFPROP's native
molar quality/phase flags separately from its legacy quality sentinel. Reports
include `Comp_<stage>_q_refprop_<in|out>` and `Comp_<stage>_phase_<in|out>`;
supercritical fluid is labelled explicitly, since its vapour quality is undefined.
These checks apply to the final solved states, so intermediate iteration states
do not prematurely terminate convergence. Previous results need fresh evaluation.

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
