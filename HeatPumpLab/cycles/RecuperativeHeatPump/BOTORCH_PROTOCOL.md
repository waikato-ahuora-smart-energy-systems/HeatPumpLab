# BoTorch versus custom BO: frozen v14 benchmark

Protocol fixed before thermodynamic search runs, 10 September 2026.

## Scope and budget

Single-objective minimisation of SEC on the same three fixed mass compositions:
Propane/Pentane 30/70, Propane/Pentane 60/40 and Propane/Isobutane 30/70.
Sink target 200 °C; other physical settings come from the same saved baseline
table as the preceding comparison. Optimise the two recuperator lifts in
[0,100] K and compressor-2 fraction in [0.1,0.9], with combined lift bounded by
123.333333333 K. The composition is held fixed during every trial.

Nine new BoTorch runs: 400 unique expensive candidate evaluations, including
failed candidates, for each composition and seeds 20260909, 20260910, 20260911.
Bootstrap evaluations count. Cheap domain rejections and duplicate proposals
do not consume expensive calls. Final winner replays are outside the budget.

The nine archived custom BO runs from `trust_region_v14_400` provide the
comparator. Before reuse, verify the complete model and custom-search source
hashes, physical settings, budgets, seeds, timeout settings, successful winner
replays and NumPy/SciPy/ctREFPROP versions. Archive byte-identical copies and
record their source paths and hashes. These are previous matched observations,
not nine new runs. No archived observations guide BoTorch's search.

All native calls use the same v14 oracle, 180-second timeout, one timeout retry,
and worker recycling after 20 candidate evaluations. One budgeted candidate
can involve multiple internal solve attempts. Successful observations have
been deterministic; failures are handled separately from observation noise.

Frozen model SHA256:
`812b3bebdaa94bd653a110f3cb2120c424bbdce01c886b4768dfb38ee1bea087`.
The new search and runner live outside `HeatPumpLab` to preserve this hash.

## Controlled implementation comparison

The intent is to evaluate the BoTorch engine with comparable modelling and
exploration policies, rather than compare the custom engine against unrelated
default model choices. The following policies remain aligned:

- Same origin [0 K,0 K,0.5] and four local bootstrap probes.
- Normalised [0,1]^3 inputs; successful SEC values standardised by their current
  mean and population standard deviation. Failed responses remain missing.
- Exact Matérn-5/2 GP with independent length scales in [0.02,3], zero mean,
  unit kernel amplitude and fixed observation variance 1e-6 in standardised
  output units. Fit length scales at the first adaptive call and every fifth
  observation, with two starting points: previous scales and [0.15,0.15,0.15].
- Existing global/local candidate pool, including explicit box faces.
- Existing maximin step every eighth adaptive call (also used with fewer than
  five successful observations).
- Existing distance-based feasibility heuristic, applied as log(feasibility)
  plus LogEI. This heuristic is not a calibrated failure classifier.

BoTorch supplies `SingleTaskGP`, GPyTorch's exact marginal likelihood fitting
via `fit_gpytorch_mll`, and analytic `LogExpectedImprovement` for q=1. It
maximises negative standardised SEC. Native GP fitting permits up to 100
iterations per initialisation, with one package fitting attempt per start;
warnings and failures are recorded. If neither start fits, use the previous
length scales. This differs from the custom direct NumPy/SciPy likelihood fit.

BoTorch also continuously optimises the acquisition with `optimize_acqf`:
12 restarts, 512 raw initial samples, maximum 100 iterations, restart batch
limit four, explicit box and linear lift constraints. Choose the best valid,
novel candidate from those results and the shared boundary/global/local pool.
Reject numerical domain violations and duplicates before the expensive call.
Acquisition optimisation errors fall back to the shared candidate pool and
are recorded. The custom comparator selects EI maxima from its candidate pool
without continuous acquisition optimisation. Thus this comparison tests both
package numerical routines and the additional acquisition optimisation, and
does not isolate which difference causes any performance change.

Use CPU double precision and one PyTorch/BLAS thread per trial, with up to
nine concurrent trials. Per-iteration PyTorch random seeds and ask/tell replay
make checkpoint reconstruction deterministic. Import BoTorch in the search
worker only, so native oracle workers do not need to import PyTorch.

Versions at setup: BoTorch 0.18.1, GPyTorch 1.15.2, PyTorch 2.14.0,
NumPy 2.4.6, SciPy 1.18.1, ctREFPROP 0.10.5. Record all versions and source
hashes in every new checkpoint and retain the locked dependency environment.

## Evaluation

Compare median and range of final incumbent COP over three seeds, choosing
incumbents by minimum SEC. Also report best designs, COP versus evaluation
count, proposal and oracle time, failures, acquisition/fit diagnostics and
fresh-winner reproducibility (SEC, COP, pressure and VHC). Replays use a new
worker and relative tolerance 1e-6 / absolute tolerance 1e-9.

Wall-time observations come from separate concurrent executions; these are
not matched serial hardware timings. Three seeds reveal variability but do
not establish statistical superiority. The experiment does not implement
multi-objective BO or change the thermodynamic model.

## References

- [BoTorch](https://botorch.org/)
- [SingleTaskGP documentation](https://botorch.readthedocs.io/en/stable/models.html#botorch.models.gp_regression.SingleTaskGP)
- [LogEI documentation](https://botorch.readthedocs.io/en/stable/acquisition.html#botorch.acquisition.analytic.LogExpectedImprovement)
- [Acquisition optimisation](https://botorch.readthedocs.io/en/stable/optim.html#botorch.optim.optimize.optimize_acqf)
- [LogEI numerical reformulation](https://arxiv.org/abs/2310.20708)
