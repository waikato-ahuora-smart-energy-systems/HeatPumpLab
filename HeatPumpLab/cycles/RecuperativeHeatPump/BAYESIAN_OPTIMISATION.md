# Production Bayesian optimisation (v18)

The public optimiser and both production runners now use a joint BoTorch LogEI
search. This promotes the experimental GP/acquisition approach to production and
extends it beyond its original three-dimensional fixed-composition benchmark.
The previous edge, split, feedback and pattern-search implementation has been
removed, including its result types, configuration and executable benchmark.
Historical reports and archived comparison data remain available.

## Domain and initialisation

The physical vector remains `[RHX1 lift, RHX2 lift, stage allocation, mass fractions...]`.
The field `comp2_power_frac` is the current model's logarithmic pressure-ratio
allocation coordinate. Free cycle coordinates are normalised to unit bounds;
fixed coordinates are omitted from the GP. A variable n-component mixture uses
n−1 independent fractions, with the last obtained by mass closure. Explicit
linear constraints enforce the simplex and total recuperator-lift allowance.
Pure-fluid and fixed-composition calculations omit composition coordinates.

The first evaluation is the supplied mass composition at zero free recuperator
lifts and stage allocation 0.5 by default. Prescribed variables stay fixed.
Up to 16 initial observations comprise feasible ±0.05 normalised-coordinate
probes, pure vertices and maximin completion. All variables remain jointly
available thereafter. There is no staged freezing of mixture recuperator lifts.

The candidate pool contains 512 base samples plus explicit boundaries. Cycle
coordinates use scrambled Sobol samples; composition uses uniform Dirichlet
samples. Coordinate boundaries, pure vertices, the last-fraction-zero face and
the zero-recuperation face receive explicit proposals. Points violating geometry
or duplicating observations are discarded before model evaluation. Boundary
coverage is finite and does not certify all lower-order optima have been found.

## Model and acquisition

Successful observations train a double-precision Matérn-5/2 GP with one length
scale per free coordinate. The response is negative standardised SEC; minimising
SEC is equivalent to maximising COP for accepted heating states. Numerical
variance is fixed at 1e-6 on this standardised scale. GP hyperparameters are fit
at each model-guided step. The old experiment's 0.02 normalised minimum length
scale is not imposed on the production model.

BoTorch LogEI is scored across the pool and compared with continuous constrained
acquisition proposals (six restarts, 128 raw samples, at most 60 acquisition
iterations). These are surrogate computations, not heat-pump evaluations.
Every eighth post-initialisation step uses maximin exploration. Insufficient
successful observations, model-fitting failures or nonfinite scores also use
maximin; acquisition-refinement failure retains the pool winner.

Failed objective penalties never train the response GP. When successful and
failed observations coexist, log-EI is augmented by
`log((d_failed + 0.02) / (d_success + d_failed + 0.04))`. These nearest-neighbour
distances discourage revisiting failed regions; this is not a calibrated
feasibility probability. Diagnostics are exposed in the returned result.

References: [BoTorch](https://arxiv.org/abs/1910.06403) and
[LogEI](https://arxiv.org/abs/2310.20708).

## Budget, output and migration

The default budget is 200 distinct model evaluations. Origin, initial samples,
failed candidates and compatible historical replays all count. Historical
preservation defaults on, admits at most three unseen compatible prior designs
as budget permits, and always evaluates current SEC before training. Set
`preserve_incumbent=False` for a cold-start study. No additional polish allowance
exists. A fully fixed problem evaluates once; a zero budget returns no winner.

Model workers retain the 180-second timeout, one retry and recycling after 20
requests. Retries use the same candidate; worker cleanup occurs in a `finally`
block. Batch cases can run in parallel, while each BO trajectory is sequential.
The selected design is the lowest-SEC successful observation; all-failure cases
have no winner. The API exposes the full observation history and diagnostics.

`OuterDesignOptimisationConfig` retains fixed/vary flags and worker controls,
and adds `seed`, `initial_evaluations`, `candidate_pool_size`,
`acquisition_restarts` and `exploration_interval`. Removed settings such as
`line_budget`, `split_budget` and `polish_budget` now raise an argument error.
Results expose `initial_x`, `best`, `evaluations` and `incumbent_replays`; the
old corner/edge/split/polish result structure has been removed.

BoTorch is a core dependency (`uv sync`); CMA-ES remains an optional benchmark
dependency. The version string is `isolated-refprop-botorch-logei-v18`.
Production filenames end in `bo_v18`, and versioned signatures prevent accidental
reuse of completed staged-search results. Experimental runners retain their own
protocols and must not be presented as this generic production configuration.

## Paper and limits of the evidence

`paper/article_outline.tex`, its bibliography and three method diagrams describe
this implementation, including joint constraints, LogEI, failure avoidance and
shared budgets. The existing v11 numerical tables are historical: the full
refrigerant ranking has not been regenerated using BO in this change.

Tests cover joint improvement on a coupled objective, fixed-seed reproduction,
arbitrary component counts, fixed-variable preservation, simplex vertices,
all-failure exploration, all-fixed and zero-budget cases, resource cleanup,
budget accounting and real REFPROP smoke execution. These checks are not a
new thermodynamic performance benchmark.

BO does not cure discontinuities or pressure-root initialisation dependence in
the inner model. Physical interpretation still requires independent-seed,
pressure-guess, small-perturbation and boundary-identity checks. A best-found COP
is neither a global-optimality certificate nor a full COP–VHC Pareto front.

Validation on 10 September 2026: all 546 tests passed, including native REFPROP
smoke tests. Ruff, whitespace checks and the offline lockfile consistency check
passed. The updated 40-page manuscript compiled successfully; its changed pages
and regenerated figures were rendered and visually checked.
