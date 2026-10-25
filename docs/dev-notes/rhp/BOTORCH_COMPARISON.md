# BoTorch LogEI versus custom BO on v14

BoTorch LogEI completed all nine new 400-evaluation trials. Against the nine matched archived custom-BO trials, median COP improved only for Propane/Pentane 60/40 (+0.0227%); it was slightly lower for Propane/Pentane 30/70 (−0.0035%) and Propane/Isobutane 30/70 (−0.0091%). All 18 saved winner replays match exactly in SEC, COP, pressure and VHC.

For this single-objective problem, retain custom BO as the default: this BoTorch configuration gave similar final COP with 32–35 times the recorded proposal overhead. Those timings came from separate concurrent experiments, and the BoTorch configuration includes 12 continuous acquisition restarts. This result does not establish an inherent package speed ratio or statistical superiority with only three seeds. BoTorch remains available as the tested backend for further work, including future multi-objective development.

BoTorch found the highest individual COP for the 60/40 pentane and 30/70 isobutane cases, but this did not translate into a consistent median advantage. It had 23 failed candidates out of 3,600, versus 21/3,600 for the archived comparator. Two individual GP fit attempts raised errors in the isobutane seed 20260909, at evaluations 286 and 311; there were no recorded acquisition fallbacks. Six implementation tests passed, including objective direction and interrupted-run resume equivalence.

400 unique candidate evaluations per run; three seeds per method and composition. All compositions are mass fractions, at a 200 °C sink target. Every run uses SciPy 1.18.1 and the same frozen v14 thermodynamic model.

3,600 new BoTorch evaluations compared with 3,600 archived custom-BO evaluations; 44 failed candidate evaluations. 18/18 winners match their saved fresh-worker replays (nine new, nine archived) (SEC, COP, pressure and VHC relative tolerance 1e-6, absolute 1e-9). Replays are outside the search budget.

Failed candidates consumed 10.9% of accumulated oracle wall time. Failure classifications: cycle infeasible/unconverged: 40, timeout: 4.

Cross-run consistency: 42 repeated-design groups, covering 150 evaluations; 0 groups with mixed success/failure and 0 with SEC disagreement above 1e-6 relative. Maximum relative SEC spread: 6.73e-11. [Consistency audit](../../../results/RHPcomposition/botorch_v14_400/consistency_audit.json).

See [protocol](BOTORCH_PROTOCOL.md), [source manifest](../../../results/RHPcomposition/botorch_v14_400/manifest.json), [full numerical summary](../../../results/RHPcomposition/botorch_v14_400/summary.json), [six implementation tests](../../../results/RHPcomposition/botorch_v14_400/search_tests.log), [validation](../../../results/RHPcomposition/botorch_v14_400/validation_summary.json), and [archived baseline provenance](../../../results/RHPcomposition/botorch_v14_400/baseline_provenance.json).

## Final COP

Incumbents are selected by minimum SEC. Ranges show the three seeds.

| Composition | Method | Median COP | COP range | Failed / 1,200 |
|---|---|---:|---:|---:|
| Propane/Pentane 30/70 | BoTorch LogEI | 2.672166 | 2.671582–2.672248 | 2 |
| Propane/Pentane 30/70 | Custom BO (archived) | 2.672259 | 2.671888–2.672279 | 3 |
| Propane/Pentane 60/40 | BoTorch LogEI | 2.537006 | 2.536866–2.537182 | 1 |
| Propane/Pentane 60/40 | Custom BO (archived) | 2.536430 | 2.536177–2.537070 | 4 |
| Propane/Isobutane 30/70 | BoTorch LogEI | 2.541491 | 2.541353–2.542173 | 20 |
| Propane/Isobutane 30/70 | Custom BO (archived) | 2.541721 | 2.541291–2.541954 | 14 |

BoTorch's median COP change relative to archived custom BO: Propane/Pentane 30/70 -0.0035%; Propane/Pentane 60/40 +0.0227%; Propane/Isobutane 30/70 -0.0091%.

![Convergence](../../../results/RHPcomposition/botorch_v14_400/convergence.png)

The lower panels show percentage COP gap to the best design observed in this experiment for each composition, not a proven global optimum. Zero gaps are plotted at 0.00001% to permit a logarithmic axis. Lines show medians and shading spans the three seeds.

## Progress and cost

Median incumbent COP at each budget. Timing is median accumulated proposal/oracle wall time per run, measured under concurrent execution, excluding final replay. The custom runs are from the preceding experiment; these are not matched serial hardware timings.

| Composition | Method | 100 | 200 | 400 | Proposal s | Oracle s | Diagnostic events by seed |
|---|---|---:|---:|---:|---:|---:|---|
| Propane/Pentane 30/70 | BoTorch LogEI | 2.669391 | 2.671252 | 2.672166 | 481.9 | 619.4 | [0, 0, 0] |
| Propane/Pentane 30/70 | Custom BO (archived) | 2.669421 | 2.671888 | 2.672259 | 15.2 | 568.8 | not recorded |
| Propane/Pentane 60/40 | BoTorch LogEI | 2.531511 | 2.536767 | 2.537006 | 476.8 | 639.4 | [0, 0, 0] |
| Propane/Pentane 60/40 | Custom BO (archived) | 2.532686 | 2.536177 | 2.536430 | 13.8 | 584.4 | not recorded |
| Propane/Isobutane 30/70 | BoTorch LogEI | 2.531621 | 2.540757 | 2.541491 | 430.2 | 1674.8 | [2, 0, 0] |
| Propane/Isobutane 30/70 | Custom BO (archived) | 2.531145 | 2.540223 | 2.541721 | 13.3 | 1453.7 | not recorded |

## Best design per method and composition

| Composition | Method | Seed | COP | SEC | RHX1 K | RHX2 K | Compressor-2 fraction |
|---|---|---:|---:|---:|---:|---:|---:|
| Propane/Pentane 30/70 | BoTorch LogEI | 20260909 | 2.67224768 | 0.59799754 | 0.000000 | 18.904689 | 0.16513192 |
| Propane/Pentane 30/70 | Custom BO (archived) | 20260910 | 2.67227858 | 0.59798649 | 0.000000 | 19.410289 | 0.16660683 |
| Propane/Pentane 60/40 | BoTorch LogEI | 20260909 | 2.53718161 | 0.65054122 | 0.000000 | 76.489351 | 0.88108210 |
| Propane/Pentane 60/40 | Custom BO (archived) | 20260910 | 2.53706966 | 0.65058860 | 0.000000 | 71.197384 | 0.80633771 |
| Propane/Isobutane 30/70 | BoTorch LogEI | 20260909 | 2.54217310 | 0.64843564 | 19.602427 | 59.874327 | 0.90000000 |
| Propane/Isobutane 30/70 | Custom BO (archived) | 20260911 | 2.54195437 | 0.64852762 | 17.044307 | 62.188283 | 0.90000000 |

## Per-seed outcomes

| Composition | Method | Seed | Final COP | Best evaluation | Failures | Replay match |
|---|---|---:|---:|---:|---:|---|
| Propane/Pentane 30/70 | BoTorch LogEI | 20260909 | 2.67224768 | 332 | 1 | True |
| Propane/Pentane 30/70 | BoTorch LogEI | 20260910 | 2.67216637 | 344 | 0 | True |
| Propane/Pentane 30/70 | BoTorch LogEI | 20260911 | 2.67158211 | 356 | 1 | True |
| Propane/Pentane 30/70 | Custom BO (archived) | 20260909 | 2.67188826 | 179 | 0 | True |
| Propane/Pentane 30/70 | Custom BO (archived) | 20260910 | 2.67227858 | 144 | 2 | True |
| Propane/Pentane 30/70 | Custom BO (archived) | 20260911 | 2.67225879 | 322 | 1 | True |
| Propane/Pentane 60/40 | BoTorch LogEI | 20260909 | 2.53718161 | 391 | 0 | True |
| Propane/Pentane 60/40 | BoTorch LogEI | 20260910 | 2.53686610 | 281 | 1 | True |
| Propane/Pentane 60/40 | BoTorch LogEI | 20260911 | 2.53700602 | 241 | 0 | True |
| Propane/Pentane 60/40 | Custom BO (archived) | 20260909 | 2.53643019 | 307 | 1 | True |
| Propane/Pentane 60/40 | Custom BO (archived) | 20260910 | 2.53706966 | 219 | 2 | True |
| Propane/Pentane 60/40 | Custom BO (archived) | 20260911 | 2.53617730 | 182 | 1 | True |
| Propane/Isobutane 30/70 | BoTorch LogEI | 20260909 | 2.54217310 | 322 | 3 | True |
| Propane/Isobutane 30/70 | BoTorch LogEI | 20260910 | 2.54149059 | 235 | 3 | True |
| Propane/Isobutane 30/70 | BoTorch LogEI | 20260911 | 2.54135263 | 170 | 14 | True |
| Propane/Isobutane 30/70 | Custom BO (archived) | 20260909 | 2.54172078 | 324 | 4 | True |
| Propane/Isobutane 30/70 | Custom BO (archived) | 20260910 | 2.54129111 | 203 | 3 | True |
| Propane/Isobutane 30/70 | Custom BO (archived) | 20260911 | 2.54195437 | 379 | 7 | True |

## BoTorch diagnostics

Recorded search diagnostic events: GP fit: 2. Warnings do not necessarily imply a failed fit or fallback. Full messages and evaluation numbers are retained in each trial checkpoint and the numerical summary. The fitting stage tries two length-scale initialisations; a GP-fit event records an exception from one attempt. Acquisition fallbacks and nonfinite acquisition values have separate event types.

## Limits of the comparison

Three seeds describe variability but do not establish statistical superiority. This is a controlled BoTorch implementation: it retains the custom BO bootstrap, Matérn GP assumptions, failed-region heuristic, boundary candidate pool and periodic maximin exploration. BoTorch supplies the GP fit, LogEI and continuous acquisition optimisation. It does not use all of BoTorch's default GP settings. The experiment does not isolate package fitting differences from the added continuous acquisition optimisation. Archived baseline observations never enter the BoTorch search.

## Reproduction

Install the locked benchmark dependencies with `uv sync --group benchmark --group dev`. The native REFPROP installation and original baseline input are also required. Run `python main_RHPbotorch_benchmark.py --n-cpu 9` from the repository root, with OPENBLAS_NUM_THREADS, OMP_NUM_THREADS and VECLIB_MAXIMUM_THREADS set to 1. Existing completed checkpoints are reused; use a separate output directory for new trials. The comparison folder must also contain the nine custom-BO checkpoints listed in baseline_provenance.json. Generate the report with `python main_RHPbotorch_report.py` and audit this frozen experiment with `python main_RHPbotorch_validate.py`.
