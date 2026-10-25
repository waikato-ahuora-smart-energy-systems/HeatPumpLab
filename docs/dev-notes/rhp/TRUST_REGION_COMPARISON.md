# COBYQA / TuRBO / BO comparison on v14

400 unique candidate evaluations per run; three seeds per method and composition. All compositions are mass fractions, at a 200 °C sink target. Every run uses SciPy 1.18.1 and the same frozen v14 thermodynamic model.

10,800 search evaluations; 49 failed candidate evaluations. 27/27 winners match in a fresh worker (SEC and COP relative tolerance 1e-6, absolute 1e-9). Replays are outside the search budget.

Failed candidates consumed 11.2% of accumulated oracle wall time. Failure classifications: cycle infeasible/unconverged: 42, timeout: 7.

Cross-run consistency: 368 repeated-design groups, covering 944 evaluations; 0 groups with mixed success/failure and 0 with SEC disagreement above 1e-6 relative. Maximum relative SEC spread: 0. [Consistency audit](../../../results/RHPcomposition/trust_region_v14_400/consistency_audit.json).

See [protocol](TRUST_REGION_PROTOCOL.md), [source manifest](../../../results/RHPcomposition/trust_region_v14_400/manifest.json), [full numerical summary](../../../results/RHPcomposition/trust_region_v14_400/summary.json), [506-test regression log](../../../results/RHPcomposition/trust_region_v14_400/pytest.log), and [three interruption/resume tests](../../../results/RHPcomposition/trust_region_v14_400/resume_tests.log).

## Interpretation

BO is the strongest default among these implementations at a 400-evaluation budget: it has the highest median COP on Propane/Pentane 30/70 and Propane/Isobutane 30/70, with narrow seed ranges on all three cases. Three seeds support this practical choice but do not establish statistical superiority.

Multistart COBYQA is useful when early improvement matters. On both Pentane mixtures it reaches near-final performance in roughly 50 evaluations. It also finds the best individual design on two compositions, but its Isobutane COP varies from 2.526484 to 2.542231 across seeds. Its final Pentane 60/40 result is identical across seeds because none of the later starts beats the first local solution.

TuRBO-1 has a 0.014% higher median COP than BO on Pentane 60/40, but its seed range is much wider (2.525356–2.537889). It has lower median COP on the other two cases. This experiment does not show a compelling reason to replace BO with this TuRBO implementation.

The winning Pentane 30/70 designs sit on, or very near, zero RHX1 lift; the best Isobutane designs sit at the upper compressor-2 fraction bound. BO's explicit boundary proposals may help in these cases, but this comparison does not isolate that effect from the other implementation differences.

All 27 fresh replays agree exactly in SEC, COP, maximum pressure and VHC. The production model hash remains unchanged. Of 10,800 candidate evaluations, 49 fail (42 infeasible/unconverged outcomes and 7 timeouts), consuming 11.2% of accumulated oracle wall time. The coarse infeasible/unconverged label does not by itself distinguish a physical violation from numerical non-convergence.

## Final COP

Incumbents are selected by minimum SEC. Ranges show the three seeds.

| Composition | Method | Median COP | COP range | Failed / 1,200 |
|---|---|---:|---:|---:|
| Propane/Pentane 30/70 | Multistart COBYQA | 2.672062 | 2.671927–2.672345 | 2 |
| Propane/Pentane 30/70 | TuRBO-1 (adapted) | 2.671305 | 2.671108–2.672242 | 0 |
| Propane/Pentane 30/70 | BO (GP-EI) | 2.672259 | 2.671888–2.672279 | 3 |
| Propane/Pentane 60/40 | Multistart COBYQA | 2.535199 | 2.535199–2.535199 | 2 |
| Propane/Pentane 60/40 | TuRBO-1 (adapted) | 2.536787 | 2.525356–2.537889 | 1 |
| Propane/Pentane 60/40 | BO (GP-EI) | 2.536430 | 2.536177–2.537070 | 4 |
| Propane/Isobutane 30/70 | Multistart COBYQA | 2.538619 | 2.526484–2.542231 | 15 |
| Propane/Isobutane 30/70 | TuRBO-1 (adapted) | 2.539817 | 2.538522–2.541340 | 8 |
| Propane/Isobutane 30/70 | BO (GP-EI) | 2.541721 | 2.541291–2.541954 | 14 |

![Convergence](../../../results/RHPcomposition/trust_region_v14_400/convergence.png)

The lower panels show percentage COP gap to the best design observed in this experiment for each composition, not a proven global optimum. Zero gaps are plotted at 0.00001% to permit a logarithmic axis. Lines show medians and shading spans the three seeds.

## Progress and cost

Median incumbent COP at each budget. Timing is median accumulated proposal/oracle wall time per run, measured under concurrent execution, excluding final replay. It is not a serial hardware benchmark.

| Composition | Method | 100 | 200 | 400 | Proposal s | Oracle s | Restarts by seed |
|---|---|---:|---:|---:|---:|---:|---|
| Propane/Pentane 30/70 | Multistart COBYQA | 2.671927 | 2.671927 | 2.672062 | 0.5 | 635.5 | [7, 5, 6] |
| Propane/Pentane 30/70 | TuRBO-1 (adapted) | 2.671031 | 2.671066 | 2.671305 | 3.1 | 614.7 | [9, 9, 8] |
| Propane/Pentane 30/70 | BO (GP-EI) | 2.669421 | 2.671888 | 2.672259 | 15.2 | 568.8 | [0, 0, 0] |
| Propane/Pentane 60/40 | Multistart COBYQA | 2.535199 | 2.535199 | 2.535199 | 0.5 | 579.9 | [6, 5, 6] |
| Propane/Pentane 60/40 | TuRBO-1 (adapted) | 2.524135 | 2.525356 | 2.536787 | 3.0 | 670.2 | [8, 8, 8] |
| Propane/Pentane 60/40 | BO (GP-EI) | 2.532686 | 2.536177 | 2.536430 | 13.8 | 584.4 | [0, 0, 0] |
| Propane/Isobutane 30/70 | Multistart COBYQA | 2.525998 | 2.526008 | 2.538619 | 0.5 | 1580.6 | [6, 6, 6] |
| Propane/Isobutane 30/70 | TuRBO-1 (adapted) | 2.534655 | 2.538522 | 2.539817 | 3.1 | 1749.0 | [8, 8, 8] |
| Propane/Isobutane 30/70 | BO (GP-EI) | 2.531145 | 2.540223 | 2.541721 | 13.3 | 1453.7 | [0, 0, 0] |

## Best design per method and composition

| Composition | Method | Seed | COP | SEC | RHX1 K | RHX2 K | Compressor-2 fraction |
|---|---|---:|---:|---:|---:|---:|---:|
| Propane/Pentane 30/70 | Multistart COBYQA | 20260909 | 2.67234535 | 0.59796262 | 0.000000 | 19.993456 | 0.16824763 |
| Propane/Pentane 30/70 | TuRBO-1 (adapted) | 20260911 | 2.67224223 | 0.59799949 | 0.008328 | 20.499410 | 0.16968418 |
| Propane/Pentane 30/70 | BO (GP-EI) | 20260910 | 2.67227858 | 0.59798649 | 0.000000 | 19.410289 | 0.16660683 |
| Propane/Pentane 60/40 | Multistart COBYQA | 20260909 | 2.53519885 | 0.65138142 | 0.029716 | 62.841549 | 0.69208510 |
| Propane/Pentane 60/40 | TuRBO-1 (adapted) | 20260909 | 2.53788887 | 0.65024204 | 2.583182 | 72.336607 | 0.84926671 |
| Propane/Pentane 60/40 | BO (GP-EI) | 20260910 | 2.53706966 | 0.65058860 | 0.000000 | 71.197384 | 0.80633771 |
| Propane/Isobutane 30/70 | Multistart COBYQA | 20260909 | 2.54223065 | 0.64841144 | 18.704872 | 60.706999 | 0.90000000 |
| Propane/Isobutane 30/70 | TuRBO-1 (adapted) | 20260909 | 2.54133974 | 0.64878623 | 22.296734 | 57.153869 | 0.89668443 |
| Propane/Isobutane 30/70 | BO (GP-EI) | 20260911 | 2.54195437 | 0.64852762 | 17.044307 | 62.188283 | 0.90000000 |

## Per-seed outcomes

| Composition | Method | Seed | Final COP | Best evaluation | Failures | Replay match |
|---|---|---:|---:|---:|---:|---|
| Propane/Pentane 30/70 | Multistart COBYQA | 20260909 | 2.67234535 | 272 | 0 | True |
| Propane/Pentane 30/70 | Multistart COBYQA | 20260910 | 2.67192743 | 51 | 1 | True |
| Propane/Pentane 30/70 | Multistart COBYQA | 20260911 | 2.67206236 | 275 | 1 | True |
| Propane/Pentane 30/70 | TuRBO-1 (adapted) | 20260909 | 2.67110843 | 159 | 0 | True |
| Propane/Pentane 30/70 | TuRBO-1 (adapted) | 20260910 | 2.67130507 | 209 | 0 | True |
| Propane/Pentane 30/70 | TuRBO-1 (adapted) | 20260911 | 2.67224223 | 275 | 0 | True |
| Propane/Pentane 30/70 | BO (GP-EI) | 20260909 | 2.67188826 | 179 | 0 | True |
| Propane/Pentane 30/70 | BO (GP-EI) | 20260910 | 2.67227858 | 144 | 2 | True |
| Propane/Pentane 30/70 | BO (GP-EI) | 20260911 | 2.67225879 | 322 | 1 | True |
| Propane/Pentane 60/40 | Multistart COBYQA | 20260909 | 2.53519885 | 57 | 0 | True |
| Propane/Pentane 60/40 | Multistart COBYQA | 20260910 | 2.53519885 | 57 | 1 | True |
| Propane/Pentane 60/40 | Multistart COBYQA | 20260911 | 2.53519885 | 57 | 1 | True |
| Propane/Pentane 60/40 | TuRBO-1 (adapted) | 20260909 | 2.53788887 | 286 | 0 | True |
| Propane/Pentane 60/40 | TuRBO-1 (adapted) | 20260910 | 2.52535586 | 135 | 1 | True |
| Propane/Pentane 60/40 | TuRBO-1 (adapted) | 20260911 | 2.53678728 | 295 | 0 | True |
| Propane/Pentane 60/40 | BO (GP-EI) | 20260909 | 2.53643019 | 307 | 1 | True |
| Propane/Pentane 60/40 | BO (GP-EI) | 20260910 | 2.53706966 | 219 | 2 | True |
| Propane/Pentane 60/40 | BO (GP-EI) | 20260911 | 2.53617730 | 182 | 1 | True |
| Propane/Isobutane 30/70 | Multistart COBYQA | 20260909 | 2.54223065 | 352 | 3 | True |
| Propane/Isobutane 30/70 | Multistart COBYQA | 20260910 | 2.53861875 | 288 | 6 | True |
| Propane/Isobutane 30/70 | Multistart COBYQA | 20260911 | 2.52648417 | 278 | 6 | True |
| Propane/Isobutane 30/70 | TuRBO-1 (adapted) | 20260909 | 2.54133974 | 86 | 4 | True |
| Propane/Isobutane 30/70 | TuRBO-1 (adapted) | 20260910 | 2.53981660 | 397 | 2 | True |
| Propane/Isobutane 30/70 | TuRBO-1 (adapted) | 20260911 | 2.53852239 | 193 | 2 | True |
| Propane/Isobutane 30/70 | BO (GP-EI) | 20260909 | 2.54172078 | 324 | 4 | True |
| Propane/Isobutane 30/70 | BO (GP-EI) | 20260910 | 2.54129111 | 203 | 3 | True |
| Propane/Isobutane 30/70 | BO (GP-EI) | 20260911 | 2.54195437 | 379 | 7 | True |

## Limits of the comparison

Three seeds show variability but do not establish statistical superiority. The TuRBO run uses an explicitly documented NumPy/SciPy adaptation: constrained candidate rejection, a failed-region heuristic and L-BFGS-B GP fitting. It is not a run of the authors' GPyTorch package. BO retains its existing initialisation, boundary candidates and GP settings. COBYQA uses a finite penalty for failed solves. Results compare these complete implementations at the fixed budget.

## Reproduce

From the repository root, using the locked environment (SciPy 1.18.1):

```bash
uv sync --locked --group dev --group benchmark
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 .venv/bin/python -u main_RHPtrust_region_benchmark.py --n-cpu 9
MPLCONFIGDIR=/tmp/heatpumplab-trust-region-mpl .venv/bin/python main_RHPtrust_region_report.py
```

The runner resumes matching checkpoints and skips completed trials. The report generator reconstructs the numerical tables and figures; this interpretation records the completed September 2026 experiment.
