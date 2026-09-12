# Isolated optimisation experiments

These directories preserve the source used for the recorded v17 experiments.
They intentionally contain older copies of the thermodynamic model so changes
to the production v18 optimiser do not change a running or resumed experiment.
Do not edit a frozen runtime while its supervisor is running: manifests and
checkpoint metadata verify source hashes.

| Runtime | Purpose |
| --- | --- |
| `mobo_v17_clean_runtime` | Original four propane-family COP–VHC_total searches, 4,000 candidates per family. |
| `mobo_v17_local128_runtime` | Continuation with 128 nearby successful observations per GP and six local regions refreshed every 48 candidates. |
| `mobo_v17_nine_screen_runtime` | Nine additional binary-fluid screens, 800 global candidates each, with four concurrent searches. |

The local-subset comparison entry points are `main_RHPlocal_subset_trial.py` and
`main_RHPlocal_subset_report.py` in the repository root. Their baseline copies
come from the paused four-family checkpoints; the trial does not import its
observations into those searches.

Source manifests, source archives, mixture setup diagnostics and the small
integration-probe record are retained for provenance. Candidate checkpoints,
live logs, figures and reports remain under the ignored `results/` directory.
The original and local-continuation supervisors require those existing local
checkpoints and manifests; the nine-family supervisor can start a fresh screen.
`reference_manifest.json` identifies the completed reference checkpoints used
for comparison. A full nine-family comparison report also requires those
reference checkpoints under `results/RHPcomposition/mobo_v17_4000/`.

Run a supervisor from the repository environment, for example:

```sh
.venv/bin/python experiments/mobo_v17_nine_screen_runtime/supervise.py
```

Each supervisor holds a run lock, records child PIDs, checkpoints completed
evaluations, and validates final frontier points in fresh native workers.
REFPROP and its fluid files must be installed locally. Cyclopentane mixtures
use alternative mixing rules in the tested REFPROP database; their predictions
are exploratory even when independent replays match.

The default `pytest` command runs the production `tests/` directory. Run tests
inside an isolated runtime in a separate Python process to avoid collisions
between its frozen modules and production modules:

```sh
cd experiments/mobo_v17_nine_screen_runtime
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ../../.venv/bin/python -m pytest test_screen.py
```

Each runtime's `PROTOCOL.md`, where present, describes its budget, controls and
validation rules. The fixed 800- or 4,000-evaluation budgets are stopping rules,
not demonstrations of global convergence or chemical stability.
