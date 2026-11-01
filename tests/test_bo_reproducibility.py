"""Same-seed BO runs must propose bit-identical designs, whatever the state of the process heap."""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import torch

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._stable_solves import (
    AddressIndependentSolves,
)

REPO_ROOT = Path(__file__).resolve().parents[1]
PROPOSALS = 170

# A closed-loop search on a synthetic objective with constraint violations and numerical failures,
# as in production. Tensors of random size, seeded per process, leave the two heaps different.
SEARCH = """
import sys
import numpy as np
import torch
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._bayesian import (
    BayesianDesignSearch, DesignDomain,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
    OuterDesignOptimisationConfig,
)

heap_seed, proposals = int(sys.argv[1]), int(sys.argv[2])
heap, held = np.random.default_rng(heap_seed), []
config = OuterDesignOptimisationConfig(seed=20261031, preserve_incumbent=False)
search = BayesianDesignSearch(DesignDomain([1 / 3] * 3, 123.0, config), config)


def evaluate(x):
    lift = (x[0] + x[1]) / 100
    gas = 0.6 + 0.8 * x[5] + 0.5 * lift - x[3]
    if abs(np.sin(37 * x[3] + 11 * x[4])) < 0.03:
        return 0.0, False, None
    margins = dict(c1_discharge_superheat=10 * gas, c1_suction_superheat=20.0, c1_pressure_ratio=2.0)
    sec = (
        0.65 - 0.08 * np.tanh(3 * lift) + 0.3 * (x[2] - 0.5) ** 2
        + 0.2 * ((x[3] - 0.2) ** 2 + (x[4] - 0.1) ** 2) + 0.01 * np.sin(9 * x[0] / 100 + 5 * x[5])
    )
    return (sec, True, margins) if gas >= 0 else (0.0, False, margins)


for _ in range(proposals):
    for _ in range(int(heap.integers(5, 200))):
        held.append(torch.empty(int(heap.integers(1, 3000)), dtype=torch.double))
        if len(held) > 300:
            held.pop(int(heap.integers(0, len(held))))
    x = search.ask()
    print(*(float(value).hex() for value in x))
    search.tell(x, *evaluate(x))
"""


def test_two_processes_propose_identical_designs_for_150_plus_steps():
    env = dict(
        os.environ,
        PYTHONPATH=str(REPO_ROOT),
        OMP_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        MKL_NUM_THREADS="1",
    )
    runs = [
        subprocess.Popen(
            [sys.executable, "-W", "ignore", "-c", SEARCH, str(heap_seed), str(PROPOSALS)],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        for heap_seed in (1, 2)
    ]
    outputs = [run.communicate(timeout=1800) for run in runs]
    for run, (_, stderr) in zip(runs, outputs):
        assert run.returncode == 0, stderr
    first, second = (stdout.splitlines() for stdout, _ in outputs)
    assert len(first) == PROPOSALS
    differing = [i for i, (a, b) in enumerate(zip(first, second)) if a != b]
    assert not differing, f"proposals differ from evaluation {differing[0] + 1}: {differing}"


def _triangular_system(rows, columns):
    rng = np.random.default_rng(rows)
    matrix = rng.standard_normal((rows, rows))
    factor = torch.linalg.cholesky(torch.tensor(matrix @ matrix.T / rows + np.eye(rows)))
    return factor, torch.tensor(rng.standard_normal((rows, columns)))


def test_padded_solves_do_not_depend_on_buffer_addresses():
    heap = np.random.default_rng(0)
    for rows in (65, 80, 128, 200):
        for columns in (1, 3):
            factor, rhs = _triangular_system(rows, columns)
            results = {"lower": set(), "upper": set(), "cholesky": set()}
            held = []
            with AddressIndependentSolves():
                for _ in range(100):
                    held.append(torch.empty(int(heap.integers(1, 800)), dtype=torch.double))
                    if len(held) > 50:
                        held.pop(int(heap.integers(0, 50)))
                    b = rhs.clone()
                    results["lower"].add(
                        torch.linalg.solve_triangular(factor, b, upper=False).numpy().tobytes()
                    )
                    results["upper"].add(
                        torch.linalg.solve_triangular(factor.mT, b, upper=True).numpy().tobytes()
                    )
                    results["cholesky"].add(torch.cholesky_solve(b, factor).numpy().tobytes())
            assert {name: len(found) for name, found in results.items()} == dict.fromkeys(results, 1)


def test_padded_solves_match_plain_solves_in_shape_value_and_gradient():
    factor, rhs = _triangular_system(70, 1)
    for columns in (rhs, torch.cat([rhs, 2 * rhs, 3 * rhs], -1), torch.cat([rhs, 2 * rhs], -1)):
        plain_b = columns.clone().requires_grad_()
        plain = torch.cholesky_solve(plain_b, factor)
        plain.square().sum().backward()
        padded_b = columns.clone().requires_grad_()
        with AddressIndependentSolves():
            padded = torch.cholesky_solve(padded_b, factor)
            padded.square().sum().backward()
            triangular = torch.linalg.solve_triangular(factor, columns, upper=False)
            method = columns.cholesky_solve(factor)
        assert padded.shape == columns.shape == triangular.shape == method.shape
        assert torch.allclose(padded, plain, rtol=1e-10)
        assert torch.allclose(padded_b.grad, plain_b.grad, rtol=1e-9)
        assert torch.allclose(factor @ triangular, columns, rtol=1e-9)
