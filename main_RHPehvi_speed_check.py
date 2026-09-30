"""Check that BoTorch's fused qLogEHVI C++ kernel loads, and time it against the pure-Python path.

Uses a problem the size of the fixed-split search: 5-D inputs, 256 training points, a few dozen front
points, 64 MC samples, candidate pools scored in chunks as in rhp_mobo_search.
"""
import time

import numpy as np
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._torch_env import ensure_ninja_on_path

ensure_ninja_on_path()  # before BoTorch first loads the fused kernel
import torch
from botorch.acquisition.multi_objective import logei
from botorch.acquisition.multi_objective.logei import qLogExpectedHypervolumeImprovement
from botorch.models import SingleTaskGP
from botorch.sampling.normal import SobolQMCNormalSampler
from botorch.utils.multi_objective.box_decompositions.non_dominated import NondominatedPartitioning


def build(seed=0):
    rng = np.random.default_rng(seed)
    X = torch.tensor(rng.random((256, 5)))
    a = rng.random(256)
    Y = torch.tensor(np.column_stack((1+a+.05*rng.standard_normal(256), 2-a+.05*rng.standard_normal(256))))
    model = SingleTaskGP(X, Y, train_Yvar=torch.full_like(Y, 1e-6)).double().eval()
    ref = torch.zeros(2, dtype=torch.double)
    return model, NondominatedPartitioning(ref_point=ref, Y=Y), ref


def run(model, partition, ref, repeats=3):
    acq = qLogExpectedHypervolumeImprovement(model=model, ref_point=ref.tolist(), partitioning=partition,
                                             sampler=SobolQMCNormalSampler(torch.Size([64]), seed=1))
    pool = torch.rand(2048, 1, 5, dtype=torch.double, generator=torch.Generator().manual_seed(2))
    best = np.inf
    for _ in range(repeats):
        started = time.perf_counter()
        with torch.no_grad():
            values = torch.cat([acq(chunk) for chunk in pool.split(64)])
        X = pool[:64].clone().requires_grad_(True)
        acq(X).sum().backward()
        best = min(best, time.perf_counter()-started)
    return best, values


def main():
    torch.set_num_threads(1)
    model, partition, ref = build()
    print(f'front points {len(partition.pareto_Y)}, hypercells {len(partition.get_hypercell_bounds()[0])}')
    fused_time, fused_values = run(model, partition, ref)
    loaded = logei._C is not None
    print('fused kernel loaded:', loaded)
    saved, logei._C = logei._C, None
    python_time, python_values = run(model, partition, ref)
    logei._C = saved
    print(f'pure Python {python_time:.3f} s, fused {fused_time:.3f} s, speed-up {python_time/fused_time:.2f}x')
    print('max |difference| in log-EHVI:', float((fused_values-python_values).abs().max()))
    if not loaded:
        print('Fused kernel not active: check that .venv/bin/ninja exists (uv sync) and a C++ compiler works '
              '(xcode-select --install). The first successful load compiles for about a minute.')


if __name__ == '__main__':
    main()
