"""Phase 3 offline benchmark: which surrogate predicts COP and VHC best near the front?

Uses the saved v22 Phase 2 runs (no REFPROP calls). For each run, successful designs
are split 80/20; each model is trained on the 80% and scored on the held-out 20%,
separately for all held-out designs, those near the observed front, and those whose
neighbourhood spans the gas-cooler pinch-regime switch. A local variant refits on the
128 nearest training designs, as the local-128 search does.

Models
  fixed_noise   current search model: BoTorch default GP, fixed noise 1e-8 (raw units)
  learned_noise same GP with a learned noise term (floor 1e-8 on the standardised scale)
  regime_split  one fixed-noise GP per pinch regime + a GP classifier for the regime;
                prediction = classifier-weighted mixture

Run: .venv/bin/python main_RHPsurrogate_benchmark.py --n-cpu 6
"""
import argparse
import glob
import json
import time
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

OUT = Path("results/RHPcomposition/phase3_surrogate_benchmark")
FLOOR = 0.025
MODELS = ("fixed_noise", "learned_noise", "regime_split")


def load_run(path):
    r = json.loads(Path(path).read_text())
    base = r["metadata"]["base_case"]
    hot_T = float(base["T_air_tar"]) + float(base["dT_GC_min_tar"]) - 0.05
    good = [x for x in r["records"] if x["success"] and x.get("COP") and x.get("VHC_MJ_per_m3")]
    X = np.array([[g["x"][0] / 100, g["x"][1] / 100, g["x"][3]] for g in good])
    Y = np.array([[g["COP"], g["VHC_MJ_per_m3"]] for g in good])
    if all(g.get("gc_pinch_active") for g in good):
        hot = np.array([g["gc_pinch_active"] == "gc1_hot" for g in good])
        source = "report"
    else:
        hot = np.array([abs(g["T_comp_discharge"] - hot_T) < 0.01 for g in good])
        source = "discharge_temperature"
    return X, Y, hot, source


def near_front(Y, tol=0.002):
    """Designs within tol (relative COP) of the best COP at equal or higher VHC."""
    order = np.argsort(-Y[:, 1])
    envelope = np.empty(len(Y))
    best = -np.inf
    for i in order:
        best = max(best, Y[i, 0])
        envelope[i] = best
    return Y[:, 0] >= (1 - tol) * envelope


def fit_gp(X, y, noise):
    import torch
    from botorch.fit import fit_gpytorch_mll
    from botorch.models import SingleTaskGP
    from gpytorch.constraints import GreaterThan
    from gpytorch.likelihoods import GaussianLikelihood
    from gpytorch.mlls import ExactMarginalLogLikelihood
    Xt = torch.tensor(X, dtype=torch.double)
    yt = torch.tensor(y, dtype=torch.double).reshape(-1, 1)
    if noise == "fixed":
        model = SingleTaskGP(Xt, yt, train_Yvar=torch.full_like(yt, 1e-8))
    else:
        model = SingleTaskGP(Xt, yt, likelihood=GaussianLikelihood(noise_constraint=GreaterThan(1e-8)))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            fit_gpytorch_mll(ExactMarginalLogLikelihood(model.likelihood, model), max_attempts=2,
                             optimizer_kwargs={"options": {"maxiter": 200}})
        except Exception:
            pass
    model.eval()
    return model


def predict(model, X, observation_noise=False):
    import torch
    with torch.no_grad():
        post = model.posterior(torch.tensor(X, dtype=torch.double), observation_noise=observation_noise)
        return post.mean.squeeze(-1).numpy(), post.variance.squeeze(-1).clamp_min(1e-18).numpy()


def lengthscales(model):
    kernel = getattr(model.covar_module, "base_kernel", model.covar_module)
    return kernel.lengthscale.detach().numpy().ravel().tolist()


def model_predict(name, Xtr, ytr, htr, Xte):
    """Return mean, variance, length scales (list of lists) for one objective."""
    if name == "fixed_noise":
        m = fit_gp(Xtr, ytr, "fixed")
        mu, var = predict(m, Xte)
        return mu, var, [lengthscales(m)]
    if name == "learned_noise":
        m = fit_gp(Xtr, ytr, "learned")
        mu, var = predict(m, Xte, observation_noise=True)
        return mu, var, [lengthscales(m)]
    # regime_split
    if htr.sum() < 8 or (~htr).sum() < 8:
        return model_predict("fixed_noise", Xtr, ytr, htr, Xte)
    from scipy.stats import norm
    a = fit_gp(Xtr[htr], ytr[htr], "fixed")
    b = fit_gp(Xtr[~htr], ytr[~htr], "fixed")
    c = fit_gp(Xtr, np.where(htr, 1.0, -1.0), "learned")
    ma, va = predict(a, Xte)
    mb, vb = predict(b, Xte)
    mc, vc = predict(c, Xte, observation_noise=True)
    p = norm.cdf(mc / np.sqrt(vc))
    mu = p * ma + (1 - p) * mb
    var = p * (va + ma**2) + (1 - p) * (vb + mb**2) - mu**2
    return mu, np.maximum(var, 1e-18), [lengthscales(a), lengthscales(b)]


def score(y, mu, var):
    rel = np.abs(mu - y) / np.abs(y)
    nlpd = 0.5 * np.log(2 * np.pi * var) + 0.5 * (y - mu) ** 2 / var
    return dict(n=int(len(y)), median_rel_err=float(np.median(rel)), p90_rel_err=float(np.quantile(rel, .9)),
                rmse_rel=float(np.sqrt(np.mean(rel**2))), mean_nlpd=float(np.mean(nlpd)))


def run_one(path, seed=20260928, local_points=15):
    import torch
    torch.set_num_threads(1)
    started = time.perf_counter()
    X, Y, hot, source = load_run(path)
    rng = np.random.default_rng(seed)
    test = rng.random(len(X)) < 0.2
    train = ~test
    front = near_front(Y)
    # Switch neighbourhood: any of the 10 nearest designs in the other regime.
    D = np.sqrt(((X[:, None] - X[None]) ** 2).sum(-1))
    nn = np.argsort(D, axis=1)[:, 1:11]
    switch = np.array([len(set(hot[[i, *nn[i]]])) > 1 for i in range(len(X))])
    groups = dict(all=test, near_front=test & front, near_switch=test & switch)
    out = dict(run=str(path), n=int(len(X)), n_test=int(test.sum()), regime_source=source,
               hot_fraction=float(hot.mean()), global_={}, local={})
    for name in MODELS:
        out["global_"][name] = {}
        for j, obj in enumerate(("COP", "VHC")):
            mu, var, ls = model_predict(name, X[train], Y[train, j], hot[train], X[test])
            full_mu = np.full(len(X), np.nan); full_var = np.full(len(X), np.nan)
            full_mu[test], full_var[test] = mu, var
            res = {g: score(Y[m, j], full_mu[m], full_var[m]) for g, m in groups.items() if m.sum()}
            flat = np.ravel(ls)
            res["lengthscales"] = ls
            res["frac_at_floor"] = float(np.mean(np.abs(flat - FLOOR) < 1e-4))
            out["global_"][name][obj] = res
    # Local-128 variant on near-front held-out designs.
    idx = np.flatnonzero(test & front)
    idx = rng.choice(idx, min(local_points, len(idx)), replace=False) if len(idx) else idx
    tr = np.flatnonzero(train)
    for name in MODELS:
        out["local"][name] = {}
        for j, obj in enumerate(("COP", "VHC")):
            ys, mus, vars_, floors = [], [], [], []
            for i in idx:
                near = tr[np.argsort(D[i, tr])[:128]]
                mu, var, ls = model_predict(name, X[near], Y[near, j], hot[near], X[[i]])
                ys.append(Y[i, j]); mus.append(mu[0]); vars_.append(var[0])
                floors.append(np.mean(np.abs(np.ravel(ls) - FLOOR) < 1e-4))
            if ys:
                res = score(np.array(ys), np.array(mus), np.array(vars_))
                res["frac_at_floor"] = float(np.mean(floors))
                out["local"][name][obj] = res
    out["seconds"] = time.perf_counter() - started
    return out


def report(results):
    lines = ["# Phase 3 surrogate benchmark", "",
             "Median relative prediction error on held-out designs (lower is better), pooled over runs "
             "as the median of per-run values. NLPD = mean negative log predictive density (lower is better).", ""]
    for scope in ("global_", "local"):
        lines += [f"## {'Global GP (all training data)' if scope == 'global_' else 'Local GP (128 nearest)'}", "",
                  "| Model | Objective | all | near front | near switch | near-front NLPD | length scales at floor |",
                  "|---|---|---:|---:|---:|---:|---:|"]
        for name in MODELS:
            for obj in ("COP", "VHC"):
                rows = [r[scope][name][obj] for r in results if obj in r[scope].get(name, {})]
                if not rows:
                    continue
                med = lambda g, k="median_rel_err": np.median([x[g][k] for x in rows if g in x]) if any(g in x for x in rows) else np.nan
                lines.append(f"| {name} | {obj} | {med('all'):.2e} | {med('near_front'):.2e} | "
                             f"{med('near_switch'):.2e} | {med('near_front', 'mean_nlpd'):.2f} | "
                             f"{np.mean([x['frac_at_floor'] for x in rows]):.0%} |" if scope == "global_" else
                             f"| {name} | {obj} | - | {np.median([x['median_rel_err'] for x in rows]):.2e} | - | "
                             f"{np.median([x['mean_nlpd'] for x in rows]):.2f} | {np.mean([x['frac_at_floor'] for x in rows]):.0%} |")
        lines.append("")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-cpu", type=int, default=6)
    parser.add_argument("--runs", default="results/RHPcomposition/phase2_superheat_*K/*/result.json")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    paths = sorted(glob.glob(args.runs))
    print(f"{len(paths)} runs", flush=True)
    with ProcessPoolExecutor(args.n_cpu) as pool:
        results = []
        for r in pool.map(run_one, paths):
            print(f"done {r['run']} in {r['seconds']:.0f}s", flush=True)
            results.append(r)
    (OUT / "results.json").write_text(json.dumps(results, indent=1))
    (OUT / "summary.md").write_text(report(results))
    print(report(results))


if __name__ == "__main__":
    main()
