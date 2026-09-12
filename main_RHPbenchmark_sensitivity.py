"""Check each method's final best design with independent inner-pressure guesses."""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path

from main_RHPoptimiser_benchmark import serialise
from main_RHPretest_report import load_trials
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective


def execute(job):
    objective = _OuterDesignObjective(refrigerant=job['refrigerant'],
                                     base_case=dict(job['base_case'], P_H_est=job['pressure_guess_MPa']),
                                     max_evals=1, exact_evaluation_timeout_seconds=180,
                                     exact_evaluation_timeout_retries=1, exact_evaluation_worker_max_tasks=1)
    try:
        record = serialise(objective.evaluate(job['x']))
    finally:
        objective.close()
    return dict(job=job, record=record)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', default='results/RHPcomposition/optimiser_v13_400')
    args = parser.parse_args()
    output = Path(args.output)
    trials = load_trials(output)
    selected = {}
    for t in trials:
        m = t['metadata']
        key = m['refrigerant'], m['method']
        if key not in selected or t['replay']['COP'] > selected[key]['replay']['COP']:
            selected[key] = t
    directory = output / 'sensitivity'
    directory.mkdir(exist_ok=True)
    jobs = []
    for t in selected.values():
        m, r = t['metadata'], t['replay']
        for guess in (3.0, 7.0):
            jobs.append(dict(refrigerant=m['refrigerant'], method=m['method'], base_case=m['base_case'],
                             seed=m['seed'], x=r['x'], default_COP=r['COP'], pressure_guess_MPa=guess,
                             model_source_sha256=m['model_source_sha256']))
    rows = []
    pending = {}
    with ProcessPoolExecutor(max_workers=3) as pool:
        for job in jobs:
            key = hashlib.sha256(json.dumps(job, sort_keys=True).encode()).hexdigest()[:16]
            path = directory / (key + '.json')
            if path.exists():
                saved = json.loads(path.read_text())
                assert saved['job'] == job
                rows.append(saved)
            else:
                pending[pool.submit(execute, job)] = path
        for f in as_completed(pending):
            row = f.result()
            pending[f].write_text(json.dumps(row, indent=2, allow_nan=False))
            rows.append(row)
            print('Pressure-guess checks', len(rows), '/', len(jobs), flush=True)
    successful = [r for r in rows if r['record']['success']]
    summary = dict(evaluations=len(rows), successes=len(successful),
                   max_absolute_COP_change_percent=max(abs(100*(r['record']['COP']/r['job']['default_COP']-1)) for r in successful),
                   model_source_sha256=jobs[0]['model_source_sha256'],
                   protocol='Best of three seeds per stochastic method; 3 and 7 MPa initial guesses; diagnostics excluded from searches')
    (directory / 'summary.json').write_text(json.dumps(summary, indent=2))
    (directory / 'results.json').write_text(json.dumps(rows, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
