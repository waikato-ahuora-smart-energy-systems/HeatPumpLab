"""Run one bounded, sequential local-Pareto GP subset comparison arm."""
import argparse
import hashlib
import json
from pathlib import Path
import resource
import sys
import time

import numpy as np
from rhp_local_pareto_search import LocalParetoSearch, front_indices, hypervolume
from main_RHPmobo import PARTNERS, metadata_for
from main_RHPoptimiser_benchmark import serialise
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT/'results/RHPcomposition/local_pareto_subset_trial'


def save(path, data):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(data, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def run(partner, subset, budget):
    subset_size = None if subset == 'full' else int(subset)
    baseline_path = OUTPUT/'baseline'/f'Propane_{partner}.json'
    baseline = json.loads(baseline_path.read_text())
    # Confirm exact thermodynamic/search source and installed libraries match the saved experiment.
    assert baseline['metadata'] == metadata_for(partner, 4000, 20260910)
    options = dict(refrigerant=baseline['metadata']['refrigerant'], base_case=baseline['metadata']['base_case'],
                   exact_evaluation_timeout_seconds=180, exact_evaluation_timeout_retries=0,
                   exact_evaluation_worker_max_tasks=1, results_dir=str(OUTPUT/'oracle'))
    path = OUTPUT/f'{partner}_{subset}.json'
    metadata = dict(partner=partner, subset=subset, budget=budget, seed=20260911,
                    baseline_sha256=hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
                    source_sha256=hashlib.sha256((ROOT/'rhp_local_pareto_search.py').read_bytes()+Path(__file__).read_bytes()).hexdigest())
    data = json.loads(path.read_text()) if path.exists() else dict(metadata=metadata, records=[], complete=False)
    assert data['metadata'] == metadata
    if data['complete']:
        print('Already complete', path, flush=True)
        return
    search = LocalParetoSearch(baseline, subset_size)
    for record in data['records']:
        search.tell(record['x'],record['COP'],record['VHC_MJ_per_m3'],record['success'])
    search.states = {('global' if k == 'global' else int(k)):v for k,v in data.get('model_states',{}).items()}
    data['anchor_evaluations'] = [i+1 for i in search.anchor_indices]
    objective = _OuterDesignObjective(max_evals=max(1,budget-len(data['records'])), **options)
    search.max_lift = objective.max_recuperator_lift
    initial_values = [[r['COP'],r['VHC_MJ_per_m3']] if r['success'] and not r.get('validation_rejected') else [0.,0.] for r in baseline['records']]
    data['baseline_hypervolume'] = hypervolume(initial_values)
    save(path,data)
    try:
        while len(data['records']) < budget:
            proposal = search.ask()
            diag = search.last_diagnostic
            started = time.perf_counter()
            evaluated = objective.evaluate(proposal)
            if 'worker_failure' in evaluated.constraint_errors:
                raise RuntimeError(evaluated.message)
            assert np.allclose(evaluated.x,proposal,rtol=0,atol=1e-7)
            record = serialise(evaluated)
            record.update(proposal_seconds=diag['proposal_seconds'], evaluation_wall_seconds=time.perf_counter()-started,
                          search_diagnostic=diag)
            search.tell(record['x'],record['COP'],record['VHC_MJ_per_m3'],record['success'])
            data['records'].append(record)
            data['model_states'] = search.states
            data['hypervolume'] = hypervolume(search.y)
            data['peak_rss_MiB'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform == 'darwin' else 1024)
            save(path,data)
            print(partner,subset,len(data['records']),f"proposal {diag['proposal_seconds']:.2f}s",
                  'success',record['success'],f"HV gain {100*(data['hypervolume']/data['baseline_hypervolume']-1):.5f}%",
                  f"peak RSS {data['peak_rss_MiB']:.0f} MiB",flush=True)
    finally:
        objective.close()
    data.setdefault('replays',{})
    replay = _OuterDesignObjective(max_evals=budget,**options)
    try:
        while True:
            values=initial_values+[[r['COP'],r['VHC_MJ_per_m3']] if r['success'] and not r.get('validation_rejected') else [0.,0.] for r in data['records']]
            pending=[i-len(initial_values) for i in front_indices(values) if i>=len(initial_values)
                     and str(i-len(initial_values)) not in data['replays']]
            if not pending:
                break
            for i in pending:
                original=data['records'][i]
                checked=serialise(replay.evaluate(original['x']))
                matches=checked['success'] and all(np.isclose(checked[k],original[k],rtol=1e-6,atol=1e-9)
                                                 for k in ('COP','VHC_MJ_per_m3','SEC','Pmax_MPa'))
                data['replays'][str(i)]=dict(result=checked,matches=bool(matches))
                if not matches:
                    original['validation_rejected']=True
                save(path,data)
    finally:
        replay.close()
    values=initial_values+[[r['COP'],r['VHC_MJ_per_m3']] if r['success'] and not r.get('validation_rejected') else [0.,0.] for r in data['records']]
    data['hypervolume']=hypervolume(values)
    data['new_front_indices']=[i-len(initial_values) for i in front_indices(values) if i>=len(initial_values)]
    data['complete']=True
    save(path,data)
    print('Completed',partner,subset,flush=True)


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partner',choices=PARTNERS,required=True)
    parser.add_argument('--subset',choices=('128','256','full'),required=True)
    parser.add_argument('--budget',type=int,default=12)
    args=parser.parse_args()
    run(args.partner,args.subset,args.budget)
