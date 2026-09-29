"""Resumable 25,000-evaluation joint temperature/COP/VHC BO, per propane pair."""
from __future__ import annotations

import argparse
import fcntl
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import time

import numpy as np

from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._objective import _OuterDesignObjective
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
    DEFAULT_BASE_CASE, DEFAULT_MAX_STAGE_PRESSURE_RATIO)
# Import the BO implementation only in parent-side functions. Spawned native
# workers import this entry point and should not load Torch/BoTorch needlessly.

ROOT = Path(__file__).resolve().parent
PARTNERS = ('Butane', 'Isobutane', 'Pentane', 'Isopentane')
RUN_VERSION = 'joint-temperature-cop-vhc-local128-v19'


def save(path, data):
    path = Path(path)
    temporary = path.with_suffix(path.suffix+'.tmp')
    with temporary.open('w') as stream:
        json.dump(data, stream, indent=2, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def read_journal(path):
    """Recover only a torn final append; malformed completed lines are errors."""
    path = Path(path)
    if not path.exists():
        return []
    data = path.read_bytes()
    if data and not data.endswith(b'\n'):
        end = data.rfind(b'\n')+1
        with path.open('r+b') as stream:
            stream.truncate(end)
        data = data[:end]
    return [json.loads(line) for line in data.splitlines()]


def append(path, record):
    with Path(path).open('a') as stream:
        stream.write(json.dumps(record, allow_nan=False)+'\n')
        stream.flush()
        os.fsync(stream.fileno())


def finite(value):
    if isinstance(value, (np.integer, np.floating)):
        value = value.item()
    return float(value) if value is not None and np.isfinite(value) else None


def exact_evaluation(partner, x, temperature, output, *, first='Propane', min_discharge_superheat_K=0.0):
    from rhp_temperature_mobo import usable
    objective = _OuterDesignObjective(max_evals=1,
        refrigerant=f'{first};{partner}|0.5;0.5|mass',
        base_case=dict(DEFAULT_BASE_CASE, T_air_tar=float(temperature)),
        exact_evaluation_timeout_seconds=180, exact_evaluation_timeout_retries=0,
        exact_evaluation_worker_max_tasks=1, results_dir=str(output/'oracle'),
        min_discharge_superheat_K=min_discharge_superheat_K)
    started = time.perf_counter()
    try:
        record = objective.evaluate(x)
        if not np.allclose(record.x, x, rtol=0, atol=1e-7):
            raise ValueError('Joint candidate was changed by objective repair')
        report = record.report
        item = dict(x=list(record.x), T_air_tar=float(temperature), success=record.success,
                    COP=finite(record.cop), SEC=finite(record.sec),
                    VHC_MJ_per_m3=finite(report.get('VHC', 0)/1000),
                    Pmax_MPa=finite(record.pmax/1000) if record.pmax is not None else None,
                    message=record.message, constraint_errors=record.constraint_errors,
                    elapsed_seconds=record.elapsed_seconds,
                    evaluation_wall_seconds=time.perf_counter()-started,
                    solve_failure_kind=report.get('solve_failure_kind'),
                    solve_error_stage=report.get('solve_error_stage'),
                    property_failure_context=report.get('property_failure_context'),
                    numerical_retries=report.get('numerical_retries', []),
                    constraint_margins=report.get('constraint_margins'),
                    gc_pinch_active=report.get('gc_pinch_active'),
                    gc_approaches=report.get('gc_approaches'))
        for field in ('T_comp_discharge', 'T_gc_out', 'Tref_evap_outlet', 'Comp_1_T_out',
                      'Comp_2_T_out', 'PR_1', 'PR_2', 'PR_3', 'Overall_PR', 'dT_gc_out_min_offset'):
            item[field] = finite(report.get(field))
        for stage in range(1, 4):
            for end in ('in', 'out'):
                field = f'Comp_{stage}_q_refprop_{end}'
                item[field] = finite(report.get(field))
                field = f'Comp_{stage}_phase_{end}'
                item[field] = report.get(field)
        if item['success']:
            if not usable(item) or not np.isclose(report['T_air_tar'], temperature, atol=1e-6, rtol=0):
                raise ValueError('Successful solve lacks valid objectives or requested temperature')
        return item
    finally:
        objective.close()


def metadata(partner, budget, seed):
    from rhp_temperature_mobo import MILESTONES, TEMPERATURES
    paths = sorted((ROOT/'HeatPumpLab').rglob('*.py'))+[ROOT/'rhp_temperature_mobo.py', ROOT/'rhp_temperature_report.py', Path(__file__).resolve()]
    return dict(run_version=RUN_VERSION, partner=partner, budget=budget, seed=seed,
        max_stage_pressure_ratio=DEFAULT_MAX_STAGE_PRESSURE_RATIO,
        compressor_gas_only=True,
        milestones=sorted(set([m for m in MILESTONES if m <= budget]+[budget])),
        temperatures_C=TEMPERATURES.tolist(), base_case=DEFAULT_BASE_CASE,
        objectives=['COP', 'VHC_total_MJ_per_m3', 'T_air_tar_C'],
        reference_point_physical=[0., 0., 95.], scales=[4., 10., 125.],
        subset_size=128, initial_evaluations=200, numerical_threads=1,
        known_objective='Temperature is exact; GP outputs are COP and VHC only',
        search='qLogEHVI q=1; all grid levels enumerated; global maximin every eighth adaptive step',
        evaluation_timeout_seconds=180, timeout_retries=0, native_tasks_per_worker=1,
        budget_accounting='Includes initial, BO, exploratory and failed candidates; final frontier replays counted separately',
        source_sha256=hashlib.sha256(b''.join(p.read_bytes() for p in paths)).hexdigest(),
        versions={name: version(name) for name in ('numpy', 'scipy', 'botorch', 'torch', 'gpytorch', 'ctREFPROP')})


def milestone_snapshot(records, previous=None):
    from rhp_temperature_mobo import TEMPERATURES, front_indices, hypervolume_3d, objectives, usable
    indices = front_indices(records)
    frontier = [dict(evaluation=i+1, **records[i]) for i in indices]
    hv = hypervolume_3d([objectives(records[i]) for i in indices])
    earlier = previous['hypervolume_normalized'] if previous else None
    per_temperature = []
    for temperature in TEMPERATURES:
        at = [r for r in records if r['T_air_tar'] == temperature]
        good = [r for r in at if usable(r)]
        per_temperature.append(dict(temperature_C=float(temperature), evaluations=len(at),
            successful=len(good), failures=sum(not r['success'] for r in at),
            frontier_points=int(sum(r['T_air_tar'] == temperature for r in frontier)),
            best_COP=max((r['COP'] for r in good), default=None),
            best_VHC_MJ_per_m3=max((r['VHC_MJ_per_m3'] for r in good), default=None)))
    return dict(evaluations=len(records), label='observed frontier, prior to final fresh-process validation',
        hypervolume_normalized=hv, hypervolume_physical=hv*4*10*125,
        hypervolume_gain_since_previous=None if earlier is None else hv-earlier,
        hypervolume_gain_percent_since_previous=None if not earlier else 100*(hv/earlier-1),
        frontier_points=len(frontier), frontier=frontier,
        failures=sum(not r['success'] for r in records),
        proposal_seconds=sum(r.get('proposal_seconds', 0) for r in records),
        evaluation_wall_seconds=sum(r.get('evaluation_wall_seconds', 0) for r in records),
        per_temperature=per_temperature)


def write_milestones(output, records, marks):
    previous, baseline, summaries = None, None, []
    for mark in marks:
        if len(records) < mark:
            break
        path = output/f'milestone_{mark:05d}.json'
        if path.exists():
            snapshot = json.loads(path.read_text())
        else:
            snapshot = milestone_snapshot(records[:mark], previous)
            save(path, snapshot)
        if not path.with_suffix('.html').exists():
            from rhp_temperature_report import render_front
            render_front(snapshot, path.with_suffix('.html'))
        if baseline is None:
            baseline = snapshot['hypervolume_normalized']
        summary = {k: v for k, v in snapshot.items() if k != 'frontier'}
        summary['hypervolume_gain_percent_since_first_mark'] = (
            100*(snapshot['hypervolume_normalized']/baseline-1) if baseline else None)
        summaries.append(summary)
        previous = snapshot
    save(output/'milestone_comparison.json', summaries)


def verify_final(partner, output, records):
    from rhp_temperature_mobo import front_indices
    path = output/'frontier_replays.jsonl'
    replays = {r['evaluation']-1: r for r in read_journal(path)}
    validated = [dict(r) for r in records]
    for index, replay in replays.items():
        if not replay['matches']:
            validated[index]['validation_rejected'] = True
    while True:
        indices = front_indices(validated)
        pending = [i for i in indices if i not in replays]
        if not pending:
            break
        for index in pending:
            if (output.parent/'PAUSE').exists() or (output/'PAUSE').exists():
                return False
            original = records[index]
            checked = exact_evaluation(partner, original['x'], original['T_air_tar'], output)
            matches = checked['success'] and all(
                checked[key] is not None and original[key] is not None
                and np.isclose(checked[key], original[key], rtol=1e-6, atol=1e-9)
                for key in ('COP', 'SEC', 'VHC_MJ_per_m3', 'Pmax_MPa'))
            replay = dict(evaluation=index+1, matches=bool(matches), record=checked)
            append(path, replay)
            replays[index] = replay
            if not matches:
                validated[index]['validation_rejected'] = True
    snapshot = milestone_snapshot(validated)
    snapshot['label'] = 'final frontier validated in fresh processes'
    snapshot['validation_evaluations'] = len(replays)
    snapshot['validation_rejections'] = sum(not r['matches'] for r in replays.values())
    save(output/'validated_frontier.json', snapshot)
    from rhp_temperature_report import render_front
    render_front(snapshot, output/'validated_frontier.html')
    return True


def run_family(partner, budget, seed, output):
    from rhp_temperature_mobo import TemperatureParetoSearch, hypervolume_3d
    output = Path(output)/f'Propane_{partner}'
    output.mkdir(parents=True, exist_ok=True)
    with (output/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        meta = metadata(partner, budget, seed)
        path = output/'metadata.json'
        if path.exists() and json.loads(path.read_text()) != meta:
            raise ValueError('Existing run settings/source differ; use a fresh output directory')
        save(path, meta)
        records = read_journal(output/'records.jsonl')
        state_path = output/'state.json'
        saved = json.loads(state_path.read_text()) if state_path.exists() else {}
        if saved.get('complete'):
            return
        # State can lag a committed record by one if interrupted during append.
        if saved and len(records) not in (saved['evaluations'], saved['evaluations']+1):
            raise ValueError('Checkpoint and journal lengths disagree')
        search = TemperatureParetoSearch(records, seed=seed, state=saved.get('search_state'))
        # Reproduce the missing ask to recover model-cache visit counts on a torn checkpoint.
        if saved and len(records) == saved['evaluations']+1:
            search = TemperatureParetoSearch(records[:-1], seed=seed, state=saved.get('search_state'))
            x, temperature = search.ask()
            if temperature != records[-1]['T_air_tar'] or not np.allclose(x, records[-1]['x'], atol=1e-10, rtol=0):
                raise ValueError('Uncheckpointed proposal does not reproduce')
            search.tell(records[-1])
        write_milestones(output, records, meta['milestones'])

        def checkpoint(status):
            save(state_path, dict(evaluations=len(records), budget=budget, status=status,
                complete=status == 'complete', search_state=search.state,
                frontier_points=len(search.front), failures=sum(not r['success'] for r in records),
                hypervolume_normalized=hypervolume_3d(np.asarray(search.y)[search.front]),
                updated_utc=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())))

        while len(records) < budget:
            if (output.parent/'PAUSE').exists() or (output/'PAUSE').exists():
                checkpoint('paused')
                return
            started = time.perf_counter()
            x, temperature = search.ask()
            proposed_seconds = time.perf_counter()-started
            item = exact_evaluation(partner, x, temperature, output)
            item.update(proposal_seconds=proposed_seconds, search_diagnostic=search.diagnostic)
            search.tell(item)
            append(output/'records.jsonl', item)
            records.append(item)
            checkpoint('running')
            if len(records) in meta['milestones']:
                write_milestones(output, records, meta['milestones'])
            if len(records) == 1 or len(records) % 10 == 0:
                print(partner, len(records), '/', budget, 'front', len(search.front),
                      'failures', sum(not r['success'] for r in records), flush=True)
        checkpoint('validating final frontier')
        complete = verify_final(partner, output, records)
        checkpoint('complete' if complete else 'paused during validation')


def main():
    import torch
    torch.set_num_interop_threads(1)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partner', required=True, choices=PARTNERS)
    parser.add_argument('--budget', type=int, default=25000)
    parser.add_argument('--seed', type=int, default=20260913)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.budget < 1 or args.budget > 25000:
        parser.error('Budget must be between 1 and 25000')
    run_family(args.partner, args.budget, args.seed, args.output)


if __name__ == '__main__':
    main()
