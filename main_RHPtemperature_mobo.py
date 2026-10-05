"""Resumable per-temperature COP-VHC BO over the 100:5:220 C grid, per propane pair (up to 25,000 evaluations)."""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
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
RUN_VERSION = 'per-temperature-cop-vhc-v26-face-coverage'
# Stop rule on the mean per-temperature hypervolume, scaled from the fixed-split runner
# (100/2e-4/3/400 at ~1,000 evaluations for one temperature).
STOP_CHECK_INTERVAL = 1000
STOP_TOLERANCE = 2e-4
STOP_BLOCKS = 3
STOP_EARLIEST = 5000
# Ridge seeding at every fourth grid temperature (100, 120, ..., 220 C), split 0.5.
SEED_TEMPERATURES = tuple(float(t) for t in np.arange(100., 221., 20.))
SEED_SPLIT = .5             # clipped into --split-bounds
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import COMP2_POWER_FRAC_BOUNDS as SPLIT_BOUNDS
RIDGE_SEED_MIN_BUDGET = 2000


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


def exact_evaluation(partner, x, temperature, output, *, first='Propane', min_discharge_superheat_K=0.0,
                     components=None, base_case=None):
    """Solve one design; ``components`` (e.g. a ternary) overrides the first/partner binary.

    ``base_case`` replaces DEFAULT_BASE_CASE (e.g. another compressor model); T_air_tar is set from
    ``temperature`` either way.
    """
    from rhp_temperature_mobo import usable
    names = list(components) if components else [first, partner]
    objective = _OuterDesignObjective(max_evals=1,
        refrigerant=';'.join(names)+'|'+';'.join([f'{1/len(names):.12g}']*len(names))+'|mass',
        base_case=dict(base_case or DEFAULT_BASE_CASE, T_air_tar=float(temperature)),
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


def metadata(partner, budget, seed, *, batch_size=4, min_discharge_superheat_K=0.0, ridge_seed=True,
             stop_rule=True, split_bounds=SPLIT_BOUNDS):
    from rhp_temperature_mobo import MILESTONES, TEMPERATURES
    paths = sorted((ROOT/'HeatPumpLab').rglob('*.py'))+[ROOT/'rhp_temperature_mobo.py', ROOT/'rhp_temperature_report.py', ROOT/'rhp_ridge.py', Path(__file__).resolve()]
    return dict(run_version=RUN_VERSION, partner=partner, budget=budget, seed=seed,
        batch_size=int(batch_size), min_discharge_superheat_K=float(min_discharge_superheat_K),
        comp2_power_frac_bounds=[float(v) for v in split_bounds],
        constraint_handling='constraint GP on margins (per region) multiplies the acquisition; '
                            'distance penalty for numerical failures only',
        ridge_seed=bool(ridge_seed),
        ridge_seed_temperatures_C=list(SEED_TEMPERATURES) if ridge_seed else None,
        ridge_seed_split=float(np.clip(SEED_SPLIT, *split_bounds)) if ridge_seed else None,
        stop_rule=dict(check_interval=STOP_CHECK_INTERVAL, relative_tolerance=STOP_TOLERANCE,
                       consecutive_stable_blocks=STOP_BLOCKS, earliest_stop=STOP_EARLIEST) if stop_rule else None,
        max_stage_pressure_ratio=DEFAULT_MAX_STAGE_PRESSURE_RATIO,
        compressor_gas_only=True,
        milestones=sorted(set([m for m in MILESTONES if m <= budget]+[budget])),
        temperatures_C=TEMPERATURES.tolist(), base_case=DEFAULT_BASE_CASE,
        objectives=['COP', 'VHC_total_MJ_per_m3'], condition='T_air_tar_C on the 100:5:220 grid',
        quality='mean over the 25 temperatures of the 2D (COP/4, VHC/10) hypervolume against the origin',
        reference_point_physical=[0., 0.], scales=[4., 10.],
        subset_size=128, initial_evaluations=200, numerical_threads=1,
        known_objective='Temperature is an exact input condition; GP outputs are COP and VHC only',
        search=f'per-temperature qLogEHVI on a shared 5-input GP; adaptive step k targets level (7k) mod 25; '
               f'batches of {int(batch_size)} with a pending-point penalty; per temperature one face-restricted '
               'region per zero-lift face plus 4 front centres; boundary candidates in every pool; '
               'maximin at the scheduled temperature every eighth adaptive step',
        evaluation_timeout_seconds=180, timeout_retries=0, native_tasks_per_worker=1,
        budget_accounting='Includes initial, BO, exploratory and failed candidates; final frontier replays counted separately',
        source_sha256=hashlib.sha256(b''.join(p.read_bytes() for p in paths)).hexdigest(),
        versions={name: version(name) for name in ('numpy', 'scipy', 'botorch', 'torch', 'gpytorch', 'ctREFPROP')})


def milestone_snapshot(records, previous=None):
    from rhp_temperature_mobo import (TEMPERATURES, mean_temperature_hypervolume, per_temperature_front_indices,
                                      per_temperature_hypervolumes, usable)
    indices = per_temperature_front_indices(records)
    frontier = [dict(evaluation=i+1, **records[i]) for i in indices]
    hv = mean_temperature_hypervolume(records)
    level_hv = per_temperature_hypervolumes(records)
    earlier = previous['hypervolume_normalized'] if previous else None
    per_temperature = []
    for temperature in TEMPERATURES:
        at = [r for r in records if r['T_air_tar'] == temperature]
        good = [r for r in at if usable(r)]
        per_temperature.append(dict(temperature_C=float(temperature), evaluations=len(at),
            successful=len(good), failures=sum(not r['success'] for r in at),
            frontier_points=int(sum(r['T_air_tar'] == temperature for r in frontier)),
            hypervolume_normalized=level_hv[float(temperature)],
            best_COP=max((r['COP'] for r in good), default=None),
            best_VHC_MJ_per_m3=max((r['VHC_MJ_per_m3'] for r in good), default=None)))
    return dict(evaluations=len(records), label='observed frontier, prior to final fresh-process validation',
        hypervolume_normalized=hv, hypervolume_physical=hv*4*10,
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


def verify_final(partner, output, records, min_discharge_superheat_K=0.0):
    from rhp_temperature_mobo import per_temperature_front_indices as front_indices
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
            checked = exact_evaluation(partner, original['x'], original['T_air_tar'], output,
                                       min_discharge_superheat_K=min_discharge_superheat_K)
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


def stability_check(records, checks):
    """Append hypervolume checks at each interval crossed; True once the front has stopped moving."""
    from rhp_temperature_mobo import mean_temperature_hypervolume
    mark = (checks[-1]['evaluations'] if checks else 0)+STOP_CHECK_INTERVAL
    while mark <= len(records):
        hv = mean_temperature_hypervolume(records[:mark])
        change = hv/checks[-1]['hypervolume_normalized']-1 if checks and checks[-1]['hypervolume_normalized'] else None
        checks.append(dict(evaluations=mark, hypervolume_normalized=hv, relative_change=change))
        mark += STOP_CHECK_INTERVAL
    recent = [c['relative_change'] for c in checks[-STOP_BLOCKS:]]
    return bool(checks and checks[-1]['evaluations'] >= STOP_EARLIEST and len(recent) == STOP_BLOCKS
                and all(c is not None and abs(c) < STOP_TOLERANCE for c in recent))


def rewrite_journal(path, records):
    temporary = Path(path).with_suffix('.jsonl.tmp')
    with temporary.open('w') as stream:
        for record in records:
            stream.write(json.dumps(record, allow_nan=False)+'\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def run_family(partner, budget, seed, output, *, batch_size=4, min_discharge_superheat_K=0.0, ridge_seed=True,
               stop_rule=True, split_bounds=SPLIT_BOUNDS):
    from rhp_temperature_mobo import TemperatureParetoSearch
    if ridge_seed and budget < RIDGE_SEED_MIN_BUDGET:
        raise ValueError(f'Ridge seeding uses about 1,000 evaluations; budget must be at least '
                         f'{RIDGE_SEED_MIN_BUDGET} (or pass ridge_seed=False / --no-ridge-seed)')
    output = Path(output)/f'Propane_{partner}'
    output.mkdir(parents=True, exist_ok=True)
    journal = output/'records.jsonl'

    def evaluate(x, temperature):
        return exact_evaluation(partner, x, temperature, output,
                                min_discharge_superheat_K=min_discharge_superheat_K)

    with (output/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        meta = metadata(partner, budget, seed, batch_size=batch_size,
                        min_discharge_superheat_K=min_discharge_superheat_K, ridge_seed=ridge_seed,
                        stop_rule=stop_rule, split_bounds=split_bounds)
        path = output/'metadata.json'
        if path.exists() and json.loads(path.read_text()) != meta:
            raise ValueError('Existing run settings/source differ; use a fresh output directory')
        save(path, meta)
        records = read_journal(journal)
        state_path = output/'state.json'
        saved = json.loads(state_path.read_text()) if state_path.exists() else {}
        if saved.get('complete'):
            return
        seeding = ridge_seed and not saved.get('ridge_seed_done')
        if seeding and len(records) > saved.get('evaluations', 0):
            # Seed records are committed one temperature at a time; drop a torn temperature and redo it.
            records = records[:saved.get('evaluations', 0)]
            rewrite_journal(journal, records)
        # State can lag the journal by up to one batch if interrupted while appending.
        lag = len(records)-saved.get('evaluations', 0) if saved else 0
        if saved and not 0 <= lag <= batch_size:
            raise ValueError('Checkpoint and journal lengths disagree')
        if lag:
            # Reproduce the interrupted batch to recover model-cache visit counts, then drop its unevaluated rest.
            search = TemperatureParetoSearch(records[:-lag], seed=seed, state=saved.get('search_state'),
                                             split_bounds=split_bounds)
            batch = search.ask_batch(min(batch_size, budget-len(records)+lag))
            for (x, temperature), record in zip(batch, records[-lag:]):
                if temperature != record['T_air_tar'] or not np.allclose(x, record['x'], atol=1e-10, rtol=0):
                    raise ValueError('Uncheckpointed proposal does not reproduce')
                search.tell(record)
            search.pending = []
        else:
            search = TemperatureParetoSearch(records, seed=seed, state=saved.get('search_state'),
                                             split_bounds=split_bounds)
        checks = list(saved.get('stability_checks', []))
        seeded = list(saved.get('ridge_seeded_temperatures', []))
        ridges = list(saved.get('ridge_seeds', []))
        seed_done = bool(saved.get('ridge_seed_done')) or not ridge_seed
        stop_reason = saved.get('stop_reason')
        write_milestones(output, records, meta['milestones'])

        def checkpoint(status):
            save(state_path, dict(evaluations=len(records), budget=budget, status=status,
                complete=status == 'complete', search_state=search.state,
                frontier_points=len(search.front), failures=sum(not r['success'] for r in records),
                hypervolume_normalized=search.hypervolume(),
                ridge_seed_done=seed_done, ridge_seeded_temperatures=seeded, ridge_seeds=ridges,
                stability_checks=checks, stop_reason=stop_reason,
                updated_utc=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())))

        def commit(item):
            search.tell(item)
            append(journal, item)
            records.append(item)

        if not seed_done:
            from rhp_ridge import seed as ridge_seed_run
            if len(records) != sum(r.get('method') == 'ridge_seed' for r in records):
                raise ValueError('Ridge seeding must run before any other evaluation')
            for temperature in SEED_TEMPERATURES:
                if temperature in seeded:
                    continue
                if (output.parent/'PAUSE').exists() or (output/'PAUSE').exists():
                    checkpoint('paused')
                    return
                started = time.perf_counter()

                def seed_evaluation(x, temperature=temperature):
                    item = evaluate(x, temperature)
                    item.update(method='ridge_seed', proposal_seconds=0., batch_size=1,
                                search_diagnostic=dict(method='ridge_seed'))
                    return item
                seed_records, seed_ridges = ridge_seed_run(seed_evaluation, (temperature-15)*2/3,
                                                           float(np.clip(SEED_SPLIT, *split_bounds)),
                                                           workers=batch_size)
                offset = len(records)
                for item in seed_records:
                    commit(item)
                ridges.extend(dict(r, T_air_tar=temperature,
                                   record=None if r['record'] is None else offset+r['record']) for r in seed_ridges)
                seeded.append(temperature)
                checkpoint('ridge seeding')
                print(f'{partner} ridge seeding {temperature:.0f} C: {len(seed_records)} evaluations, '
                      f'{len(seed_ridges)} ridges ({time.perf_counter()-started:.0f} s)', flush=True)
            seed_done = True
            search.offset = len(records)
            checkpoint('running')

        with ThreadPoolExecutor(max(1, batch_size)) as pool:
            while len(records) < budget and not stop_reason:
                if (output.parent/'PAUSE').exists() or (output/'PAUSE').exists():
                    checkpoint('paused')
                    return
                q = min(batch_size, budget-len(records))
                started = time.perf_counter()
                batch = search.ask_batch(q)
                proposed_seconds = (time.perf_counter()-started)/q
                items = list(pool.map(lambda job: evaluate(*job), batch))
                before = len(records)
                for item, diagnostic in zip(items, search.batch_diagnostics):
                    item.update(proposal_seconds=proposed_seconds, search_diagnostic=diagnostic, batch_size=q)
                    commit(item)
                search.pending = []
                if stop_rule and stability_check(records, checks):
                    stop_reason = 'stable_hypervolume'
                checkpoint('running')
                if any(before < m <= len(records) for m in meta['milestones']):
                    write_milestones(output, records, meta['milestones'])
                if before == 0 or len(records)//100 > before//100:
                    print(partner, len(records), '/', budget, 'front', len(search.front),
                          'failures', sum(not r['success'] for r in records), flush=True)
        if stop_reason and len(records) not in meta['milestones']:
            # Record the stopping point as a final milestone snapshot.
            write_milestones(output, records, sorted(set([m for m in meta['milestones'] if m < len(records)]
                                                         + [len(records)])))
        checkpoint('validating final frontier')
        complete = verify_final(partner, output, records, min_discharge_superheat_K)
        checkpoint('complete' if complete else 'paused during validation')


def main():
    import torch
    torch.set_num_interop_threads(1)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partner', required=True, choices=PARTNERS)
    parser.add_argument('--budget', type=int, default=25000)
    parser.add_argument('--seed', type=int, default=20260913)
    parser.add_argument('--output', required=True)
    parser.add_argument('--batch-size', type=int, default=4,
                        help='Proposals evaluated in parallel per BO step (pending-point penalty)')
    parser.add_argument('--min-discharge-superheat', type=float, default=0.0,
                        help='Required compressor discharge superheat margin, K')
    parser.add_argument('--ridge-seed', action=argparse.BooleanOptionalAction, default=True,
                        help='Seed with pinch-balance ridge solves at 100, 120, ..., 220 C (rhp_ridge.py)')
    parser.add_argument('--stop-rule', action=argparse.BooleanOptionalAction, default=True,
                        help=f'Stop once 3D hypervolume changes < {STOP_TOLERANCE:g} over {STOP_BLOCKS} blocks of '
                             f'{STOP_CHECK_INTERVAL} evaluations (not before {STOP_EARLIEST})')
    parser.add_argument('--split-bounds', type=float, nargs=2, default=list(SPLIT_BOUNDS), metavar=('LOW', 'HIGH'),
                        help=f'Narrower bounds on the compressor-2 power fraction, within the global '
                             f'{SPLIT_BOUNDS[0]}-{SPLIT_BOUNDS[1]} (default)')
    args = parser.parse_args()
    if not SPLIT_BOUNDS[0] <= args.split_bounds[0] < args.split_bounds[1] <= SPLIT_BOUNDS[1]:
        parser.error(f'Split bounds must satisfy {SPLIT_BOUNDS[0]} <= LOW < HIGH <= {SPLIT_BOUNDS[1]}')
    if args.budget < 1 or args.budget > 25000:
        parser.error('Budget must be between 1 and 25000')
    if args.batch_size < 1:
        parser.error('Batch size must be at least 1')
    try:
        run_family(args.partner, args.budget, args.seed, args.output, batch_size=args.batch_size,
                   min_discharge_superheat_K=args.min_discharge_superheat, ridge_seed=args.ridge_seed,
                   stop_rule=args.stop_rule, split_bounds=tuple(args.split_bounds))
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == '__main__':
    main()
