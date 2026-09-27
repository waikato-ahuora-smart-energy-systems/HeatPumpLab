"""Fixed 200 C Pareto BO: assess every 500, stop when stable, cap at 5,000."""
import argparse
import fcntl
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from main_RHPtemperature_mobo import exact_evaluation, save
from rhp_frontier_stability import assess, frontier, render

ROOT = Path(__file__).resolve().parent
PARTNERS = ('Propane', 'Butane', 'Isobutane', 'Pentane', 'Isopentane', 'Neopentane', 'CO2')


def settings(partner, budget, seed, first='Propane', fixed_comp_frac=None):
    from importlib.metadata import version
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import (
        DEFAULT_BASE_CASE, DEFAULT_MAX_STAGE_PRESSURE_RATIO)
    paths = sorted(ROOT.glob('*.py'))+sorted((ROOT/'HeatPumpLab').rglob('*.py'))
    digest = hashlib.sha256()
    for path in paths:
        digest.update(str(path.relative_to(ROOT)).encode())
        digest.update(path.read_bytes())
    return dict(partner=partner, first=first, components=[first, partner],
                budget=budget, seed=seed, temperature_C=200,
                max_stage_pressure_ratio=DEFAULT_MAX_STAGE_PRESSURE_RATIO,
                compressor_gas_only=True, fixed_comp_frac=fixed_comp_frac,
                search_dimensions=4 if fixed_comp_frac is None else 3,
                base_case=dict(DEFAULT_BASE_CASE, T_air_tar=200.), check_interval=500,
                relative_tolerance=.001, consecutive_stable_blocks=3, earliest_stop=2000,
                minimum_frontier_points=10, minimum_recent_success_fraction=.8,
                curve_check='101 VHC thresholds fixed at first check; maximum relative COP change; max-VHC endpoint change',
                reference_point=[0., 0.], subset_size=128, local_start_minimum=256,
                global_exploration_every=8, numerical_threads=1,
                budget_accounting='All search candidates including failures; fresh final replays counted separately',
                source_sha256=digest.hexdigest(),
                versions={k: version(k) for k in ('numpy', 'scipy', 'botorch', 'torch', 'gpytorch', 'ctREFPROP')})


def pause_requested(output):
    return (output/'PAUSE').exists() or (output.parent/'PAUSE').exists()


def restore(result):
    from rhp_local_pareto_search import restore_search
    from rhp_mobo_search import ParetoBOSearch
    if result.get('local128_phase'):
        return restore_search(result)
    search = ParetoBOSearch((200.-15)*2/3, result['metadata']['seed'],
                           fixed_comp_frac=result['metadata'].get('fixed_comp_frac'))
    for r in result['records']:
        search.tell(r['x'], r['COP'], r['VHC_MJ_per_m3'], r['success'])
    search.kernel_state = result.get('kernel_state')
    search.mean_state = result.get('mean_state')
    return search, None


def validate(result, output):
    """Replay every exposed frontier point, preserving the search journal."""
    replays = result.setdefault('frontier_replays', {})
    while True:
        records = [dict(r, validation_rejected=(str(i+1) in replays and not replays[str(i+1)]['matches']))
                   for i, r in enumerate(result['records'])]
        pending = [r for r in frontier(records) if str(r['evaluation']) not in replays]
        if not pending:
            snapshot = assess(records)
            snapshot['label'] = 'Final frontier validated in fresh processes'
            snapshot['validation_evaluations'] = len(replays)
            snapshot['validation_rejections'] = sum(not r['matches'] for r in replays.values())
            save(output/'validated_frontier.json', snapshot)
            render(snapshot, output/'validated_frontier.html', result['metadata']['partner'],
                   first=result['metadata']['first'])
            return True
        for original in pending:
            if pause_requested(output):
                return False
            checked = exact_evaluation(result['metadata']['partner'], original['x'], 200., output,
                                       first=result['metadata']['first'])
            matches = checked['success'] and all(
                checked.get(k) is not None and original.get(k) is not None
                and np.isclose(checked[k], original[k], rtol=1e-6, atol=1e-9)
                for k in ('COP', 'VHC_MJ_per_m3', 'SEC', 'Pmax_MPa'))
            replays[str(original['evaluation'])] = dict(matches=bool(matches), record=checked)
            save(output/'result.json', result)


def run_family(partner, budget, seed, output, first='Propane', fixed_comp_frac=None):
    from rhp_local_pareto_search import choose_anchors, restore_search, REFRESH_EVERY
    if first not in PARTNERS or partner not in PARTNERS or first == partner:
        raise ValueError('Choose two distinct components from the paper primary pool')
    if fixed_comp_frac is not None and (not np.isfinite(fixed_comp_frac) or not .1 <= fixed_comp_frac <= .9):
        raise ValueError('fixed_comp_frac must be finite and between 0.1 and 0.9')
    output = Path(output)/f'{first}_{partner}'
    output.mkdir(parents=True, exist_ok=True)
    with (output/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        meta = settings(partner, budget, seed, first, fixed_comp_frac)
        path = output/'result.json'
        result = json.loads(path.read_text()) if path.exists() else dict(metadata=meta, records=[], checks=[])
        if result['metadata'] != meta:
            raise ValueError('Source/settings mismatch; use a fresh output directory')
        if result.get('complete'):
            return
        search, epoch = restore(result)

        def checkpoint(status):
            if result.get('local128_phase'):
                result['local128_state'] = dict(epoch_start=epoch,
                    model_states={str(k): v for k, v in search.states.items()})
            else:
                result['kernel_state'], result['mean_state'] = search.kernel_state, search.mean_state
            result.update(status=status, complete=status == 'complete',
                          updated_utc=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()))
            save(path, result)
            points = frontier(result['records'])
            vhc = np.array([p['VHC_MJ_per_m3'] for p in points])
            hv = float(np.dot(np.diff(np.r_[0., vhc]), [p['COP'] for p in points])) if points else 0.
            save(output/'state.json', dict(status=status, complete=result['complete'],
                evaluations=len(result['records']), budget=budget,
                failures=sum(not r['success'] for r in result['records']),
                frontier_points=len(points), hypervolume=hv,
                stop_reason=result.get('stop_reason'), last_check=result['checks'][-1] if result['checks'] else None,
                updated_utc=result['updated_utc']))

        def check_progress():
            n = len(result['records'])
            if n and n % 500 == 0 and (not result['checks'] or result['checks'][-1]['evaluations'] != n):
                snapshot = assess(result['records'], result['checks'][-1] if result['checks'] else None)
                save(output/f'check_{n:05d}.json', snapshot)
                render(snapshot, output/f'check_{n:05d}.html', partner, first=first)
                result['checks'].append({k: v for k, v in snapshot.items() if k != 'frontier'})
                if snapshot['stop_for_stability']:
                    result['stop_reason'] = 'stable_curve'
                checkpoint('running')

        # Recover a check if interruption followed the evaluation checkpoint.
        check_progress()
        while len(result['records']) < budget and not result.get('stop_reason'):
            if pause_requested(output):
                checkpoint('paused')
                return
            n = len(result['records'])
            if not result.get('local128_phase') and n >= 256:
                try:
                    choose_anchors(result['records'])
                except ValueError:
                    pass
                else:
                    result['local128_phase'] = dict(start_evaluations=n, seed=seed)
                    search, epoch = restore_search(result)
            elif result.get('local128_phase'):
                start = result['local128_phase']['start_evaluations']
                next_epoch = start+((n-start)//REFRESH_EVERY)*REFRESH_EVERY
                if next_epoch != epoch:
                    search, epoch = restore_search(result)
            started = time.perf_counter()
            # Retain whole-domain exploration after switching to local GP fits.
            if n >= 32 and (n-32) % 8 == 7:
                pool, distances = search._candidates()
                x = search.decode(pool[np.argmax(distances)])
                method = 'global_maximin'
            else:
                x = search.ask()
                method = 'local128' if epoch is not None else 'global_bo'
            if fixed_comp_frac is not None and x[2] != fixed_comp_frac:
                raise ValueError('Proposal changed the fixed compressor split')
            proposal_seconds = time.perf_counter()-started
            item = exact_evaluation(partner, x, 200., output, first=first)
            item.update(proposal_seconds=proposal_seconds, method=method)
            search.tell(item['x'], item['COP'], item['VHC_MJ_per_m3'], item['success'])
            result['records'].append(item)
            checkpoint('running')
            check_progress()
            if len(result['records']) == 1 or len(result['records']) % 10 == 0:
                print(f'{first}/{partner}', len(result['records']), '/', budget, method, flush=True)
        result.setdefault('stop_reason', 'evaluation_cap')
        checkpoint('validating')
        complete = validate(result, output)
        checkpoint('complete' if complete else 'paused during validation')


def main():
    import torch
    torch.set_num_interop_threads(1)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--partner', choices=PARTNERS, required=True)
    parser.add_argument('--first', choices=PARTNERS, default='Propane')
    parser.add_argument('--budget', type=int, default=5000)
    parser.add_argument('--seed', type=int, default=20260915)
    parser.add_argument('--output', required=True)
    parser.add_argument('--comp-frac', type=float, default=None,
                        help='Fix compressor stage allocation and remove it from BO variables')
    args = parser.parse_args()
    if not 1 <= args.budget <= 5000:
        parser.error('Budget must be between 1 and 5000')
    if args.first == args.partner:
        parser.error('The two mixture components must differ')
    if args.comp_frac is not None and (not np.isfinite(args.comp_frac) or not .1 <= args.comp_frac <= .9):
        parser.error('comp-frac must be finite and between 0.1 and 0.9')
    run_family(args.partner, args.budget, args.seed, args.output, first=args.first,
               fixed_comp_frac=args.comp_frac)


if __name__ == '__main__':
    main()
