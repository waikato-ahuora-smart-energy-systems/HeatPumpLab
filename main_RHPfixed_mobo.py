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


def settings(partner, budget, seed, first='Propane', fixed_comp_frac=None, min_discharge_superheat_K=0.0,
             batch_size=4, ridge_seed=False):
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
                min_discharge_superheat_K=float(min_discharge_superheat_K),
                constraint_model='GP on discharge/suction superheat and PR margins; distance penalty for numerical failures only',
                search_dimensions=4 if fixed_comp_frac is None else 3,
                base_case=dict(DEFAULT_BASE_CASE, T_air_tar=200.), check_interval=100,
                relative_tolerance=2e-4, consecutive_stable_blocks=3, earliest_stop=400,
                stop_rule='hypervolume only; curve and extent changes reported, not required',
                batch_size=int(batch_size),
                ridge_seed=bool(ridge_seed),
                ridge_seed_design=('pinch-balance ridges A[k+1]=A[k] at propane fractions 0, 1/3, 2/3, 1 and '
                                   'theta 0, 0.5, 1; 6-point coarse S sweep then Illinois secant; one refinement line per theta '
                                   'in the widest COP-VHC gap between solved compositions') if ridge_seed else None,
                search_coverage='half of global BO steps restricted to a zero-lift face; one local region per face',
                minimum_frontier_points=10, minimum_recent_success_fraction=.8,
                curve_check='101 VHC thresholds fixed at first check; maximum relative COP change; max-VHC endpoint change (report only)',
                reference_point=[0., 0.], subset_size=128, local_start_minimum=256,
                global_exploration_every=8, numerical_threads=1,
                budget_accounting='All search candidates including failures; fresh final replays counted separately',
                source_sha256=digest.hexdigest(),
                versions={k: version(k) for k in ('numpy', 'scipy', 'botorch', 'torch', 'gpytorch', 'ctREFPROP')})


def pause_requested(output):
    return (output/'PAUSE').exists() or (output.parent/'PAUSE').exists()


def restore(result):
    from rhp_local_pareto_search import max_lift_for, restore_search
    from rhp_mobo_search import ParetoBOSearch
    if result.get('local128_phase'):
        return restore_search(result)
    search = ParetoBOSearch(max_lift_for(result['metadata']), result['metadata']['seed'],
                           fixed_comp_frac=result['metadata'].get('fixed_comp_frac'))
    for r in result['records']:
        search.tell(r['x'], r['COP'], r['VHC_MJ_per_m3'], r['success'], r.get('constraint_margins'))
    search.offset = seed_count(result)
    search.kernel_state = result.get('kernel_state')
    search.mean_state = result.get('mean_state')
    search.constraint_states = dict(result.get('constraint_states') or {})
    return search, None


def seed_count(result):
    n = 0
    for r in result['records']:
        if r.get('method') != 'ridge_seed':
            break
        n += 1
    return n


def local_start(result):
    """Local-128 search starts 256 evaluations in, or 128 after ridge seeding."""
    n = seed_count(result)
    return n+128 if n else 256


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
                                       first=result['metadata']['first'],
                                       min_discharge_superheat_K=result['metadata'].get('min_discharge_superheat_K', 0.0))
            matches = checked['success'] and all(
                checked.get(k) is not None and original.get(k) is not None
                and np.isclose(checked[k], original[k], rtol=1e-6, atol=1e-9)
                for k in ('COP', 'VHC_MJ_per_m3', 'SEC', 'Pmax_MPa'))
            replays[str(original['evaluation'])] = dict(matches=bool(matches), record=checked)
            save(output/'result.json', result)


RIDGE_SEED_MIN_BUDGET = 200


def run_family(partner, budget, seed, output, first='Propane', fixed_comp_frac=None, min_discharge_superheat_K=0.0,
               batch_size=4, ridge_seed=True):
    from concurrent.futures import ThreadPoolExecutor
    from rhp_local_pareto_search import choose_anchors, restore_search, REFRESH_EVERY
    if first not in PARTNERS or partner not in PARTNERS or first == partner:
        raise ValueError('Choose two distinct components from the paper primary pool')
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import COMP2_POWER_FRAC_BOUNDS as SPLIT
    if fixed_comp_frac is not None and (not np.isfinite(fixed_comp_frac) or not SPLIT[0] <= fixed_comp_frac <= SPLIT[1]):
        raise ValueError(f'fixed_comp_frac must be finite and between {SPLIT[0]} and {SPLIT[1]}')
    if ridge_seed and budget < RIDGE_SEED_MIN_BUDGET:
        raise ValueError(f'Ridge seeding uses about 130 evaluations; budget must be at least '
                         f'{RIDGE_SEED_MIN_BUDGET} (or pass ridge_seed=False / --no-ridge-seed)')
    output = Path(output)/f'{first}_{partner}'
    output.mkdir(parents=True, exist_ok=True)
    with (output/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        meta = settings(partner, budget, seed, first, fixed_comp_frac, min_discharge_superheat_K, batch_size,
                        ridge_seed)
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
                    model_states={str(k): v for k, v in search.states.items()},
                    constraint_states={str(k): v for k, v in search.constraint_states.items()})
            else:
                result['kernel_state'], result['mean_state'] = search.kernel_state, search.mean_state
                result['constraint_states'] = search.constraint_states
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
            interval = meta['check_interval']
            if n and n % interval == 0 and (not result['checks'] or result['checks'][-1]['evaluations'] != n):
                snapshot = assess(result['records'], result['checks'][-1] if result['checks'] else None,
                                  tolerance=meta['relative_tolerance'], minimum_evaluations=meta['earliest_stop'],
                                  hv_only=True, required_blocks=meta['consecutive_stable_blocks'])
                save(output/f'check_{n:05d}.json', snapshot)
                if n % 500 == 0:
                    render(snapshot, output/f'check_{n:05d}.html', partner, first=first)
                result['checks'].append({k: v for k, v in snapshot.items() if k != 'frontier'})
                if snapshot['stop_for_stability']:
                    result['stop_reason'] = 'stable_curve'
                checkpoint('running')

        if ridge_seed and not result.get('ridge_seed_done'):
            # Seeding runs as one unit; an interrupted seed phase is simply repeated.
            from rhp_local_pareto_search import max_lift_for
            from rhp_ridge import seed as ridge_seed_run
            if result['records']:
                raise ValueError('Ridge seeding must run before any other evaluation')

            def seed_evaluation(x):
                item = exact_evaluation(partner, x, 200., output, first=first,
                                        min_discharge_superheat_K=min_discharge_superheat_K)
                item.update(method='ridge_seed', proposal_seconds=0., batch_size=1)
                return item
            started = time.perf_counter()
            seed_records, ridges = ridge_seed_run(
                seed_evaluation, max_lift_for(meta), fixed_comp_frac if fixed_comp_frac is not None else .5,
                workers=batch_size)
            for item in seed_records:
                search.tell(item['x'], item['COP'], item['VHC_MJ_per_m3'], item['success'],
                            item.get('constraint_margins'))
                result['records'].append(item)
            search.offset = len(seed_records)
            result.update(ridge_seed_done=True, ridge_seeds=ridges,
                          ridge_seed_seconds=time.perf_counter()-started)
            print(f'{first}/{partner} ridge seeding: {len(seed_records)} evaluations, {len(ridges)} ridges', flush=True)
            checkpoint('running')
        # Recover a check if interruption followed the evaluation checkpoint.
        check_progress()
        while len(result['records']) < budget and not result.get('stop_reason'):
            if pause_requested(output):
                checkpoint('paused')
                return
            n = len(result['records'])
            if not result.get('local128_phase') and n >= local_start(result):
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
            # Propose a batch; slots keep their serial order for reproducible records.
            # A batch never crosses the local-phase switch or a local refresh epoch.
            limit = budget-n
            if not result.get('local128_phase') and n < local_start(result):
                limit = min(limit, local_start(result)-n)
            elif result.get('local128_phase'):
                limit = min(limit, REFRESH_EVERY-((n-result['local128_phase']['start_evaluations']) % REFRESH_EVERY))
            batch = []
            for m in range(n, n+min(batch_size, limit)):
                # Retain whole-domain exploration after switching to local GP fits.
                g = m-seed_count(result)
                if g >= 32 and (g-32) % 8 == 7:
                    pool, distances = search._candidates()
                    x = search.decode(pool[np.argmax(distances)])
                    method = 'global_maximin'
                    timing = {}
                else:
                    ask_started = time.perf_counter()
                    x = search.ask()
                    method = 'local128' if epoch is not None else 'global_bo'
                    diag = getattr(search, 'last_diagnostic', {}) if epoch is not None else {}
                    timing = {k: diag[k] for k in ('fit_seconds', 'constraint_seconds', 'pool_score_seconds',
                                                   'acquisition_opt_seconds') if k in diag}
                    timing['ask_seconds'] = time.perf_counter()-ask_started
                if fixed_comp_frac is not None and x[2] != fixed_comp_frac:
                    raise ValueError('Proposal changed the fixed compressor split')
                search.pending.append(search.encode(x).tolist())
                batch.append((x, method, timing))
            proposal_seconds = (time.perf_counter()-started)/len(batch)
            with ThreadPoolExecutor(len(batch)) as pool_:
                items = list(pool_.map(lambda b: exact_evaluation(
                    partner, b[0], 200., output, first=first,
                    min_discharge_superheat_K=min_discharge_superheat_K), batch))
            for (x, method, timing), item in zip(batch, items):
                item.update(proposal_seconds=proposal_seconds, method=method, batch_size=len(batch),
                            proposal_timing=timing)
                search.tell(item['x'], item['COP'], item['VHC_MJ_per_m3'], item['success'],
                            item.get('constraint_margins'))
                result['records'].append(item)
                checkpoint('running')
                check_progress()
                if len(result['records']) % 20 == 0:
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
    parser.add_argument('--ridge-seed', action=argparse.BooleanOptionalAction, default=True,
                        help='Seed the search with pinch-balance ridge solves (see rhp_ridge.py); '
                             'on by default, --no-ridge-seed to disable')
    parser.add_argument('--batch-size', type=int, default=4,
                        help='Candidates proposed and evaluated in parallel per step')
    parser.add_argument('--min-discharge-superheat', type=float, default=0.0,
                        help='Required compressor discharge superheat in K (supercritical discharge exempt)')
    args = parser.parse_args()
    if not 1 <= args.budget <= 5000:
        parser.error('Budget must be between 1 and 5000')
    if args.first == args.partner:
        parser.error('The two mixture components must differ')
    from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation._types import COMP2_POWER_FRAC_BOUNDS as SPLIT
    if args.comp_frac is not None and (not np.isfinite(args.comp_frac) or not SPLIT[0] <= args.comp_frac <= SPLIT[1]):
        parser.error(f'comp-frac must be finite and between {SPLIT[0]} and {SPLIT[1]}')
    run_family(args.partner, args.budget, args.seed, args.output, first=args.first,
               fixed_comp_frac=args.comp_frac, min_discharge_superheat_K=args.min_discharge_superheat,
               batch_size=args.batch_size, ridge_seed=args.ridge_seed)


if __name__ == '__main__':
    main()
