"""Supervise the isolated Pareto experiment and publish progress and final reports."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import zipfile

from main_RHPmobo import PARTNERS, metadata_for
from main_RHPmobo_report import build_report, hypervolume

ROOT = Path(__file__).resolve().parent
OUTPUT = ROOT.parents[1]/'results/RHPcomposition/mobo_v17_4000'


def write_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temp.replace(path)


def progress(status):
    families = []
    for partner in PARTNERS:
        data = json.loads((OUTPUT/f'Propane_{partner}.json').read_text())
        records = data['records']
        good = [r for r in records if r['success']]
        families.append(dict(partner=partner, evaluations=len(records), budget=4000,
                             complete=data.get('complete', False),
                             best_COP=max((r['COP'] for r in good), default=None),
                             hypervolume=hypervolume(records), failed=len(records)-len(good)))
    stamp = time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())
    write_json(OUTPUT/'progress.json', dict(status=status, updated_utc=stamp, families=families))
    lines = ['# Clean COP–VHC_total Pareto search', '', f'Status: **{status}**. Updated {stamp}.', '',
             '4,000 candidate evaluations per family, starting from zero. Isolated v17 thermodynamic model. '
             'No previous results or GP states were imported. Interim best values await final validation.', '',
             '| Propane partner | Evaluations | Best observed COP | Failures |', '|---|---:|---:|---:|']
    for row in families:
        cop = f"{row['best_COP']:.7f}" if row['best_COP'] is not None else '—'
        lines.append(f"| {row['partner']} | {row['evaluations']:,} / 4,000 | {cop} | {row['failed']} |")
    if status == 'complete':
        lines += ['', '[Validated final report](REPORT.md)']
    temp = OUTPUT/'STATUS.tmp'
    temp.write_text('\n'.join(lines)+'\n')
    temp.replace(OUTPUT/'STATUS.md')


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with (OUTPUT/'run.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        manifest = json.loads((OUTPUT/'source_manifest.json').read_text())
        for name, digest in manifest['files'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest() == digest, name
        for partner in PARTNERS:
            path = OUTPUT/f'Propane_{partner}.json'
            if not path.exists():
                write_json(path, dict(metadata=metadata_for(partner,4000,20260910), records=[],
                                     frontier_indices=[], complete=False, clean_start=True))
        state = dict(status='running', pid=os.getpid(), process_group=os.getpgrp(),
                     started_utc=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()),
                     runtime=str(ROOT), output=str(OUTPUT), clean_start=True, concurrent_searches=1)
        write_json(OUTPUT/'runner_state.json', state)
        progress('running')
        for partner in PARTNERS:
            if json.loads((OUTPUT/f'Propane_{partner}.json').read_text()).get('complete'):
                continue
            # A fresh process per family releases its GP allocations before
            # loading the next family, as well as limiting concurrency to one.
            state['active_partner'] = partner
            with (OUTPUT/'run.log').open('a') as log:
                child = subprocess.Popen([sys.executable, '-u', str(ROOT/'main_RHPmobo.py'),
                                          '--budget','4000','--workers','1','--partner',partner,
                                          '--output',str(OUTPUT)],
                                         cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
                state['search_pid'] = child.pid
                write_json(OUTPUT/'runner_state.json', state)
                while child.poll() is None:
                    progress(f'running: Propane/{partner} — one search at a time')
                    try:
                        child.wait(timeout=30)
                    except subprocess.TimeoutExpired:
                        pass
            state['search_exit_code'] = child.returncode
            if child.returncode:
                state['status'] = 'failed'
                write_json(OUTPUT/'runner_state.json', state)
                progress('failed — inspect run.log')
                raise RuntimeError(f'{partner} search exited with status {child.returncode}')
        state.pop('active_partner', None)
        try:
            state['status'] = 'building report'
            write_json(OUTPUT/'runner_state.json', state)
            progress('building report')
            build_report(OUTPUT)
            state['status'] = 'complete'
            write_json(OUTPUT/'runner_state.json', state)
            progress('complete')
            files = [p for p in OUTPUT.iterdir() if p.is_file() and p.suffix not in ('.zip','.tmp','.lock')
                     and p.name not in ('evidence_manifest.json','supervisor.log')]
            write_json(OUTPUT/'evidence_manifest.json', dict(files={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}))
            with zipfile.ZipFile(OUTPUT/'evidence_snapshot.zip','w',zipfile.ZIP_DEFLATED) as archive:
                for p in files:
                    archive.write(p,p.name)
            print('Completed validated reports and evidence archive', flush=True)
        except Exception:
            state['status'] = 'report failed'
            write_json(OUTPUT/'runner_state.json', state)
            progress('report failed — inspect supervisor.log')
            raise


if __name__ == '__main__':
    main()
