"""Run four isolated local-128 continuations and validate their final fronts."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parents[1]
FROZEN = PROJECT/'experiments/mobo_v17_clean_runtime'
OUTPUT = PROJECT/'results/RHPcomposition/mobo_v17_4000'
PARTNERS = ('Butane', 'Isopentane', 'Pentane', 'Isobutane')


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    temporary.replace(path)


def progress(state):
    rows = []
    for partner in PARTNERS:
        data = json.loads((OUTPUT/f'Propane_{partner}.json').read_text())
        good = [r for r in data['records'] if r['success'] and not r.get('validation_rejected')]
        rows.append(dict(partner=partner, evaluations=len(data['records']), budget=4000,
                         complete=data.get('complete', False), best_COP=max(r['COP'] for r in good),
                         failed=sum(not r['success'] for r in data['records'])))
    stamp = time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime())
    save(OUTPUT/'progress.json', dict(status=state['status'], updated_utc=stamp, families=rows,
                                    concurrent_searches=4, subset_size=128))
    lines = ['# COP–VHC_total Pareto search', '', f"Status: **{state['status']}**. Updated {stamp}.", '',
             'Four parallel searches; one numerical thread per search. Nearest 128 successful observations per GP. '
             'Six local boxes refreshed every 48 new candidates. Existing v17 observations retained; budget 4,000 per family.', '',
             '| Propane partner | Evaluations | Best observed COP | Failures |', '|---|---:|---:|---:|']
    for row in rows:
        lines.append(f"| {row['partner']} | {row['evaluations']:,} / 4,000 | {row['best_COP']:.7f} | {row['failed']} |")
    lines += ['', '[Continuation protocol](LOCAL128_PROTOCOL.md)']
    if state['status'] == 'complete':
        lines += ['', '[Validated final report](REPORT.md)']
    (OUTPUT/'STATUS.tmp').write_text('\n'.join(lines)+'\n')
    (OUTPUT/'STATUS.tmp').replace(OUTPUT/'STATUS.md')


def main():
    with (OUTPUT/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for manifest_path, base in ((OUTPUT/'source_manifest.json', FROZEN),
                                    (ROOT/'source_manifest.json', ROOT)):
            for name, digest in json.loads(manifest_path.read_text())['files'].items():
                assert hashlib.sha256((base/name).read_bytes()).hexdigest() == digest, name
        archive = OUTPUT/'before_local128'
        if not archive.exists():
            previous = json.loads((OUTPUT/'runner_state.json').read_text())
            assert previous['status'] == 'paused', previous
            archive.mkdir()
            for name in [f'Propane_{p}.json' for p in PARTNERS]+['runner_state.json', 'progress.json', 'STATUS.md']:
                shutil.copy2(OUTPUT/name, archive/name)
            save(archive/'manifest.json', dict(files={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                                                    for p in archive.iterdir() if p.is_file()}))
        shutil.copy2(ROOT/'PROTOCOL.md', OUTPUT/'LOCAL128_PROTOCOL.md')
        shutil.copy2(ROOT/'source_manifest.json', OUTPUT/'local128_source_manifest.json')
        state = dict(status='running', pid=os.getpid(), process_group=os.getpgrp(),
                     started_utc=time.strftime('%Y-%m-%d %H:%M:%S UTC', time.gmtime()),
                     runtime=str(ROOT), output=str(OUTPUT), concurrent_searches=4,
                     numerical_threads_per_search=1, subset_size=128, search_pids={}, exit_codes={})
        children = {}
        env = dict(os.environ, OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
                   VECLIB_MAXIMUM_THREADS='1', NUMEXPR_NUM_THREADS='1', OMP_THREAD_LIMIT='1',
                   MPLCONFIGDIR='/tmp/heatpumplab-botorch-mpl',
                   TORCH_EXTENSIONS_DIR='/tmp/heatpumplab-torch-extensions')
        save(OUTPUT/'runner_state.json', state)
        try:
            for partner in PARTNERS:
                if json.loads((OUTPUT/f'Propane_{partner}.json').read_text()).get('complete'):
                    continue
                with (OUTPUT/f'local128_{partner}.log').open('a') as log:
                    children[partner] = subprocess.Popen(
                        [sys.executable, '-u', str(ROOT/'run_local.py'), '--partner', partner, '--output', str(OUTPUT)],
                        cwd=PROJECT, env=env, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)
                state['search_pids'][partner] = children[partner].pid
                save(OUTPUT/'runner_state.json', state)
            while any(child.poll() is None for child in children.values()):
                state['exit_codes'] = {p: c.returncode for p, c in children.items() if c.poll() is not None}
                if any(state['exit_codes'].values()):
                    state['status'] = 'running with failed family'
                save(OUTPUT/'runner_state.json', state)
                progress(state)
                time.sleep(20)
            state['exit_codes'] = {p: c.returncode for p, c in children.items()}
            if any(state['exit_codes'].values()):
                raise RuntimeError(state['exit_codes'])
            state['status'] = 'building report'
            save(OUTPUT/'runner_state.json', state)
            progress(state)
            subprocess.run([sys.executable, str(FROZEN/'main_RHPmobo_report.py'), '--output', str(OUTPUT)],
                           cwd=FROZEN, env=env, check=True)
            # The frozen reporter describes the original global phase; add the continuation method.
            report = OUTPUT/'REPORT.md'
            content = report.read_text()
            start = content.index('BoTorch qLogEHVI uses separate GP outputs')
            end = content.index('## Validation and limits', start)
            content = content[:start] + (ROOT/'PROTOCOL.md').read_text().split('\n', 2)[2] + '\n\n' + content[end:]
            start = content.index('Resume the saved checkpoints with ')
            content = content[:start] + f'Resume this continuation with `{sys.executable} {ROOT / "supervise.py"}`.\n'
            report.write_text(content)
            state['status'] = 'complete'
        except BaseException as exc:
            state['status'] = 'failed'
            state['error'] = repr(exc)
            # Stop child searches if supervisor setup/reporting fails unexpectedly.
            for child in children.values():
                if child.poll() is None:
                    child.terminate()
            raise
        finally:
            save(OUTPUT/'runner_state.json', state)
            progress(state)


if __name__ == '__main__':
    main()
