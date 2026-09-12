"""Maintained entry points must not import archived experiment runtimes."""

import os
from pathlib import Path
import subprocess
import sys
import textwrap


def test_search_and_reporting_imports_use_maintained_source():
    root = Path(__file__).resolve().parents[1]
    script = textwrap.dedent("""
        from pathlib import Path
        import sys

        root = Path.cwd().resolve()
        import rhp_local_pareto_search
        import rhp_mobo_search
        import main_RHPmobo
        import main_RHPmobo_snapshot
        import main_RHPlocal_subset_trial
        import main_RHPlocal_subset_report
        import HeatPumpLab

        for module in (rhp_local_pareto_search, rhp_mobo_search, main_RHPmobo,
                       main_RHPmobo_snapshot, main_RHPlocal_subset_trial,
                       main_RHPlocal_subset_report):
            assert Path(module.__file__).resolve().parent == root, module.__file__
        assert Path(HeatPumpLab.__file__).resolve().is_relative_to(root / 'HeatPumpLab')
        assert rhp_local_pareto_search.restore_search.__module__ == 'rhp_local_pareto_search'
        archive = root / 'experiments'
        for entry in sys.path:
            assert not Path(entry or '.').resolve().is_relative_to(archive), entry
        for name, module in tuple(sys.modules.items()):
            path = getattr(module, '__file__', None)
            if path:
                assert not Path(path).resolve().is_relative_to(archive), (name, path)
    """)
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1", MPLBACKEND="Agg")
    env.pop("PYTHONPATH", None)
    result = subprocess.run(
        [sys.executable, "-c", script], cwd=root, env=env,
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr
