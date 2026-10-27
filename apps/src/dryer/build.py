"""Assemble the Advanced dryer heat pumps app (apps/dryer/index.html, CoolProp in the browser).

  uv run python apps/src/dryer/build.py [--out DIR] [--check]

index.template.html is the page with each <script> body replaced by an ``{{include ...}}`` marker. The build
puts the named files back in order, with no separators:

  sghp_model.js      the SGHP cycle model (saturation dome and solvers) the dryer models build on
  dryer_model.js     the RHP-A / SLTHP cycle models (``DRYER``)
  dryer_cascades.js  the multi-loop cycles: TT cascade, n-level cascade, HP + MVR, sCO2 Brayton (``DRYCASC``)
  ammonia_water.js   the IAPWS G4-01 ammonia-water properties (``AMMWAT``)
  dryer_more.js      flash tank / economiser, dual-pressure, Stirling and ammonia-water hybrid cycles (``DRYMORE``)
  cp_preset.js       promises the page waits on until CoolProp is ready
  @coolprop          the CoolProp Emscripten loader, generated from ../coolprop/coolprop.js
  adapter.js         CoolProp loading and the CoolPropFluid property adapter
  dryer_study.js     the REFPROP design-study views
  dryer_live.js      the live-design controls and results
  dryer_cases.js     saved cases, case links and the comparison
  boot.js            starts the page

coolprop.wasm is copied next to index.html. The study results (dashboard_data.json, cycles_data.json,
state_5C_1bar.json) already live in apps/dryer and are not rebuilt. --check rebuilds into a temporary folder and
fails if the result differs from apps/dryer/index.html.
"""
from __future__ import annotations

import argparse
import re
import shutil
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
CP_DIR = HERE.parent / 'coolprop'
OUT = HERE.parent.parent / 'dryer'
MARKER = re.compile(r'\{\{include ([^}]+)\}\}')


def coolprop_module(coolprop_js: Path) -> str:
    """The body of the page's <script type="module">: the CoolProp loader with the stack helpers exposed (so the
    page can undo the stack left behind when a CoolProp C++ exception reaches JavaScript) and the factory handed
    to the page."""
    loader = coolprop_js.read_text(encoding='utf-8').replace('</script', '<\\/script')
    anchor = 'Module["decrementExceptionRefcount"]=decrementExceptionRefcount;'
    if anchor not in loader:
        raise SystemExit(f'{coolprop_js}: cannot find {anchor!r} to expose the stack helpers')
    loader = loader.replace(anchor, anchor + 'Module["sghpStackSave"]=stackSave;Module["sghpStackRestore"]=stackRestore;', 1)
    return '\n' + re.sub(r'export\s+default\s+Module\s*;?', '', loader) + '\nwindow.__cpFactoryReady(Module);\n'


def render(coolprop_js: Path = CP_DIR / 'coolprop.js') -> str:
    template = (HERE / 'index.template.html').read_text(encoding='utf-8')

    def body(m: re.Match) -> str:
        return ''.join(coolprop_module(coolprop_js) if name == '@coolprop' else (HERE / name).read_text(encoding='utf-8')
                       for name in m.group(1).split())

    return MARKER.sub(body, template)


def build(out_dir: Path = OUT, coolprop_js: Path = CP_DIR / 'coolprop.js') -> Path:
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / 'index.html'
    path.write_text(render(coolprop_js), encoding='utf-8')
    wasm = coolprop_js.parent / 'coolprop.wasm'
    if wasm.resolve() != (out_dir / 'coolprop.wasm').resolve():
        shutil.copyfile(wasm, out_dir / 'coolprop.wasm')
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', type=Path, default=OUT)
    ap.add_argument('--check', action='store_true', help='fail if a rebuild differs from apps/dryer/index.html')
    a = ap.parse_args()
    if a.check:
        with tempfile.TemporaryDirectory() as tmp:
            same = build(Path(tmp)).read_bytes() == (OUT / 'index.html').read_bytes()
        print('apps/dryer/index.html is up to date' if same else 'apps/dryer/index.html differs from its sources')
        sys.exit(0 if same else 1)
    print('->', build(a.out))


if __name__ == '__main__':
    main()
