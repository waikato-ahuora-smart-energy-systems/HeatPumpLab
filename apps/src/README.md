# Web app sources

The pages in `apps/steam/` and `apps/dryer/` are built from the files here. Edit the sources, rebuild, and commit
the sources and the rebuilt page together. CI rebuilds both apps and fails if a committed page differs from what
its sources produce. GitHub Pages publishes `apps/` without this `src/` folder.

```sh
uv run python apps/src/steam/build.py     # -> apps/steam/index.html and coolprop.wasm
uv run python apps/src/dryer/build.py     # -> apps/dryer/index.html and coolprop.wasm
```

| Folder | Contents |
|---|---|
| `coolprop/` | CoolProp 8 compiled to WebAssembly: the Emscripten loader (`coolprop.js`) and `coolprop.wasm`. Both builds inline the loader and copy the `.wasm` next to the page. CoolProp is MIT-licensed; see the apps' `THIRD_PARTY_NOTICES.md`. |
| `steam/` | Steam-generating heat pump: the cycle model (`sghp_model.js`, a port of `ConventionalHeatPump`), the CoolProp adapter and controls (`sghp_live.js`), saved cases and comparison (`sghp_cases.js`), the page body (`body.html`) and the styles and charts it shares with the Python dashboard (`dashboard_template.html`). `build.py --mock` builds a layout-only page with a test fluid. |
| `dryer/` | Advanced dryer heat pumps: `index.template.html` is the page with each script replaced by an `{{include ...}}` marker, and `build.py` puts the scripts back in order. The scripts are the dryer's own copy of the SGHP model and CoolProp adapter (`sghp_model.js`, `adapter.js`), the cycle models (`dryer_model.js`, `dryer_cascades.js`), the study views (`dryer_study.js`), the live controls (`dryer_live.js`) and saved cases (`dryer_cases.js`). The study results (`*.json`) live in `apps/dryer/` and are not rebuilt. |

The JavaScript models are checked against the Python ones by the tests that load these files with Node
(`tests/test_dryer_live.py`, `tests/test_sghp_*.py` and others); they run in CI.
