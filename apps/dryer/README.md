# Advanced dryer heat pumps (web app)

An interactive page for designing and comparing heat pumps that heat drying air. Its "Live
design" mode solves the recuperative heat pump (RHP-A), split-stream heat pump (SLTHP) and
other cycles in the browser. The models are ported from HeatPumpLab, plus published cycles
(Kong et al., Hoang et al., Dalal et al.). Fluid properties come from CoolProp compiled to
WebAssembly. Its "Study" mode shows the results of the REFPROP design study.

The page runs entirely in the browser. It has no server-side code and sends none of your inputs anywhere; the only external
requests are for the web fonts.
Saved cases stay in the browser's local storage, and a case link carries the design inputs in
the URL.

## Files

| File | Contents |
|---|---|
| `index.html` | The page: layout, the cycle models in JavaScript and the CoolProp Emscripten loader |
| `coolprop.wasm` | CoolProp 8 compiled to WebAssembly |
| `dashboard_data.json`, `cycles_data.json`, `state_5C_1bar.json` | REFPROP design study results (results only, no REFPROP code or data files) |
| `THIRD_PARTY_NOTICES.md` | Licences for CoolProp and the web fonts |

`index.html` is built from the sources in [`apps/src/dryer`](../src/dryer) (`uv run python apps/src/dryer/build.py`); edit those, not the page. See [apps/src](../src/README.md).

## Running locally

The page fetches the `.wasm` and `.json` files, so open it through a static web server rather
than as a `file://` URL. From the repository root:

```sh
python3 -m http.server 8000 --directory apps/dryer
```

Then open <http://localhost:8000>. Any static server that serves `.wasm` as
`application/wasm` works. GitHub Pages works the same way.

## Licence

Apache License 2.0; see [LICENSE](../../LICENSE). Third-party components keep their own
licences; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
