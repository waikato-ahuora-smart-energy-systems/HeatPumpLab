# Steam-generating heat pump (web app)

An interactive rating model for a single-stage heat pump that raises steam from low-grade waste
heat. The cycle has an internal heat exchanger and a hot-water subcooler, and the condenser
raises steam. Set the source, steam and hot-water temperatures, the refrigerant (pure or up to
three-component mixtures), the compressor (three reciprocating R600a efficiency maps, Recip-S,
Recip-M and Recip-L) and its speed, and the cycle re-solves in the browser. Cases can be saved,
compared side by side and swept.

The model is a JavaScript port of HeatPumpLab's `ConventionalHeatPump` cycle, checked against
the Python model. Fluid properties come from CoolProp compiled to WebAssembly.

The page runs entirely in the browser. It has no server-side code and sends none of your inputs anywhere; the only external
requests are for the web fonts. Saved cases stay in the browser's local storage, and a case link
carries the design inputs in the URL.

## Files

| File | Contents |
|---|---|
| `index.html` | The page: layout, the cycle model in JavaScript and the CoolProp Emscripten loader |
| `coolprop.wasm` | CoolProp 8 compiled to WebAssembly |
| `THIRD_PARTY_NOTICES.md` | Licences for CoolProp and the web fonts |

## Running locally

The page fetches `coolprop.wasm`, so open it through a static web server rather than as a
`file://` URL. From the repository root:

```sh
python3 -m http.server 8000 --directory apps/steam
```

Then open <http://localhost:8000>. Any static server that serves `.wasm` as
`application/wasm` works. GitHub Pages works the same way.

## Licence

Apache License 2.0; see [LICENSE](../../LICENSE). Third-party components keep their own
licences; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
