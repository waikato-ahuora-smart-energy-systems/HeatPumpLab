# The Heat Pump Lab

[![Tests](https://github.com/waikato-ahuora-smart-energy-systems/HeatPumpLab/actions/workflows/tests.yml/badge.svg)](https://github.com/waikato-ahuora-smart-energy-systems/HeatPumpLab/actions/workflows/tests.yml)
[![Release](https://img.shields.io/github/v/release/waikato-ahuora-smart-energy-systems/HeatPumpLab)](https://github.com/waikato-ahuora-smart-energy-systems/HeatPumpLab/releases/latest)
[![PyPI](https://img.shields.io/pypi/v/heatpumplab)](https://pypi.org/project/heatpumplab/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23604060.svg)](https://doi.org/10.5281/zenodo.23604060)
[![Licence: Apache-2.0](https://img.shields.io/badge/licence-Apache--2.0-blue.svg)](LICENSE)
[![Python 3.12 | 3.13 | 3.14](https://img.shields.io/badge/python-3.12%20%7C%203.13%20%7C%203.14-blue.svg)](https://www.python.org/)
[![Properties: CoolProp | REFPROP](https://img.shields.io/badge/properties-CoolProp%20%7C%20REFPROP-informational.svg)](http://www.coolprop.org)
[![Web apps](https://img.shields.io/badge/web%20apps-live-brightgreen.svg)](https://waikato-ahuora-smart-energy-systems.github.io/HeatPumpLab/)

The Heat Pump Lab, from the Ahuora Centre at the University of Waikato,
is a Python package for building and analysing high-temperature heat pump cycles. It solves
steady-state cycle designs with real-fluid properties, including zeotropic refrigerant mixtures
and transcritical and supercritical operation, and reports COP, pressures, compressor duties and heat-exchanger
approaches.

The package contains seven cycle models. Each is a module under `HeatPumpLab.cycles`:

| Cycle | Module | What it models | Reference |
| --- | --- | --- | --- |
| Recuperative heat pump (RHP) | `RecuperativeHeatPump` | One or two compressors with up to three recuperative heat exchangers (RHX) placed between the gas coolers, for high-temperature glide heating. Includes Bayesian optimisation of the design. | [1] |
| Split-stream heat pump (SLTHP) | `SplitHeatPump` | Two-stage cycle with RHX1 in parallel with the low-temperature gas cooler. | [1] |
| Steam-generating heat pump | `ConventionalHeatPump` | Single-stage rating model with compressor maps, heat losses and steam raising. | [2] |
| Transcritical–transcritical CO₂ cascade | `TTCascadeHeatPump` | The bottom cycle gives its heat to the top cycle's suction gas in a cascade gas cooler; the top cycle heats the air in two gas coolers, with optional secondary air. IHXs in both cycles, optional expanders. | [3–5] |
| n-level cascade | `CascadeHeatPump` | Two to four levels heating the air in series through each level's subcooler and its share of condenser heat. Optional IHX per level and expanders. | [6] |
| Heat pump with MVR steam | `MVRHeatPump` | A vapour-compression cycle heats the air and raises low-pressure steam, which a chain of mechanical vapour recompression fans, with water injection and splits to the process, lifts to heat the air to the target. | [2] |
| sCO₂ reverse Brayton heat pump | `BraytonHeatPump` | Supercritical CO₂ cycle with split compression, a preheater, an IHX and a turbine expander, architectures 1 and 2. | [7] |

## Web apps

Two of the models also run in the browser, with CoolProp compiled to WebAssembly. Nothing is
installed and no inputs leave the browser:

- **Advanced dryer heat pumps** ([`apps/dryer`](apps/dryer)): the air-drying cycles (all but
  the steam-generating heat pump), solved live, plus the results of the REFPROP design study.
- **Steam-generating heat pump** ([`apps/steam`](apps/steam)): the single-stage cycle raising
  steam, with compressor maps, saved cases, side-by-side comparison and sweeps.

A landing page ([`apps/index.html`](apps/index.html)) links the two, and their sources and build scripts are
in [`apps/src`](apps/src/README.md). The apps are published at
<https://waikato-ahuora-smart-energy-systems.github.io/HeatPumpLab/>. To run them locally:

```sh
python3 -m http.server 8000 --directory apps
```

Then open <http://localhost:8000>. The Python models reproduce the apps' results.

## Requirements

- Python 3.12 to 3.14, and [uv](https://docs.astral.sh/uv/) (recommended) or pip.
- Fluid properties come from [CoolProp](http://www.coolprop.org) (free, installed with the
  package) by default, for every cycle. Mixtures use phase-imposed flashes with cached bubble and
  dew points, plus a phase-envelope fallback near the critical point, so zeotropic blends solve
  robustly and quickly.
- **NIST REFPROP 10** (optional) can be used instead. REFPROP is commercial software and is not
  included; you need your own licence. HeatPumpLab looks for it at:
  - macOS: `/Library/REFPROP/` (`librefprop.dylib`, `librefprop-source.dylib`, `librefprop-sink.dylib`)
  - Windows: `C:/Program Files (x86)/REFPROP/` (`REFPROP.DLL`, `REFPROP-source.DLL`, `REFPROP-sink.DLL`)

  The three library files are copies of the same REFPROP library, so the refrigerant, source and
  sink fluids can be loaded at the same time.

## Installation

From [PyPI](https://pypi.org/project/heatpumplab/):

```sh
pip install heatpumplab                 # or: uv add heatpumplab
pip install "heatpumplab[refprop]"      # optional: the REFPROP backend (needs a REFPROP installation)
```

From source, to run the tests, the research scripts or the latest development version:

```sh
git clone https://github.com/waikato-ahuora-smart-energy-systems/HeatPumpLab.git
cd HeatPumpLab
uv sync                      # core dependencies (including CoolProp), and the package itself
uv sync --extra refprop      # optional: the REFPROP backend (needs a REFPROP installation)
uv sync --group dev          # optional: pytest and ruff
```

To use REFPROP, set `HEATPUMPLAB_PROPERTY_BACKEND=refprop`, or pass `property_backend="refprop"`
to the steam-generating heat pump's configuration (`ConventionalHPConfig`). The two backends
agree closely; mixture results can differ slightly where their binary interaction parameters
differ, and absolute enthalpies and entropies of some mixtures differ by a constant (reference
state), which does not affect cycle results. The fluid property overview table
(`get_fluid_prop_overview_table`) uses REFPROP-only outputs such as ODP and GWP.

Run the examples below with `python your_script.py`, or `uv run python your_script.py` in a source checkout.

## Quick start

### Recuperative heat pump

```python
from HeatPumpLab.cycles.RecuperativeHeatPump import RecuperativeHP

hp = RecuperativeHP(
    refrigerant="Propane;Isopentane|0.2;0.8|mass",  # components; fractions; basis
    T_source_in=40.0,       # heat source (water) inlet, °C
    Tref_evap_inlet=10.0,   # refrigerant entering the evaporator, °C
    T_air_inlet=15.0,       # sink (air) inlet, °C
    T_air_tar=200.0,        # sink outlet target, °C
    dT_GC_min_tar=20.0,     # gas-cooler minimum approach, K
    dTmin_rhx=5.0,          # recuperator minimum approach, K
    n_comp=2,               # compression stages
    dT_rhx1_cold=60.0,      # RHX1 lift on the stage-1 suction, K
    dT_rhx2_cold=0.0,       # RHX2 lift before stage 2, K
    comp2_power_frac=0.5,   # pressure-ratio allocation to stage 2
    sink="Air",
    source="water",
).solve()

print(hp.converged, hp.COP)              # True 2.59
print(hp.get_report()["Pmax"])           # maximum pressure, kPa
```

### Design optimisation

`optimise_recuperative_heat_pump` searches the RHX lifts, stage allocation and, optionally, the
mixture composition with Bayesian optimisation (BoTorch). Solves run in worker processes, so
call it from a script with a main guard:

```python
from HeatPumpLab.cycles.RecuperativeHeatPump.optimisation import (
    OuterDesignOptimisationConfig, optimise_recuperative_heat_pump,
)

if __name__ == "__main__":
    result = optimise_recuperative_heat_pump(
        refrigerant="Propane;Butane|0.3;0.7|mass",
        config=OuterDesignOptimisationConfig(
            max_exact_evaluations=200, seed=20260910,
            vary_composition=False, preserve_incumbent=False,
        ),
    )
    print(result.best.cop, result.best.x)
```

A 200-evaluation search takes roughly ten to fifteen minutes on a laptop. See the
[optimisation protocol](HeatPumpLab/cycles/RecuperativeHeatPump/BAYESIAN_OPTIMISATION.md) for
the settings and their limits.

### Split-stream heat pump

```python
from HeatPumpLab.cycles.SplitHeatPump import SplitHeatPump, SplitHPConfig

res = SplitHeatPump(SplitHPConfig(
    refrigerant="Propane;Isopentane|0.2;0.8|mass",
    dT_rhx1_cold=80.0,          # RHX1 lift, K
)).solve()
print(res.success, res.COP, res.P_H)    # gas-cooler pressure, MPa
```

### Steam-generating heat pump

```python
from HeatPumpLab.cycles.ConventionalHeatPump import ConventionalHeatPump, ConventionalHPConfig

res = ConventionalHeatPump(ConventionalHPConfig()).solve()       # n-butane, Recip-M compressor map
res = ConventionalHeatPump(ConventionalHPConfig(T_evap=40.0, T_cond=110.0)).solve()
print(res.COP, res.report["Q_heat_kW"], res.report["T_discharge"])

# REFPROP instead of CoolProp (needs `uv sync --extra refprop` and a REFPROP installation)
res = ConventionalHeatPump(ConventionalHPConfig(property_backend="refprop")).solve()
```

The defaults model n-butane on the Recip-M reciprocating compressor map, raising 120 °C steam from a 55 °C
heat source, with hot water from the subcooler. Compressors can be described by efficiency maps,
efficiency against pressure ratio, or EN 12900 polynomials (`compressor_curves.py`).

## Refrigerants

Fluids use REFPROP-style names (`Propane`, `Butane`, `Isopentane`, ...), which are mapped to
CoolProp's names automatically. A mixture is written as `"Fluid1;Fluid2|x1;x2|basis"` with basis `mass`
or `mol`, for example `"Propane;Isopentane|0.2;0.8|mass"`. Ternary mixtures use the same form.

## Repository layout

```
HeatPumpLab/
  components/      compressors, heat exchangers, valves, streams and solver relationships
  cycles/          one module per cycle (see the table above), and shared helpers in common/
  utils/           property packages (CoolProp, REFPROP), caching, parallel solving, plotting
apps/              web apps: landing page, dryer/ and steam/ (published with GitHub Pages)
  src/             the apps' JavaScript sources and build scripts
docs/dev-notes/    development notes (recuperative heat pump solver and optimiser)
tests/             test suite
scripts/           research search scripts for the recuperative heat pump (Pareto, MOBO, trust region)
```

## Tests

```sh
uv sync --group dev
uv run pytest
```

Tests that check the REFPROP backend skip automatically when it is not installed.

## Further documentation

- [ABOUT.md](HeatPumpLab/cycles/RecuperativeHeatPump/ABOUT.md): the recuperative heat pump cycle variants
- [BAYESIAN_OPTIMISATION.md](HeatPumpLab/cycles/RecuperativeHeatPump/BAYESIAN_OPTIMISATION.md): the production optimiser
- [`docs/dev-notes/rhp/`](docs/dev-notes/rhp/README.md): development notes on the recuperative heat pump solver
  and optimiser, including the optimiser benchmarks
  ([OPTIMISER_COMPARISON.md](docs/dev-notes/rhp/OPTIMISER_COMPARISON.md),
  [ORIGIN_OPTIMISER_COMPARISON.md](docs/dev-notes/rhp/ORIGIN_OPTIMISER_COMPARISON.md)) and the solver fixes by
  evaluation version
- [`apps/src/`](apps/src/README.md): the web apps' sources and how to rebuild them

## References

1. T.G. Walmsley, K.-M. Adamson, C. Arpagaus, F. Schlosser, Recuperative heat pump cycle concept
   for high-temperature glide process heating applications, *Applied Thermal Engineering* 301
   (2026) 131660. https://doi.org/10.1016/j.applthermaleng.2026.131660
2. T.G. Walmsley, B. Ong, K. Hall, B. Lincoln, F. Schlosser, J. Walden, R. Padullés, D. Hoang,
   M. Walmsley, C. Arpagaus, Steam cogeneration heat pumps for the food, dairy and meat industries,
   in: 5th High-Temperature Heat Pump Symposium 2026, Book of Presentations, Copenhagen, Denmark,
   February 2026, pp. 270–275.
3. L. Kong, New high temperature transcritical CO₂ heat pump concept for spray drying, PhD thesis,
   The University of Waikato (2026). https://hdl.handle.net/10289/18056 (the TT cascade model
   follows cycle TT-4, with the additional gas cooler)
4. L. Kong, T.G. Walmsley, D.K. Hoang, F. Schlosser, Q. Chen, J.K. Carson, D.J. Cleland,
   Transcritical-transcritical cascade CO₂ heat pump cycles for high-temperature heating: A
   numerical evaluation, *Applied Thermal Engineering* 238 (2024) 122005.
   https://doi.org/10.1016/j.applthermaleng.2023.122005
5. L. Kong, F. Schlosser, S. Kloeppel, J.K. Carson, D.J. Cleland, T.G. Walmsley,
   Transcritical–transcritical cascade CO₂ heat pump with expansion work recovery: A
   thermodynamic analysis, *Energies* 19 (2026) 3767. https://doi.org/10.3390/en19163767
6. D.K. Hoang, T.G. Walmsley, D.J. Cleland, Q. Chen, J.K. Carson, Thermo-economic investigation
   and multi objective optimization of cascade high temperature heat pump using low global warming
   refrigerants, *Applied Thermal Engineering* 279 (2025) 127901.
   https://doi.org/10.1016/j.applthermaleng.2025.127901
7. V. Dalal, S. Klöppel, C. Rose, T.G. Walmsley, Evaluating sCO₂ heat pump architectures for
   large-scale industrial air drying above 150 °C, European sCO₂ Conference (Euro sCO₂ 2027),
   under review.

## Citing

If you use HeatPumpLab in published work, please cite reference [1], and the reference for any
other cycle you use. To cite the software itself, use its Zenodo DOI,
[10.5281/zenodo.23604060](https://doi.org/10.5281/zenodo.23604060), which always resolves to the
latest release. GitHub's "Cite this repository" button gives both, from `CITATION.cff`.

## Contributing

Bug reports, questions and contributions are welcome; see [CONTRIBUTING.md](CONTRIBUTING.md).

## Licence

Copyright 2021–2026 The University of Waikato.

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE).

## Contact

Tim Walmsley, Ahuora Centre, University of Waikato: tim.walmsley@waikato.ac.nz
