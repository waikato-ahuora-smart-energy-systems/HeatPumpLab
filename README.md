# The Heat Pump Lab

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
| Conventional heat pump | `ConventionalHeatPump` | Single-stage rating model with compressor maps, heat losses and steam raising (the steam-generating heat pump). | — |
| Transcritical–transcritical CO₂ cascade | `TTCascadeHeatPump` | The bottom cycle gives its heat to the top cycle's suction gas in a cascade gas cooler; the top cycle heats the air in two gas coolers, with optional secondary air. IHXs in both cycles, optional expanders. | [2] |
| n-level cascade | `CascadeHeatPump` | Two to four levels heating the air in series through each level's subcooler and its share of condenser heat. Optional IHX per level and expanders. | [3] |
| Heat pump with MVR steam | `MVRHeatPump` | A vapour-compression cycle heats the air and raises low-pressure steam, which a chain of mechanical vapour recompression fans, with water injection and splits to the process, lifts to heat the air to the target. | — |
| sCO₂ reverse Brayton heat pump | `BraytonHeatPump` | Supercritical CO₂ cycle with split compression, a preheater, an IHX and a turbine expander, architectures 1 and 2. | [4] |

The air-drying cycles (all but the conventional one) are also available in the Advanced dryer heat
pumps web app, and the conventional cycle in the Steam-generating heat pump web app, from the Ahuora
Centre; the Python models reproduce the apps' results.

## Background

Internal heat exchangers (IHXs) are widely used in vapour-compression and transcritical heat
pumps. They transfer heat from the warm refrigerant just before expansion to the suction gas.
Pinch analysis suggests a heat pump should absorb heat below the pinch and reject upgraded heat
above it. Applied to the cycle itself, this means recuperators placed between the gas coolers
can raise the compressor discharge temperature to the level a high-temperature process needs,
at a lower pressure. To distinguish them from a conventional IHX, the model calls these
recuperative heat exchangers (RHX). The concept and its results are described in:

> T.G. Walmsley, K.-M. Adamson, C. Arpagaus, F. Schlosser, Recuperative heat pump cycle concept
> for high-temperature glide process heating applications, *Applied Thermal Engineering* 301
> (2026) 131660. https://doi.org/10.1016/j.applthermaleng.2026.131660

See [ABOUT.md](HeatPumpLab/cycles/RecuperativeHeatPump/ABOUT.md) for the cycle variants.

## Requirements

- Python 3.12 or later, and [uv](https://docs.astral.sh/uv/) (recommended) or pip.
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

```sh
git clone https://github.com/waikato-ahuora-smart-energy-systems/HeatPumpLab.git
cd HeatPumpLab
uv sync                      # core dependencies (including CoolProp), and the package itself
uv sync --extra refprop      # optional: the REFPROP backend (needs a REFPROP installation)
uv sync --group dev          # optional: pytest and ruff
```

To use REFPROP, set `HEATPUMPLAB_PROPERTY_BACKEND=refprop`, or pass `property_backend="refprop"`
to the conventional cycle's configuration. The two backends agree closely; mixture results can
differ slightly where their binary interaction parameters differ, and absolute enthalpies and
entropies of some mixtures differ by a constant (reference state), which does not affect cycle
results. The fluid property overview table (`get_fluid_prop_overview_table`) uses REFPROP-only
outputs such as ODP and GWP.

Run the examples below with `uv run python your_script.py` from the repository root.

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

### Conventional heat pump

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
tests/             test suite
rhp_*.py           research search scripts for the recuperative heat pump (Pareto, MOBO, trust region)
```

## Tests

```sh
uv sync --group dev
uv run pytest
```

Tests that check the REFPROP backend skip automatically when it is not installed.

## Further documentation

Detailed notes on the recuperative heat pump solver and optimiser are in
[`HeatPumpLab/cycles/RecuperativeHeatPump/`](HeatPumpLab/cycles/RecuperativeHeatPump/), including:

- [ABOUT.md](HeatPumpLab/cycles/RecuperativeHeatPump/ABOUT.md): the cycle variants
- [BAYESIAN_OPTIMISATION.md](HeatPumpLab/cycles/RecuperativeHeatPump/BAYESIAN_OPTIMISATION.md): the production optimiser
- [OPTIMISER_COMPARISON.md](HeatPumpLab/cycles/RecuperativeHeatPump/OPTIMISER_COMPARISON.md) and
  [ORIGIN_OPTIMISER_COMPARISON.md](HeatPumpLab/cycles/RecuperativeHeatPump/ORIGIN_OPTIMISER_COMPARISON.md): optimiser benchmarks
- `V14_SOLVER_FIXES.md`, `V15_SOLVER_FIXES.md` and the other investigation notes: solver development history

## References

1. T.G. Walmsley, K.-M. Adamson, C. Arpagaus, F. Schlosser, Recuperative heat pump cycle concept
   for high-temperature glide process heating applications, *Applied Thermal Engineering* 301
   (2026) 131660. https://doi.org/10.1016/j.applthermaleng.2026.131660
2. Kong (2026), PhD thesis, University of Waikato: cycle TT-4 with the additional gas cooler.
   [full reference to be added]
3. Hoang et al. [full reference to be added]
4. V. Dalal, S. Klöppel, C. Rose, T.G. Walmsley, Evaluating sCO₂ heat pump architectures for
   large-scale industrial air drying above 150 °C, European sCO₂ Conference (Euro sCO₂ 2027).
   [details to be confirmed]

## Citing

If you use HeatPumpLab in published work, please cite reference [1], and the reference for any
other cycle you use.

## Licence

Copyright 2021–2026 The University of Waikato.

Licensed under the Apache License, Version 2.0. See [LICENSE](LICENSE).

## Contact

Tim Walmsley, Ahuora Centre, University of Waikato: tim.walmsley@waikato.ac.nz
