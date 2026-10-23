# v17 internal PH recovery

The validated nearby-state fallback is now implemented inside `PropertyPackage._flash`. It recovers **18 of the 25 cycles still failing under v15**, while retaining all 19 earlier recoveries: **37 of the original 44 failed candidates now solve**. Seven candidates still return explicit PH error 226.

## Internal call behaviour

Both `get_std_properties("PH", ...)` and `get_std_properties("HP", ...)` reach this path before unit conversion and cache insertion. For a native binary-mixture error 226, the existing saturation-bounded PT recovery runs first. If that fails, the new fallback searches nearby native PH states at fixed pressure or fixed enthalpy. No optimiser wrapper or diagnostic monkeypatch is required.

The search uses at most **24 native PH probes**: seven enthalpy distances in each direction, from 10⁻⁵ to 5 kJ/kg, and five relative pressure distances in each direction, from 10⁻⁶ to 0.01. Probe calls bypass `_flash` to prevent recursive recovery. Only native success with finite, ordered two-phase states is considered; endpoints must also pass independent phase validation, with the existing local precision correction when needed.

Validated states on opposite sides of the requested point form a bracket. Unequal-distance pairs are supported; no extrapolation is used. Temperature, phase densities, compositions, quality and thermodynamic properties are interpolated, and bulk density is reconstructed from phase specific volumes. A candidate is accepted only after the existing pressure, fugacity, material, density and enthalpy checks pass. Otherwise it seeds the local binary equilibrium correction at the requested pressure and enthalpy, with unchanged close-state guards and tolerances. At most 74 bracket corrections are attempted, in addition to endpoint precision corrections.

The fallback never uses property values from a failed native flash. Exhausted recovery propagates the original error and its property context. Failed results are not cached. Successes enter the ordinary property cache, and `refprop_diagnostics` records whether the nearby search succeeded or exhausted its brackets. Profiling records the search under `refprop_ph_nearby_recovery`.

Successful native flashes, pure-fluid recovery and unrelated error codes retain their existing paths. This new fallback is limited to actual two-component mixtures. The optimiser version is now `isolated-refprop-pattern-search-v17`, so obsolete completed evaluations and versioned composition replays are not reused.

## Validation

The production implementation was tested in **82 fresh-process evaluations**, using the same requests as the frozen v15 validation. Each complete evaluation had a 200-second diagnostic deadline.

| Group | Evaluations | Successful | Timeouts |
|---|---:|---:|---:|
| Original failed candidates | 44 | 37 | 0 |
| Previous benchmark winners | 18 | 18 | 0 |
| Sampled successful controls | 12 | 12 | 0 |
| Near-zero probes | 8 | 8 | 0 |

All **57 previously successful evaluations remain successful**. The 18 winners, 12 controls and eight near-zero probes have exactly unchanged SEC, COP and discharge pressure in this replay. Among the 19 earlier recovered candidates, the largest relative COP change is 2.94 × 10⁻⁷; all checked comparisons remain below 10⁻⁶ relative tolerance.

The requested total-VHC definition is retained: heat duty divided by the sum of active compressor suction volume flows. Every successful replay verifies that formula directly. Historical v15 VHC values used only the first suction flow, so regression comparisons reconstruct that historical denominator instead of comparing unlike definitions. That reconstructed comparison also remains within 10⁻⁶ relative tolerance.

The full test suite passes: **564 passed**, with one existing PyTorch deprecation warning. Tests exercise internal PH/HP routing, cache insertion after successful recovery, rejection of failed native output, unequal brackets, finite search bounds, independent validation, unsupported inputs and four real REFPROP regressions. A profiling-test fixture was also updated to supply the compressor suction flow required by the total-VHC definition. Ruff and whitespace checks pass for the changes.

## Remaining failures and cost

Seven complete cycles still exhaust validated recovery and return PH error 226. Recovering a representative property state does not guarantee that every later property call in its cycle will recover. These failures remain numerical failures rather than a claim of physical infeasibility.

Recovery increases effort on difficult cycles. The slowest failed replay took about **164.7 seconds**, and the slowest successful replay about **113.0 seconds**, under concurrent instrumented execution. All stayed below the diagnostic deadline; the normal isolated-objective timeout remains in force. No new 400-evaluation optimiser comparison was run.

## Implementation and evidence

- [Internal properties implementation](../../utils/properties_package.py)
- [Fallback regression tests](../../../tests/test_ph_nearby_recovery.py)
- [Numerical summary](../../../results/RHPcomposition/v17_internal_ph_fallback/summary.json)
- [All original candidate outcomes](../../../results/RHPcomposition/v17_internal_ph_fallback/candidate_outcomes.md)
- [Verified full test run](../../../results/RHPcomposition/v17_internal_ph_fallback/full_tests_verified.log)
- [Frozen source manifest](../../../results/RHPcomposition/v17_internal_ph_fallback/source_manifest.json)
- [Replay audit script](../../../main_RHPv17_report.py)

Raw requests, results, traces, status records and a source snapshot are retained alongside the summary. The final thermodynamic source SHA-256 is `505a78158f12993ced8924499bb1710984bd33bc0ef82842628877bed43200c8`. Rebuild the numerical audit with `python main_RHPv17_report.py`. The earlier `full_tests.log` records the superseded test run before correcting the VHC profiling fixture; `full_tests_verified.log` is the final passing run.
