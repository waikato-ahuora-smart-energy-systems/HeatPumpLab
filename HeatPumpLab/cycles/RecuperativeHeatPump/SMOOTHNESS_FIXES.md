# Smoothness fixes in evaluation version 13

The production solver now initializes REFPROP mixture phase envelopes, refines
single-phase PH calculations, and interpolates the mixture discharge-pressure
target. These changes address the near-zero RHX2 failure region and the six small
COP steps found by the random smoothness audit.

## Property calculations

`HeatPumpLab/utils/properties_package.py` owns fluid/composition preparation for
each shared native library. Switching a package between the default, source and
sink libraries restores its fluid and phase envelope. Switching compositions
rebuilds the envelope; a failed setup invalidates the previous envelope too.
Packages sharing a configured library reuse its setup. Native REFPROP calls still
require the existing process isolation when solving concurrently.

The high-level spline flag is enabled. Low-level flashes retain flag zero, so
their phase identification and molar units are unchanged. Positive native errors,
uncertain-phase warning -235, and nonfinite/invalid standard states are rejected
before entering the property cache. Native warnings, errors, precision-refinement
declines and successful recoveries remain available in each package's
`refprop_diagnostics` dictionary.

A successful single-phase mixture PH flash is checked against PT at the returned
temperature. Up to four bounded Newton corrections use the PT heat capacity to
match the requested enthalpy. A two-phase result, warning, failure, or correction
greater than 0.01 K cancels refinement and retains the validated original flash.
The precision policy is included in the property cache key. This corrects the
microscopic RHX2 temperature-direction failure without relaxing exchanger checks.

Stricter error handling also exposed a spurious solid-state PH root for ordinary
liquid water. A failed pure-liquid PH calculation can now recover using PT root
finding strictly between checked melting/model and boiling limits. The final
enthalpy residual must meet its tolerance. Actual solid-state inputs remain
rejected. This recovery does not force a phase or reuse the failed state's
temperature, density or enthalpy.

## Pressure search

`solver/_discharge_pressure.py` interpolates the lower boundary of the existing
accepted pinch band with Brent's method in logarithmic pressure. For the usual
20 K target and 0.05 K band, it targets 19.95001 K; the small interior margin
supports fresh final-state validation. The original final pinch validation,
lowest-pressure branch selection and fallback searches remain in place.
Interpolation failures are counted in `hp.pressure_interpolation_fallbacks`.

This removes the coarse pressure grid responsible for the observed local COP
steps. Accepted-endpoint and unbracketed fallback paths are retained, so this is
not a claim of global differentiability across all feasible/infeasible boundaries.

## Regression changes

- Added coverage for library and composition switches, failed setup invalidation,
  rejected/non-cached native errors, preservation of warnings, PH precision,
  genuine two-phase states, liquid-water recovery, solid-state rejection, and
  near-zero RHX2 cycle solves.
- Added an analytic pressure-continuity regression and interpolation-fallback
  coverage.
- Replaced the old Isopentane/Neopentane expected-temperature-cross regression
  with 501-point exchanger-profile checks. The old spike came from REFPROP error
  248 and an enthalpy error of approximately 334 kJ/kg.
- Corrected two standalone exchanger fixtures' air flow rates so their intended
  temperature limits are physically attainable with actual butane properties.
  They previously relied on the wrong fluid remaining loaded in a shared library.
- Evaluation signatures and default output filenames advance to v13. Completed
  v11 and v12 objective results are not reused.

## Remaining native failure surfaced by validation

For Propane/Isobutane 30/70 mass, RHX2 = 28.908922096239987 K and compressor-2
fraction = 0.7116845176396405, the samples at RHX1 = 74.79577898468112 K and
74.80577898468113 K encounter REFPROP two-phase PH error 226 while computing a
temperature profile. The failed PH calculation at 3.946390399004728 MPa and
34257.983530203266 J/mol returned an incorrect enthalpy. Alternative PQ probes
also showed convergence problems/inconsistency; they were not adopted as a
production recovery.

These two evaluations now report a property-solver failure rather than silently
using the returned state. Their physical feasibility is not established by that
failure. They are distinct from the six pressure steps and the near-zero RHX2 gap
that this change fixes.

The reproducible production replay is `paper/analysis/validate_v13_smoothness.py`.
Its inputs, individual results, source hash and summary are saved in
`results/RHPcomposition/smoothness_v13_validation`. The summary explicitly lists
native failures separately from successful physical states.

## Completed validation

- Full test suite: **486 passed**. Ruff and `git diff --check` passed.
- Fresh-process production replay: **312 successes from 314 evaluations**.
  The two rejections are the error-226 points described above.
- All **14 RHX2 gap probes** succeeded, including RHX2 = 0 and 1e-8 K.
- All **12 endpoints of the six narrow brackets** succeeded. The largest COP
  difference across those roughly 9.765625e-6 K intervals is **0.000002406%**,
  down from observed steps of roughly 0.0086–0.0167% before the fixes.
- All successful replays passed the existing final gas-cooler pinch check.
  The production Python source hash was unchanged throughout the final replay.
