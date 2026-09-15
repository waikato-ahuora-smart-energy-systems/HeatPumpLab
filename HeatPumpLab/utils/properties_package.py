"""Wrapper around REFPROP for thermodynamic property evaluation."""

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from collections.abc import Callable
from math import isfinite, log
from numbers import Real
from sys import platform
import numpy as np
from scipy.optimize import root, root_scalar
from scipy.special import expit, logit
from .property_cache import PropertyCacheKey, active_property_call_cache
from .solve_profile import profile_section


class RefpropError(ValueError):
    """A REFPROP operation failed or returned an unusable state."""


def phase_from_refprop_quality(quality) -> str:
    """Interpret ABFLSH's molar quality, including its single-phase flags.

    REFPROP uses q < 0 for liquid, q > 1 for superheated vapour,
    998 for vapour with undefined quality and 999 for supercritical fluid.
    Missing values and the wrapper's historical sentinel are not phase evidence.
    """
    if not isinstance(quality, Real) or not isfinite(quality):
        return "unknown"
    if quality == 999:
        return "supercritical"
    if 1 < quality <= 998:
        return "superheated_vapour"
    if quality == 1:
        return "saturated_vapour"
    if 0 < quality < 1:
        return "two_phase"
    if quality == 0:
        return "saturated_liquid"
    if -998 <= quality < 0:
        return "liquid"
    return "unknown"


COMPRESSOR_GAS_PHASES = frozenset(
    {"saturated_vapour", "superheated_vapour", "supercritical"}
)


def _get_refprop_paths(platform_name: str) -> tuple[str, str, str, str]:
    """Return REFPROP library paths and root path for the current platform."""
    if platform_name == "darwin":
        refprop_root = r"/Library/REFPROP/"
        return (
            r"/Library/REFPROP/librefprop.dylib",
            r"/Library/REFPROP/librefprop-source.dylib",
            r"/Library/REFPROP/librefprop-sink.dylib",
            refprop_root,
        )

    refprop_root = "C:/Program Files (x86)/REFPROP"
    return (
        "C:/Program Files (x86)/REFPROP/REFPROP.DLL",
        "C:/Program Files (x86)/REFPROP/REFPROP-source.DLL",
        "C:/Program Files (x86)/REFPROP/REFPROP-sink.DLL",
        refprop_root,
    )


def _load_refprop_libraries(
    platform_name: str,
    library_factory=REFPROPFunctionLibrary,
) -> tuple[REFPROPFunctionLibrary, REFPROPFunctionLibrary, REFPROPFunctionLibrary]:
    """Load REFPROP libraries for default, source, and sink property contexts."""
    ref_path, src_path, snk_path, refprop_root = _get_refprop_paths(platform_name)
    rp_ref = library_factory(ref_path)
    rp_src = library_factory(src_path)
    rp_snk = library_factory(snk_path)
    for rp in rp_ref, rp_snk, rp_src:
        rp.SETPATHdll(refprop_root)
    return rp_ref, rp_src, rp_snk


def _normalise_cache_value(value):
    """Return a stable, hashable representation for property cache keys."""
    if isinstance(value, Real):
        return float(value)
    if isinstance(value, np.ndarray):
        return tuple(_normalise_cache_value(item) for item in value.tolist())
    if isinstance(value, list | tuple):
        return tuple(_normalise_cache_value(item) for item in value)
    return value


def _copy_property_result(result):
    """Return a caller-owned copy of a cached property result."""
    if isinstance(result, dict):
        return dict(result)
    if isinstance(result, list):
        return list(result)
    if isinstance(result, tuple):
        return tuple(result)
    return result


def _get_cached_property_result(key: PropertyCacheKey, get_result: Callable):
    """Return a cached property result or calculate and cache a miss."""
    cache = active_property_call_cache()
    if cache is None:
        return get_result()

    hit, result = cache.lookup(key)
    if hit:
        return _copy_property_result(result)

    result = get_result()
    cache.store(key, _copy_property_result(result))
    return _copy_property_result(result)


class PropertyPackage:
    """Encapsulate REFPROP library access for a configured working fluid."""

    RP: REFPROPFunctionLibrary
    RP_ref: REFPROPFunctionLibrary
    RP_src: REFPROPFunctionLibrary
    RP_snk: REFPROPFunctionLibrary
    RP_ref, RP_src, RP_snk = _load_refprop_libraries(platform)

    # REFPROPdll's spline flag, NOT ABFLSHdll's phase/units flag.
    iFlag: int = 1
    polish_ph: bool = True
    _composition_ready: bool = False
    # Native contexts are shared by streams. Switching fluids or composition
    # invalidates the phase envelope for every package using that library.
    _library_contexts: dict[object, tuple[str, tuple[float, ...] | None]] = {}
    UNITS: int = RP_ref.GETENUMdll(0, "SI with C").iEnum
    iMass: int = 0
    P_crit: float
    P_max: float
    T_crit: float
    D_crit: float
    molar_mass: float
    fld_name: str
    fld_composition: list[float]
    composition_basis: str
    PP_fld_names: list[str]
    n_components: int = 1
    TC_to_K: float = 273.15
    DEFAULT_VAL: float = -9999990.0

    def __init__(self, fld_name: str = None, pp_instance: str = "") -> None:
        """Initialise an instance bound to a particular fluid definition."""
        (
            self.fld_name,
            self.fld_composition,
            self.n_components,
            self.composition_basis,
        ) = _get_fluid_info(fld_name)
        self.iMass = 0
        self._set_RP_lib(pp_instance)
        self.fld_composition = self._get_mole_fractions(self.fld_composition)
        self._composition_ready = True
        self._ensure_refprop_context(self.fld_composition)
        self.molar_mass, self.P_crit, self.T_crit, self.D_crit = (
            self._get_critical_properties()
        )
        # LIMITSdll reports pressure in kPa; HeatPumpLab stores pressure in MPa.
        self.P_max = self.RP.LIMITSdll("EOS", self.fld_composition).Pmax / 1000
        self.PP_fld_names = self._get_database_fluid_name()

    def get_specified_properties(
        self, propsOut: str, propsIn: str, prop1: float, prop2: float
    ):
        """Return selected properties for a given input pair."""
        key = self._get_cache_key(
            "specified", propsOut, propsIn, float(prop1), float(prop2)
        )
        return _get_cached_property_result(
            key,
            lambda: self._get_specified_properties_uncached(
                propsOut, propsIn, prop1, prop2
            ),
        )

    def get_std_properties(
        self, propsIn_str: str, prop1: float, prop2: float, comp: list[float]
    ):
        """Convenience wrapper that converts to/from REFPROP base units."""
        key = self._get_cache_key(
            "std",
            propsIn_str,
            float(prop1),
            float(prop2),
            tuple(float(value) for value in comp),
        )
        return _get_cached_property_result(
            key,
            lambda: self._get_std_properties_uncached(
                propsIn_str, prop1, prop2, comp
            ),
        )

    def _get_cache_key(self, kind: str, *values) -> PropertyCacheKey:
        """Return a cache key for a configured REFPROP property request."""
        return (
            kind,
            id(self.RP),
            self.fld_name,
            tuple(float(value) for value in self.fld_composition),
            self.UNITS,
            self.iMass,
            self.iFlag,
            self.polish_ph,
            float(self.molar_mass),
            float(self.TC_to_K),
            float(self.DEFAULT_VAL),
            *[_normalise_cache_value(value) for value in values],
        )

    def _get_specified_properties_uncached(
        self, propsOut: str, propsIn: str, prop1: float, prop2: float
    ):
        """Return selected properties without consulting the active cache."""
        self._ensure_refprop_context(self.fld_composition)
        with profile_section("refprop_specified"):
            result = self.RP.REFPROPdll(
                self.fld_name,
                propsIn,
                propsOut,
                self.UNITS,
                self.iMass,
                self.iFlag,
                float(prop1),
                float(prop2),
                self.fld_composition,
            )
        self._check_result(result, "REFPROPdll")
        if not all(isfinite(value) for value in result.Output[:len(propsOut.split(";"))]):
            raise RefpropError(f"REFPROPdll returned nonfinite {propsOut}.")
        return result.Output

    def _get_std_properties_uncached(
        self, propsIn_str: str, prop1: float, prop2: float, comp: list[float]
    ):
        """Convenience wrapper without consulting the active cache."""
        self._ensure_refprop_context(comp)
        molar_props = self._convert_props_from_SI_with_C(propsIn_str, prop1, prop2)
        try:
            refprop_result = self._flash(propsIn_str, *molar_props, comp)
        except RefpropError as exc:
            exc.property_context = {"fluid": self.fld_name, "input_pair": propsIn_str,
                                    "inputs": [float(prop1), float(prop2)],
                                    "composition": list(comp),
                                    "units": "model units: MPa, Celsius, kJ/kg"}
            raise
        if (self.polish_ph and propsIn_str in ("PH", "HP")
                and refprop_result.ierr == 0
                and sum(value > 0 for value in comp) > 1):
            pressure, enthalpy = molar_props if propsIn_str == "PH" else molar_props[::-1]
            if not 0 <= refprop_result.q <= 1:
                refprop_result = self._polish_ph(refprop_result, pressure, enthalpy, comp)
            elif (len(comp) == 2 and 0 < refprop_result.q < 1
                  and abs(refprop_result.h-enthalpy) > 1e-6):
                refined = self._refine_binary_ph(refprop_result, pressure, enthalpy, comp)
                if refined is not None:
                    refprop_result = refined
        return self._convert_props_to_SI_with_C(refprop_result)

    def _flash(self, pair, first, second, comp):
        """Check every native flash before its result can reach the cache."""
        with profile_section("refprop_std"):
            result = self.RP.ABFLSHdll(pair, first, second, comp, 0)
        try:
            self._check_result(result, f"ABFLSHdll({pair})")
        except RefpropError:
            if (pair in ("PH", "HP") and result.ierr in (16, 248)
                    and sum(value > 0 for value in comp) == 1):
                pressure, enthalpy = (first, second) if pair == "PH" else (second, first)
                recovered = self._recover_liquid_ph(pressure, enthalpy, comp)
                if recovered is not None:
                    return recovered
            if (pair in ("PH", "HP") and result.ierr == 226
                    and sum(value > 0 for value in comp) > 1):
                pressure, enthalpy = (first, second) if pair == "PH" else (second, first)
                recovered = self._recover_two_phase_ph(pressure, enthalpy, comp)
                if recovered is None:
                    recovered = self._recover_nearby_binary_ph(pressure, enthalpy, comp)
                if recovered is not None:
                    return recovered
            raise
        if (not all(isfinite(getattr(result, name)) for name in ("P", "T", "D", "h", "s", "q", "Cp"))
                or min(result.P, result.T, result.D) <= 0):
            raise RefpropError(f"ABFLSHdll({pair}) returned an invalid state.")
        return result

    def _recover_two_phase_ph(self, pressure, enthalpy, comp):
        """Invert PT within saturation limits; independently verify equilibrium.

        Native success alone is insufficient near a phase boundary. Evaluate
        each phase separately and reject any inconsistent pressure, fugacity,
        material balance or enthalpy. Never use the failed PH output.
        """
        try:
            bubble = self._flash("PQ", pressure, 0., comp)
            dew = self._flash("PQ", pressure, 1., comp)
            if (bubble.ierr or dew.ierr or not bubble.T < dew.T
                    or not bubble.h < enthalpy < dew.h):
                return None
            with profile_section("refprop_ph_recovery"):
                root = root_scalar(
                    lambda temperature: self._flash("PT", pressure, temperature, comp).h - enthalpy,
                    method="brentq", bracket=(bubble.T, dew.T), xtol=1e-10, maxiter=50,
                )
                state = self._flash("PT", pressure, root.root, comp)
            if (not root.converged or state.ierr or not 0. < state.q < 1.
                    or not state.Dl > state.Dv > 0.
                    or abs(state.h - enthalpy) > 1e-6):
                return None
            if not self._valid_two_phase_ph(state, pressure, enthalpy, comp):
                self._record_diagnostic("PH recovery rejected", "Two-phase balance checks failed")
                return None
            self._record_diagnostic("PH recovery succeeded", "Saturation-bounded PT; independently checked phases")
            return state
        except (ValueError, RuntimeError, ZeroDivisionError, OverflowError):
            return None

    def _recover_nearby_binary_ph(self, pressure, enthalpy, comp):
        """Bracket a failed binary PH flash using checked nearby native states.

        At most 24 direct native PH calls and 74 bracket corrections are tried.
        Direct calls avoid recursively invoking recovery. No failed flash output
        is used; only independently validated two-phase endpoints may seed an
        interpolation or local equilibrium correction at the requested P and H.
        """
        if (len(comp) != 2 or not all(isfinite(x) and 0 < x < 1 for x in comp)
                or not isfinite(pressure) or pressure <= 0 or not isfinite(enthalpy)
                or not isfinite(self.molar_mass) or self.molar_mass <= 0):
            return None
        with profile_section("refprop_ph_nearby_recovery"):
            for axis, steps in (("h", (1e-5, 1e-4, .001, .01, .1, 1., 5.)),
                                ("p", (1e-6, 1e-5, 1e-4, .001, .01))):
                sides = {-1: [], 1: []}
                for step in steps:
                    for sign in (-1, 1):
                        p = pressure*(1+sign*step) if axis == "p" else pressure
                        h = enthalpy+sign*step*self.molar_mass if axis == "h" else enthalpy
                        try:
                            with profile_section("refprop_std"):
                                state = self.RP.ABFLSHdll("PH", p, h, comp, 0)
                            if (state.ierr != 0 or not 0 < state.q < 1
                                    or not all(isfinite(getattr(state, field)) for field in
                                               ("P", "T", "D", "Dl", "Dv", "h", "e", "s", "Cp"))
                                    or min(state.P, state.T, state.D, state.Dv) <= 0
                                    or state.Dl <= state.Dv):
                                continue
                            corrected = self._refine_binary_ph(state, p, h, comp)
                            state = corrected if corrected is not None else state
                            if not self._valid_two_phase_ph(state, p, h, comp):
                                continue
                        except (ValueError, RuntimeError, ZeroDivisionError, OverflowError):
                            continue
                        sides[sign].append((step, state))
                        for other_step, other in sides[-sign]:
                            low, high = (state, other) if sign == -1 else (other, state)
                            low_step = step if sign == -1 else other_step
                            try:
                                candidate = self._interpolate_ph_states(low, high, low_step/(step+other_step))
                                method = "interpolation"
                                if not self._valid_two_phase_ph(candidate, pressure, enthalpy, comp):
                                    candidate = self._refine_binary_ph(candidate, pressure, enthalpy, comp)
                                    method = "interpolation and equilibrium refinement"
                                if candidate is not None:
                                    self._record_diagnostic("PH nearby recovery succeeded",
                                                            f"Validated {axis} bracket; {method}")
                                    return candidate
                            except (ValueError, RuntimeError, ZeroDivisionError, OverflowError):
                                continue
        self._record_diagnostic("PH nearby recovery exhausted", "No validated bracket correction")
        return None

    @staticmethod
    def _interpolate_ph_states(low, high, weight):
        """Interpolate inside a checked bracket; reconstruct bulk specific volume."""
        fields = {field: (1-weight)*getattr(low, field)+weight*getattr(high, field)
                  for field in ("T", "P", "Dl", "Dv", "q", "h", "s", "e")}
        for field in ("x", "y"):
            fields[field] = [(1-weight)*a+weight*b for a, b in zip(getattr(low, field), getattr(high, field))]
        fields["D"] = 1/((1-fields["q"])/fields["Dl"]+fields["q"]/fields["Dv"])
        return low._replace(**fields)

    def _valid_two_phase_ph(self, state, pressure, enthalpy, comp):
        """Independently check a recovered or refined two-phase equilibrium."""
        if state.ierr or not 0 < state.q < 1 or not state.Dl > state.Dv > 0:
            return False
        liquid = self.RP.THERMdll(state.T, state.Dl, state.x)
        vapor = self.RP.THERMdll(state.T, state.Dv, state.y)
        fl = self.RP.FUGCOFdll(state.T, state.Dl, state.x)
        fv = self.RP.FUGCOFdll(state.T, state.Dv, state.y)
        self._check_result(fl, "PH recovery liquid fugacity")
        self._check_result(fv, "PH recovery vapor fugacity")
        if fl.ierr or fv.ierr:
            return False
        errors = [abs(state.h-enthalpy)/1e-6,
                  abs((1-state.q)*liquid.h + state.q*vapor.h - enthalpy)/1e-6,
                  abs(liquid.P/pressure-1)/1e-7, abs(vapor.P/pressure-1)/1e-7,
                  abs(state.D*((1-state.q)/state.Dl + state.q/state.Dv)-1)/1e-8,
                  abs(sum(state.x)-1)/1e-8, abs(sum(state.y)-1)/1e-8]
        for i, fraction in enumerate(comp):
            errors.append(abs((1-state.q)*state.x[i]+state.q*state.y[i]-fraction)/1e-8)
            if fraction > 0:
                if min(state.x[i], state.y[i], fl.phi[i], fv.phi[i]) <= 0:
                    return False
                errors.append(abs(log(state.x[i]*fl.phi[i]/(state.y[i]*fv.phi[i])))/1e-7)
        return all(isfinite(error) and error <= 1 for error in errors)

    def _refine_binary_ph(self, original, pressure, enthalpy, comp):
        """Correct a successful PH flash or a bracket of validated native states.

        This never seeds from a failed native flash. Log densities and logistic
        fractions retain positive densities and bounded compositions/quality.
        A candidate must stay near its initial state and pass exactly
        the same independent checks used by the existing PH recovery path.
        """
        try:
            def unpack(values):
                temperature = float(values[0])
                dl, dv = np.exp(values[1:3])
                x, y, q = expit(values[3:])
                return temperature, float(dl), float(dv), [float(x), float(1-x)], [float(y), float(1-y)], float(q)

            def residual(values):
                t, dl, dv, x, y, q = unpack(values)
                liquid, vapor = self.RP.THERMdll(t, dl, x), self.RP.THERMdll(t, dv, y)
                fl, fv = self.RP.FUGCOFdll(t, dl, x), self.RP.FUGCOFdll(t, dv, y)
                self._check_result(fl, "PH refinement liquid fugacity")
                self._check_result(fv, "PH refinement vapor fugacity")
                return [(liquid.P-pressure)/pressure, (vapor.P-pressure)/pressure,
                        ((1-q)*liquid.h+q*vapor.h-enthalpy)/max(abs(enthalpy), 1.),
                        (1-q)*x[0]+q*y[0]-comp[0],
                        log(x[0]*fl.phi[0]/(y[0]*fv.phi[0])),
                        log(x[1]*fl.phi[1]/(y[1]*fv.phi[1]))]

            start = [original.T, log(original.Dl), log(original.Dv),
                     logit(original.x[0]), logit(original.y[0]), logit(original.q)]
            solved = root(residual, start, method="hybr", options={"xtol": 1e-10, "maxfev": 60})
            if not solved.success or not np.isfinite(solved.x).all():
                return None
            t, dl, dv, x, y, q = unpack(solved.x)
            if (abs(t-original.T) > .05 or max(abs(dl/original.Dl-1), abs(dv/original.Dv-1)) > .001
                    or max(abs(q-original.q), abs(x[0]-original.x[0]), abs(y[0]-original.y[0])) > .001):
                return None
            liquid, vapor = self.RP.THERMdll(t, dl, x), self.RP.THERMdll(t, dv, y)
            candidate = original._replace(T=t, P=pressure, Dl=dl, Dv=dv, x=x, y=y, q=q,
                                          D=1/((1-q)/dl+q/dv), h=(1-q)*liquid.h+q*vapor.h,
                                          e=(1-q)*liquid.e+q*vapor.e, s=(1-q)*liquid.s+q*vapor.s)
            if self._valid_two_phase_ph(candidate, pressure, enthalpy, comp):
                self._record_diagnostic("PH refinement succeeded", "Local binary equilibrium; independently checked phases")
                return candidate
        except (ValueError, RuntimeError, ZeroDivisionError, OverflowError, AttributeError):
            pass
        return None

    def _recover_liquid_ph(self, pressure, enthalpy, comp):
        """Recover a pure liquid between checked melting and boiling limits.

        This also handles REFPROP's spurious solid-state root for liquid water.
        No failed PH output is used as a temperature, density or phase guess.
        """
        try:
            bubble = self._flash("PQ", pressure, 0., comp)
            if bubble.ierr or enthalpy >= bubble.h:
                return None
            limits = self.RP.LIMITSdll("EOS", comp)
            melting = self.RP.MELTPdll(pressure, comp)
            self._check_result(melting, "MELTPdll")
            # -503 selects a reported melting-line solution (not a failed
            # flash). Any other unavailable/uncertain boundary is unusable.
            if melting.ierr not in (0, -503) or not isfinite(melting.T) or melting.T <= 0:
                return None
            bounds = (max(limits.Tmin, melting.T) + 1e-4, bubble.T - 1e-4)
            if not all(isfinite(value) for value in bounds) or bounds[0] >= bounds[1]:
                return None

            def residual(value):
                candidate = self._flash("PT", pressure, value, comp)
                if candidate.ierr:
                    raise RefpropError("PH recovery probe has a native warning.")
                return candidate.h - enthalpy

            with profile_section("refprop_ph_recovery"):
                root = root_scalar(residual, method="brentq", bracket=bounds,
                                   xtol=1e-9, maxiter=50)
                candidate = self._flash("PT", pressure, root.root, comp)
            if (root.converged and not candidate.ierr
                    and abs(candidate.h - enthalpy) <= max(1e-8, abs(enthalpy) * 1e-10)):
                self._record_diagnostic("PH recovery succeeded", "Bracketed liquid PT flash")
                return candidate
        except (ValueError, RuntimeError):
            pass
        self._record_diagnostic("PH recovery failed", "No validated saturation-bounded root")
        return None

    def _check_result(self, result, operation):
        """Retain native diagnostics; errors and uncertain phases are unusable."""
        # ctREFPROP's SETFLUIDSdll returns only an integer status.
        code = result if isinstance(result, int) else result.ierr
        message = f"Fluid setup for {self.fld_name}" if isinstance(result, int) else result.herr
        if code:
            self._record_diagnostic(f"{operation}:{code}", message)
            if code > 0 or code == -235:
                raise RefpropError(f"{operation} [{code}]: {message}")

    def _record_diagnostic(self, key, message):
        if not hasattr(self, "refprop_diagnostics"):
            self.refprop_diagnostics = {}
        entry = self.refprop_diagnostics.setdefault(key, {"count": 0, "message": message})
        entry["count"] += 1

    def _polish_ph(self, original, pressure, enthalpy, comp):
        """Remove PH/PT round-off in a known single phase, without forcing it."""
        temperature = original.T
        reason = "iteration limit"
        for _ in range(4):
            try:
                result = self._flash("PT", pressure, temperature, comp)
            except RefpropError:
                reason = "PT flash failed"
                break
            if result.ierr or 0 <= result.q <= 1 or result.Cp <= 0:
                reason = "PT phase or heat capacity unsuitable"
                break
            error = result.h - enthalpy
            if abs(error) <= max(1e-8, abs(enthalpy) * 1e-13):
                return result
            step = error / result.Cp
            if not isfinite(step) or abs(step) > .01:
                reason = "temperature correction too large"
                break
            temperature -= step
        # This is precision refinement, not recovery of a failed original flash.
        # _flash already validated the original; never replace it with a failed
        # PT probe or apply a single-phase correction to a two-phase state.
        self._record_diagnostic("PH polish declined: " + reason, reason)
        return original

    def _ensure_refprop_context(self, comp):
        """Restore the loaded fluid and its composition-specific phase envelope."""
        self._set_database_fluid()
        composition = tuple(float(value) for value in comp)
        if self._library_contexts[self.RP][1] == composition:
            return
        # An unsuccessful native spline build can destroy the previous
        # envelope too. Never advertise that previous composition as ready.
        self._library_contexts[self.RP] = (self.fld_name, None)
        if sum(value > 0 for value in composition) > 1:
            with profile_section("refprop_saturation_setup"):
                result = self.RP.SATSPLNdll(comp)
            self._check_result(result, "SATSPLNdll")
            if result.ierr:
                raise RefpropError(f"SATSPLNdll [{result.ierr}]: {result.herr}")
        self._library_contexts[self.RP] = (self.fld_name, composition)

    def _convert_props_to_SI_with_C(self, molar_props):
        """Convert molar REFPROP outputs to SI with Celsius offsets."""
        return {
            "P": molar_props.P / 1000,
            "h": molar_props.h / self.molar_mass,
            "q": molar_props.q
            if (molar_props.q >= 0.0 and molar_props.q <= 1.0)
            else self.DEFAULT_VAL,
            # Preserve phase evidence independently of the legacy quality field.
            "q_refprop": molar_props.q,
            "phase": phase_from_refprop_quality(molar_props.q),
            "T": molar_props.T - self.TC_to_K,
            "s": molar_props.s / self.molar_mass,
            "Cp": molar_props.Cp / self.molar_mass
            if molar_props.Cp > 0.0
            else self.DEFAULT_VAL,
            "D": molar_props.D * self.molar_mass,
        }

    def _convert_props_from_SI_with_C(
        self, propsIn_str: str, prop1: float, prop2: float
    ):
        """Convert SI-with-Celsius inputs to REFPROP molar units."""
        if propsIn_str == "PH" or propsIn_str == "PS":
            molar_props = prop1 * 1000, prop2 * self.molar_mass
        elif propsIn_str == "HP" or propsIn_str == "SP":
            molar_props = prop1 * self.molar_mass, prop2 * 1000
        elif propsIn_str == "PT":
            molar_props = prop1 * 1000, prop2 + self.TC_to_K
        elif propsIn_str == "TP":
            molar_props = prop1 + self.TC_to_K, prop2 * 1000
        elif propsIn_str == "PQ":
            molar_props = prop1 * 1000, prop2
        elif propsIn_str == "QP":
            molar_props = prop1, prop2 * 1000
        elif propsIn_str == "TQ":
            molar_props = prop1 + self.TC_to_K, prop2
        elif propsIn_str == "QT":
            molar_props = prop1, prop2 + self.TC_to_K
        elif propsIn_str == "TS":
            molar_props = prop1 + self.TC_to_K, prop2 * self.molar_mass
        elif propsIn_str == "ST":
            molar_props = prop1 * self.molar_mass, prop2 + self.TC_to_K
        elif propsIn_str == "TD":
            molar_props = prop1 + self.TC_to_K, prop2 / self.molar_mass
        elif propsIn_str == "DT":
            molar_props = prop1 / self.molar_mass, prop2 + self.TC_to_K
        elif propsIn_str == "PD":
            molar_props = prop1 * 1000, prop2 / self.molar_mass
        elif propsIn_str == "DP":
            molar_props = prop1 / self.molar_mass, prop2 * 1000
        else:
            raise Exception("Unrecognised fluid property input combination.")
        return molar_props

    def _set_RP_lib(self, pp_instance):
        """Select the REFPROP instance based on source/sink/default context."""
        if pp_instance.lower() == "src" or pp_instance.lower() == "source":
            self.RP = self.RP_src
        elif pp_instance.lower() == "snk" or pp_instance.lower() == "sink":
            self.RP = self.RP_snk
        else:
            self.RP = self.RP_ref
        self._set_database_fluid()
        if self._composition_ready:
            self._ensure_refprop_context(self.fld_composition)

    def _set_database_fluid(self) -> None:
        """Register the current fluid composition with REFPROP."""
        context = self._library_contexts.get(self.RP)
        if context is not None and context[0] == self.fld_name:
            return
        # Invalidate before calling native setup, including its failure path.
        self._library_contexts.pop(self.RP, None)
        result = self.RP.SETFLUIDSdll(self.fld_name)
        self._check_result(result, "SETFLUIDSdll")
        self._library_contexts[self.RP] = (self.fld_name, None)

    def _check_correct_fluid(self) -> bool:
        """Return True if the library is still loaded with the expected fluid."""
        return self.PP_fld_names == self._get_database_fluid_name()

    def _get_database_fluid_name(self) -> list[str]:
        """Return the fluid names currently loaded in the REFPROP instance."""
        return [self.RP.NAMEdll(i + 1)[0] for i in range(self.n_components)]

    def _get_critical_properties(self):
        """Fetch molar mass and critical properties for the configured fluid."""
        self._ensure_refprop_context(self.fld_composition)
        result = self.RP.REFPROPdll(
            self.fld_name,
            "PH",
            "M;PC;TC;DC",
            self.UNITS,
            self.iMass,
            self.iFlag,
            1.0,
            500,
            self.fld_composition,
        )
        self._check_result(result, "REFPROPdll(critical)")
        values = result.Output[0:4]
        if not all(isfinite(value) for value in values):
            raise RefpropError("REFPROPdll returned nonfinite critical properties.")
        return values

    def _get_mole_fractions(self, composition):
        """Return mole fractions for REFPROP low-level flash routines."""
        composition = np.array(composition, dtype=float)
        if self.composition_basis != "mass":
            return composition.tolist()

        molar_masses = np.array(
            [
                self.RP.WMOLdll(self._get_component_fraction_vector(i))
                for i in range(self.n_components)
            ],
            dtype=float,
        )
        mole_amounts = composition / molar_masses
        return (mole_amounts / mole_amounts.sum()).tolist()

    def _get_component_fraction_vector(self, index: int) -> list[float]:
        """Return a pure-component composition vector for the loaded fluid."""
        return [
            1.0 if component == index else 0.0
            for component in range(self.n_components)
        ]

    def get_component_molar_masses(self) -> tuple[float, ...]:
        """Return component molar masses in configured fluid order."""
        self._ensure_refprop_context(self.fld_composition)
        return tuple(
            float(self.RP.WMOLdll(self._get_component_fraction_vector(index)))
            for index in range(self.n_components)
        )


def _get_fluid_info(fluid: str) -> tuple[str, list[float], int, str]:
    """Parse a REFPROP fluid composition string."""
    _fld_full: list[str]
    _composition: np.ndarray[float]
    composition_basis: str = "mol"

    _fld_full = fluid.split("|")
    fld_name_ls = _fld_full[0]
    if len(_fld_full) == 1:
        fld_name: str = fld_name_ls
        fld_composition = [1.0]
        n_components = 1
    else:
        fld_name = "*".join([f + ".FLD" for f in fld_name_ls.split(";")])
        _composition = np.array(_fld_full[1].split(";"), dtype=float)
        if np.any(_composition < 0.0):
            raise ValueError("Fluid composition cannot contain negative fractions.")
        if _composition.sum() <= 0.0:
            raise ValueError("Fluid composition must contain at least one positive fraction.")
        fld_composition = (_composition / _composition.sum()).tolist()
        if len(_fld_full) == 3:
            composition_basis = _fld_full[2].strip().lower()
            if composition_basis not in ("mol", "mass"):
                raise ValueError("Fluid composition basis must be 'mol' or 'mass'.")
        n_components = 1 + fld_name.count("*")

    if len(fld_composition) != n_components:
        raise ValueError("Fluid composition length does not match component count.")

    return fld_name, fld_composition, n_components, composition_basis
