"""CoolProp (v8) property backend that can stand in for REFPROP's PropertyPackage.

``CoolPropPackage`` takes the same fluid strings as ``PropertyPackage`` ("Butane",
"Propane;Isopentane|0.2;0.8|mass") and returns the same dictionaries from ``get_std_properties``:
P (MPa), T (degC), h (kJ/kg), s (kJ/kg/K), D (kg/m3), Cp (kJ/kg/K), q, q_refprop and phase, with
REFPROP's phase labels. Failures raise ``RefpropError`` so existing error handling applies.

Backend: CoolProp HEOS (the same reference equations of state as REFPROP for the common
refrigerants). Reference state: IIR for each component (h = 200 kJ/kg, s = 1 kJ/kg/K, saturated
liquid at 0 C), REFPROP's default for the hydrocarbon refrigerants; water and air keep CoolProp's
default. Absolute h and s may still differ from REFPROP where its fluid file uses another default;
differences of state (and every cycle result) do not depend on the reference.

Mixtures use CoolProp's HEOS mixture model with its own binary interaction parameters, which can
differ from REFPROP's (and are estimated for pairs without fitted data); two-phase flashes for
mixtures are less robust than REFPROP's. Compare against REFPROP before relying on mixture results.

Select the backend with ``make_property_package(fluid, backend="coolprop")`` or the environment
variable ``HEATPUMPLAB_PROPERTY_BACKEND=coolprop`` (see ``property_backend``).
"""
from __future__ import annotations

import math
import os

from .properties_package import RefpropError
from .property_cache import active_property_call_cache

try:
    import CoolProp.CoolProp as _CP
except ImportError:  # pragma: no cover - optional dependency
    _CP = None

DEFAULT_VAL = -9999990.0
ENV_VAR = "HEATPUMPLAB_PROPERTY_BACKEND"

# REFPROP fluid-file names used in HeatPumpLab -> CoolProp names (others are passed through,
# CoolProp also accepts many aliases and R-numbers).
NAME_MAP = {
    "butane": "n-Butane", "nbutane": "n-Butane", "n-butane": "n-Butane",
    "isobutane": "IsoButane", "ibutane": "IsoButane",
    "pentane": "n-Pentane", "npentane": "n-Pentane", "n-pentane": "n-Pentane",
    "isopentane": "Isopentane", "ipentane": "Isopentane",
    "neopentane": "Neopentane", "hexane": "n-Hexane", "nhexane": "n-Hexane", "isohexane": "Isohexane",
    "cyclopentane": "Cyclopentane", "cyclohexane": "Cyclohexane", "propane": "Propane",
    "propylen": "Propylene", "propylene": "Propylene", "water": "Water", "air": "Air",
    "co2": "CarbonDioxide", "ammonia": "Ammonia", "nh3": "Ammonia", "dme": "DimethylEther",
    "r1234zee": "R1234ze(E)", "r1234ze(e)": "R1234ze(E)", "r1233zde": "R1233zd(E)",
    "r1233zd(e)": "R1233zd(E)", "r1234yf": "R1234yf", "r1336mzzz": "R1336mzz(Z)",
    "r245fa": "R245fa", "r134a": "R134a", "r600": "n-Butane", "r600a": "IsoButane",
    "r290": "Propane", "r601": "n-Pentane", "r601a": "Isopentane", "r718": "Water", "r744": "CarbonDioxide",
}
_KEEP_DEFAULT_REFERENCE = {"Water", "Air", "CarbonDioxide", "Nitrogen", "Oxygen", "Argon"}

GAS_LIKE = frozenset({"saturated_vapour", "superheated_vapour", "supercritical"})


def coolprop_available() -> bool:
    return _CP is not None


def coolprop_version() -> str | None:
    if _CP is None:
        return None
    try:
        return _CP.get_global_param_string("version")
    except Exception:  # pragma: no cover
        return getattr(__import__("CoolProp"), "__version__", None)


def parse_fluid(fluid: str) -> tuple[list[str], list[float], str]:
    """HeatPumpLab fluid string -> (CoolProp names, fractions, basis 'mol' | 'mass')."""
    parts = fluid.split("|")
    names = [n.strip() for n in parts[0].split(";")]
    if len(parts) == 1:
        if len(names) != 1:
            raise ValueError("a mixture needs fractions: 'A;B|0.3;0.7|mass'")
        fractions = [1.0]
    else:
        fractions = [float(v) for v in parts[1].split(";")]
    basis = parts[2].strip().lower() if len(parts) == 3 else "mol"
    if basis not in ("mol", "mass"):
        raise ValueError("Fluid composition basis must be 'mol' or 'mass'.")
    if len(fractions) != len(names):
        raise ValueError("Fluid composition length does not match component count.")
    if any(v < 0 for v in fractions) or sum(fractions) <= 0:
        raise ValueError("Fluid composition must be non-negative with a positive sum.")
    total = sum(fractions)
    return [coolprop_name(n) for n in names], [v / total for v in fractions], basis


def coolprop_name(name: str) -> str:
    key = name.strip()
    if key.upper().endswith(".FLD"):
        key = key[:-4]
    return NAME_MAP.get(key.lower().replace(" ", ""), key)


def _namespaces():
    """Where CoolProp may define its input-pair and phase constants (v6/7: CoolProp.CoolProp and
    CoolProp.constants; v8 nanobind: possibly the package or enum classes on these modules)."""
    spaces = [_CP]
    import importlib
    for mod in ("CoolProp", "CoolProp.constants", "CoolProp._constants"):
        try:
            spaces.append(importlib.import_module(mod))
        except Exception:
            pass
    out = []
    for sp in spaces:
        if sp is None or sp in out:
            continue
        out.append(sp)
        for attr in dir(sp):                     # enum classes, e.g. input_pairs / phases
            if attr.startswith("_"):
                continue
            obj = getattr(sp, attr, None)
            if isinstance(obj, type) and obj not in out:
                out.append(obj)
    return out


def _const(name: str, index_function: str | None = None, index_name: str | None = None):
    """A CoolProp constant by name, searched across the namespaces above; then the index function."""
    spaces = _namespaces()
    for sp in spaces:
        value = getattr(sp, name, None)
        if value is not None and not callable(value):
            return value
    if index_function:
        for sp in spaces:
            fn = getattr(sp, index_function, None)
            if callable(fn):
                try:
                    return fn(index_name or name)
                except Exception:
                    pass
    raise RefpropError(f"CoolProp constant '{name}' not found (searched {len(spaces)} namespaces)")


def _as_int(value):
    try:
        return int(value)
    except Exception:
        return value


class CoolPropPackage:
    """CoolProp stand-in for ``PropertyPackage`` (same units, keys and phase labels)."""

    DEFAULT_VAL = DEFAULT_VAL
    TC_to_K = 273.15

    def __init__(self, fld_name: str, pp_instance: str = "", backend: str = "HEOS",
                 reference_state: str | None = "IIR"):
        if _CP is None:
            raise RefpropError("CoolProp is not installed (pip install 'coolprop>=8')")
        self.fld_name = fld_name
        self.pp_instance = pp_instance
        self.names, fractions, self.composition_basis = parse_fluid(fld_name)
        self.n_components = len(self.names)
        self.reference_notes: list[str] = []
        if reference_state:
            for n in self.names:
                if n in _KEEP_DEFAULT_REFERENCE:
                    continue
                try:
                    setter = getattr(_CP, "set_reference_state", None) or getattr(_CP, "set_reference_stateS")
                    setter(n, reference_state)
                except Exception as exc:  # e.g. IIR undefined for a fluid with Tc < 0 C
                    self.reference_notes.append(f"{n}: kept default reference ({exc})")
        try:
            self.AS = _CP.AbstractState(backend, "&".join(self.names))
            if self.n_components > 1:
                if self.composition_basis == "mass":
                    self.AS.set_mass_fractions(fractions)
                else:
                    self.AS.set_mole_fractions(fractions)
        except Exception as exc:
            raise RefpropError(f"CoolProp could not set up '{fld_name}': {exc}") from exc
        self.fld_composition = list(self.AS.get_mole_fractions()) if self.n_components > 1 else [1.0]
        self.molar_mass = float(self.AS.molar_mass()) * 1000.0      # g/mol
        self.P_crit, self.T_crit, self.D_crit = self._critical_point()
        try:
            self.P_max = float(self.AS.pmax()) / 1e6
        except Exception:
            self.P_max = float("nan")
        self._pairs = {
            "PT": (_const("PT_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e6, b + 273.15)),
            "TP": (_const("PT_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e6, a + 273.15)),
            "PH": (_const("HmassP_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e3, a * 1e6)),
            "HP": (_const("HmassP_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e3, b * 1e6)),
            "PS": (_const("PSmass_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e6, b * 1e3)),
            "SP": (_const("PSmass_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e6, a * 1e3)),
            "PQ": (_const("PQ_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e6, b)),
            "QP": (_const("PQ_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e6, a)),
            "TQ": (_const("QT_INPUTS", "get_input_pair_index"), lambda a, b: (b, a + 273.15)),
            "QT": (_const("QT_INPUTS", "get_input_pair_index"), lambda a, b: (a, b + 273.15)),
            "TD": (_const("DmassT_INPUTS", "get_input_pair_index"), lambda a, b: (b, a + 273.15)),
            "DT": (_const("DmassT_INPUTS", "get_input_pair_index"), lambda a, b: (a, b + 273.15)),
            "PD": (_const("DmassP_INPUTS", "get_input_pair_index"), lambda a, b: (b, a * 1e6)),
            "DP": (_const("DmassP_INPUTS", "get_input_pair_index"), lambda a, b: (a, b * 1e6)),
            "TS": (_const("SmassT_INPUTS", "get_input_pair_index"), lambda a, b: (b * 1e3, a + 273.15)),
            "ST": (_const("SmassT_INPUTS", "get_input_pair_index"), lambda a, b: (a * 1e3, b + 273.15)),
        }
        self._phase = {}
        for cname, label in (("iphase_liquid", "liquid"), ("iphase_gas", "superheated_vapour"),
                             ("iphase_supercritical", "supercritical"),
                             ("iphase_supercritical_gas", "superheated_vapour"),
                             ("iphase_supercritical_liquid", "liquid"), ("iphase_twophase", "two_phase")):
            try:
                self._phase[_as_int(_const(cname, "get_phase_index", cname[1:]))] = label
            except RefpropError:
                pass                          # fall back to quality/temperature below

    def _sat_T(self, P_Pa: float) -> float:
        """Dew temperature (degC) at P, for phase inference when CoolProp's phase index is unavailable."""
        if P_Pa >= self.P_crit * 1e6:
            return -273.15
        st = _CP.AbstractState(self.AS.backend_name(), "&".join(self.names)) if self.n_components == 1 else None
        try:
            probe = st or self.AS
            probe.update(self._pairs["PQ"][0], P_Pa, 1.0)
            return probe.T() - 273.15
        except Exception:
            return float("inf")

    # ---- critical point ------------------------------------------------------------

    def _critical_point(self):
        try:
            return (float(self.AS.p_critical()) / 1e6, float(self.AS.T_critical()) - 273.15,
                    float(self.AS.rhomass_critical()))
        except Exception:
            pass
        try:   # mixtures: CoolProp's critical point search
            pts = self.AS.all_critical_points()
            if pts:
                pt = max(pts, key=lambda c: c.T)
                return pt.p / 1e6, pt.T - 273.15, pt.rhomolar * self.molar_mass / 1000.0
        except Exception:
            pass
        return float("nan"), float("nan"), float("nan")

    # ---- properties ----------------------------------------------------------------

    def get_std_properties(self, propsIn_str: str, prop1: float, prop2: float, comp=None) -> dict:
        """PropertyPackage-compatible state dict (MPa, degC, kJ/kg, kJ/kg/K, kg/m3)."""
        cache = active_property_call_cache()
        if cache is None:
            return self._std(propsIn_str, prop1, prop2)
        key = ("coolprop", id(self), self.fld_name, propsIn_str, float(prop1), float(prop2))
        hit, result = cache.lookup(key)
        if hit:
            return dict(result)
        result = self._std(propsIn_str, prop1, prop2)
        cache.store(key, dict(result))
        return dict(result)

    def _std(self, pair: str, a: float, b: float) -> dict:
        if pair not in self._pairs:
            raise RefpropError(f"Unsupported CoolProp input pair '{pair}'")
        index, convert = self._pairs[pair]
        v1, v2 = convert(float(a), float(b))
        try:
            self.AS.update(index, v1, v2)
            T, P = self.AS.T(), self.AS.p()
            h, s, D = self.AS.hmass(), self.AS.smass(), self.AS.rhomass()
            ph = self.AS.phase()
            try:
                Q = self.AS.Q()
            except Exception:
                Q = float("nan")
            try:
                cp = self.AS.cpmass() / 1000.0
            except Exception:
                cp = DEFAULT_VAL
        except Exception as exc:
            err = RefpropError(f"CoolProp {pair} flash failed for {self.fld_name} at ({a}, {b}): {exc}")
            err.property_context = {"fluid": self.fld_name, "input_pair": pair, "inputs": [a, b],
                                    "backend": "CoolProp"}
            raise err from exc
        if not all(math.isfinite(x) for x in (T, P, h, s, D)):
            raise RefpropError(f"CoolProp {pair} returned a nonfinite state for {self.fld_name}")
        phase = self._phase.get(_as_int(ph), "unknown")
        if phase == "unknown":               # phase constants unavailable: infer from quality
            phase = "two_phase" if math.isfinite(Q) and 0.0 <= Q <= 1.0 else (
                "superheated_vapour" if math.isfinite(self.T_crit) and T - 273.15 >= self._sat_T(P) else "liquid")
        if pair in ("PQ", "QP", "TQ", "QT"):
            q_in = b if pair in ("PQ", "TQ") else a
            Q = q_in
        if phase == "two_phase" or (math.isfinite(Q) and 0.0 <= Q <= 1.0 and pair in ("PQ", "QP", "TQ", "QT")):
            if Q <= 1e-12:
                phase, Q = "saturated_liquid", 0.0
            elif Q >= 1 - 1e-12:
                phase, Q = "saturated_vapour", 1.0
            else:
                phase = "two_phase"
        q = Q if phase in ("two_phase", "saturated_liquid", "saturated_vapour") else DEFAULT_VAL
        q_refprop = {"liquid": -1.0, "superheated_vapour": 998.0, "supercritical": 999.0}.get(phase, q)
        return {"P": P / 1e6, "h": h / 1000.0, "q": q, "q_refprop": q_refprop, "phase": phase,
                "T": T - 273.15, "s": s / 1000.0, "Cp": cp if cp > 0 else DEFAULT_VAL, "D": D}

    def get_component_molar_masses(self) -> tuple[float, ...]:
        out = []
        for n in self.names:
            out.append(float(_CP.PropsSI("molar_mass", n)) * 1000.0)
        return tuple(out)


# ----------------------------------------------------------------------------------------
# Backend selection
# ----------------------------------------------------------------------------------------


def property_backend(backend: str | None = None) -> str:
    """'refprop' or 'coolprop': the argument, else $HEATPUMPLAB_PROPERTY_BACKEND, else 'refprop'."""
    value = (backend or os.environ.get(ENV_VAR) or "refprop").strip().lower()
    if value not in ("refprop", "coolprop"):
        raise ValueError(f"property backend must be 'refprop' or 'coolprop', not '{value}'")
    return value


def make_property_package(fluid: str, pp_instance: str = "", backend: str | None = None):
    """PropertyPackage (REFPROP) or CoolPropPackage for ``fluid``."""
    if property_backend(backend) == "coolprop":
        return CoolPropPackage(fluid, pp_instance)
    from .properties_package import PropertyPackage
    return PropertyPackage(fluid, pp_instance)
