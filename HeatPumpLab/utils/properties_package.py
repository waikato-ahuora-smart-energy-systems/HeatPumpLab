"""Wrapper around REFPROP for thermodynamic property evaluation."""

from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
from collections.abc import Callable
from numbers import Real
from sys import platform
import numpy as np
from .property_cache import PropertyCacheKey, active_property_call_cache
from .solve_profile import profile_section


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

    iFlag: int = 0
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
            float(self.molar_mass),
            float(self.TC_to_K),
            float(self.DEFAULT_VAL),
            *[_normalise_cache_value(value) for value in values],
        )

    def _get_specified_properties_uncached(
        self, propsOut: str, propsIn: str, prop1: float, prop2: float
    ):
        """Return selected properties without consulting the active cache."""
        with profile_section("refprop_specified"):
            return self.RP.REFPROPdll(
                self.fld_name,
                propsIn,
                propsOut,
                self.UNITS,
                self.iMass,
                self.iFlag,
                float(prop1),
                float(prop2),
                self.fld_composition,
            ).Output

    def _get_std_properties_uncached(
        self, propsIn_str: str, prop1: float, prop2: float, comp: list[float]
    ):
        """Convenience wrapper without consulting the active cache."""
        molar_props = self._convert_props_from_SI_with_C(propsIn_str, prop1, prop2)
        with profile_section("refprop_std"):
            refprop_result = self.RP.ABFLSHdll(
                propsIn_str, molar_props[0], molar_props[1], comp, 0
            )
        return self._convert_props_to_SI_with_C(refprop_result)

    def _convert_props_to_SI_with_C(self, molar_props):
        """Convert molar REFPROP outputs to SI with Celsius offsets."""
        return {
            "P": molar_props.P / 1000,
            "h": molar_props.h / self.molar_mass,
            "q": molar_props.q
            if (molar_props.q >= 0.0 and molar_props.q <= 1.0)
            else self.DEFAULT_VAL,
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

    def _set_database_fluid(self) -> None:
        """Register the current fluid composition with REFPROP."""
        self.RP.SETFLUIDSdll(self.fld_name)

    def _check_correct_fluid(self) -> bool:
        """Return True if the library is still loaded with the expected fluid."""
        return self.PP_fld_names == self._get_database_fluid_name()

    def _get_database_fluid_name(self) -> list[str]:
        """Return the fluid names currently loaded in the REFPROP instance."""
        return [self.RP.NAMEdll(i + 1)[0] for i in range(self.n_components)]

    def _get_critical_properties(self):
        """Fetch molar mass and critical properties for the configured fluid."""
        return self.RP.REFPROPdll(
            self.fld_name,
            "PH",
            "M;PC;TC;DC",
            self.UNITS,
            self.iMass,
            self.iFlag,
            1.0,
            500,
            self.fld_composition,
        ).Output[0:4]

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
