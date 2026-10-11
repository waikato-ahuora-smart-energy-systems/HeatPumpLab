"""Compressor curves for the conventional heat pump model.

A curve turns the compressor inlet conditions into a discharge state and, when the curve knows the
compressor size, a mass flow and power. Three families are provided:

``EfficiencyMap``
    Isentropic and volumetric efficiency as functions of the dew-point evaporating and condensing
    temperatures (to, tc), e.g. the bicubic fits of manufacturer data in ``compressor_data.py``
    (manufacturer Recip-S, Recip-M, Recip-L on R600a). Swept volume gives mass flow and power. With a
    ``density_reference`` the map is read at the reference-fluid point with the same saturated-vapour
    densities (see ``density_basis.py``), so it transfers to other refrigerants and mixtures.

``EfficiencyCurve``
    Isentropic efficiency as a function of pressure ratio (a constant, a polynomial in PR, any
    callable, or HeatPumpLab's Pareto-front correlations), optionally with a swept volume and a
    volumetric efficiency so that mass flow and power follow. Without a swept volume the cycle is
    reported per kg of refrigerant (or scaled to a heat duty, see ``ConventionalHPConfig``).

``MapCurve`` / ``EN12900Curve``
    Manufacturer maps: refrigerant mass flow and power as functions of the dew-point evaporating
    and condensing temperatures (degC), rated at a fixed suction superheat or suction temperature.
    ``EN12900Curve`` takes the standard 10-coefficient polynomials (EN 12900 / AHRI 540)::

        X = c1 + c2*to + c3*tc + c4*to^2 + c5*to*tc + c6*tc^2
               + c7*to^3 + c8*tc*to^2 + c9*to*tc^2 + c10*tc^3

    Off-rating suction superheat is corrected by scaling mass flow with suction density
    (constant volumetric efficiency); power is left as rated by default, or scaled with
    mass flow x isentropic enthalpy rise (``power_correction="isentropic_work"``).

Units: temperatures degC, pressures MPa, enthalpy kJ/kg, density kg/m3, mass flow kg/s, power kW.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

from .density_basis import SaturatedVapourDensity, equivalent_point


@dataclass
class CompressorInlet:
    """Everything a curve may need, prepared by the cycle model."""

    suction: dict            # suction state: T, P, h, s, D (property-package dict)
    P_discharge: float       # MPa
    h_isentropic: float      # kJ/kg, discharge enthalpy at suction entropy
    T_evap_dew: float        # degC, dew point at suction pressure (manufacturer-map basis)
    T_cond_dew: float        # degC, dew point at discharge pressure
    rated_suction: Callable[[float], dict] | None = None   # T (degC) -> state at suction pressure
    rho_evap_dew: float | None = None   # kg/m3, saturated vapour density at the suction pressure
    rho_cond_dew: float | None = None   # kg/m3, saturated vapour density at the discharge pressure

    @property
    def pressure_ratio(self) -> float:
        return self.P_discharge / self.suction["P"]


@dataclass
class CompressorPoint:
    """Curve output. Specific values are per kg of refrigerant."""

    h_discharge: float              # kJ/kg, enthalpy leaving the compressor
    w_shaft: float                  # kJ/kg shaft work
    w_electric: float               # kJ/kg electrical input
    eta_is: float                   # isentropic efficiency on the shaft basis (implied for maps)
    m_dot: float | None = None      # kg/s, None when the curve has no size
    eta_vol: float | None = None
    q_loss: float = 0.0             # kJ/kg lost from the compressor shell to ambient
    warnings: list = field(default_factory=list)
    to_map: float | None = None     # degC, evaporating temperature the map was read at
    tc_map: float | None = None     # degC, condensing temperature the map was read at

    @property
    def W_electric(self) -> float | None:
        return None if self.m_dot is None else self.m_dot * self.w_electric

    @property
    def W_shaft(self) -> float | None:
        return None if self.m_dot is None else self.m_dot * self.w_shaft


# ----------------------------------------------------------------------------------------
# Envelopes and polynomial helpers
# ----------------------------------------------------------------------------------------


def in_envelope(to: float, tc: float, envelope) -> bool:
    """True if (to, tc) lies inside ``envelope`` (boundary included).

    ``envelope`` is either a rectangle ((to_min, to_max), (tc_min, tc_max)) or a polygon given as a
    sequence of three or more (to, tc) vertices.
    """
    pts = [tuple(map(float, p)) for p in envelope]
    if len(pts) == 2:
        (to_lo, to_hi), (tc_lo, tc_hi) = pts
        return to_lo <= to <= to_hi and tc_lo <= tc <= tc_hi
    inside = False
    n = len(pts)
    for i in range(n):
        (x1, y1), (x2, y2) = pts[i], pts[(i + 1) % n]
        # on an edge counts as inside
        cross = (x2 - x1) * (tc - y1) - (y2 - y1) * (to - x1)
        if abs(cross) < 1e-9 and min(x1, x2) - 1e-9 <= to <= max(x1, x2) + 1e-9 \
                and min(y1, y2) - 1e-9 <= tc <= max(y1, y2) + 1e-9:
            return True
        if (y1 > tc) != (y2 > tc):
            x_cross = x1 + (tc - y1) * (x2 - x1) / (y2 - y1)
            if to < x_cross:
                inside = not inside
    return inside


def bicubic_polynomial(coeffs: Sequence[float]) -> Callable[[float, float], float]:
    """eta(to, tc) = sum_i sum_j a[4i + j] * tc**j * to**i for i, j = 0..3 (16 coefficients).

    This is the ordering of the R600a compressor workbook: coefficients 0-3 are the to**0 terms
    (polynomial in tc), 4-7 the to**1 terms, and so on.
    """
    a = [float(v) for v in coeffs]
    if len(a) != 16:
        raise ValueError("bicubic efficiency polynomials have 16 coefficients")

    def f(to: float, tc: float) -> float:
        total = 0.0
        for i in range(4):
            row = a[4 * i] + tc * (a[4 * i + 1] + tc * (a[4 * i + 2] + tc * a[4 * i + 3]))
            total += row * to ** i
        return total
    return f


# ----------------------------------------------------------------------------------------
# Efficiency maps in (to, tc)
# ----------------------------------------------------------------------------------------


@dataclass
class EfficiencyMap:
    """Isentropic and volumetric efficiency maps in dew-point to, tc (degC).

    eta_is(to, tc): isentropic efficiency on ``efficiency_basis``:
        "electrical": eta_is = m (h_is - h_suc) / P_electrical  (overall, includes the motor)
        "shaft":      eta_is = m (h_is - h_suc) / P_shaft
    eta_vol(to, tc): m / (rho_suction * displacement).
    displacement_m3h: swept volume flow at the rated frequency (m3/h); None reports per kg only.
    frequency_Hz / rated_frequency_Hz: supply frequency and the frequency the displacement and maps
        refer to; the swept volume scales with frequency_Hz / rated_frequency_Hz. The efficiency maps
        are kept unchanged with speed (an assumption: they are rated-speed data).
    rated_rpm: shaft speed at the rated frequency (4-pole motor: 1450 rpm at 50 Hz).
    motor_eff: converts between shaft and electrical power.
    heat_loss_fraction: share of the input power (electrical basis) or shaft power (shaft basis)
        that leaves the shell instead of heating the gas. 0 suits suction-gas-cooled semi-hermetics.
    envelope: application envelope for warnings (see ``in_envelope``).
    eta_is_bounds / eta_vol_bounds: guard rails; values outside are clipped with a warning
        (the fits are not physical far outside the envelope).
    density_reference: saturated-vapour density curve of the fluid the maps were fitted for. When
        set, (to, tc) in the maps and the envelope are reference-fluid temperatures matched on density
        (``density_basis.equivalent_point``): tc on the saturated vapour at the condensing pressure,
        to on the saturated vapour at the evaporating pressure or, with ``suction_reference``, on
        the actual suction density. None reads the maps at the working fluid's own (to, tc).
    suction_reference: the reference fluid's suction density against to at the rated suction
        superheat (R600a at to + 15 K); captures the effect of the actual superheat.
    """

    eta_is: Callable[[float, float], float]
    eta_vol: Callable[[float, float], float] | None = None
    displacement_m3h: float | None = None
    efficiency_basis: str = "electrical"
    frequency_Hz: float = 50.0
    rated_frequency_Hz: float = 50.0
    rated_rpm: float = 1450.0
    motor_eff: float = 0.93
    heat_loss_fraction: float = 0.0
    envelope: object = None
    name: str = ""
    eta_is_bounds: tuple[float, float] = (0.2, 1.0)
    eta_vol_bounds: tuple[float, float] = (0.0, 1.1)
    density_reference: SaturatedVapourDensity | None = None
    suction_reference: SaturatedVapourDensity | None = None

    @property
    def rpm(self) -> float:
        return self.rated_rpm * self.frequency_Hz / self.rated_frequency_Hz

    @property
    def swept_volume_m3h(self) -> float | None:
        """Swept volume at frequency_Hz."""
        if self.displacement_m3h is None:
            return None
        return self.displacement_m3h * self.frequency_Hz / self.rated_frequency_Hz

    def map_point(self, to: float, tc: float, rho_evap_dew: float | None = None,
                  rho_cond_dew: float | None = None, rho_suction: float | None = None) -> tuple[float, float, list]:
        """(to, tc) at which the maps are read, and any warnings, for a working-fluid point."""
        if self.density_reference is None:
            return to, tc, []
        if rho_cond_dew is None or (rho_evap_dew is None and self.suction_reference is None):
            raise ValueError("a density-basis map needs the saturated vapour densities at the "
                             "evaporating and condensing pressures (rho_evap_dew, rho_cond_dew)")
        if self.suction_reference is not None and rho_suction is None:
            raise ValueError("a suction-density map needs the suction density (rho_suction)")
        return equivalent_point(self.density_reference, rho_evap_dew, rho_cond_dew,
                                self.suction_reference, rho_suction)

    def efficiencies(self, to: float, tc: float, rho_evap_dew: float | None = None,
                     rho_cond_dew: float | None = None, rho_suction: float | None = None) -> dict:
        """Map values for a working-fluid point: eta_is, eta_vol (unclipped), the map point and envelope."""
        tm, tcm, warnings = self.map_point(to, tc, rho_evap_dew, rho_cond_dew, rho_suction)
        return dict(eta_is=self.eta_is(tm, tcm), eta_vol=None if self.eta_vol is None else self.eta_vol(tm, tcm),
                    to_map=tm, tc_map=tcm, warnings=warnings,
                    in_envelope=self.envelope is None or in_envelope(tm, tcm, self.envelope))

    def fluid_grid(self, to_axis: Sequence[float], tc_axis: Sequence[float],
                   curve: SaturatedVapourDensity | None = None,
                   suction_curve: SaturatedVapourDensity | None = None) -> dict:
        """The maps recomputed on a working fluid's (to, tc) axes, rows by tc.

        curve: the working fluid's dew-point vapour density (``density_basis.fluid_curve``);
        suction_curve: its suction density against to at the superheat being drawn (needed with a
        ``suction_reference``). to_map / tc_map give the map temperature of each column / row (each
        depends on one axis only); cells beyond the fluid tables are None.
        """
        if self.density_reference is None:
            to_map, tc_map = list(to_axis), list(tc_axis)
        else:
            if curve is None:
                raise ValueError("a density-basis map needs the working fluid's vapour density curve")
            from .density_basis import to_reference
            if self.suction_reference is not None:
                if suction_curve is None:
                    raise ValueError("a suction-density map needs the working fluid's suction density curve")
                to_map = [to_reference(suction_curve, self.suction_reference, t)[0] for t in to_axis]
            else:
                to_map = [to_reference(curve, self.density_reference, t)[0] for t in to_axis]
            tc_map = [to_reference(curve, self.density_reference, t)[0] for t in tc_axis]
        cell = lambda f, a, b: None if a is None or b is None or f is None else f(a, b)
        return dict(
            to=list(to_axis), tc=list(tc_axis), to_map=to_map, tc_map=tc_map,
            eta_is=[[cell(self.eta_is, a, b) for a in to_map] for b in tc_map],
            eta_vol=[[cell(self.eta_vol, a, b) for a in to_map] for b in tc_map],
            in_env=[[None if a is None or b is None else (self.envelope is None or in_envelope(a, b, self.envelope))
                     for a in to_map] for b in tc_map],
        )

    def evaluate(self, inlet: CompressorInlet) -> CompressorPoint:
        to, tc, warnings = self.map_point(inlet.T_evap_dew, inlet.T_cond_dew, inlet.rho_evap_dew, inlet.rho_cond_dew,
                                          inlet.suction.get("D"))
        if self.envelope is not None and not in_envelope(to, tc, self.envelope):
            where = (f"to {to:.1f} C, tc {tc:.1f} C" if self.density_reference is None else
                     f"{self.density_reference.fluid}-equivalent to {to:.1f} C, tc {tc:.1f} C")
            warnings.append(f"outside the {self.name or 'compressor'} envelope ({where})")
        eta = _clip(self.eta_is(to, tc), self.eta_is_bounds, "eta_is", warnings)
        h_in = inlet.suction["h"]
        dh_is = inlet.h_isentropic - h_in
        if self.efficiency_basis == "electrical":
            w_el = dh_is / eta
            w_sh = w_el * self.motor_eff
            to_gas = (1 - self.heat_loss_fraction) * w_el
        elif self.efficiency_basis == "shaft":
            w_sh = dh_is / eta
            w_el = w_sh / self.motor_eff
            to_gas = (1 - self.heat_loss_fraction) * w_sh
        else:
            raise ValueError("efficiency_basis must be 'electrical' or 'shaft'")
        m_dot = eta_vol = None
        if self.eta_vol is not None:
            eta_vol = _clip(self.eta_vol(to, tc), self.eta_vol_bounds, "eta_vol", warnings)
        if self.displacement_m3h is not None:
            if eta_vol is None:
                raise ValueError("a displacement needs a volumetric efficiency map")
            m_dot = eta_vol * self.swept_volume_m3h / 3600.0 * inlet.suction["D"]
        return CompressorPoint(h_discharge=h_in + to_gas, w_shaft=w_sh, w_electric=w_el,
                               eta_is=dh_is / w_sh, m_dot=m_dot, eta_vol=eta_vol,
                               q_loss=(w_el if self.efficiency_basis == "electrical" else w_sh) - to_gas,
                               warnings=warnings, to_map=to, tc_map=tc)


def _clip(value: float, bounds: tuple[float, float], name: str, warnings: list) -> float:
    lo, hi = bounds
    if not lo <= value <= hi:
        warnings.append(f"{name} {value:.3f} clipped to {lo}-{hi}")
        return min(max(value, lo), hi)
    return value


# ----------------------------------------------------------------------------------------
# Efficiency-based curves
# ----------------------------------------------------------------------------------------


def brendel_volumetric_efficiency(pressure_ratio: float) -> float:
    """Brendel et al. correlation, as used elsewhere in HeatPumpLab."""
    return 1 - 0.08244 * pow(max(pressure_ratio, 1.0) - 1, 0.72773)


def _as_function(spec, name: str) -> Callable[[float], float]:
    if callable(spec):
        return spec
    if isinstance(spec, (int, float)):
        value = float(spec)
        return lambda pr: value
    if isinstance(spec, str):
        key = spec.lower()
        if key in ("pareto", "pareto_v2", "pareto_v1"):
            from ...components.compressor import pareto_front_efficiency
            version = 1 if key == "pareto_v1" else 2
            return lambda pr: pareto_front_efficiency(pr, "isentropic", version=version)
        if key == "brendel":
            return brendel_volumetric_efficiency
        raise ValueError(f"unknown {name} correlation '{spec}'")
    if isinstance(spec, Sequence):
        coeffs = [float(c) for c in spec]
        return lambda pr: sum(c * pr ** i for i, c in enumerate(coeffs))
    raise TypeError(f"{name} must be a number, a coefficient sequence, a callable or a name")


@dataclass
class EfficiencyCurve:
    """Isentropic (and optionally volumetric) efficiency versus pressure ratio.

    eta_is:  float, coefficients [c0, c1, c2, ...] of a polynomial in PR, callable(PR), or
             "pareto_v2" / "pareto_v1" (HeatPumpLab compressor Pareto fronts, shaft basis).
    eta_vol: as above, or "brendel"; only used when ``displacement_m3h`` is given.
    displacement_m3h: swept volume flow (m3/h). Gives mass flow = eta_vol * V * rho_suction.
    motor_eff: electrical-to-shaft efficiency.
    heat_loss_fraction: share of shaft work lost from the shell instead of entering the gas.
    """

    eta_is: object = "pareto_v2"
    motor_eff: float = 0.93
    displacement_m3h: float | None = None
    eta_vol: object = "brendel"
    heat_loss_fraction: float = 0.0
    eta_is_bounds: tuple[float, float] = (0.2, 1.0)

    def evaluate(self, inlet: CompressorInlet) -> CompressorPoint:
        pr = inlet.pressure_ratio
        eta = _as_function(self.eta_is, "eta_is")(pr)
        lo, hi = self.eta_is_bounds
        warnings = []
        if not lo <= eta <= hi:
            warnings.append(f"eta_is {eta:.3f} at PR {pr:.2f} clipped to {lo}-{hi}")
            eta = min(max(eta, lo), hi)
        h_in = inlet.suction["h"]
        w = (inlet.h_isentropic - h_in) / eta
        q_loss = self.heat_loss_fraction * w
        m_dot = eta_vol = None
        if self.displacement_m3h is not None:
            eta_vol = _as_function(self.eta_vol, "eta_vol")(pr)
            if eta_vol <= 0:
                warnings.append(f"volumetric efficiency {eta_vol:.3f} <= 0 at PR {pr:.2f}")
            m_dot = max(eta_vol, 0.0) * self.displacement_m3h / 3600.0 * inlet.suction["D"]
        return CompressorPoint(h_discharge=h_in + w - q_loss, w_shaft=w, w_electric=w / self.motor_eff,
                               eta_is=eta, m_dot=m_dot, eta_vol=eta_vol, q_loss=q_loss,
                               warnings=warnings)


# ----------------------------------------------------------------------------------------
# Manufacturer maps
# ----------------------------------------------------------------------------------------


_MASS_FLOW_UNITS = {"kg/s": 1.0, "kg/h": 1 / 3600.0, "g/s": 1e-3}
_POWER_UNITS = {"kW": 1.0, "W": 1e-3}


@dataclass
class MapCurve:
    """Mass flow and power maps in dew-point evaporating/condensing temperature (degC).

    mass_flow(to, tc), power(to, tc): callables in the units given by ``mass_flow_unit`` and
    ``power_unit``, rated at suction superheat ``rated_superheat`` (K above the evaporating dew
    point) or, if set, a fixed suction temperature ``rated_suction_T`` (degC).
    power_basis: "electrical" (semi-hermetic / hermetic data) or "shaft" (open drive);
        ``motor_eff`` converts between them.
    heat_loss_fraction: share of the input power that leaves the shell instead of the gas.
    envelope: optional application envelope, see ``in_envelope``; outside it a warning is raised.
    """

    mass_flow: Callable[[float, float], float]
    power: Callable[[float, float], float]
    mass_flow_unit: str = "kg/h"
    power_unit: str = "W"
    rated_superheat: float | None = 10.0
    rated_suction_T: float | None = None
    power_basis: str = "electrical"
    motor_eff: float = 0.93
    heat_loss_fraction: float = 0.0
    correct_mass_flow: bool = True
    power_correction: str = "none"      # "none" | "isentropic_work"
    envelope: tuple | None = None

    def _rated_suction_temperature(self, inlet: CompressorInlet) -> float | None:
        if self.rated_suction_T is not None:
            return self.rated_suction_T
        if self.rated_superheat is not None:
            return inlet.T_evap_dew + self.rated_superheat
        return None

    def evaluate(self, inlet: CompressorInlet) -> CompressorPoint:
        to, tc = inlet.T_evap_dew, inlet.T_cond_dew
        warnings = []
        if self.envelope is not None and not in_envelope(to, tc, self.envelope):
            warnings.append(f"outside the compressor envelope (to {to:.1f} C, tc {tc:.1f} C)")
        m_rated = self.mass_flow(to, tc) * _MASS_FLOW_UNITS[self.mass_flow_unit]
        W_rated = self.power(to, tc) * _POWER_UNITS[self.power_unit]
        if m_rated <= 0 or W_rated <= 0:
            raise ValueError(f"compressor map gives non-positive mass flow or power at to {to:.1f}, tc {tc:.1f} C")
        suction = inlet.suction
        m_dot, W = m_rated, W_rated
        T_rated = self._rated_suction_temperature(inlet)
        needs_rated = (self.correct_mass_flow or self.power_correction != "none") and T_rated is not None \
            and abs(T_rated - suction["T"]) > 1e-9
        if needs_rated:
            if inlet.rated_suction is None:
                raise ValueError("rated suction state needed for the superheat correction")
            rated = inlet.rated_suction(T_rated)
            if self.correct_mass_flow:
                m_dot = m_rated * suction["D"] / rated["D"]
            if self.power_correction == "isentropic_work":
                dh_actual = inlet.h_isentropic - suction["h"]
                dh_rated = rated["h_isentropic"] - rated["h"]
                W = W_rated * (m_dot / m_rated) * (dh_actual / dh_rated)
            elif self.power_correction != "none":
                raise ValueError(f"unknown power_correction '{self.power_correction}'")
        if self.power_basis == "electrical":
            W_el, W_sh = W, W * self.motor_eff
            to_gas = (1 - self.heat_loss_fraction) * W_el
        elif self.power_basis == "shaft":
            W_sh, W_el = W, W / self.motor_eff
            to_gas = (1 - self.heat_loss_fraction) * W_sh
        else:
            raise ValueError("power_basis must be 'electrical' or 'shaft'")
        w_gas = to_gas / m_dot
        eta_is = (inlet.h_isentropic - suction["h"]) / (W_sh / m_dot)
        return CompressorPoint(h_discharge=suction["h"] + w_gas, w_shaft=W_sh / m_dot,
                               w_electric=W_el / m_dot, eta_is=eta_is, m_dot=m_dot,
                               q_loss=(W_el - to_gas) / m_dot if self.power_basis == "electrical"
                               else (W_sh - to_gas) / m_dot,
                               warnings=warnings)


def en12900_polynomial(coeffs: Sequence[float]) -> Callable[[float, float], float]:
    """X(to, tc) for the 10 EN 12900 coefficients c1..c10."""
    c = [float(v) for v in coeffs]
    if len(c) != 10:
        raise ValueError("EN 12900 polynomials have 10 coefficients")

    def f(to: float, tc: float) -> float:
        return (c[0] + c[1] * to + c[2] * tc + c[3] * to ** 2 + c[4] * to * tc + c[5] * tc ** 2
                + c[6] * to ** 3 + c[7] * tc * to ** 2 + c[8] * to * tc ** 2 + c[9] * tc ** 3)
    return f


class EN12900Curve(MapCurve):
    """MapCurve from EN 12900 10-coefficient polynomials (default units W and kg/h)."""

    def __init__(self, power_coeffs: Sequence[float], mass_flow_coeffs: Sequence[float], **kwargs):
        super().__init__(mass_flow=en12900_polynomial(mass_flow_coeffs),
                         power=en12900_polynomial(power_coeffs), **kwargs)
        self.power_coeffs = tuple(power_coeffs)
        self.mass_flow_coeffs = tuple(mass_flow_coeffs)

    def __repr__(self):
        return (f"EN12900Curve(power_coeffs={self.power_coeffs}, mass_flow_coeffs={self.mass_flow_coeffs}, "
                f"rated_superheat={self.rated_superheat}, rated_suction_T={self.rated_suction_T})")


def fit_en12900(to: Sequence[float], tc: Sequence[float], values: Sequence[float]) -> list[float]:
    """Least-squares EN 12900 coefficients from tabulated data (e.g. a manufacturer table)."""
    import numpy as np
    to, tc, y = (np.asarray(v, dtype=float) for v in (to, tc, values))
    A = np.column_stack([np.ones_like(to), to, tc, to ** 2, to * tc, tc ** 2,
                         to ** 3, tc * to ** 2, to * tc ** 2, tc ** 3])
    coeffs, *_ = np.linalg.lstsq(A, y, rcond=None)
    return coeffs.tolist()


