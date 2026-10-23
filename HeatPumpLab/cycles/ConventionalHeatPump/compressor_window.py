"""Compressor efficiency beyond the application window, and the mechanical limits.

The window's edges are of two kinds:

ratio edges
    The pressure (or density) ratio is too high or too low and efficiency falls away fast. For the R600a
    window: the high-ratio edge (6, 120) -> (-10, 97) sits at a density ratio of about 19 (condensing saturated
    vapour over the 15 K-superheat suction), where the volumetric efficiency is about 0.5; the low-ratio edge
    (27, 42) -> (48, 70) at a pressure ratio of about 1.5-1.7.
mechanical edges
    Seals, materials and the motor: minimum and maximum evaporating, maximum and minimum condensing. They say
    nothing about efficiency.

Outside the window the efficiency maps are read at the nearest window point (``nearest_on_envelope``, continuous
across the edge, so no rise or jump from extrapolating the polynomial fits), then fall off with the ratio:

- above the high-ratio line (density ratio r above its limit r_hi), clearance re-expansion:
  eta_vol = 1 - (1 - eta_vol,q)(r - 1)/(r_hi - 1), and eta_is = eta_is,q (eta_vol / eta_vol,q)^beta, so both
  reach zero together; beta = 0.4 gives about 25 % lower eta_is at 1.5 r_hi (a best guess, no data past the
  window);
- below the low-ratio line, valve losses take over as the lift vanishes: with x = r - 1,
  eta_is = eta_is,q [x / (x + delta)] / [x_lo / (x_lo + delta)], delta = 0.3; eta_vol is held (at most 1).

The ratio lines follow the ratio edges (piecewise linear in to, flat beyond their ends), so the falloff starts
exactly at the window and the result is continuous. The density ratio is the working fluid's own (condensing
saturated vapour, or the critical isochore when supercritical, over the suction), which on the density basis is
the map point's R600a ratio; on the saturated and temperature bases the map point's R600a ratio is used, so the
falloff starts at the window there too (``EfficiencyMap.extension_ratio``).

The mechanical limits are checked on the real cycle and flagged, not enforced (``limit_flags``): high-side
pressure, suction pressure, discharge temperature and electrical work per swept volume (the shaft-torque load).
The same numbers are in local_scripts/sghp_live/sghp_model.js.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .density_basis import R600A_SAT_VAPOUR, R600A_SUCTION_15K

# R600A_ENVELOPE edge i -> i + 1 (compressor_data.R600A_ENVELOPE)
R600A_EDGE_KINDS = ("mech", "ratio_low", "mech", "mech", "ratio_high", "mech")
R600A_EDGE_NAMES = ("minimum condensing", "low-ratio", "maximum evaporating", "maximum condensing",
                    "high-ratio", "minimum evaporating")

# Mechanical limits from R600a over the window at 15 K suction superheat with no shell loss: high-side pressure
# (dew at 120 C), suction pressure (dew at -10 C), discharge temperature and electrical work per swept volume,
# each the window's maximum (minimum for the suction pressure), rounded up.
R600A_LIMITS = {
    "Recip-S": dict(P_dis_max_bar=28.4, P_suc_min_bar=1.08, T_dis_max=140.0, w_swept_max_kJm3=1085.0),
    "Recip-M": dict(P_dis_max_bar=28.4, P_suc_min_bar=1.08, T_dis_max=135.0, w_swept_max_kJm3=1015.0),
    "Recip-L": dict(P_dis_max_bar=28.4, P_suc_min_bar=1.08, T_dis_max=137.0, w_swept_max_kJm3=1015.0),
}


def map_density_ratio(to: float, tc: float) -> float:
    """R600a density ratio at a map point: condensing saturated vapour over the 15 K-superheat suction."""
    def cl(t, ref):
        lo, hi = ref.T_range
        return min(max(t, lo), hi)
    return R600A_SAT_VAPOUR.rho_dew(cl(tc, R600A_SAT_VAPOUR)) / R600A_SUCTION_15K.rho_dew(cl(to, R600A_SUCTION_15K))


def ratio_line(envelope, kinds, kind: str, n: int = 16) -> list[tuple[float, float]]:
    """(to, density ratio) along the edges of one kind, by to."""
    pts = []
    m = len(envelope)
    for i, a in enumerate(envelope):
        if kinds[i] != kind:
            continue
        b = envelope[(i + 1) % m]
        for k in range(n + 1):
            f = k / n
            to, tc = a[0] + f * (b[0] - a[0]), a[1] + f * (b[1] - a[1])
            pts.append((to, map_density_ratio(to, tc)))
    return sorted(pts, key=lambda p: p[0])


def line_at(line, x: float) -> float:
    if x <= line[0][0]:
        return line[0][1]
    if x >= line[-1][0]:
        return line[-1][1]
    i = 0
    while line[i + 1][0] < x:
        i += 1
    (x0, y0), (x1, y1) = line[i], line[i + 1]
    return y1 if x1 == x0 else y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def nearest_on_envelope(to: float, tc: float, envelope) -> dict:
    """Nearest point of the window's boundary to (to, tc), and the edge it lies on."""
    best = None
    m = len(envelope)
    for i, a in enumerate(envelope):
        b = envelope[(i + 1) % m]
        dx, dy = b[0] - a[0], b[1] - a[1]
        t = min(max(((to - a[0]) * dx + (tc - a[1]) * dy) / (dx * dx + dy * dy), 0.0), 1.0)
        x, y = a[0] + t * dx, a[1] + t * dy
        d = math.hypot(to - x, tc - y)
        if best is None or d < best["d"] - 1e-12:
            best = dict(to=x, tc=y, d=d, edge=i, t=t)
    return best


@dataclass
class EnvelopeExtension:
    """How efficiency continues beyond a polygonal window (see the module docstring)."""

    envelope: tuple
    edge_kinds: tuple = R600A_EDGE_KINDS
    edge_names: tuple | None = None   # None: the R600a names for its kinds, else "edge i"
    beta: float = 0.4
    delta_low: float = 0.3
    eta_vol_floor: float = 0.02
    eta_is_floor: float = 0.05
    _lines: dict = field(default=None, init=False, repr=False)

    def __post_init__(self):
        if len(self.edge_kinds) != len(self.envelope):
            raise ValueError(f"edge_kinds needs one kind per envelope edge ({len(self.envelope)})")
        if self.edge_names is None:
            self.edge_names = (R600A_EDGE_NAMES if tuple(self.edge_kinds) == R600A_EDGE_KINDS
                               else tuple(f"edge {i}" for i in range(len(self.envelope))))
        elif len(self.edge_names) != len(self.envelope):
            raise ValueError(f"edge_names needs one name per envelope edge ({len(self.envelope)})")
        self._lines = dict(high=ratio_line(self.envelope, self.edge_kinds, "ratio_high"),
                           low=ratio_line(self.envelope, self.edge_kinds, "ratio_low"))

    def ratio_limits(self, to: float) -> tuple[float, float]:
        """(low, high) density-ratio limits at map evaporating temperature to; no edge of a kind, no limit."""
        low, high = self._lines["low"], self._lines["high"]
        return line_at(low, to) if low else -math.inf, line_at(high, to) if high else math.inf

    def efficiencies(self, eta_is, eta_vol, to: float, tc: float, r: float) -> dict:
        """eta_is, eta_vol at map point (to, tc) with density ratio r: the maps inside, the extension outside."""
        from .compressor_curves import in_envelope
        if in_envelope(to, tc, self.envelope):
            return dict(inside=True, eta_is=eta_is(to, tc), eta_vol=None if eta_vol is None else eta_vol(to, tc),
                        read=(to, tc), ratio=r, ratio_beyond=0.0)
        q = nearest_on_envelope(to, tc, self.envelope)
        eis = eta_is(q["to"], q["tc"])
        ev = None if eta_vol is None else eta_vol(q["to"], q["tc"])
        out = dict(inside=False, read=(q["to"], q["tc"]), edge=q["edge"], edge_kind=self.edge_kinds[q["edge"]],
                   edge_name=self.edge_names[q["edge"]], eta_is_edge=eis, eta_vol_edge=ev, ratio=r, ratio_beyond=0.0)
        r_lo, r_hi = self.ratio_limits(to)
        if r > r_hi:
            evq = ev if ev is not None else 0.5          # no volumetric map: taken as 0.5 at the edge
            ev_new = max(1 - (1 - evq) * (r - 1) / (r_hi - 1), self.eta_vol_floor)
            eis = eis * (ev_new / evq) ** self.beta
            if ev is not None:
                ev = ev_new
            out.update(ratio_limit=r_hi, ratio_beyond=r / r_hi - 1)
        elif r < r_lo:
            x, xl, d = max(r - 1, 0.0), r_lo - 1, self.delta_low
            eis = eis * (x / (x + d)) / (xl / (xl + d))
            if ev is not None:
                ev = min(ev, 1.0)
            out.update(ratio_limit=r_lo, ratio_beyond=r / r_lo - 1)
        out.update(eta_is=max(eis, self.eta_is_floor), eta_vol=ev)
        return out


def limit_detail(f: dict, rev_m3: float | None = None, shaft_fraction: float = 1.0) -> str:
    """A sentence for an exceeded limit: what it is, by how much, and why it matters.

    shaft_fraction: shaft over electrical work, turning the (electrical) work per swept volume into torque."""
    d = 2 if f["unit"] == "bar" else 0
    n = lambda x: f"{x:.{d}f}"
    pct = f"{f['over_pct']:.1f} %"
    k = f["key"]
    if k == "P_dis":
        return (f"High-side pressure {n(f['value'])} bar is {pct} beyond the operating window ({n(f['limit'])} bar): "
                f"check with the supplier that the casing, gaskets and valves are rated for it.")
    if k == "T_dis":
        return (f"Discharge temperature {n(f['value'])} °C is {f['value'] - f['limit']:.0f} K beyond the operating window "
                f"({n(f['limit'])} °C): check with the supplier that the oil and discharge valves are suitable.")
    if k == "w_swept":
        what = (f"Shaft torque ≈ {f['value'] * shaft_fraction * rev_m3 / (2 * math.pi) * 1000:.0f} N·m ({n(f['value'])} kJ/m³ swept)"
                if rev_m3 else f"Shaft load {n(f['value'])} kJ/m³ swept")
        return (f"{what} is {pct} beyond the operating window ({n(f['limit'])} kJ/m³): check with the supplier that the "
                f"motor, crank and bearings are suitable.")
    if k == "P_suc":
        return (f"Suction pressure {n(f['value'])} bar is {pct} below the operating window ({n(f['limit'])} bar): check "
                f"with the supplier for air ingress and motor cooling.")
    return f"{f['label']} {n(f['value'])} {f['unit']} is beyond the {n(f['limit'])} {f['unit']} limit."


def limit_flags(limits: dict | None, P_dis_bar: float, T_dis: float, w_swept_kJm3: float | None,
                P_suc_bar: float, rev_m3: float | None = None, shaft_fraction: float = 1.0) -> list[dict]:
    """Mechanical limits exceeded on the real cycle: [{key, label, value, limit, unit, kind, over_pct, detail}].
    rev_m3: swept volume per revolution (m3), for the shaft torque in the shaft-load detail; shaft_fraction:
    shaft over electrical work (w_swept is electrical work per swept volume)."""
    if not limits:
        return []
    checks = (("P_dis", "High-side pressure", P_dis_bar, limits.get("P_dis_max_bar"), "bar", 1),
              ("T_dis", "Discharge temperature", T_dis, limits.get("T_dis_max"), "°C", 1),
              ("w_swept", "Work per swept volume (shaft load)", w_swept_kJm3, limits.get("w_swept_max_kJm3"), "kJ/m³", 1),
              ("P_suc", "Suction pressure", P_suc_bar, limits.get("P_suc_min_bar"), "bar", -1))
    out = [dict(key=k, label=lab, value=v, limit=lim, unit=u, kind="max" if s > 0 else "min", over_pct=100 * abs(v / lim - 1))
           for k, lab, v, lim, u, s in checks
           if lim is not None and v is not None and math.isfinite(v) and s * (v - lim) > 0]
    for f in out:
        f["detail"] = limit_detail(f, rev_m3, shaft_fraction)
    return out
