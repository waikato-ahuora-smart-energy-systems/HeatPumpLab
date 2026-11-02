"""Screening estimate of the refrigerant charge of a solved cycle, as a range, without a heat exchanger design.

Charge is density times internal volume (Dawoud, Mugnini and Arteconi, Energy Conversion and Management 372 (2027)
122328, Eq. 20). The cycle fixes the duty, the temperature driving force and the refrigerant state along each exchanger.
The area follows from an assumed overall coefficient U, and the volume from an assumed internal volume per unit of that
area, both by exchanger service in ``GEOMETRY``. Each exchanger is split into equal-duty segments:

    dA = dQ / (U dT),  dm = dA (V/A) rho

The specific charge (g per kW of heating) does not depend on the capacity. The range comes from:

    low      compact geometry, homogeneous two-phase density
    central  typical geometry, Zivi (1964) void fraction
    high     bulky geometry, Zivi void fraction

The system figures divide the exchanger charge by the share the exchangers hold of the total, for the compressor oil,
piping and receivers. ``EXCHANGER_SHARE`` spans what that paper reports: more than two-thirds in low-charge systems, and
43.5 % in one of its optimised designs. Those are small propane units with PAG oil, so at other scales and with other
fluids and oils the system figures are an allowance, not a prediction.

GEOMETRY holds engineering assumptions, not fitted or validated data. The Table 6 geometry of that paper (a plate
condenser of about 0.6 mm V/A and a 7 mm fin-and-tube evaporator of about 0.09 mm) falls inside the compact to typical
range. Pass ``geometry`` to replace them. A zeotrope's circulating composition can differ from its charged composition,
which this does not model: densities use the circulating composition throughout. Saturated densities at that overall
composition stand in for the local liquid and vapour phases in the void fraction.

Units: U in W/m2/K, V/A in m, duties in kW (or kJ per kg of basis flow), specific charge in g/kW.
"""

from dataclasses import asdict, dataclass, field

import numpy as np

CASES = ("compact", "typical", "bulky")
# case: service: (U, V/A). Air: U on the finned area. Recuperator: V/A of each refrigerant side.
GEOMETRY = {
    "compact": {"air": (60.0, 0.05e-3), "liquid": (2500.0, 0.5e-3), "recuperator": (800.0, 0.5e-3)},
    "typical": {"air": (40.0, 0.15e-3), "liquid": (1500.0, 1.0e-3), "recuperator": (500.0, 1.0e-3)},
    "bulky": {"air": (25.0, 0.30e-3), "liquid": (800.0, 3.0e-3), "recuperator": (250.0, 3.0e-3)},
}
EXCHANGER_SHARE = {"low": 2.0 / 3.0, "central": 0.55, "high": 0.435}
# bound: (geometry case, density index in a ChargeSection's rho pairs: 0 homogeneous, 1 void fraction)
BOUNDS = {"low": ("compact", 0), "central": ("typical", 1), "high": ("bulky", 1)}
N_SEGMENTS = 24


def void_fraction_density(rho: float, rho_l: float, rho_v: float) -> float:
    """Zivi two-phase density from the homogeneous density and the saturated liquid and vapour densities."""
    x = min(max((1.0 / rho - 1.0 / rho_l) / (1.0 / rho_v - 1.0 / rho_l), 0.0), 1.0)
    if x == 0.0:
        return rho_l
    if x == 1.0:
        return rho_v
    alpha = 1.0 / (1.0 + (1.0 - x) / x * (rho_v / rho_l) ** (2.0 / 3.0))
    return alpha * rho_v + (1.0 - alpha) * rho_l


@dataclass
class ChargeSection:
    """One exchanger: equal-duty segments, each with a driving force (K) and, per refrigerant side, a (homogeneous,
    void fraction) density pair (kg/m3)."""

    name: str
    service: str          # a key of the GEOMETRY cases: "air", "liquid" or "recuperator"
    duty: float
    dT: list
    rho: list


@dataclass
class ChargeEstimate:
    """Specific charge, g per kW of heating. ``exchanger_*`` holds the exchangers only; ``system_*`` adds the allowance
    for the oil, piping and receivers. ``by_exchanger`` gives each exchanger's g/kW and, for the typical geometry, its
    area and internal volume per kW of heating."""

    exchanger_low: float
    exchanger_central: float
    exchanger_high: float
    system_low: float
    system_central: float
    system_high: float
    by_exchanger: dict = field(default_factory=dict)
    warnings: list = field(default_factory=list)

    def system_kg(self, Q_kW: float) -> tuple:
        """The system (low, central, high) charge in kg at a heating capacity in kW."""
        return tuple(g * Q_kW / 1000.0 for g in (self.system_low, self.system_central, self.system_high))

    def to_dict(self):
        return asdict(self)


def estimate_charge(sections: list, Q_heat: float, *, mixture: bool = False, geometry: dict | None = None) -> ChargeEstimate:
    """The range for exchangers ``sections``, per ``Q_heat`` of heating (the same units as the sections' duties)."""
    if not Q_heat > 0.0:
        raise ValueError("charge estimate needs a positive heating duty")
    geometry = GEOMETRY if geometry is None else geometry
    totals, by = dict.fromkeys(BOUNDS, 0.0), {}
    for sec in sections:
        dT = np.asarray(sec.dT, dtype=float)
        if dT.size == 0 or dT.min() <= 0.0:
            raise ValueError(f"{sec.name}: charge estimate needs a positive driving force")
        dQ = sec.duty / dT.size
        row = {}
        for bound, (case, density) in BOUNDS.items():
            U, v_per_a = geometry[case][sec.service]
            dA = 1000.0 * dQ / (U * dT)                    # m2 per (kW of the duty's basis)
            kg = sum(float(np.sum(dA * v_per_a * np.array([p[density] for p in side]))) for side in sec.rho)
            row[bound] = 1000.0 * kg / Q_heat
            totals[bound] += row[bound]
        U, v_per_a = geometry["typical"][sec.service]
        area = float(np.sum(1000.0 * dQ / (U * dT))) / Q_heat
        row.update(area_m2_per_kW=area, volume_L_per_kW=area * v_per_a * len(sec.rho) * 1000.0)
        by[sec.name] = row
    warnings = ["screening range from assumed U and V/A, not a design or filling quantity"]
    if mixture:
        warnings.append("zeotrope: the composition shift between charged and circulating refrigerant is not modelled")
    return ChargeEstimate(*(totals[b] for b in BOUNDS), *(totals[b] / EXCHANGER_SHARE[b] for b in BOUNDS),
                          by_exchanger=by, warnings=warnings)


def _pair(st: dict, saturated) -> tuple:
    """(homogeneous, void fraction) density of a state; ``saturated(q)`` gives the saturated density at its pressure."""
    rho = st["D"]
    if st.get("phase") != "two_phase":
        return rho, rho
    return rho, void_fraction_density(rho, saturated(0.0), saturated(1.0))


def _stream_pair(stream, memo: dict) -> tuple:
    def saturated(q):
        key = (stream.pressure, q)
        if key not in memo:
            memo[key] = stream.PP.get_std_properties("PQ", stream.pressure, q, stream.fld_composition)["D"]
        return memo[key]
    return _pair({"D": stream.density, "phase": stream.phase}, saturated)


def _mean_pair(a: tuple, b: tuple) -> tuple:
    return 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1])


def stream_section(hx, refrigerant: str, n: int = N_SEGMENTS) -> ChargeSection:
    """A section from a solved component exchanger (``BaseHeatExchanger``): its profiles run from each side's inlet,
    linear in enthalpy, so sink segment i faces source segment n - 1 - i in counterflow."""
    data, memo = hx.get_graph_data(n), {}
    source, sink = data["source"], data["sink"]
    src_ref, snk_ref = hx.source_in.fluid == refrigerant, hx.sink_in.fluid == refrigerant
    if src_ref and snk_ref:
        service = "recuperator"
    else:
        secondary = hx.sink_in if src_ref else hx.source_in
        service = "air" if secondary.density < 100.0 else "liquid"
    dT, rho_src, rho_snk = [], [], []
    for i in range(n):
        a, b = source[n - 1 - i], source[n - i]
        c, d = sink[i], sink[i + 1]
        dT.append(0.5 * ((a.temperature - d.temperature) + (b.temperature - c.temperature)))
        if src_ref:
            rho_src.append(_mean_pair(_stream_pair(a, memo), _stream_pair(b, memo)))
        if snk_ref:
            rho_snk.append(_mean_pair(_stream_pair(c, memo), _stream_pair(d, memo)))
    return ChargeSection(hx.name, service, hx.duty, dT, [side for side in (rho_src, rho_snk) if side])


def loop_section(fl, exchanger: dict, path: dict, n: int = N_SEGMENTS) -> ChargeSection:
    """A section from a ``LoopResult`` exchanger, whose hot and cold [q, T] profiles run from the cold end in
    counterflow, and the refrigerant's isobaric path through it ([T, P bar, h, s] points)."""
    hot, cold = np.asarray(exchanger["hot"], dtype=float), np.asarray(exchanger["cold"], dtype=float)
    duty = float(hot[-1, 0])
    service = "air" if "air" in str(exchanger.get("other")) else "liquid"
    P, hs = path["pts"][0][1] / 10.0, [p[2] for p in path["pts"]]
    h_cold, h_hot = min(hs), max(hs)
    memo = {}

    def saturated(q):
        return memo.setdefault(q, fl.PQ(P, q)["D"])
    f = (np.arange(n) + 0.5) / n
    dT = np.interp(f * duty, hot[:, 0], hot[:, 1]) - np.interp(f * duty, cold[:, 0], cold[:, 1])
    rho = [_pair(fl.PH(P, h_cold + x * (h_hot - h_cold)), saturated) for x in f]
    return ChargeSection(exchanger["name"], service, duty, dT.tolist(), [rho])


def estimate_rhp_charge(hp, *, geometry: dict | None = None) -> ChargeEstimate:
    """The range for a solved ``RecuperativeHP``."""
    if not hp.converged:
        raise ValueError("charge estimate needs a converged design")
    active = {id(hx) for hx in list(hp.gas_coolers) + list(hp.evaporators)}
    sections = []
    for attr in ("gas_cooler1", "gas_cooler2", "gas_cooler3", "gas_cooler4", "evaporator", "evap_sup", "rhx1", "rhx2",
                 "rhx3", "ihx"):
        hx = getattr(hp, attr, None)
        if hx is None or hx.bypassed or hx.duty <= hx.tol or (attr[:3] in ("gas", "eva") and id(hx) not in active):
            continue
        section = stream_section(hx, hp.refrigerant)
        section.name = attr
        sections.append(section)
    return estimate_charge(sections, sum(gc.duty for gc in hp.gas_coolers), mixture=";" in hp.refrigerant,
                           geometry=geometry)


def estimate_loop_charge(fl, result, *, geometry: dict | None = None) -> ChargeEstimate:
    """The range for a successful single-fluid ``LoopResult``. Each exchanger takes the path of the same name, as a gas
    cooler where the design is supercritical."""
    if not result.success:
        raise ValueError("charge estimate needs a successful design")
    paths = {p["name"]: p for loop in result.paths.values() for p in loop}
    sections = []
    for ex in result.exchangers:
        name = ex["name"].removesuffix(" to air")
        sections.append(loop_section(fl, ex, paths.get(name) or paths[f"{name} / gas cooler"]))
    return estimate_charge(sections, result.q_air, mixture=";" in fl.fluid, geometry=geometry)
