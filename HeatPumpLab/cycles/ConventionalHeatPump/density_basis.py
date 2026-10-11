"""Density basis for compressor efficiency maps.

A compressor map fitted for one refrigerant (the reference, here R600a) is applied to another by
matching densities rather than temperatures. An operating point with dew-point evaporating and
condensing temperatures (to, tc) on the working fluid is evaluated at the reference-fluid point
(to*, tc*) with

    rho_ref,suction(to*) = rho_fluid,suction        (actual compressor suction state)
    rho''_ref(tc*)       = rho''_fluid(P_cond)      (saturated vapour at the condensing pressure)

where rho_ref,suction(to*) is the R600a density at its dew pressure for to* and the rated suction
superheat (15 K: the map is taken to have been rated at to + 15 K suction). The suction match
carries the secondary effect of superheat: more superheat, a thinner suction gas, a lower to*.
For R600a at 15 K superheat to* = to and tc* = tc, so the original map is recovered exactly.
``equivalent_point`` without a suction reference matches the saturated vapour at the evaporating
pressure instead (no superheat effect).

Each reference curve is a table interpolated with a monotone cubic (PCHIP) in ln rho, so the
Python model and the browser model use identical numbers. The R600a tables come from CoolProp
8.0.0 (HEOS, IsoButane: Buecker and Wagner 2006, the equation REFPROP also uses). Interpolation
error: saturated vapour < 4e-7 up to 100 C, < 2e-5 up to 125 C, < 2e-3 within 0.5 K of the
critical point; suction at 15 K superheat < 1e-6 throughout. Densities outside a table are
clamped to its ends with a warning.
"""
from __future__ import annotations

import bisect
import math
from dataclasses import dataclass, field

# R600a (isobutane) saturated vapour density, kg/m3, against temperature, degC.
R600A_SAT_VAPOUR_T = (
    -50.0, -49.0, -48.0, -47.0, -46.0, -45.0, -44.0, -43.0, -42.0, -41.0,
    -40.0, -39.0, -38.0, -37.0, -36.0, -35.0, -34.0, -33.0, -32.0, -31.0,
    -30.0, -29.0, -28.0, -27.0, -26.0, -25.0, -24.0, -23.0, -22.0, -21.0,
    -20.0, -19.0, -18.0, -17.0, -16.0, -15.0, -14.0, -13.0, -12.0, -11.0,
    -10.0, -9.0, -8.0, -7.0, -6.0, -5.0, -4.0, -3.0, -2.0, -1.0,
    0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
    10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0,
    20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,
    30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0,
    40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0, 49.0,
    50.0, 51.0, 52.0, 53.0, 54.0, 55.0, 56.0, 57.0, 58.0, 59.0,
    60.0, 61.0, 62.0, 63.0, 64.0, 65.0, 66.0, 67.0, 68.0, 69.0,
    70.0, 71.0, 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0,
    80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0, 88.0, 89.0,
    90.0, 91.0, 92.0, 93.0, 94.0, 95.0, 96.0, 97.0, 98.0, 99.0,
    100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0,
    110.0, 111.0, 112.0, 113.0, 114.0, 115.0, 116.0, 117.0, 118.0, 119.0,
    120.0, 121.0, 122.0, 123.0, 124.0, 125.0, 125.5, 126.0, 126.5, 127.0,
    127.5, 128.0, 128.5, 129.0, 129.5, 130.0, 130.5, 131.0, 131.5, 132.0,
    132.25, 132.5, 132.75, 133.0, 133.25, 133.5, 133.75, 134.0, 134.25, 134.5,
)
R600A_SAT_VAPOUR_RHO = (
    0.53214474, 0.56059848, 0.5902529, 0.62114383, 0.65330766, 0.68678135, 0.72160239, 0.75780885,
    0.79543935, 0.83453306, 0.87512968, 0.9172695, 0.96099334, 1.0063426, 1.0533591, 1.1020854,
    1.1525646, 1.2048401, 1.2589562, 1.3149575, 1.3728892, 1.4327973, 1.4947279, 1.5587282,
    1.6248456, 1.6931281, 1.7636244, 1.8363838, 1.9114561, 1.9888916, 2.0687415, 2.1510572,
    2.2358911, 2.323296, 2.4133254, 2.5060334, 2.6014747, 2.6997047, 2.8007796, 2.9047561,
    3.0116916, 3.1216443, 3.2346731, 3.3508375, 3.4701978, 3.5928153, 3.7187517, 3.8480697,
    3.9808328, 4.1171053, 4.2569524, 4.4004401, 4.5476352, 4.6986057, 4.8534203, 5.0121488,
    5.1748618, 5.3416312, 5.5125297, 5.6876312, 5.8670106, 6.050744, 6.2389089, 6.4315835,
    6.6288477, 6.8307824, 7.0374701, 7.2489944, 7.4654404, 7.6868946, 7.9134452, 8.1451817,
    8.3821953, 8.6245789, 8.8724271, 9.1258362, 9.3849044, 9.6497318, 9.9204203, 10.197074,
    10.479799, 10.768704, 11.0639, 11.365498, 11.673616, 11.98837, 12.309882, 12.638274,
    12.973673, 13.316208, 13.666011, 14.023218, 14.387967, 14.7604, 15.140664, 15.528908,
    15.925286, 16.329955, 16.743078, 17.164821, 17.595355, 18.034858, 18.483509, 18.941496,
    19.409011, 19.886253, 20.373426, 20.870741, 21.378416, 21.896675, 22.425752, 22.965885,
    23.517324, 24.080325, 24.655155, 25.242088, 25.841411, 26.453418, 27.078418, 27.716728,
    28.36868, 29.034617, 29.714896, 30.409889, 31.119982, 31.84558, 32.587101, 33.344984,
    34.119686, 34.911684, 35.721477, 36.549588, 37.396563, 38.262976, 39.149427, 40.056548,
    40.985003, 41.935488, 42.90874, 43.905533, 44.926686, 45.973065, 47.045585, 48.145217,
    49.27299, 50.430001, 51.617414, 52.836473, 54.088505, 55.374931, 56.697274, 58.057169,
    59.456379, 60.896802, 62.380493, 63.909679, 65.486782, 67.114439, 68.795537, 70.533241,
    72.331035, 74.192771, 76.122721, 78.125648, 80.206886, 82.372443, 84.629123, 86.984684,
    89.448035, 92.029492, 94.741105, 97.597099, 100.61446, 103.81377, 107.22031, 110.86578,
    112.79036, 114.79079, 116.87418, 119.0488, 121.3244, 123.71257, 126.22725, 128.88551,
    131.70853, 134.72319, 137.96442, 141.47894, 145.33162, 149.61695, 151.96461, 154.48162,
    157.20256, 160.17507, 163.46775, 167.18532, 171.50123, 176.73987, 183.64597, 194.90209,
)
# R600a density, kg/m3, at its dew pressure for T and T + 15 K (the rated suction state), same T.
R600A_RATED_SUPERHEAT_K = 15.0
R600A_SUCTION_15K_RHO = (
    0.49750999, 0.52421361, 0.55204793, 0.58104671, 0.61124424, 0.64267531, 0.67537522, 0.70937977,
    0.74472529, 0.78144859, 0.81958697, 0.85917827, 0.90026077, 0.9428733, 0.98705513, 1.0328461,
    1.0802864, 1.1294168, 1.1802787, 1.2329136, 1.2873639, 1.3436723, 1.401882, 1.4620366,
    1.5241803, 1.5883578, 1.6546143, 1.7229953, 1.7935471, 1.8663162, 1.9413498, 2.0186955,
    2.0984015, 2.1805164, 2.2650893, 2.35217, 2.4418086, 2.5340559, 2.6289632, 2.7265821,
    2.8269652, 2.9301652, 3.0362356, 3.1452305, 3.2572045, 3.3722128, 3.4903111, 3.6115559,
    3.7360042, 3.8637136, 3.9947424, 4.1291496, 4.2669948, 4.4083383, 4.553241, 4.7017646,
    4.8539717, 5.0099253, 5.1696894, 5.3333287, 5.5009086, 5.6724955, 5.8481566, 6.0279598,
    6.2119741, 6.4002692, 6.5929159, 6.7899858, 6.9915515, 7.1976867, 7.4084661, 7.6239654,
    7.8442615, 8.0694322, 8.2995567, 8.5347152, 8.7749893, 9.0204617, 9.2712165, 9.5273391,
    9.7889162, 10.056036, 10.328788, 10.607264, 10.891556, 11.181759, 11.477968, 11.780281,
    12.088798, 12.40362, 12.724851, 13.052594, 13.386958, 13.728051, 14.075984, 14.430872,
    14.79283, 15.161976, 15.538431, 15.922318, 16.313761, 16.712891, 17.119837, 17.534734,
    17.957719, 18.388932, 18.828515, 19.276615, 19.733382, 20.198969, 20.673533, 21.157234,
    21.650237, 22.152709, 22.664824, 23.186759, 23.718693, 24.260814, 24.813311, 25.37638,
    25.95022, 26.535039, 27.131046, 27.738458, 28.357499, 28.988396, 29.631385, 30.286707,
    30.954611, 31.635351, 32.329191, 33.036402, 33.75726, 34.492053, 35.241076, 36.004633,
    36.783037, 37.57661, 38.385686, 39.210608, 40.05173, 40.909416, 41.784045, 42.676005,
    43.585698, 44.51354, 45.459958, 46.425396, 47.410312, 48.415179, 49.440486, 50.486739,
    51.554461, 52.644194, 53.756497, 54.891951, 56.051153, 57.234726, 58.443311, 59.677574,
    60.938203, 62.225911, 63.541437, 64.885546, 66.25903, 67.662711, 69.097441, 70.564102,
    72.06361, 73.596919, 75.165016, 76.768932, 78.409742, 80.088569, 81.806592, 83.565055,
    84.459862, 85.365282, 86.281496, 87.208691, 88.147064, 89.096824, 90.05819, 91.031395,
    92.016693, 93.014358, 94.024694, 95.048041, 96.08479, 97.135402, 97.666078, 98.20044,
    98.738587, 99.280631, 99.826706, 100.37698, 100.93168, 101.49111, 102.05575, 102.62656,
)
R600A_T_CRIT = 134.66      # degC
R600A_RHO_CRIT = 225.5     # kg/m3


def _pchip_end(h0, h1, d0, d1):
    v = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
    if v * d0 <= 0:
        return 0.0
    if d0 * d1 <= 0 and abs(v) > abs(3 * d0):
        return 3 * d0
    return v


@dataclass(frozen=True)
class SaturatedVapourDensity:
    """A density curve rho(T) rising with T (saturated vapour, or vapour at a fixed superheat above
    the dew point), monotone-cubic in ln rho, with its inverse. ``rho_dew`` / ``T_dew`` are indexed
    by the dew-point (evaporating or condensing) temperature."""

    T: tuple
    rho: tuple
    fluid: str = "R600a"
    _y: tuple = field(init=False, repr=False)
    _m: tuple = field(init=False, repr=False)

    def __post_init__(self):
        T, rho = self.T, self.rho
        if len(T) != len(rho) or len(T) < 3:
            raise ValueError("T and rho need the same length, at least 3")
        if any(b <= a for a, b in zip(T, T[1:])) or any(b <= a for a, b in zip(rho, rho[1:])):
            raise ValueError("T and rho must both increase strictly")
        y = [math.log(r) for r in rho]
        h = [b - a for a, b in zip(T, T[1:])]
        d = [(y[i + 1] - y[i]) / h[i] for i in range(len(h))]
        m = [0.0] * len(T)
        for i in range(1, len(T) - 1):
            if d[i - 1] * d[i] > 0:
                w1, w2 = 2 * h[i] + h[i - 1], h[i] + 2 * h[i - 1]
                m[i] = (w1 + w2) / (w1 / d[i - 1] + w2 / d[i])
        m[0] = _pchip_end(h[0], h[1], d[0], d[1])
        m[-1] = _pchip_end(h[-1], h[-2], d[-1], d[-2])
        object.__setattr__(self, "_y", tuple(y))
        object.__setattr__(self, "_m", tuple(m))

    @property
    def T_range(self) -> tuple[float, float]:
        return self.T[0], self.T[-1]

    @property
    def rho_range(self) -> tuple[float, float]:
        return self.rho[0], self.rho[-1]

    def rho_dew(self, T: float) -> float:
        """Saturated vapour density (kg/m3) at T (degC), within the table range."""
        lo, hi = self.T_range
        if not lo <= T <= hi:
            raise ValueError(f"{self.fluid} table covers {lo} to {hi} C, not {T}")
        return math.exp(self._ln_rho(T))

    def _ln_rho(self, T: float) -> float:
        n = len(self.T)
        i = min(max(bisect.bisect_right(self.T, T) - 1, 0), n - 2)
        h = self.T[i + 1] - self.T[i]
        s = (T - self.T[i]) / h
        h00, h10 = 2 * s ** 3 - 3 * s ** 2 + 1, s ** 3 - 2 * s ** 2 + s
        h01, h11 = -2 * s ** 3 + 3 * s ** 2, s ** 3 - s ** 2
        return h00 * self._y[i] + h10 * h * self._m[i] + h01 * self._y[i + 1] + h11 * h * self._m[i + 1]

    def T_dew(self, rho: float) -> tuple[float, bool]:
        """Temperature (degC) at which the reference vapour has density rho, and whether rho was
        outside the table (then the nearest end is returned)."""
        if not rho > 0:
            raise ValueError(f"density must be positive, not {rho}")
        if rho <= self.rho[0]:
            return self.T[0], rho < self.rho[0]
        if rho >= self.rho[-1]:
            return self.T[-1], rho > self.rho[-1]
        target = math.log(rho)
        i = bisect.bisect_right(self.rho, rho) - 1
        if self.rho[i] == rho:
            return self.T[i], False
        a, b = self.T[i], self.T[i + 1]          # ln rho(a) < target < ln rho(b); monotone between
        for _ in range(100):
            mid = 0.5 * (a + b)
            if self._ln_rho(mid) <= target:
                a = mid
            else:
                b = mid
            if b - a < 1e-11:
                break
        return 0.5 * (a + b), False


R600A_SAT_VAPOUR = SaturatedVapourDensity(R600A_SAT_VAPOUR_T, R600A_SAT_VAPOUR_RHO, "R600a")
R600A_SUCTION_15K = SaturatedVapourDensity(R600A_SAT_VAPOUR_T, R600A_SUCTION_15K_RHO, "R600a at 15 K suction superheat")


def equivalent_point(reference: SaturatedVapourDensity, rho_evap: float | None, rho_cond: float,
                     suction_reference: SaturatedVapourDensity | None = None, rho_suction: float | None = None):
    """Reference-fluid (to*, tc*) plus warnings for clamping.

    tc* from the saturated-vapour density at the condensing pressure (on ``reference``). to* from
    the actual suction density on ``suction_reference`` when given, else from the saturated-vapour
    density at the evaporating pressure on ``reference``.
    """
    if suction_reference is not None:
        if rho_suction is None:
            raise ValueError("a suction-density basis needs the suction density")
        evap_ref, rho_e, evap_name = suction_reference, rho_suction, "suction"
    else:
        evap_ref, rho_e, evap_name = reference, rho_evap, "evaporating vapour"
    to_eq, lo_clamped = evap_ref.T_dew(rho_e)
    tc_eq, hi_clamped = reference.T_dew(rho_cond)
    warnings = []
    for name, ref, clamped, rho, T in ((evap_name, evap_ref, lo_clamped, rho_e, to_eq),
                                       ("condensing vapour", reference, hi_clamped, rho_cond, tc_eq)):
        if clamped:
            warnings.append(f"{name} density {rho:.4g} kg/m3 is beyond the {ref.fluid} range "
                            f"({ref.rho[0]:.4g} to {ref.rho[-1]:.4g} kg/m3): map read at {T:g} C")
    return to_eq, tc_eq, warnings


# ---------------------------------------------------------------- maps in a working fluid's (to, tc)


def fluid_curve(rho_of_T, T_lo: float, T_hi: float, step: float = 1.0, fluid: str = "") -> SaturatedVapourDensity:
    """Tabulate a working fluid's dew-point vapour density rho_of_T(T) from T_lo up to T_hi (degC).

    Stops early where a flash fails or the density stops rising (near a critical point), so the
    table is always invertible. Used to draw a density-basis map in the fluid's own (to, tc).
    """
    T, rho = [], []
    n = int(math.floor((T_hi - T_lo) / step + 1e-9))
    for k in range(n + 1):
        t = T_lo + k * step
        try:
            r = float(rho_of_T(t))
        except Exception:
            break
        if not r > 0 or not math.isfinite(r) or (rho and r <= rho[-1]):
            break
        T.append(t)
        rho.append(r)
    if len(T) < 3:
        raise ValueError(f"could not tabulate the saturated vapour density of {fluid or 'the fluid'}")
    return SaturatedVapourDensity(tuple(T), tuple(rho), fluid)


def to_reference(curve: SaturatedVapourDensity, reference: SaturatedVapourDensity, T: float):
    """Reference-fluid temperature with the working fluid's vapour density at T, or None outside
    the working-fluid table; the flag says whether the reference table clamped."""
    lo, hi = curve.T_range
    if not lo <= T <= hi:
        return None, False
    return reference.T_dew(curve.rho_dew(T))


def from_reference(curve: SaturatedVapourDensity, reference: SaturatedVapourDensity, T_ref: float):
    """Working-fluid temperature with the reference fluid's vapour density at T_ref, or None when
    that density is outside the working-fluid table."""
    T, clamped = curve.T_dew(reference.rho_dew(T_ref))
    return None if clamped else T


def envelope_in_fluid(envelope, reference: SaturatedVapourDensity, curve: SaturatedVapourDensity,
                      step: float = 1.0, evap_reference: SaturatedVapourDensity | None = None,
                      evap_curve: SaturatedVapourDensity | None = None) -> list:
    """A reference-fluid (to, tc) envelope redrawn in the working fluid's (to, tc). Edges are
    sampled every ``step`` K (they are not straight in the new coordinates); points whose density
    is outside the working-fluid table are dropped. The evaporating axis uses ``evap_reference`` /
    ``evap_curve`` (suction densities) when given, else ``reference`` / ``curve``."""
    er, ec = (evap_reference, evap_curve) if evap_reference is not None else (reference, curve)
    pts = []
    for (x1, y1), (x2, y2) in zip(envelope, list(envelope[1:]) + [envelope[0]]):
        n = max(1, int(math.ceil(max(abs(x2 - x1), abs(y2 - y1)) / step)))
        for k in range(n):
            a = x1 + (x2 - x1) * k / n
            b = y1 + (y2 - y1) * k / n
            fa, fb = from_reference(ec, er, a), from_reference(curve, reference, b)
            if fa is not None and fb is not None:
                pts.append((fa, fb))
    return pts
