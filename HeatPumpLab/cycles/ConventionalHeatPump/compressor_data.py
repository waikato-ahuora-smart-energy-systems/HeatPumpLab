"""Compressor efficiency maps from manufacturer data.

Three semi-hermetic reciprocating compressors on R600a (isobutane), named by size (Recip-S, Recip-M,
Recip-L), fitted as bicubic polynomials
in the dew-point evaporating and condensing temperatures (to, tc in degC); see
``compressor_curves.bicubic_polynomial`` for the coefficient order. Source: workbook
"R600a_compressor_parameters_2025-03-24" (volumetric and isentropic efficiency, coefficients 0-15),
with an estimated application envelope shared by the three compressors (read from its figure).

Basis and assumptions:
- The isentropic efficiency is on the electrical-input basis, m (h_is - h_suc) / P_el
  (confirmed by the user).
- Fitted for R600a, taken as rated at 15 K suction superheat. For any refrigerant (n-butane
  included) the fits are read on a density basis by default (``density_basis.py``): to* is the R600a
  evaporating temperature whose suction density at 15 K superheat equals the actual suction density,
  and tc* the R600a condensing temperature with the working fluid's saturated-vapour density at the
  condensing pressure. For R600a at 15 K superheat this is the original map; other superheats move
  to* (the secondary density effect). basis="saturated" matches the saturated vapour at the
  evaporating pressure instead (no superheat effect); basis="temperature" reads the fits at the
  working fluid's own (to, tc).
- The swept volume is not in the workbook. At 50 Hz (1450 rpm): Recip-S 13.42 m3/h, Recip-M
  56.25 m3/h and Recip-L 151.6 m3/h, from the manufacturer's published data.
"""
from __future__ import annotations

from .compressor_curves import EfficiencyMap, bicubic_polynomial
from .compressor_window import R600A_EDGE_KINDS, R600A_LIMITS, EnvelopeExtension
from .density_basis import R600A_SAT_VAPOUR, R600A_SUCTION_15K

# Estimated application envelope (to, tc) degC, R600a.
R600A_ENVELOPE = ((-10.0, 30.0), (27.0, 42.0), (48.0, 70.0), (48.0, 120.0), (6.0, 120.0), (-10.0, 97.0))

# Beyond the window (compressor_window): edge kinds, efficiency falloff past the ratio edges.
R600A_EXTENSION = EnvelopeExtension(R600A_ENVELOPE, R600A_EDGE_KINDS)

# Swept volume at 50 Hz (m3/h), where known.
DISPLACEMENT_50HZ_M3H = {"Recip-S": 13.42, "Recip-M": 56.25, "Recip-L": 151.6}

R600A_COEFFICIENTS = {
    "Recip-S": dict(
        volumetric=(
            1.12984011409572, -0.01000517868509311, 9.959408735935163e-05, -5.429511239235132e-07,
            0.003567133668441818, 9.479581042021442e-05, -9.565274823906224e-07, 7.170717115542697e-09,
            -0.0001232257031541566, 1.028269483251338e-06, -2.282557219145039e-08, 7.981615275880793e-11,
            2.550844912615862e-06, -6.231819262761143e-08, 8.149157481312754e-10, -3.246732661005363e-12,
        ),
        isentropic=(
            0.4116836162879849, 0.005039822811869186, -2.631708318422278e-05, -2.046577662134925e-08,
            -0.003314441060494978, 0.0001339324043644695, -8.108428787902513e-07, 2.761876674705595e-09,
            -0.0008637056504366176, 3.207454510381633e-05, -4.196854713876601e-07, 1.713676254812856e-09,
            4.34014145040436e-06, -3.48778095942875e-07, 5.507361160638747e-09, -2.42783380539514e-11,
        ),
    ),
    "Recip-M": dict(
        volumetric=(
            1.096611343131835, -0.00726227358913289, 5.309047868640223e-05, -3.096676393991497e-07,
            0.006733364512289938, -3.458999501344552e-05, 1.169742469456288e-06, -5.714608237170313e-09,
            -0.0002840846510202909, 6.967940930769206e-06, -9.7576674515825e-08, 4.673105300527226e-10,
            5.100400187957733e-06, -1.539976618200062e-07, 1.904689442801216e-09, -8.23309880717932e-12,
        ),
        isentropic=(
            0.2820090966150441, 0.0120891735264767, -0.0001194130317419781, 3.921449567398475e-07,
            0.001020597758718241, -0.0001552582461411895, 3.819292261440514e-06, -2.125798153420446e-08,
            -0.001021385474733074, 4.213732121164092e-05, -5.832829781115806e-07, 2.561191347258478e-09,
            8.812924593464052e-06, -5.206597937866331e-07, 7.768970794352923e-09, -3.446962012422449e-11,
        ),
    ),
    "Recip-L": dict(
        volumetric=(
            0.9456698645061833, -0.0009241333917138055, -1.981682414911695e-05, -2.974086664626571e-08,
            0.006439330476049707, -0.0001232201170140139, 2.390487484843521e-06, -8.830117868237834e-09,
            -0.0001643829586845571, 5.825566412853439e-06, -9.913209075311744e-08, 4.150372919033212e-10,
            2.535938980404551e-06, -9.615602095713242e-08, 1.436561992150825e-09, -6.097111927442945e-12,
        ),
        isentropic=(
            0.365732694101855, 0.008380402635408853, -7.40105429843498e-05, 1.995263626032545e-07,
            -0.0006020506345257311, -2.396568897230092e-05, 1.303701754995347e-06, -6.770570118114042e-09,
            -0.0009300664051148016, 3.286808772249003e-05, -3.98272290261303e-07, 1.544452281550714e-09,
            1.954486640420136e-06, -2.661964600946924e-07, 4.30163369026706e-09, -1.872496284512809e-11,
        ),
    ),
}


def recip_r600a(model: str = "Recip-M", displacement_m3h: float | None | str = "50Hz",
                 basis: str = "density", **kwargs) -> EfficiencyMap:
    """Efficiency map for a compressor in the R600a workbook ("Recip-S", "Recip-M", "Recip-L").

    displacement_m3h: swept volume (m3/h); "50Hz" uses DISPLACEMENT_50HZ_M3H; None gives a
    size-less map (results per kg of refrigerant).
    basis: "density" (default) reads the R600a fits at the R600a point with the working fluid's
    suction density (R600a at 15 K superheat) and condensing saturated-vapour density;
    "saturated" matches saturated-vapour densities on both sides; "temperature" reads them at the
    working fluid's own (to, tc).
    kwargs override EfficiencyMap fields (e.g. efficiency_basis, motor_eff, heat_loss_fraction).
    """
    if model not in R600A_COEFFICIENTS:
        raise KeyError(f"unknown compressor '{model}'; available: {sorted(R600A_COEFFICIENTS)}")
    c = R600A_COEFFICIENTS[model]
    if displacement_m3h == "50Hz":
        displacement_m3h = DISPLACEMENT_50HZ_M3H.get(model)
    if basis not in ("density", "saturated", "temperature"):
        raise ValueError("basis must be 'density', 'saturated' or 'temperature'")
    options = dict(envelope=R600A_ENVELOPE, name=model, extension=R600A_EXTENSION, limits=dict(R600A_LIMITS[model]),
                   density_reference=None if basis == "temperature" else R600A_SAT_VAPOUR,
                   suction_reference=R600A_SUCTION_15K if basis == "density" else None) | kwargs
    if "envelope" in kwargs and "extension" not in kwargs:
        # the extension's edge kinds belong to the R600a window; another envelope has none unless given
        options["extension"] = None
    return EfficiencyMap(eta_is=bicubic_polynomial(c["isentropic"]),
                         eta_vol=bicubic_polynomial(c["volumetric"]),
                         displacement_m3h=displacement_m3h, **options)


def recip_models() -> list[str]:
    """Compressors in the R600a workbook, smallest first."""
    return sorted(R600A_COEFFICIENTS, key=lambda m: DISPLACEMENT_50HZ_M3H[m])
