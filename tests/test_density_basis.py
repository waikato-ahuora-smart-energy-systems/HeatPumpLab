"""Density basis for the R600a efficiency maps (HeatPumpLab/cycles/ConventionalHeatPump/density_basis.py).

The maps are fitted for R600a in its dew-point (to, tc), taken as rated at 15 K suction superheat.
Another refrigerant (or another superheat) is read at the R600a point (to*, tc*) with

    rho_R600a(P_dew(to*), to* + 15 K) = actual suction density of the working fluid
    rho''_R600a(tc*)                  = rho''_fluid at the condensing pressure

These tests prove the recomputed maps are right:

- both R600a density tables and their inverses are exact round trips and match the equation of
  state (CoolProp / REFPROP, where installed);
- for R600a at 15 K superheat the recomputed map IS the original map, in the map and in the cycle;
- for synthetic fluids with a known density relation the equivalent point is known analytically
  (curves shifted by 12 K read the map 12 K lower, exactly), for points, grids and the envelope;
- in the cycle, the map point satisfies both density equalities, every efficiency is the fit
  evaluated there, and the envelope is checked there;
- superheat acts only through the suction density: more superheat, lower to*, same tc*;
- the browser model (sghp_model.js) gives the same numbers as Python.
"""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

from HeatPumpLab.cycles.ConventionalHeatPump import (
    ConventionalHeatPump,
    ConventionalHPConfig,
    R600A_COEFFICIENTS,
    R600A_ENVELOPE,
    R600A_SAT_VAPOUR,
    R600A_SUCTION_15K,
    SaturatedVapourDensity,
    bicubic_polynomial,
    recip_r600a,
    envelope_in_fluid,
    equivalent_point,
    fluid_curve,
    in_envelope,
)
from HeatPumpLab.cycles.ConventionalHeatPump.density_basis import (
    R600A_SAT_VAPOUR_RHO,
    R600A_SAT_VAPOUR_T,
    R600A_SUCTION_15K_RHO,
)
from test_conventional_heat_pump import ToyFluid, model

ROOT = Path(__file__).resolve().parents[1]
R, S = R600A_SAT_VAPOUR, R600A_SUCTION_15K           # condensing (saturated) and suction (15 K) references
ETA_IS = bicubic_polynomial(R600A_COEFFICIENTS["Recip-M"]["isentropic"])
ETA_VOL = bicubic_polynomial(R600A_COEFFICIENTS["Recip-M"]["volumetric"])
GRID = [(to, tc) for to in (-10.0, 0.0, 10.0, 20.0, 30.0, 40.0, 48.0) for tc in (40.0, 60.0, 80.0, 100.0, 120.0)]


def r600a_point(to, tc):
    """Densities of R600a itself at (to, tc) with the rated 15 K suction superheat."""
    return dict(rho_evap_dew=R.rho_dew(to), rho_cond_dew=R.rho_dew(tc), rho_suction=S.rho_dew(to))


# ---------------------------------------------------------------- the R600a curves


def test_tables_are_monotone_and_span_the_map():
    T = R600A_SAT_VAPOUR_T
    for curve, rho in ((R, R600A_SAT_VAPOUR_RHO), (S, R600A_SUCTION_15K_RHO)):
        assert curve.T == T and curve.rho == rho
        assert all(b > a for a, b in zip(T, T[1:])) and all(b > a for a, b in zip(rho, rho[1:]))
    assert T[0] <= -10.0 and T[-1] >= 120.0 and T[-1] < 134.66   # the envelope, below the critical point


def test_suction_table_is_thinner_than_saturated_vapour():
    for T in R600A_SAT_VAPOUR_T:                                  # 15 K superheat at the same pressure
        assert S.rho_dew(T) < R.rho_dew(T)
    assert R.rho_dew(50.0) == pytest.approx(17.595355, rel=1e-7) and S.rho_dew(50.0) == pytest.approx(16.313761, rel=1e-7)


def test_inverse_round_trips_everywhere():
    for curve in (R, S):
        lo, hi = curve.T_range
        worst = 0.0
        for k in range(4001):
            T = lo + (hi - lo) * k / 4000
            T_back, clamped = curve.T_dew(curve.rho_dew(T))
            assert not clamped
            worst = max(worst, abs(T_back - T))
        assert worst < 1e-9
        for T, rho in zip(curve.T, curve.rho):                   # nodes invert exactly
            assert curve.T_dew(rho) == (T, False)


def test_interpolant_is_monotone_between_nodes():
    for curve in (R, S):
        lo, hi = curve.T_range
        prev = 0.0
        for k in range(20001):
            r = curve.rho_dew(lo + (hi - lo) * k / 20000)
            assert r > prev
            prev = r


def test_out_of_range_densities_clamp_with_a_warning():
    assert R.T_dew(R.rho_range[0] * 0.5) == (R.T_range[0], True)
    assert R.T_dew(R.rho_range[1] * 1.5) == (R.T_range[1], True)
    to, tc, warnings = equivalent_point(R, None, R.rho_range[1] * 1.5, S, S.rho_range[0] * 0.5)
    assert (to, tc) == R.T_range and len(warnings) == 2
    assert "suction" in warnings[0] and "15 K" in warnings[0] and "condensing" in warnings[1]
    with pytest.raises(ValueError):
        R.T_dew(0.0)
    with pytest.raises(ValueError):
        R.rho_dew(R.T_range[1] + 1.0)
    with pytest.raises(ValueError):
        equivalent_point(R, None, 50.0, S, None)                  # suction basis needs the suction density


def test_bad_tables_are_rejected():
    with pytest.raises(ValueError):
        SaturatedVapourDensity((0.0, 1.0, 1.0), (1.0, 2.0, 3.0))
    with pytest.raises(ValueError):
        SaturatedVapourDensity((0.0, 1.0, 2.0), (1.0, 3.0, 2.0))


# ---------------------------------------------------------------- recomputed maps at a point


def expected_map(comp, to, tc):
    """The original maps inside the window; outside, the extension at the R600a density ratio."""
    if in_envelope(to, tc, R600A_ENVELOPE):
        return ETA_IS(to, tc), ETA_VOL(to, tc)
    w = comp.extension.efficiencies(ETA_IS, ETA_VOL, to, tc, R.rho_dew(tc) / S.rho_dew(to))
    return w["eta_is"], w["eta_vol"]


def test_r600a_at_15k_recovers_the_original_map():
    comp = recip_r600a("Recip-M")
    assert comp.density_reference is R and comp.suction_reference is S
    for to, tc in GRID:
        e = comp.efficiencies(to, tc, **r600a_point(to, tc))
        assert e["to_map"] == pytest.approx(to, abs=1e-9) and e["tc_map"] == pytest.approx(tc, abs=1e-9)
        assert (e["eta_is"], e["eta_vol"]) == pytest.approx(expected_map(comp, to, tc), abs=1e-10)
        assert e["in_envelope"] == in_envelope(to, tc, R600A_ENVELOPE)
        assert not e["warnings"]


def test_saturated_basis_also_recovers_the_original_map():
    comp = recip_r600a("Recip-M", basis="saturated")
    assert comp.density_reference is R and comp.suction_reference is None
    for to, tc in GRID:
        e = comp.efficiencies(to, tc, R.rho_dew(to), R.rho_dew(tc))
        assert (e["to_map"], e["tc_map"]) == (pytest.approx(to, abs=1e-9), pytest.approx(tc, abs=1e-9))


def test_shifted_fluid_reads_the_shifted_map():
    """Suction and condensing densities equal to R600a's 12 K lower: the map is read 12 K lower."""
    comp = recip_r600a("Recip-M")
    lo, hi = R.T_range
    for to, tc in GRID:
        if not (lo <= to - 12.0 and tc - 12.0 <= hi):
            continue
        e = comp.efficiencies(to, tc, rho_evap_dew=None, rho_cond_dew=R.rho_dew(tc - 12.0), rho_suction=S.rho_dew(to - 12.0))
        assert e["to_map"] == pytest.approx(to - 12.0, abs=1e-9)
        assert e["tc_map"] == pytest.approx(tc - 12.0, abs=1e-9)
        assert (e["eta_is"], e["eta_vol"]) == pytest.approx(expected_map(comp, to - 12.0, tc - 12.0), abs=1e-10)
        assert e["in_envelope"] == in_envelope(to - 12.0, tc - 12.0, R600A_ENVELOPE)


def test_equivalent_point_matches_densities():
    # suction gas 30 % thinner than R600a's rated suction state: to* < to; tc* from the condensing side only
    for to, tc in GRID:
        rho_s, rho_c = 0.7 * S.rho_dew(to), 0.7 * R.rho_dew(tc)
        to_m, tc_m, warnings = equivalent_point(R, None, rho_c, S, rho_s)
        if warnings:
            assert to_m == S.T_range[0]
            continue
        assert S.rho_dew(to_m) == pytest.approx(rho_s, rel=1e-10)
        assert R.rho_dew(tc_m) == pytest.approx(rho_c, rel=1e-10)
        assert to_m < to and tc_m < tc


def test_temperature_basis_keeps_the_old_reading():
    comp = recip_r600a("Recip-M", basis="temperature")
    assert comp.density_reference is None and comp.suction_reference is None
    e = comp.efficiencies(20.0, 90.0)                             # no densities needed
    assert (e["to_map"], e["tc_map"]) == (20.0, 90.0)
    assert e["eta_is"] == pytest.approx(ETA_IS(20.0, 90.0))
    with pytest.raises(ValueError):
        recip_r600a("Recip-M", basis="pressure")


def test_density_map_needs_densities():
    with pytest.raises(ValueError, match="saturated vapour densities"):
        recip_r600a("Recip-M").efficiencies(20.0, 90.0)
    with pytest.raises(ValueError, match="suction density"):
        recip_r600a("Recip-M").efficiencies(20.0, 90.0, R.rho_dew(20.0), R.rho_dew(90.0))


# ---------------------------------------------------------------- maps redrawn in a working fluid's (to, tc)


def shifted_curve(base, dT):
    lo, hi = base.T_range
    return fluid_curve(lambda T: base.rho_dew(T - dT) if lo <= T - dT <= hi else float("nan"), -38.0, 160.0, 1.0, "shifted")


def test_fluid_curve_stops_where_it_cannot_invert():
    c = fluid_curve(lambda T: 10.0 + T if T < 20 else 5.0, 0.0, 40.0, 1.0)      # stops rising at 20 C
    assert c.T_range == (0.0, 19.0)
    c = fluid_curve(lambda T: 10.0 + T if T < 30 else 1 / 0, 0.0, 40.0, 2.0)   # flash failure
    assert c.T_range == (0.0, 28.0)
    with pytest.raises(ValueError):
        fluid_curve(lambda T: 1 / 0, 0.0, 40.0)


def test_fluid_grid_on_r600a_at_15k_is_the_original_map():
    sat = fluid_curve(R.rho_dew, -50.0, 134.0, 2.0, "R600a")
    suc = fluid_curve(S.rho_dew, -50.0, 134.0, 2.0, "R600a 15 K")
    to_ax, tc_ax = [-10.0 + 2 * i for i in range(30)], [30.0 + 2 * i for i in range(46)]
    g = recip_r600a("Recip-M").fluid_grid(to_ax, tc_ax, sat, suc)
    assert g["to_map"] == pytest.approx(to_ax, abs=1e-9) and g["tc_map"] == pytest.approx(tc_ax, abs=1e-9)
    for j, tc in enumerate(tc_ax):
        for i, to in enumerate(to_ax):
            inside = in_envelope(to, tc, R600A_ENVELOPE)
            assert g["in_env"][j][i] == inside
            if inside:          # outside, the cells show the extension beyond the window (test_compressor_window)
                assert g["eta_is"][j][i] == pytest.approx(ETA_IS(to, tc), abs=1e-10)
                assert g["eta_vol"][j][i] == pytest.approx(ETA_VOL(to, tc), abs=1e-10)


def test_fluid_grid_on_a_shifted_fluid_is_the_shifted_map():
    sat, suc = shifted_curve(R, 12.0), shifted_curve(S, 12.0)
    to_ax, tc_ax = [-30.0 + 2 * i for i in range(50)], [20.0 + 2 * i for i in range(70)]
    g = recip_r600a("Recip-M").fluid_grid(to_ax, tc_ax, sat, suc)
    for j, tc in enumerate(tc_ax):
        for i, to in enumerate(to_ax):
            if not (R.T_range[0] <= to - 12.0 and tc - 12.0 <= R.T_range[1]) or tc > sat.T_range[1]:
                continue
            assert g["to_map"][i] == pytest.approx(to - 12.0, abs=1e-9)
            assert g["in_env"][j][i] == in_envelope(to - 12.0, tc - 12.0, R600A_ENVELOPE)
            if g["in_env"][j][i]:          # outside: the extension beyond the window
                assert g["eta_is"][j][i] == pytest.approx(ETA_IS(to - 12.0, tc - 12.0), abs=1e-10)
    assert g["tc_map"][-1] is None and g["eta_is"][-1][0] is None                # beyond the fluid table
    with pytest.raises(ValueError):
        recip_r600a("Recip-M").fluid_grid(to_ax, tc_ax, sat)                  # needs the suction curve


def test_fluid_grid_moves_only_the_evaporating_axis_with_superheat():
    sat = fluid_curve(R.rho_dew, -50.0, 134.0, 2.0)
    thin = fluid_curve(lambda T: 0.9 * S.rho_dew(T), -50.0, 134.0, 2.0)        # more superheat than rated
    to_ax, tc_ax = [0.0 + 2 * i for i in range(20)], [40.0 + 2 * i for i in range(30)]
    comp = recip_r600a("Recip-M")
    rated = comp.fluid_grid(to_ax, tc_ax, sat, fluid_curve(S.rho_dew, -50.0, 134.0, 2.0))
    more = comp.fluid_grid(to_ax, tc_ax, sat, thin)
    assert more["tc_map"] == rated["tc_map"]
    assert all(m < r for m, r in zip(more["to_map"], rated["to_map"]))


def test_envelope_moves_with_the_density():
    env = envelope_in_fluid(R600A_ENVELOPE, R, shifted_curve(R, 12.0), 1.0, S, shifted_curve(S, 12.0))
    same = envelope_in_fluid(R600A_ENVELOPE, R, fluid_curve(R.rho_dew, -50.0, 134.0, 1.0), 1.0,
                             S, fluid_curve(S.rho_dew, -50.0, 134.0, 1.0))
    assert len(env) == len(same) > 6 * 10
    for (a, b), (x, y) in zip(env, same):
        assert (a, b) == (pytest.approx(x + 12.0, abs=1e-8), pytest.approx(y + 12.0, abs=1e-8))
    for x, y in R600A_ENVELOPE:                                                  # corners are kept
        assert any(abs(a - x) < 1e-8 and abs(b - y) < 1e-8 for a, b in same)
    # a suction-only shift moves the envelope along the evaporating axis only
    evap_only = envelope_in_fluid(R600A_ENVELOPE, R, fluid_curve(R.rho_dew, -50.0, 134.0, 1.0), 1.0,
                                  S, shifted_curve(S, 5.0))
    assert len(evap_only) == len(same)
    for (a, b), (x, y) in zip(evap_only, same):
        assert (a, b) == (pytest.approx(x + 5.0, abs=1e-8), pytest.approx(y, abs=1e-8))


@pytest.mark.skipif(not (ROOT / "local_scripts" / "main_SGHP_dashboard.py").exists(),
                    reason="local_scripts/main_SGHP_dashboard.py not present (gitignored)")
def test_dashboard_map_reads_the_fluid_densities():
    import sys
    sys.path.insert(0, str(ROOT / "local_scripts"))
    import main_SGHP_dashboard as dash
    toy = ToyFluid()
    for sh in (15.0, 25.0):
        m = dash.fluid_maps(recip_r600a("Recip-M"), toy, "toy n-butane", dT_superheat=sh)
        e = m["eff"]
        assert m["mapSH"] == sh and len(e["to"]) > 10 and len(e["tc"]) > 10 and len(m["envelope"]) > 20
        for i in range(0, len(e["to"]), 7):
            for j in range(0, len(e["tc"]), 9):
                to, tc, a, b = e["to"][i], e["tc"][j], e["to_map"][i], e["tc_map"][j]
                if a is None or b is None:
                    continue
                dew = toy.ref_TQ(to, 1.0)
                assert S.rho_dew(a) == pytest.approx(toy.ref_PT(dew["P"], to + sh)["D"], rel=1e-4)   # rounded to 4 dp
                assert R.rho_dew(b) == pytest.approx(toy.ref_TQ(tc, 1.0)["D"], rel=1e-4)
                if e["in_env"][j][i]:
                    assert e["eta_is"][j][i] == pytest.approx(ETA_IS(a, b), abs=2e-4)


# ---------------------------------------------------------------- in the cycle (toy n-butane)


def test_cycle_reads_the_map_at_the_equivalent_point():
    toy = ToyFluid()
    comp = recip_r600a("Recip-M", displacement_m3h=60.0)
    res = model(toy, T_evap=20.0, T_cond=90.0, compressor=comp).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    to_m, tc_m = r["map_to"], r["map_tc"]
    # the defining equalities: actual suction density, condensing saturated-vapour density
    assert S.rho_dew(to_m) == pytest.approx(st["suction"]["D"], rel=1e-10)
    assert R.rho_dew(tc_m) == pytest.approx(st["cond_dew"]["D"], rel=1e-10)
    assert st["cond_dew"]["D"] == pytest.approx(toy.ref_PQ(res.P_H, 1.0)["D"])
    # every efficiency is the R600a fit at that point
    assert r["eta_is_overall"] == pytest.approx(ETA_IS(to_m, tc_m))
    assert r["eta_vol"] == pytest.approx(ETA_VOL(to_m, tc_m))
    assert r["m_dot_kgs"] == pytest.approx(ETA_VOL(to_m, tc_m) * 60.0 / 3600.0 * st["suction"]["D"])
    assert r["in_envelope"] == in_envelope(to_m, tc_m, R600A_ENVELOPE)
    # and it differs from reading the fits at the working fluid's own temperatures
    old = model(toy, T_evap=20.0, T_cond=90.0,
                compressor=recip_r600a("Recip-M", displacement_m3h=60.0, basis="temperature")).solve()
    assert (old.report["map_to"], old.report["map_tc"]) == (pytest.approx(20.0), pytest.approx(90.0))
    assert abs(old.COP - res.COP) > 1e-3


def test_superheat_moves_the_map_through_suction_density():
    toy = ToyFluid()
    pts = {}
    for basis in ("density", "saturated"):
        for sh in (5.0, 15.0, 25.0):
            res = model(toy, T_evap=20.0, T_cond=90.0, dT_superheat=sh,
                        compressor=recip_r600a("Recip-M", displacement_m3h=60.0, basis=basis)).solve()
            assert res.success, res.message
            pts[basis, sh] = (res.report["map_to"], res.report["map_tc"], res.states["suction"]["D"])
    d5, d15, d25 = (pts["density", sh] for sh in (5.0, 15.0, 25.0))
    assert d5[2] > d15[2] > d25[2]                                # hotter suction gas is thinner ...
    assert d5[0] > d15[0] > d25[0]                                # ... so the map is read at a lower to
    assert d5[1] == d15[1] == d25[1]                              # condensing side unaffected
    for sh in (5.0, 15.0, 25.0):                                  # the saturated basis ignores superheat
        assert pts["saturated", sh][:2] == pytest.approx(pts["saturated", 15.0][:2], abs=1e-12)


def test_envelope_is_checked_at_the_equivalent_point():
    toy = ToyFluid()
    comp = recip_r600a("Recip-M", displacement_m3h=None)
    for to, tc in ((40.0, 45.0), (20.0, 90.0), (0.0, 60.0)):
        res = model(toy, T_evap=to, T_cond=tc, dT_superheat=0.0, compressor=comp).solve()   # no IHX needed
        assert res.success, res.message
        inside = in_envelope(res.report["map_to"], res.report["map_tc"], R600A_ENVELOPE)
        assert res.report["in_envelope"] == inside
        assert any("R600a-equivalent" in w for w in res.warnings) == (not inside)


# ---------------------------------------------------------------- against the equation of state


def _check_table(curve, rho_eos, tol100=5e-6, tol125=5e-5):
    worst = {"<=100": 0.0, "<=125": 0.0}
    T = -45.0
    while T <= 125.0:
        k = "<=100" if T <= 100.0 else "<=125"
        worst[k] = max(worst[k], abs(curve.rho_dew(T) / rho_eos(T) - 1.0))
        T += 0.37
    assert worst["<=100"] < tol100 and worst["<=125"] < tol125, worst


def test_tables_match_coolprop():
    CP = pytest.importorskip("CoolProp.CoolProp")
    _check_table(R, lambda T: CP.PropsSI("D", "T", T + 273.15, "Q", 1, "IsoButane"))

    def suction(T):
        P = CP.PropsSI("P", "T", T + 273.15, "Q", 1, "IsoButane")
        return CP.PropsSI("D", "P", P, "T", T + 15.0 + 273.15, "IsoButane")
    _check_table(S, suction)


def test_tables_match_refprop(hpl):
    from HeatPumpLab.cycles.ConventionalHeatPump import RefpropProperties
    iso = RefpropProperties("Isobutane")
    _check_table(R, lambda T: iso.ref_TQ(T, 1.0)["D"])

    def suction(T):
        return iso.ref_PT(iso.ref_TQ(T, 1.0)["P"], T + 15.0)["D"]
    _check_table(S, suction)


def test_r600a_cycle_identity_refprop(hpl):
    """On R600a at 15 K suction superheat the density basis gives the original map and cycle."""
    base = dict(refrigerant="Isobutane", T_evap=20.0, T_cond=90.0, dT_subcool=10.0, dT_superheat=15.0,
                loss_UA_kWK=0.015)
    dens = ConventionalHeatPump(ConventionalHPConfig(**base, compressor=recip_r600a("Recip-M"))).solve()
    temp = ConventionalHeatPump(ConventionalHPConfig(**base, compressor=recip_r600a("Recip-M", basis="temperature"))).solve()
    assert dens.success and temp.success, (dens.message, temp.message)
    assert dens.report["map_to"] == pytest.approx(20.0, abs=2e-3)       # tables vs REFPROP: < 2e-5 in density
    assert dens.report["map_tc"] == pytest.approx(90.0, abs=2e-3)
    for k in ("COP", "Q_heat_kW", "W_electric_kW", "m_dot_kgs", "eta_vol", "eta_is_overall"):
        assert dens.report[k] == pytest.approx(temp.report[k], rel=1e-4), k
    hot = ConventionalHeatPump(ConventionalHPConfig(**(base | dict(dT_superheat=25.0)),
                                                    compressor=recip_r600a("Recip-M"))).solve()
    assert hot.success and hot.report["map_to"] < 19.0 and hot.report["map_tc"] == pytest.approx(90.0, abs=2e-3)


def test_butane_design_point_map_refprop(hpl):
    res = ConventionalHeatPump(ConventionalHPConfig()).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    # n-butane at 50 / 125 C, 15 K suction superheat: the R600a map is read at about 35.9 / 108.6 C
    assert r["map_to"] == pytest.approx(35.87, abs=0.1) and r["map_tc"] == pytest.approx(108.59, abs=0.1)
    assert S.rho_dew(r["map_to"]) == pytest.approx(st["suction"]["D"], rel=1e-9)
    assert R.rho_dew(r["map_tc"]) == pytest.approx(st["cond_dew"]["D"], rel=1e-9)
    assert r["in_envelope"]
    assert r["eta_vol"] == pytest.approx(ETA_VOL(r["map_to"], r["map_tc"]))
    assert r["eta_is_overall"] == pytest.approx(ETA_IS(r["map_to"], r["map_tc"]))


# ---------------------------------------------------------------- browser model parity


def _node(script):
    node = shutil.which("node")
    js = ROOT / "apps" / "src" / "steam" / "sghp_model.js"
    if node is None or not js.exists():
        pytest.skip("node or sghp_model.js not available")
    return json.loads(subprocess.run([node, "-e", f"const S = require({json.dumps(str(js))});" + script],
                                     capture_output=True, text=True, check=True).stdout)


def test_browser_model_matches_python():
    cases = [(R.rho_dew(c), S.rho_dew(t)) for t, c in GRID] + [(47.0, 3.1), (250.0, 0.3), (98.76, 12.345)]
    out = _node(f"const cases = {json.dumps(cases)};"
                "console.log(JSON.stringify({pts: cases.map(([c, s]) => S.equivalentPoint(S.R600A_SAT_VAPOUR, null, c, S.R600A_SUCTION_15K, s).slice(0, 2)),"
                " T: S.R600A_SAT_VAPOUR.T, rho: S.R600A_SAT_VAPOUR.rho, rhoS: S.R600A_SUCTION_15K.rho,"
                " mp: S.mapPoint(S.compressorMap(), 50, 125, 17.0, 90.0, 15.5).slice(0, 2),"
                " mpS: S.mapPoint(S.compressorMap('Recip-M', {map_basis: 'saturated'}), 50, 125, 17.0, 90.0, 15.5).slice(0, 2),"
                " mpT: S.mapPoint(S.compressorMap('Recip-M', {map_basis: 'temperature'}), 50, 125, 17.0, 90.0).slice(0, 2)}));")
    assert out["T"] == list(R600A_SAT_VAPOUR_T) and out["rho"] == list(R600A_SAT_VAPOUR_RHO)
    assert out["rhoS"] == list(R600A_SUCTION_15K_RHO)
    for (c, s), (to_js, tc_js) in zip(cases, out["pts"]):
        to_py, tc_py, _ = equivalent_point(R, None, c, S, s)
        assert to_js == pytest.approx(to_py, abs=1e-9) and tc_js == pytest.approx(tc_py, abs=1e-9)
    assert out["mp"] == pytest.approx(list(equivalent_point(R, 17.0, 90.0, S, 15.5)[:2]), abs=1e-9)
    assert out["mpS"] == pytest.approx(list(equivalent_point(R, 17.0, 90.0)[:2]), abs=1e-9)
    assert out["mpT"] == [50, 125]


def test_browser_fluid_grid_matches_python():
    to_ax, tc_ax = [-20.0 + 4 * i for i in range(20)], [30.0 + 4 * i for i in range(32)]
    out = _node("const R = S.R600A_SAT_VAPOUR, Q = S.R600A_SUCTION_15K;"
                "const c = S.fluidCurve(t => (t - 7.5 >= R.T_range[0] && t - 7.5 <= R.T_range[1]) ? 0.9 * R.rho_dew(t - 7.5) : NaN, -40, 160, 1, 'x');"
                "const s = S.fluidCurve(t => (t - 9 >= Q.T_range[0] && t - 9 <= Q.T_range[1]) ? 0.85 * Q.rho_dew(t - 9) : NaN, -40, 160, 1, 'x');"
                f"const g = S.fluidGrid(S.compressorMap(), {json.dumps(to_ax)}, {json.dumps(tc_ax)}, c, s);"
                "console.log(JSON.stringify({g, env: S.envelopeInFluid(S.ENVELOPE, R, c, 2, Q, s)}));")
    lo, hi = R.T_range
    sat = fluid_curve(lambda t: 0.9 * R.rho_dew(t - 7.5) if lo <= t - 7.5 <= hi else float("nan"), -40.0, 160.0, 1.0, "x")
    suc = fluid_curve(lambda t: 0.85 * S.rho_dew(t - 9) if lo <= t - 9 <= hi else float("nan"), -40.0, 160.0, 1.0, "x")
    py = recip_r600a("Recip-M").fluid_grid(to_ax, tc_ax, sat, suc)
    for k in ("to_map", "tc_map"):
        assert [v is None for v in out["g"][k]] == [v is None for v in py[k]]
        assert [v for v in out["g"][k] if v is not None] == pytest.approx([v for v in py[k] if v is not None], abs=1e-9)
    for k in ("eta_is", "eta_vol"):
        for rj, rp in zip(out["g"][k], py[k]):
            assert [v for v in rj if v is not None] == pytest.approx([v for v in rp if v is not None], abs=1e-10)
    assert out["g"]["in_env"] == py["in_env"]
    env_py = envelope_in_fluid(R600A_ENVELOPE, R, sat, 2.0, S, suc)
    assert len(out["env"]) == len(env_py)
    for (a, b, *_), (x, y) in zip(out["env"], env_py):          # the browser also tags each point's edge
        assert a == pytest.approx(x, abs=1e-9) and b == pytest.approx(y, abs=1e-9)
