"""Compressor efficiency beyond the application window (compressor_window.py, and sghp_model.js).

Ratio edges (efficiency falls away) and mechanical edges (seals, materials, motor: flagged on the cycle,
efficiency unchanged); outside the window the maps are read at the nearest window point, then fall off with the
density ratio past the ratio lines.
"""
import json
import math
import shutil
import subprocess
from pathlib import Path

import pytest

from HeatPumpLab.cycles.ConventionalHeatPump import R600A_ENVELOPE, recip_r600a, in_envelope
from HeatPumpLab.cycles.ConventionalHeatPump.compressor_window import (
    R600A_EDGE_KINDS, R600A_LIMITS, limit_flags, map_density_ratio, nearest_on_envelope)
from HeatPumpLab.cycles.ConventionalHeatPump.compressor_data import R600A_EXTENSION as X

ROOT = Path(__file__).resolve().parents[1]
COMP = recip_r600a("Recip-M", basis="temperature")


def eff(to, tc, r=None):
    return X.efficiencies(COMP.eta_is, COMP.eta_vol, to, tc, map_density_ratio(to, tc) if r is None else r)


def test_edges_are_ratio_or_mechanical():
    # the high-ratio edge sits near a density ratio of 19, the low-ratio edge near a pressure ratio of 1.5-1.7
    hi = [map_density_ratio(*p) for p in ((6, 120), (-2, 108.5), (-10, 97))]
    assert all(18.5 < r < 19.6 for r in hi)
    lo = [map_density_ratio(*p) for p in ((27, 42), (48, 70))]
    assert all(1.55 < r < 1.85 for r in lo)
    assert R600A_EDGE_KINDS.count("mech") == 4


def test_inside_is_the_map():
    for to, tc in ((20, 80), (-10, 97), (48, 120), (27, 42), (0, 60)):
        e = eff(to, tc)
        assert e["inside"] and e["eta_is"] == COMP.eta_is(to, tc) and e["eta_vol"] == COMP.eta_vol(to, tc)


def test_continuous_across_every_edge():
    n = len(R600A_ENVELOPE)
    for i, a in enumerate(R600A_ENVELOPE):
        b = R600A_ENVELOPE[(i + 1) % n]
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        nx, ny = b[1] - a[1], -(b[0] - a[0])                       # outward normal (clockwise polygon)
        L = math.hypot(nx, ny)
        if in_envelope(mx + nx / L, my + ny / L, R600A_ENVELOPE):
            nx, ny = -nx, -ny
        edge = eff(mx, my)
        out = eff(mx + 1e-4 * nx / L, my + 1e-4 * ny / L)
        assert not out["inside"] and out["edge"] == i
        assert out["eta_is"] == pytest.approx(edge["eta_is"], abs=1e-4)
        assert out["eta_vol"] == pytest.approx(edge["eta_vol"], abs=1e-4)


def test_mechanical_edges_hold_efficiency():
    # the bicubic fit rises above tc = 120 C (0.714 -> 0.788 at to = 48); the extension holds the edge value
    assert COMP.eta_is(48, 130) > COMP.eta_is(48, 120) + 0.05
    e = eff(48, 130)
    assert e["edge_kind"] == "mech" and e["ratio_beyond"] == 0
    assert e["eta_is"] == pytest.approx(COMP.eta_is(48, 120)) and e["eta_vol"] == pytest.approx(COMP.eta_vol(48, 120))
    e = eff(60, 100)                                              # beyond the maximum evaporating edge
    assert e["read"] == pytest.approx((48, 100)) and e["eta_is"] == pytest.approx(COMP.eta_is(48, 100))


def test_high_ratio_falloff():
    q = nearest_on_envelope(-20, 100, R600A_ENVELOPE)
    assert R600A_EDGE_KINDS[q["edge"]] == "ratio_high"
    r_hi = X.ratio_limits(-20)[1]
    ev_q, eis_q = COMP.eta_vol(q["to"], q["tc"]), COMP.eta_is(q["to"], q["tc"])
    for k in (1.1, 1.5, 1.8):
        e = eff(-20, 100, k * r_hi)
        ev = 1 - (1 - ev_q) * (k * r_hi - 1) / (r_hi - 1)            # clearance re-expansion
        assert e["eta_vol"] == pytest.approx(ev) and e["eta_is"] == pytest.approx(eis_q * (ev / ev_q) ** 0.4)
        assert e["ratio_beyond"] == pytest.approx(k - 1)
    assert 0.68 < eff(-20, 100, 1.5 * r_hi)["eta_is"] / eis_q < 0.82     # about 25 % lower at 1.5 r_hi
    far = eff(-20, 100, 4 * r_hi)                                    # past zero flow: floors
    assert far["eta_vol"] == X.eta_vol_floor and far["eta_is"] < 0.2
    # a corner: beyond the maximum condensing edge and the high-ratio line together
    c = eff(6, 125)
    assert c["edge_kind"] == "mech" and c["ratio_beyond"] > 0 and c["eta_is"] < COMP.eta_is(6, 120)


def test_low_ratio_falloff():
    to, tc = 55, 60
    q = nearest_on_envelope(to, tc, R600A_ENVELOPE)
    assert R600A_EDGE_KINDS[q["edge"]] == "ratio_low"
    r_lo = X.ratio_limits(to)[0]
    eis_q = COMP.eta_is(q["to"], q["tc"])
    for r in (1.5, 1.2, 1.05):
        e = eff(to, tc, r)
        x, xl = r - 1, r_lo - 1
        assert e["eta_is"] == pytest.approx(max(eis_q * (x / (x + 0.3)) / (xl / (xl + 0.3)), X.eta_is_floor))
        assert e["eta_vol"] <= 1.0
    assert eff(to, tc, 1.0)["eta_is"] == X.eta_is_floor


def test_limit_flags():
    lim = R600A_LIMITS["Recip-M"]
    assert limit_flags(lim, 24.0, 120.0, 900.0, 5.0) == []
    f = limit_flags(lim, 65.0, 177.0, 1200.0, 0.9)
    assert [x["key"] for x in f] == ["P_dis", "T_dis", "w_swept", "P_suc"]
    assert f[0]["limit"] == 28.4 and f[-1]["kind"] == "min"
    w = limit_flags(lim, 24.0, 120.0, 1051.0, 5.0, 56.25 / (1450 * 60))[0]
    assert w["key"] == "w_swept" and w["over_pct"] == pytest.approx(3.55, abs=0.01)
    assert w["detail"] == ("Shaft torque ≈ 108 N·m (1051 kJ/m³ swept) is 3.5 % beyond the operating window (1015 kJ/m³): "
                           "check with the supplier that the motor, crank and bearings are suitable.")


def test_browser_extension_matches_python():
    node = shutil.which("node")
    js = ROOT / "apps" / "src" / "steam" / "sghp_model.js"
    if node is None or not js.exists():
        pytest.skip("node or sghp_model.js not available")
    pts = [[to, tc] for to in range(-30, 71, 7) for tc in range(20, 161, 9)]
    script = ("const S = require(%s); const pts = %s; const c = S.compressorMap('Recip-M', {map_basis: 'temperature'});"
              "const out = []; for (const [to, tc] of pts) for (const r of [S.mapDensityRatio(to, tc), 1.3]) {"
              " const e = S.extendedEfficiencies(c, to, tc, r); out.push([e.eta_is, e.eta_vol, e.inside, e.ratio_beyond]); }"
              "console.log(JSON.stringify({out, flags: S.limitFlags(S.LIMITS['Recip-M'], {P_dis_bar: 65, T_dis: 177, w_swept_kJm3: 1200, P_suc_bar: 0.9, rev_m3: 56.25 / (1450 * 60)})}));"
              % (json.dumps(str(js)), json.dumps(pts)))
    o = json.loads(subprocess.run([node, "-e", script], capture_output=True, text=True, check=True).stdout)
    k = 0
    for to, tc in pts:
        for r in (map_density_ratio(to, tc), 1.3):
            e = eff(to, tc, r)
            js_e = o["out"][k]; k += 1
            assert js_e[2] == e["inside"]
            assert js_e[0] == pytest.approx(e["eta_is"], abs=1e-12) and js_e[1] == pytest.approx(e["eta_vol"], abs=1e-12)
            assert js_e[3] == pytest.approx(e["ratio_beyond"], abs=1e-12)
    py = limit_flags(R600A_LIMITS["Recip-M"], 65, 177, 1200, 0.9, 56.25 / (1450 * 60))
    assert [f["detail"] for f in o["flags"]] == [f["detail"] for f in py]
    for a, b in zip(o["flags"], py):
        assert {k: v for k, v in a.items() if k != "over_pct"} == {k: v for k, v in b.items() if k != "over_pct"}
        assert a["over_pct"] == pytest.approx(b["over_pct"], rel=1e-12)


def test_shaft_torque_uses_shaft_work():
    # w_swept is electrical work per swept volume; the torque takes the shaft share of it
    w = limit_flags(R600A_LIMITS["Recip-M"], 24.0, 120.0, 1051.0, 5.0, 56.25 / (1450 * 60), shaft_fraction=0.93)[0]
    assert w["detail"].startswith("Shaft torque ≈ 101 N·m (1051 kJ/m³ swept)")


@pytest.mark.parametrize("basis", ["saturated", "temperature"])
def test_extension_continuous_at_the_edge_whatever_the_superheat(basis):
    """Beyond the high-ratio edge (to -2, tc 108.5) the falloff starts at the window on every basis, even with a
    suction density far from R600a's at 15 K superheat."""
    from HeatPumpLab.cycles.ConventionalHeatPump.compressor_curves import CompressorInlet
    from HeatPumpLab.cycles.ConventionalHeatPump.density_basis import R600A_SAT_VAPOUR as R, R600A_SUCTION_15K as S
    comp = recip_r600a("Recip-M", basis=basis)

    def point(tc):
        inlet = CompressorInlet(suction=dict(T=28.0, P=0.15, h=600.0, s=2.5, D=0.8 * S.rho_dew(-2.0)), P_discharge=2.0,
                                h_isentropic=700.0, T_evap_dew=-2.0, T_cond_dew=tc,
                                rho_evap_dew=R.rho_dew(-2.0), rho_cond_dew=R.rho_dew(tc))
        return comp.evaluate(inlet)

    a, b = point(108.5 - 1e-4), point(108.5 + 1e-4)
    assert a.window["inside"] and not b.window["inside"]
    assert b.window["eta_is"] == pytest.approx(a.window["eta_is"], abs=1e-4)
    assert b.eta_vol == pytest.approx(a.eta_vol, abs=1e-4)


def test_efficiencies_and_grid_use_the_extension():
    comp = recip_r600a("Recip-M", basis="temperature")
    # the bicubic fit rises above tc = 120 C; the extension holds the edge value
    assert comp.efficiencies(48.0, 130.0)["eta_is"] == pytest.approx(comp.eta_is(48.0, 120.0), abs=1e-9)
    g = comp.fluid_grid([-2.0, 48.0], [112.0, 130.0])
    for i, to in enumerate(g["to"]):
        for j, tc in enumerate(g["tc"]):
            e = comp.efficiencies(to, tc)
            assert (g["eta_is"][j][i], g["eta_vol"][j][i]) == pytest.approx((e["eta_is"], e["eta_vol"]), abs=1e-12)


def test_envelope_override_drops_the_r600a_extension():
    from HeatPumpLab.cycles.ConventionalHeatPump.compressor_curves import EfficiencyMap
    assert recip_r600a("Recip-M", envelope=None).extension is None
    box = ((-20.0, 20.0), (60.0, 20.0), (60.0, 140.0), (-20.0, 140.0))
    assert recip_r600a("Recip-M", envelope=box).extension is None
    assert recip_r600a("Recip-M").extension is X
    with pytest.raises(ValueError, match="envelope"):
        EfficiencyMap(eta_is=COMP.eta_is, envelope=box, extension=X)
    assert EfficiencyMap(eta_is=COMP.eta_is).eta_is_bounds == (0.2, 1.0)


def test_extension_needs_a_kind_per_edge():
    from HeatPumpLab.cycles.ConventionalHeatPump.compressor_window import EnvelopeExtension
    box = ((-20.0, 20.0), (60.0, 20.0), (60.0, 140.0), (-20.0, 140.0))
    with pytest.raises(ValueError, match="edge_kinds"):
        EnvelopeExtension(box)
    e = EnvelopeExtension(box, ("mech",) * 4)
    assert e.edge_names == ("edge 0", "edge 1", "edge 2", "edge 3")
    assert e.efficiencies(COMP.eta_is, COMP.eta_vol, 70.0, 80.0, 5.0)["edge_name"] == "edge 1"
