"""The JavaScript ammonia-water properties (apps/src/dryer/ammonia_water.js) against the IAPWS G4-01 verification values
and against the Python module (HeatPumpLab/utils/ammonia_water.py, through teqp): the same model and algorithms, so the
same numbers to the solvers' tolerances."""
import json
import shutil
import subprocess
from pathlib import Path

import pytest

JS = Path(__file__).resolve().parents[1] / "apps" / "src" / "dryer" / "ammonia_water.js"

# IAPWS Guideline G4-01 (2001), section 8 (as in test_ammonia_water.py)
TABLE6 = [(0.1, 600, 35, 32.1221333), (0.1, 600, 4, 12.7721090), (0.5, 500, 32, 21.3208159),
          (0.5, 500, 1, 3.6423080), (0.9, 400, 30, 22.2830797), (0.9, 400, 0.5, 1.5499708)]
TABLE7 = [(0.2, 300, 0.040710, 0.9360, 51.941, 0.01640), (0.4, 400, 2.5545, 0.9363, 43.318, 0.8608),
          (0.6, 500, 16.698, 0.7844, 25.459, 8.86)]
TABLE8 = [(0.2, 300, 0.00437062, 0.010672, 55.16434, 0.00175506), (0.4, 400, 0.394694, 0.051541, 50.83187, 0.122658),
          (0.6, 500, 6.52607, 0.22135, 39.93714, 2.00730)]


def _node(body: str, timeout=300):
    """Run ``body`` in node with the module as ``A`` (model as ``M``); its last output line is JSON."""
    node = shutil.which("node")
    if node is None:
        pytest.skip("node not available")
    code = f"const A = require({json.dumps(str(JS))}); const M = A.model;\n{body}"
    res = subprocess.run([node, "-e", code], capture_output=True, text=True, timeout=timeout)
    assert res.returncode == 0, res.stderr[-1500:]
    return json.loads(res.stdout.strip().splitlines()[-1])


@pytest.fixture(scope="module")
def py():
    pytest.importorskip("teqp")
    from HeatPumpLab.utils.ammonia_water import AmmoniaWater
    return AmmoniaWater()


def test_iapws_verification_values():
    o = _node(f"""
const t6 = {json.dumps(TABLE6)}, t7 = {json.dumps(TABLE7)}, t8 = {json.dumps(TABLE8)};
console.log(JSON.stringify({{
  p6: t6.map(([x, T, rho]) => M.pressure(T, rho * 1000, x) / 1e6),
  b7: t7.map(([x, T]) => {{ const [p, y] = M.bubble_p(T, x); return [p / 1e6, y, M.density(T, p, x, 'liquid') / 1000,
                                                                  M.density(T, p, y, 'vapour') / 1000]; }}),
  d8: t8.map(([y, T]) => {{ const [p, x] = M.dew_p(T, y); return [p / 1e6, x, M.density(T, p, x, 'liquid') / 1000,
                                                               M.density(T, p, y, 'vapour') / 1000]; }}),
}}));""")
    for (x, T, rho, p), pj in zip(TABLE6, o["p6"]):
        assert pj == pytest.approx(p, abs=1e-7)
    for (x, T, p, y, rhoL, rhoV), (pj, yj, rLj, rVj) in zip(TABLE7, o["b7"]):
        assert pj == pytest.approx(p, rel=5e-5) and yj == pytest.approx(y, abs=5e-5)
        assert rLj == pytest.approx(rhoL, rel=5e-5) and rVj == pytest.approx(rhoV, rel=1e-3)
    for (y, T, p, x, rhoL, rhoV), (pj, xj, rLj, rVj) in zip(TABLE8, o["d8"]):
        assert pj == pytest.approx(p, rel=5e-6) and xj == pytest.approx(x, abs=5e-6)
        assert rLj == pytest.approx(rhoL, rel=5e-6) and rVj == pytest.approx(rhoV, rel=5e-6)


def test_helmholtz_matches_teqp(py):
    """Pressure, enthalpy, entropy and ln(phi) at (T, rho, x), liquid, vapour and supercritical, including water's
    critical region (the IAPWS-95 Gaussian and non-analytic terms) and both composition limits."""
    pts = [(600, 35000, 0.1), (500, 1000, 0.5), (400, 500, 0.9), (650, 17000, 0.02), (640, 15000, 0.05),
           (380, 30, 0.99), (550, 300, 1e-10), (330, 55000, 1e-6), (300, 51941, 0.2), (500, 25459, 0.6), (420, 44051.58302989, 0.3069605662707), (420, 795.90274936872, 0.8613935031725),
           (600, 8000, 0.4), (360, 20, 1 - 1e-10)]
    o = _node(f"""console.log(JSON.stringify({json.dumps(pts)}.map(([T, r, x]) =>
                  [M.pressure(T, r, x), ...M.h_s(T, r, x), ...M.ln_phi(T, r, x)])));""")
    for (T, r, x), js in zip(pts, o):
        ref = [py.pressure(T, r, x), *py.h_s(T, r, x), *py.ln_phi(T, r, x)]
        # p = rho R T (1 + delta ar_delta) cancels in a compressed liquid: compare it on the scale of rho R T
        assert js[0] == pytest.approx(ref[0], rel=1e-10, abs=1e-12 * r * 8.314471 * T), (T, r, x)
        for a, b in zip(js[1:3], ref[1:3]):
            assert a == pytest.approx(b, rel=1e-10, abs=1e-10), (T, r, x)
        # ln(phi) holds -ln Z: rounding in Z = 1 + delta ar_delta is amplified by 1/Z in a compressed liquid
        Z = ref[0] / (r * 8.314471 * T)
        for a, b in zip(js[3:], ref[3:]):
            assert a == pytest.approx(b, rel=1e-10, abs=1e-13 / Z), (T, r, x)


def test_saturation_matches_python(py):
    """Bubble and dew temperatures and pressures: the same fixed points, reached from different starting values (the
    JS starts from the pure-fluid ancillaries, Python from CoolProp tables), so agreement to the solver tolerances."""
    TP = [(0.3e6, 0.3), (1e6, 0.05), (2.5e6, 0.5), (5e6, 0.9), (1.5e6, 0.7)]
    PT = [(300, 0.2), (360, 0.5), (400, 0.4), (450, 0.8), (500, 0.6)]
    o = _node(f"""const TP = {json.dumps(TP)}, PT = {json.dumps(PT)};
console.log(JSON.stringify({{ bT: TP.map(([p, x]) => M.bubble_T(p, x)), dT: TP.map(([p, x]) => M.dew_T(p, x)),
  bp: PT.map(([T, x]) => M.bubble_p(T, x)), dp: PT.map(([T, x]) => M.dew_p(T, x)),
  eq: [[420, 2.5e6], [380, 1e6], [330, 0.3e6]].map(([T, p]) => M.equilibrium_TP(T, p)) }}));""")
    for (p, x), (Tb, yb), (Td, xd) in zip(TP, o["bT"], o["dT"]):
        rb, rd = py.bubble_T(p, x), py.dew_T(p, x)
        assert Tb == pytest.approx(rb[0], rel=1e-9) and yb == pytest.approx(rb[1], abs=1e-9)
        assert Td == pytest.approx(rd[0], rel=1e-9) and xd == pytest.approx(rd[1], abs=1e-9)
    for (T, x), (pb, yb), (pd, xd) in zip(PT, o["bp"], o["dp"]):
        rb, rd = py.bubble_p(T, x), py.dew_p(T, x)
        # the pressure substitution stops at |dp/p| < 1e-10 per step and converges linearly: ~1e-8 apart
        assert pb == pytest.approx(rb[0], rel=1e-7) and yb == pytest.approx(rb[1], abs=1e-7)
        assert pd == pytest.approx(rd[0], rel=1e-7) and xd == pytest.approx(rd[1], abs=1e-7)
    for (T, p), js in zip([(420, 2.5e6), (380, 1e6), (330, 0.3e6)], o["eq"]):
        assert js == pytest.approx(list(py.equilibrium_TP(T, p)), rel=1e-9)


CASES = [(0.3, 0.3), (0.3, 5.0), (0.5, 1.5), (0.5, 5.0), (0.9, 0.3), (0.9, 2.5)]   # (w, P MPa)


def test_stream_states_match_python(py):
    """PT (liquid, two-phase, vapour), PH and PS at those states' h and s, PQ and TQ, for each (w, P) case.
    Temperatures within 1e-6 K and enthalpies within 1e-6 kJ/kg: the flashes stop at xtol = 1e-10 K, and the
    two-phase equilibria agree to ~1e-12 in composition."""
    from HeatPumpLab.utils.ammonia_water import AmmoniaWaterStream
    o = _node(f"""const out = {json.dumps(CASES)}.map(([w, P]) => {{
  const s = A.stream(w), [Tb, Td] = s.bubble_dew(P), Ts = [Tb - 293.15, (Tb + Td) / 2 - 273.15, Td - 253.15];
  const pt = Ts.map(T => s.PT(P, T));
  return {{ bd: s.bubble_dew(P), Ts, pt, ph: pt.map(st => s.PH(P, st.h + 0.5)), ps: pt.map(st => s.PS(P, st.s - 1e-3)),
    pq: [0, 0.4, 1].map(q => s.PQ(P, q)), tq: [0, 0.5, 1].map(q => s.TQ(Ts[1], q)) }};
}});
console.log(JSON.stringify(out));""")
    keys = ("T", "P", "h", "s", "D", "q", "w_liq", "w_vap")
    tol = dict(T=1e-6, P=1e-9, h=1e-6, s=1e-8, D=None, q=1e-8, w_liq=1e-9, w_vap=1e-9)

    def same(js, ref, what):
        assert js["phase"] == ref["phase"], what
        for k in keys:
            a, b = js[k], ref[k]
            if b is None:
                assert a is None, (what, k)
            elif k == "D":
                assert a == pytest.approx(b, rel=1e-8), (what, k)
            else:
                assert a == pytest.approx(b, abs=tol[k]), (what, k)

    for (w, P), r in zip(CASES, o):
        s = AmmoniaWaterStream(w, py)
        assert r["bd"] == pytest.approx(list(s.bubble_dew(P)), rel=1e-9)
        for i, T in enumerate(r["Ts"]):
            ref = s.PT(P, T)
            same(r["pt"][i], ref, ("PT", w, P, T))
            same(r["ph"][i], s.PH(P, ref["h"] + 0.5), ("PH", w, P, T))
            same(r["ps"][i], s.PS(P, ref["s"] - 1e-3), ("PS", w, P, T))
        assert [st["phase"] for st in r["pt"]] == ["liquid", "two_phase", "superheated_vapour"]
        for q, js in zip((0, 0.4, 1), r["pq"]):
            same(js, s.PQ(P, q), ("PQ", w, P, q))
        for q, js in zip((0, 0.5, 1), r["tq"]):
            ref = s.TQ(r["Ts"][1], q)
            assert js["P"] == pytest.approx(ref["P"], rel=1e-7), ("TQ", w, P, q)
            assert js["T"] == pytest.approx(ref["T"], abs=1e-6) and js["h"] == pytest.approx(ref["h"], abs=1e-5)
            assert js["q"] == pytest.approx(ref["q"], abs=1e-8) and js["phase"] == ref["phase"]


def test_fast_isobar_table():
    """T_ph_fast follows the PH flash within 0.02 K (as in the Python tests), and matches Python's table."""
    o = _node("""const s = A.stream(0.5), lo = s.PT(2.5, 40).h, hi = s.PT(2.5, 240).h, out = [];
for (let k = 1; k < 20; k++) { const h = lo + (hi - lo) * k / 20; out.push([h, s.T_ph_fast(2.5, h), s.PH(2.5, h).T]); }
console.log(JSON.stringify(out));""")
    for h, fast, flash in o:
        assert fast == pytest.approx(flash, abs=0.02)
    pytest.importorskip("teqp")
    from HeatPumpLab.utils.ammonia_water import AmmoniaWaterStream
    s = AmmoniaWaterStream(0.5)
    for h, fast, _ in o:
        assert fast == pytest.approx(s.T_ph_fast(2.5, h), abs=1e-6)


def test_interface():
    o = _node("""const s = A.stream(0.5); let err = null;
try { s.PH(2.5, 1e5); } catch (e) { err = e instanceof Error ? e.message : 'not an Error'; }
const st = s.PT(1.0, 30.0);
console.log(JSON.stringify({ name: s.name, fluid: s.fluid, backend: s.backend, abk: A.backend, crit: [s.P_crit, s.T_crit].map(Number.isNaN),
  same: s.fast_view() === s, memo: s.PT(1.0, 30.0) === st, gas: [s.is_gas(st), s.is_gas(s.PQ(1.0, 1))],
  tph: s.Tph(1.0, st.h), err, mf: [A.mole_fraction(0.5), A.mass_fraction(A.mole_fraction(0.3))], M: [A.M_NH3, A.M_H2O] }));""")
    assert o["name"] == o["fluid"] == "ammonia-water (50.0 % NH3 by mass)"
    assert o["backend"] == o["abk"] == "IAPWS G4-01 (JavaScript)"
    assert o["crit"] == [True, True] and o["same"] and o["memo"] and o["gas"] == [False, True]
    assert o["tph"] == pytest.approx(30.0, abs=1e-6)
    assert o["err"] and "above the vapour range" in o["err"]
    assert o["mf"][0] == pytest.approx(0.5140532624, rel=1e-9) and o["mf"][1] == pytest.approx(0.3, rel=1e-14)
    assert o["M"] == [17.03026e-3, 18.015268e-3]
