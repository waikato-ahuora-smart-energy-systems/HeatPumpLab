"""Source, steam and hot-water specification of the SGHP cycle (ConventionalHeatPump).

The heat source (isothermal) sets the evaporator outlet dTmin_source below it; the steam sets the
condenser outlet (bubble point) dTmin_steam above it; the subcooler heats water from T_water_in to
T_water_out and cools the liquid as far as the minimum approach dTmin_water allows. These tests
check the pressures, the subcooler pinch (cold end and, for an unusual liquid, inside the
exchanger), the water flow, the failure messages, the equivalence with the direct T_evap / T_cond /
dT_subcool specification, and the browser model against Python.
"""
import json
import shutil
import subprocess
import tempfile
from pathlib import Path

import CoolProp.CoolProp as CP
import pytest

from HeatPumpLab.cycles.ConventionalHeatPump import ConventionalHeatPump, ConventionalHPConfig, EfficiencyCurve
from test_conventional_heat_pump import ToyFluid

ROOT = Path(__file__).resolve().parents[1]
WATER = dict(T_source=15.0, T_steam=75.0, T_water_in=30.0, T_water_out=60.0)   # evap 10 C, bubble 80 C


def water_P(T_out):
    """Hot water pressure, Pa: 1 bar above saturation at the outlet temperature."""
    return CP.PropsSI("P", "T", T_out + 273.15, "Q", 0, "Water") + 1e5


def water_dh(T_in, T_out):
    """Hot water enthalpy rise, kJ/kg, independent of the model."""
    P = water_P(T_out)
    return (CP.PropsSI("H", "T", T_out + 273.15, "P", P, "Water") - CP.PropsSI("H", "T", T_in + 273.15, "P", P, "Water")) / 1000


def steam_model(props=None, **kw):
    cfg = ConventionalHPConfig(**({"refrigerant": "toy", "dT_superheat": 5.0, "loss_UA_kWK": None,
                                   "compressor": EfficiencyCurve(eta_is=0.7, motor_eff=0.9), "Q_heat_kW": 50.0}
                                  | WATER | kw))
    return ConventionalHeatPump(cfg, props=props or ToyFluid())


def test_source_and_steam_set_the_pressures():
    toy = ToyFluid()
    res = steam_model(toy).solve()
    assert res.success, res.message
    r = res.report
    assert res.P_L == pytest.approx(toy.Psat(10.0)) and res.P_H == pytest.approx(toy.Psat(80.0))
    assert (r["T_source"], r["T_steam"]) == (15.0, 75.0)
    assert r["T_evap"] == pytest.approx(10.0) and r["T_evap_dew"] == pytest.approx(10.0)
    assert r["T_cond"] == pytest.approx(80.0) and r["T_cond_bubble"] == pytest.approx(80.0)
    assert r["COP_carnot"] == pytest.approx((80.0 + 273.15) / 70.0)


def test_hot_water_pinch_is_at_the_cold_end():
    res = steam_model().solve()
    r, st = res.report, res.states
    assert st["subcool_out"]["T"] == pytest.approx(35.0)                 # water in + 5 K
    assert r["dT_subcool"] == pytest.approx(45.0)
    assert r["subcooler_min_approach_K"] == pytest.approx(5.0)
    assert (r["T_water_in"], r["T_water_out"]) == (30.0, 60.0)
    assert r["m_water_kgs"] == pytest.approx(r["Q_subcool_kW"] / water_dh(30.0, 60.0))
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-9)


def test_hot_water_properties_come_from_coolprop():
    """The flow uses real enthalpies, at 1 bar above saturation at the outlet temperature, whatever the refrigerant."""
    r = steam_model(T_steam=75.0, T_water_in=30.0, T_water_out=60.0).solve().report
    assert r["P_water_bar"] == pytest.approx(water_P(60.0) / 1e5, rel=1e-6)
    assert r["h_water_out_kJkg"] - r["h_water_in_kJkg"] == pytest.approx(water_dh(30.0, 60.0), rel=1e-6)
    assert r["m_water_kgs"] == pytest.approx(r["Q_subcool_kW"] / water_dh(30.0, 60.0), rel=1e-6)
    # hot water that a fixed cp of 4.19 kJ/kg K would get about 4 % wrong
    hot = ConventionalHeatPump(ConventionalHPConfig(refrigerant="n-Pentane", property_backend="coolprop", T_source=150.0, T_steam=180.0,
                                                    T_water_in=160.0, T_water_out=170.0)).solve()
    assert hot.success, hot.message
    assert hot.report["P_water_bar"] == pytest.approx(water_P(170.0) / 1e5, rel=1e-6)
    cp_mean = water_dh(160.0, 170.0) / 10.0
    assert 4.3 < cp_mean < 4.45 and hot.report["m_water_kgs"] == pytest.approx(hot.report["Q_subcool_kW"] / (cp_mean * 10.0), rel=1e-6)
    assert steam_model(T_evap=10.0, T_cond=80.0, dT_subcool=45.0).solve().report["P_water_bar"] is None


def test_same_cycle_as_the_direct_specification():
    via_water = steam_model().solve()
    direct = steam_model(T_evap=10.0, T_cond=80.0, dT_subcool=45.0).solve()
    assert via_water.success and direct.success
    for k in ("COP", "q_heat_kJkg", "q_cond_kJkg", "q_subcool_kJkg", "w_electric_kJkg", "m_dot_kgs"):
        assert via_water.report[k] == pytest.approx(direct.report[k], rel=1e-10), k
    assert direct.report["T_source"] is None and direct.report["T_steam"] is None
    assert direct.report["m_water_kgs"] is None                          # no water specified


def test_ihx_after_the_subcooler_cools_the_liquid_further():
    res = steam_model(dT_superheat=15.0, T_water_in=40.0, T_water_out=70.0).solve()
    assert res.success, res.message
    st, r = res.states, res.report
    assert st["subcool_out"]["T"] == pytest.approx(45.0)                  # water pinch at the cold end
    assert st["ihx_out"]["T"] < st["subcool_out"]["T"]                    # then the IHX
    assert r["q_subcool_kJkg"] == pytest.approx(st["cond_out"]["h"] - st["subcool_out"]["h"])


def test_condensing_above_the_critical_point_is_reported():
    toy = ToyFluid()                                                       # T_crit 151.98 C
    res = steam_model(toy, T_steam=150.0).solve()                          # condensing at 155 C
    assert not res.success and "at or above the critical temperature" in res.message and "(152.0 C)" in res.message


def test_water_failures_are_reported():
    hot = steam_model(T_water_out=78.0).solve()
    assert not hot.success and "hot water outlet 78 C is less than 5 K below" in hot.message
    back = steam_model(T_water_in=60.0, T_water_out=50.0).solve()
    assert not back.success and "outlet must be above its inlet" in back.message


class SaggingLiquid(ToyFluid):
    """Toy fluid whose liquid heat capacity falls with temperature, so the subcooler pinch moves
    inside the exchanger instead of the cold end."""

    k = 0.012

    def _liq(self, T_C, P):
        st = super()._liq(T_C, P)
        st["h"] -= self.k * (T_C + 273.15 - self.T0) ** 2
        return st


def test_interior_pinch_raises_the_subcooler_outlet():
    fluid = SaggingLiquid()
    res = steam_model(fluid, T_water_out=70.0).solve()
    assert res.success, res.message
    st, r = res.states, res.report
    T_out = st["subcool_out"]["T"]
    assert T_out > 35.0 + 0.05                                             # not the cold-end value
    assert r["subcooler_min_approach_K"] == pytest.approx(5.0, abs=1e-4)
    # independent dense check of the approach along the exchanger (water enthalpy linear in duty)
    h_in, h_out, P = st["cond_out"]["h"], st["subcool_out"]["h"], res.P_H
    Pw = water_P(70.0)
    hw0, hw1 = (CP.PropsSI("H", "T", T + 273.15, "P", Pw, "Water") for T in (30.0, 70.0))
    worst = min(fluid.ref_PH(P, h_out + (h_in - h_out) * i / 400)["T"]
                - (CP.PropsSI("T", "H", hw0 + (hw1 - hw0) * i / 400, "P", Pw, "Water") - 273.15) for i in range(401))
    assert worst == pytest.approx(5.0, abs=2e-3) and worst > 5.0 - 2e-3


def test_rating_map_over_source_and_steam():
    hp = steam_model()
    recs = hp.rating_map([10.0, 15.0], [70.0, 75.0])                    # default keys: T_source, T_steam
    assert len(recs) == 4 and all(r["success"] for r in recs)
    by = {(r["T_source"], r["T_steam"]): r for r in recs}
    assert by[15.0, 70.0]["T_cond_bubble"] == pytest.approx(75.0) and by[10.0, 70.0]["T_evap_dew"] == pytest.approx(5.0)
    assert by[15.0, 70.0]["COP"] > by[10.0, 75.0]["COP"]
    with pytest.raises(ValueError):
        hp.rating_map([10.0], [70.0], low_key="T_cold")


# ---------------------------------------------------------------- browser model parity

JS_TOY = r"""
const R = 8.314462618 / 58.12, cpv = 1.25, cpl = 2.40, L0 = 385.0, T0 = 273.15, P0 = 0.1035, rho_l = 580.0;
const Psat = Tc => { const T = Tc + 273.15, L = L0 + (cpv - cpl) * (T - T0); return P0 * Math.exp((T * (L0 / T0 + (cpv - cpl) * Math.log(T / T0)) - L) / (R * T)); };
const Tsat = P => S.brent(T => Math.log(Psat(T) / P), -150, 250, 1e-12);
const liq = (Tc, P) => { const T = Tc + 273.15; return { T: Tc, P, h: cpl * (T - T0), s: cpl * Math.log(T / T0), D: rho_l, phase: "liquid", q: -1 }; };
const vap = (Tc, P) => { const T = Tc + 273.15; return { T: Tc, P, h: L0 + cpv * (T - T0), s: L0 / T0 + cpv * Math.log(T / T0) - R * Math.log(P / P0), D: P * 1000 / (R * T), phase: "superheated_vapour", q: -1 }; };
const mix = (P, q) => { const Ts = Tsat(P), l = liq(Ts, P), v = vap(Ts, P);
  return { T: Ts, P, q, h: l.h + q * (v.h - l.h), s: l.s + q * (v.s - l.s), D: 1 / (1 / l.D + q * (1 / v.D - 1 / l.D)),
           phase: q === 0 ? "saturated_liquid" : q === 1 ? "saturated_vapour" : "two_phase" }; };
const from = (P, k, val) => { const l = mix(P, 0), v = mix(P, 1);
  if (val < l[k]) return liq(S.brent(T => liq(T, P)[k] - val, -200, l.T + 1e-9, 1e-12), P);
  if (val > v[k]) return vap(S.brent(T => vap(T, P)[k] - val, v.T - 1e-9, 600, 1e-12), P);
  return mix(P, (val - l[k]) / (v[k] - l[k])); };
const toy = { P_crit: 3.796, T_crit: 151.98,
  PT: (P, T) => (T > Tsat(P) ? vap(T, P) : liq(T, P)), PH: (P, h) => from(P, "h", h), PS: (P, s) => from(P, "s", s),
  PQ: (P, q) => mix(P, q), TQ: (T, q) => mix(Psat(T), q) };
"""


def node_toy(script):
    """Run a script with S (the model) and toy (the toy fluid, with the page's CoolProp water as toy.water)."""
    node, live, cp = shutil.which("node"), ROOT / "apps" / "src" / "steam", ROOT / "apps" / "src" / "coolprop"
    if node is None or not (cp / "coolprop.js").exists():
        pytest.skip("node or apps/src/coolprop/coolprop.js not available")
    with tempfile.TemporaryDirectory() as d:
        tmp = Path(d)
        shutil.copy(cp / "coolprop.js", tmp / "coolprop.mjs")
        shutil.copy(cp / "coolprop.wasm", tmp / "coolprop.wasm")
        (tmp / "t.mjs").write_text(f"""
import F from './coolprop.mjs'; import fs from 'fs'; import {{ createRequire }} from 'module';
const S = createRequire(import.meta.url)({json.dumps(str(live / "sghp_model.js"))});
const src = fs.readFileSync({json.dumps(str(live / "sghp_live.js"))}, 'utf8');
const cut = (a, b) => src.slice(src.indexOf(a), src.indexOf(b));
const code = cut('function cpCall', '/* ---------------------------------------------------------------- property adapter')
  + cut('const LABEL', '/* ---------------------------------------------------------------- state */');
const CP = await F({{ locateFile: p => new URL(p, import.meta.url).pathname }});
const Fluid = new Function('CP', 'SGHP', code + '\\nreturn CoolPropFluid;')(CP, S);
""" + JS_TOY + "toy.water = new Fluid(CP, [['Water', 100]]);\n" + script)
        out = subprocess.run([node, str(tmp / "t.mjs")], capture_output=True, text=True, cwd=tmp)
    assert out.returncode == 0, out.stderr[-800:]
    return json.loads(out.stdout)


def test_browser_model_matches_python():
    cases = [dict(),
             dict(dT_superheat=15.0, T_water_in=40.0, T_water_out=70.0), dict(T_water_out=78.0), dict(T_steam=150.0),
             dict(loss_UA_cond_kWK=0.02, loss_UA_subcool_kWK=0.01, loss_UA_evap_kWK=0.03), dict(preheat_subcool_share=0.5),
             dict(preheat_subcool_share=1.0), dict(preheat_subcool_share=1.0, loss_UA_cond_kWK=0.02),
             dict(dP_K_cond=0.05, dP_K_subcool=0.5, dP_K_evap=0.02), dict(dP_K_cond=0.05, preheat_subcool_share=1.0),
             dict(second_steam=True, T_steam2=68.0), dict(second_steam=True, T_steam2=62.0, preheat_subcool_share=1.0),
             dict(second_steam=True, T_steam2=68.0, loss_UA_cond_kWK=0.02, loss_UA_subcool_kWK=0.01, loss_UA_steam2_kWK=0.015),
             dict(second_steam=True, T_steam2=68.0, dP_K_cond=0.05, dP_K_subcool=0.5, dP_K_evap=0.02),
             dict(second_steam=True, T_steam2=58.0), dict(second_steam=True, T_steam2=78.0)]
    # the JS model always reads the R600a map; give Python the same map, size-less in neither
    script = ("const cases = %s; const base = %s;" % (json.dumps(cases), json.dumps(WATER | {"dT_superheat": 5.0})) +
              "const comp = S.compressorMap('Recip-M', {map_basis: 'temperature'});"
              "console.log(JSON.stringify(cases.map(c => { const r = S.solve(toy, Object.assign({}, base, c), comp, 0.0);"
              " return r.success ? {ok: true, r: r.report} : {ok: false, msg: r.message}; })));")
    out = node_toy(script)
    from HeatPumpLab.cycles.ConventionalHeatPump import recip_r600a
    for c, js_res in zip(cases, out):
        py = steam_model(**c, Q_heat_kW=None, loss_UA_kWK=0.0,
                         compressor=recip_r600a("Recip-M", basis="temperature")).solve()
        assert py.success == js_res["ok"], (c, py.message, js_res)
        if not py.success:
            assert py.message.split(" (")[0] == js_res["msg"].split(" (")[0]
            continue
        for k in ("COP", "Q_heat_kW", "Q_subcool_kW", "m_water_kgs", "T_subcool_out", "dT_subcool",
                  "subcooler_min_approach_K", "T_evap_dew", "T_cond_bubble", "rho_suction_kgm3", "m_steam_kgs",
                  "steam_bleed_fraction", "P_steam_bar", "V_suction_m3h", "V_discharge_m3h", "Q_cond_kW", "Q_source_kW",
                  "Q_loss_cond_kW", "Q_loss_subcool_kW", "Q_loss_evap_kW", "VHC", "T_cond_out", "dT_subcool_cond", "steam_preheat_fraction", "steam_gen_min_approach_K",
                  "preheat_subcool_share", "PR", "P_discharge_bar", "P_evap_in_bar", "dP_cond_kPa", "dP_subcool_kPa",
                  "dP_evap_kPa", "T_cond_dew", "W_electric_kW", "density_ratio", "density_ratio_beyond", "to_read",
                  "tc_read", "w_swept_kJm3", "second_steam", "T_steam2", "T_steam2_out", "Q_steam2_kW", "m_steam2_kgs",
                  "steam2_min_approach_K", "q_steam2_kJkg", "m_hot_water_net_kgs", "energy_balance_residual_kJkg",
                  "P_water_bar", "h_water_in_kJkg", "h_water_out_kJkg"):
            if py.report[k] is None:                                  # second steam off: both models report None
                assert js_res["r"][k] is None, (c, k)
                continue
            assert js_res["r"][k] == pytest.approx(py.report[k], rel=1e-8, abs=1e-6 if k == "energy_balance_residual_kJkg" else 1e-9), (c, k)
        assert js_res["r"]["limit_flags"] == py.report["limit_flags"] and js_res["r"]["window_edge"] == py.report["window_edge"]


# ---------------------------------------------------------------- steam raised from bled hot water

def test_water_saturation_table():
    from HeatPumpLab.cycles.ConventionalHeatPump import water_saturation
    s100, s150 = water_saturation(100.0), water_saturation(150.0)          # IAPWS-95 reference values
    assert s100["P"] == pytest.approx(0.101418, rel=1e-5) and s100["h_f"] == pytest.approx(419.17, abs=0.01)
    assert s100["h_g"] == pytest.approx(2675.57, abs=0.02) and s100["v_f"] == pytest.approx(0.0010435, rel=1e-4)
    assert s150["P"] == pytest.approx(0.47616, rel=1e-4) and s150["h_f"] == pytest.approx(632.18, abs=0.02)
    assert s150["h_g"] == pytest.approx(2745.9, abs=0.1)
    mid = water_saturation(100.5)                                            # linear between grid points
    assert water_saturation(100.0)["h_g"] < mid["h_g"] < water_saturation(101.0)["h_g"]
    with pytest.raises(ValueError):
        water_saturation(250.0)


def test_steam_flow_from_bled_hot_water():
    from HeatPumpLab.cycles.ConventionalHeatPump import steam_from_feedwater, water_saturation
    sf = steam_from_feedwater(32.0, 120.0, 95.0)
    s120, s95 = water_saturation(120.0), water_saturation(95.0)
    h_feed = s95["h_f"] + s95["v_f"] * (s120["P"] - s95["P"]) * 1000
    assert sf["dh_kJkg"] == pytest.approx(s120["h_g"] - h_feed) and sf["m_kgs"] == pytest.approx(32.0 / sf["dh_kJkg"])
    assert sf["m_kgs"] * 3600 == pytest.approx(49.92, abs=0.01)             # about 50 kg/h of steam from 32 kW
    assert h_feed - s95["h_f"] == pytest.approx(0.115, abs=0.01)            # pump work to the steam pressure, small
    assert sf["preheat_kJkg"] == pytest.approx(s120["h_f"] - h_feed)         # feedwater to boiling, before evaporation
    with pytest.raises(ValueError, match="not below the steam temperature"):
        steam_from_feedwater(10.0, 95.0, 95.0)
    # in the cycle report: steam from Q_cond, bled from the hot water
    res = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=EfficiencyCurve(eta_is=0.7, motor_eff=0.9, displacement_m3h=56.25)).solve()
    r = res.report
    assert res.success and r["steam_feed_T"] == 60.0
    assert r["m_steam_kgs"] == pytest.approx(steam_from_feedwater(r["Q_cond_kW"], 75.0, 60.0)["m_kgs"])
    assert r["steam_bleed_fraction"] == pytest.approx(r["m_steam_kgs"] / r["m_water_kgs"])
    assert r["m_hot_water_net_kgs"] == pytest.approx(r["m_water_kgs"] - r["m_steam_kgs"])
    st = res.states                                                          # volume flows at the compressor
    assert r["V_suction_m3h"] == pytest.approx(r["m_dot_kgs"] / st["suction"]["D"] * 3600)
    assert r["V_discharge_m3h"] == pytest.approx(r["m_dot_kgs"] / st["discharge"]["D"] * 3600) and r["V_discharge_m3h"] < r["V_suction_m3h"]
    direct = steam_model(T_evap=10.0, T_cond=80.0, dT_subcool=45.0).solve()
    assert direct.report["m_steam_kgs"] is None                              # no steam temperature specified


# ---------------------------------------------------------------- exchanger heat losses

def test_exchanger_heat_losses():
    comp = EfficiencyCurve(eta_is=0.7, motor_eff=0.9, displacement_m3h=56.25)
    base = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp).solve()
    UA = dict(loss_UA_cond_kWK=0.02, loss_UA_subcool_kWK=0.01, loss_UA_evap_kWK=0.03)
    res = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp, **UA).solve()
    assert base.success and res.success
    r, b, st = res.report, base.report, res.states
    # refrigerant side unchanged; the losses use the mean refrigerant temperature in each exchanger (ambient 15 C)
    assert r["m_dot_kgs"] == pytest.approx(b["m_dot_kgs"]) and r["Q_cond_refrigerant_kW"] == pytest.approx(b["Q_cond_kW"])
    assert r["Q_loss_cond_kW"] == pytest.approx(0.02 * ((st["discharge"]["T"] + st["cond_out"]["T"]) / 2 - 15.0))
    assert r["Q_loss_subcool_kW"] == pytest.approx(0.01 * ((st["cond_out"]["T"] + st["subcool_out"]["T"]) / 2 - 15.0))
    assert r["Q_loss_evap_kW"] == pytest.approx(0.03 * ((st["evap_in"]["T"] + st["evap_out"]["T"]) / 2 - 15.0))
    # delivered heat is net of the condenser and subcooler losses; the source covers the evaporator loss
    assert r["Q_cond_kW"] == pytest.approx(b["Q_cond_kW"] - r["Q_loss_cond_kW"])
    assert r["Q_subcool_kW"] == pytest.approx(b["Q_subcool_kW"] - r["Q_loss_subcool_kW"])
    assert r["Q_heat_kW"] == pytest.approx(r["Q_cond_kW"] + r["Q_subcool_kW"])
    assert r["Q_source_kW"] == pytest.approx(r["Q_evap_kW"] + r["Q_loss_evap_kW"])
    assert r["COP"] == pytest.approx(r["Q_heat_kW"] / r["W_electric_kW"]) and res.COP == pytest.approx(r["COP"])
    assert res.COP < base.COP
    # overall balance: electrical + source = delivered + all losses to ambient
    losses = r["Q_compressor_loss_kW"] + r["Q_loss_cond_kW"] + r["Q_loss_subcool_kW"] + r["Q_loss_evap_kW"]
    assert r["W_electric_kW"] + r["Q_source_kW"] == pytest.approx(r["Q_heat_kW"] + losses, rel=1e-9)
    # water and steam follow the delivered duties
    assert r["m_water_kgs"] == pytest.approx(r["Q_subcool_kW"] / water_dh(30.0, 60.0))
    assert r["m_steam_kgs"] < b["m_steam_kgs"]
    # sizing to a heat output covers the losses
    sized = steam_model(Q_heat_kW=50.0, loss_UA_kWK=None, **UA).solve()
    assert sized.success and sized.report["Q_heat_kW"] == pytest.approx(50.0)


# ---------------------------------------------------------------- subcooling at the condenser exit

def test_condenser_exit_subcooling_makes_more_steam():
    """preheat_subcool_share: the share of the feedwater preheat (bled hot water heated to boiling) taken from
    subcooling the refrigerant at the condenser exit; that heat goes to the steam instead of the hot water."""
    comp = EfficiencyCurve(eta_is=0.7, motor_eff=0.9, displacement_m3h=56.25)
    base = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp).solve()
    assert base.success and base.report["dT_subcool_cond"] == pytest.approx(0.0, abs=1e-9)
    assert base.report["steam_gen_min_approach_K"] == pytest.approx(5.0)
    for share in (0.5, 1.0):
        sub = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp, preheat_subcool_share=share).solve()
        assert sub.success, sub.message
        b, r = base.report, sub.report
        assert sub.P_H == pytest.approx(base.P_H) and r["T_cond_bubble"] == pytest.approx(80.0)    # same pressure
        assert r["preheat_subcool_share"] == share and r["dT_subcool_cond"] > 0
        # the subcooling heat is that share of the steam's preheat duty
        Q_sub = r["Q_cond_refrigerant_kW"] - b["Q_cond_refrigerant_kW"]
        preheat = r["m_steam_kgs"] * r["steam_preheat_fraction"] * r["steam_dh_kJkg"]
        assert Q_sub == pytest.approx(share * preheat, rel=1e-6)
        assert r["Q_cond_kW"] > b["Q_cond_kW"] and r["m_steam_kgs"] > b["m_steam_kgs"]             # more steam
        assert r["Q_subcool_kW"] < b["Q_subcool_kW"]                                               # less hot water
        assert r["Q_heat_kW"] == pytest.approx(b["Q_heat_kW"])     # same refrigerant cycle: the heat moves to the steam
        assert r["steam_gen_min_approach_K"] >= 5.0 - 1e-6
    # with a condenser loss the share still holds (on the steam actually made)
    lossy = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp, preheat_subcool_share=1.0, loss_UA_cond_kWK=0.02).solve()
    lb = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp, loss_UA_cond_kWK=0.02).solve()
    r = lossy.report
    Q_sub = r["Q_cond_refrigerant_kW"] - lb.report["Q_cond_refrigerant_kW"]
    assert Q_sub == pytest.approx(r["m_steam_kgs"] * r["steam_preheat_fraction"] * r["steam_dh_kJkg"], rel=1e-4)


# ---------------------------------------------------------------- second steam level from the subcooling

SIZED = dict(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=EfficiencyCurve(eta_is=0.7, motor_eff=0.9, displacement_m3h=56.25))


def test_second_steam_off_changes_nothing():
    base = steam_model(**SIZED).solve()
    off = steam_model(**SIZED, second_steam=False, T_steam2=68.0, loss_UA_steam2_kWK=0.5).solve()
    assert base.report == off.report and base.states == off.states and "steam2_out" not in off.states
    assert off.report["second_steam"] is False and off.report["m_steam2_kgs"] is None


def test_second_steam_comes_out_of_the_hot_water():
    """Same refrigerant cycle: the second steam takes all the subcooling its approach allows, the hot water the rest."""
    from HeatPumpLab.cycles.ConventionalHeatPump import steam_from_feedwater
    base = steam_model(**SIZED).solve()
    res = steam_model(**SIZED, second_steam=True, T_steam2=68.0).solve()
    assert res.success, res.message
    b, r, st = base.report, res.report, res.states
    for k in ("COP", "Q_heat_kW", "Q_cond_kW", "m_steam_kgs", "T_subcool_out", "dT_subcool", "m_dot_kgs"):
        assert r[k] == pytest.approx(b[k], rel=1e-10), k
    assert r["Q_steam2_kW"] + r["Q_subcool_kW"] == pytest.approx(b["Q_subcool_kW"]) and r["Q_steam2_kW"] > 0
    assert r["Q_heat_kW"] == pytest.approx(r["Q_cond_kW"] + r["Q_steam2_kW"] + r["Q_subcool_kW"])
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-9)
    # the pinch is where boiling starts: the refrigerant there is dTmin_steam above the second steam
    sf = steam_from_feedwater(r["Q_steam2_kW"], 68.0, 60.0)
    fp, q = sf["preheat_kJkg"] / sf["dh_kJkg"], st["cond_out"]["h"] - st["steam2_out"]["h"]
    assert ToyFluid().ref_PH(res.P_H, st["steam2_out"]["h"] + fp * q)["T"] == pytest.approx(73.0, abs=1e-6)
    assert r["steam2_min_approach_K"] == pytest.approx(5.0, abs=1e-6) and r["T_steam2_out"] < 73.0
    assert r["subcooler_min_approach_K"] == pytest.approx(5.0)
    # water follows the smaller duty, and both steam levels bleed their feedwater from it
    assert r["m_water_kgs"] == pytest.approx(r["Q_subcool_kW"] / water_dh(30.0, 60.0))
    assert r["m_steam2_kgs"] == pytest.approx(sf["m_kgs"]) and r["P_steam2_bar"] == pytest.approx(sf["P_steam_MPa"] * 10)
    assert r["m_hot_water_net_kgs"] == pytest.approx(r["m_water_kgs"] - r["m_steam_kgs"] - r["m_steam2_kgs"])
    assert r["steam_bleed_fraction"] == pytest.approx((r["m_steam_kgs"] + r["m_steam2_kgs"]) / r["m_water_kgs"])


def test_second_steam_after_preheat_share_and_with_losses():
    # the condenser-exit subcooling comes first, so the second generator starts colder and raises less
    plain = steam_model(**SIZED, second_steam=True, T_steam2=68.0).solve()
    pre = steam_model(**SIZED, second_steam=True, T_steam2=68.0, preheat_subcool_share=1.0).solve()
    assert pre.success, pre.message
    assert pre.report["m_steam_kgs"] > plain.report["m_steam_kgs"] and pre.report["Q_steam2_kW"] < plain.report["Q_steam2_kW"]
    assert pre.report["steam_gen_min_approach_K"] >= 5.0 - 1e-6
    # the loss comes off the second steam duty, at the mean refrigerant temperature in that generator
    UA = dict(loss_UA_steam2_kWK=0.002, loss_UA_subcool_kWK=0.01)
    res = steam_model(**SIZED, second_steam=True, T_steam2=68.0, **UA).solve()
    r, st = res.report, res.states
    assert r["Q_loss_steam2_kW"] == pytest.approx(0.002 * ((st["cond_out"]["T"] + st["steam2_out"]["T"]) / 2 - 15.0))
    assert r["Q_loss_subcool_kW"] == pytest.approx(0.01 * ((st["steam2_out"]["T"] + st["subcool_out"]["T"]) / 2 - 15.0))
    assert r["Q_steam2_kW"] == pytest.approx(plain.report["Q_steam2_kW"] - r["Q_loss_steam2_kW"])
    assert r["Q_heat_kW"] == pytest.approx(r["Q_cond_kW"] + r["Q_steam2_kW"] + r["Q_subcool_kW"])
    # sizing to a heat output covers that loss too
    sized = steam_model(second_steam=True, T_steam2=68.0, **UA).solve()
    assert sized.success and sized.report["Q_heat_kW"] == pytest.approx(50.0)


def test_second_steam_limits():
    below = steam_model(second_steam=True, T_steam2=60.0).solve()              # not above its feedwater
    assert not below.success and "not above the hot water outlet" in below.message
    tight = steam_model(second_steam=True, T_steam2=75.0).solve()              # no room below the condenser exit
    assert not tight.success and "below the refrigerant leaving the condenser" in tight.message
    # needs the steam and hot-water specification: ignored, with a warning, on a direct specification
    direct = steam_model(second_steam=True, T_steam2=68.0, T_evap=10.0, T_cond=80.0, dT_subcool=45.0).solve()
    assert direct.success and direct.report["second_steam"] is False
    assert any("second steam level" in w for w in direct.warnings)


def test_subcooled_evaporator_inlet_is_flagged():
    """A heat source hotter than the liquid: with no IHX to catch it the valve outlet stays liquid. It solves, but the
    source then also heats the liquid and the COP is not a heat-pump COP, so it is flagged."""
    res = steam_model(T_source=45.0, dT_superheat=0.0, compressor=EfficiencyCurve(eta_is=0.4, motor_eff=0.9)).solve()
    assert res.success and res.report["evap_in_subcooled"] is True
    assert any("subcooled" in w and "inflates the COP" in w for w in res.warnings)
    assert steam_model().solve().report["evap_in_subcooled"] is False


def test_evaporating_above_the_critical_temperature_says_so():
    res = ConventionalHeatPump(ConventionalHPConfig(refrigerant="Propane", property_backend="coolprop", T_source=150.0, T_steam=180.0)).solve()
    assert not res.success and "at or above the critical temperature" in res.message


# ---------------------------------------------------------------- refrigerant pressure drops

def test_pressure_drops_with_pinned_exits():
    """dP [kPa] = K V^1.8 (V in m3/h at the mean specific volume of each exchanger's inlet and outlet). The exits
    stay pinned by the approaches, so the compressor lifts from P_L to P_H + dP_cond and the evaporator inlet sits
    dP_evap above P_L."""
    comp = EfficiencyCurve(eta_is=0.7, motor_eff=0.9, displacement_m3h=56.25)
    base = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp).solve()
    K = dict(dP_K_cond=0.05, dP_K_subcool=0.5, dP_K_evap=0.02)
    res = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp, **K).solve()
    assert base.success and res.success, res.message
    b, r, st = base.report, res.report, res.states
    assert b["dP_cond_kPa"] == 0 and b["P_discharge_bar"] == pytest.approx(b["P_cond_bar"])
    # converged: each drop matches its coefficient at the solved flow
    assert r["dP_cond_kPa"] == pytest.approx(0.05 * r["Vmean_cond_m3h"] ** 1.8, rel=1e-6)
    assert r["dP_subcool_kPa"] == pytest.approx(0.5 * r["Vmean_subcool_m3h"] ** 1.8, rel=1e-6)
    assert r["dP_evap_kPa"] == pytest.approx(0.02 * r["Vmean_evap_m3h"] ** 1.8, rel=1e-6)
    assert all(r[k] > 0 for k in ("dP_cond_kPa", "dP_subcool_kPa", "dP_evap_kPa"))
    v = lambda s: 1 / s["D"]
    m = r["m_dot_kgs"]
    assert r["Vmean_cond_m3h"] == pytest.approx(m * (v(st["discharge"]) + v(st["cond_out"])) / 2 * 3600)
    # exits pinned: same condensing and evaporating pressures; the compressor sees the drops
    assert res.P_H == pytest.approx(base.P_H) and res.P_L == pytest.approx(base.P_L)
    assert r["T_cond_bubble"] == pytest.approx(80.0) and r["T_evap_dew"] == pytest.approx(10.0)
    assert r["P_discharge_bar"] == pytest.approx(r["P_cond_bar"] + r["dP_cond_kPa"] / 100)
    assert r["P_evap_in_bar"] == pytest.approx(r["P_evap_bar"] + r["dP_evap_kPa"] / 100)
    assert st["discharge"]["P"] == pytest.approx(r["P_discharge_bar"] / 10)
    assert r["PR"] == pytest.approx(r["P_discharge_bar"] / r["P_evap_bar"]) and r["PR"] > b["PR"]
    assert res.COP < base.COP and r["W_electric_kW"] > b["W_electric_kW"]


def test_preheat_share_with_loss_when_sized_by_heat_output():
    """Sized by Q_heat_kW the mass flow is found after the cycle; the share still holds with a condenser loss."""
    UA = dict(loss_UA_cond_kWK=0.02)
    for share in (0.5, 1.0):
        lossy = steam_model(Q_heat_kW=30.0, preheat_subcool_share=share, **UA).solve()
        base = steam_model(Q_heat_kW=30.0, **UA).solve()
        assert lossy.success and base.success and lossy.report["sized_by"] == "Q_heat_kW"
        r = lossy.report
        q_sub = r["Q_cond_refrigerant_kW"] / r["m_dot_kgs"] - base.report["Q_cond_refrigerant_kW"] / base.report["m_dot_kgs"]
        preheat = r["m_steam_kgs"] * r["steam_preheat_fraction"] * r["steam_dh_kJkg"]
        assert r["m_dot_kgs"] * q_sub == pytest.approx(share * preheat, rel=1e-6)


def test_preheat_share_needs_the_hot_water_subcooler():
    res = steam_model(dT_subcool=5.0, preheat_subcool_share=0.5).solve()
    assert res.success and res.report["preheat_subcool_share"] == 0.0
    assert res.report["dT_subcool_cond"] == pytest.approx(0.0, abs=1e-9)
    assert any("preheat from subcooling needs" in w for w in res.warnings)


def test_subcooler_loss_larger_than_its_duty_fails():
    comp = EfficiencyCurve(eta_is=0.7, motor_eff=0.9, displacement_m3h=56.25)
    res = steam_model(Q_heat_kW=None, loss_UA_kWK=0.0, compressor=comp, loss_UA_subcool_kWK=0.5).solve()
    assert not res.success and "no hot water" in res.message
