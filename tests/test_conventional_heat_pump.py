"""Conventional single-stage heat pump (evaporator, compressor, condenser, subcooler, valve).

Most tests use a small self-consistent toy pure fluid (ideal-gas vapour, incompressible liquid,
saturation from equal Gibbs energy, roughly n-butane-like) so the cycle logic runs without REFPROP.
The REFPROP-backed tests at the end use the shared ``hpl`` fixture and skip without REFPROP.
"""
import math

import pytest
from scipy.optimize import brentq

from HeatPumpLab.cycles.ConventionalHeatPump import (
    ConventionalHeatPump,
    ConventionalHPConfig,
    EfficiencyCurve,
    EN12900Curve,
    MapCurve,
    R600A_ENVELOPE,
    bicubic_polynomial,
    recip_r600a,
    en12900_polynomial,
    fit_en12900,
    in_envelope,
)
from HeatPumpLab.cycles.ConventionalHeatPump.compressor_data import R600A_COEFFICIENTS


# ---------------------------------------------------------------- toy fluid


class ToyFluid:
    """Pure fluid; T in degC, P in MPa, h kJ/kg, s kJ/kg/K, D kg/m3 (HeatPumpLab dict keys)."""

    R = 8.314462618 / 58.12
    cpv, cpl, L0, T0, P0 = 1.25, 2.40, 385.0, 273.15, 0.1035
    rho_l = 580.0
    P_crit, T_crit = 3.796, 151.98

    def Psat(self, T_C):
        T = T_C + 273.15
        L = self.L0 + (self.cpv - self.cpl) * (T - self.T0)
        lnp = (T * (self.L0 / self.T0 + (self.cpv - self.cpl) * math.log(T / self.T0)) - L) / (self.R * T)
        return self.P0 * math.exp(lnp)

    def Tsat(self, P):
        return brentq(lambda T: math.log(self.Psat(T) / P), -150.0, 250.0, xtol=1e-12)

    def _liq(self, T_C, P):
        T = T_C + 273.15
        return dict(T=T_C, P=P, h=self.cpl * (T - self.T0), s=self.cpl * math.log(T / self.T0), D=self.rho_l)

    def _vap(self, T_C, P):
        T = T_C + 273.15
        return dict(T=T_C, P=P, h=self.L0 + self.cpv * (T - self.T0),
                    s=self.L0 / self.T0 + self.cpv * math.log(T / self.T0) - self.R * math.log(P / self.P0),
                    D=P * 1000.0 / (self.R * T))

    def _mix(self, P, q):
        Ts = self.Tsat(P)
        l, v = self._liq(Ts, P), self._vap(Ts, P)
        out = {k: l[k] + q * (v[k] - l[k]) for k in ("h", "s")}
        out.update(T=Ts, P=P, q=q, D=1.0 / (1 / l["D"] + q * (1 / v["D"] - 1 / l["D"])),
                   phase="saturated_liquid" if q == 0 else "saturated_vapour" if q == 1 else "two_phase")
        return out

    def ref_PQ(self, P, q):
        return self._mix(P, q)

    def ref_TQ(self, T, q):
        return self._mix(self.Psat(T), q)

    def ref_PT(self, P, T):
        Ts = self.Tsat(P)
        if T > Ts:
            return self._vap(T, P) | dict(phase="superheated_vapour", q=-9999990.0)
        return self._liq(T, P) | dict(phase="liquid", q=-9999990.0)

    def _from(self, P, key, value):
        sat_l, sat_v = self._mix(P, 0.0), self._mix(P, 1.0)
        if value < sat_l[key]:
            T = brentq(lambda T: self._liq(T, P)[key] - value, -200.0, sat_l["T"] + 1e-9)
            return self._liq(T, P) | dict(phase="liquid", q=-9999990.0)
        if value > sat_v[key]:
            T = brentq(lambda T: self._vap(T, P)[key] - value, sat_v["T"] - 1e-9, 600.0, xtol=1e-12)
            return self._vap(T, P) | dict(phase="superheated_vapour", q=-9999990.0)
        return self._mix(P, (value - sat_l[key]) / (sat_v[key] - sat_l[key]))

    def ref_PH(self, P, h):
        return self._from(P, "h", h)

    def ref_PS(self, P, s):
        return self._from(P, "s", s)


@pytest.fixture
def toy():
    return ToyFluid()


def model(toy, **kw):
    cfg = ConventionalHPConfig(**({"refrigerant": "toy", "T_evap": 10.0, "T_cond": 80.0,
                                   "dT_subcool": 5.0, "dT_superheat": 5.0, "dT_discharge_superheat": None,
                                   "loss_UA_kWK": None,
                                   "compressor": EfficiencyCurve(eta_is=0.7, motor_eff=0.9)} | kw))
    return ConventionalHeatPump(cfg, props=toy)


# ---------------------------------------------------------------- toy fluid sanity


def test_toy_fluid_is_consistent(toy):
    assert toy.Psat(-0.5) == pytest.approx(0.1013, rel=0.03)      # n-butane normal boiling point
    for T in (0.0, 40.0, 90.0):
        P = toy.Psat(T)
        assert toy.Tsat(P) == pytest.approx(T, abs=1e-8)
        st = toy.ref_PH(P, toy.ref_PQ(P, 0.3)["h"])
        assert st["q"] == pytest.approx(0.3)


# ---------------------------------------------------------------- cycle on the toy fluid


def test_cycle_states_and_energy_balance(toy):
    res = model(toy).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    assert res.P_L == pytest.approx(toy.Psat(10.0))
    assert res.P_H == pytest.approx(toy.Psat(80.0))
    assert st["suction"]["T"] == pytest.approx(15.0)
    assert st["cond_out"]["T"] == pytest.approx(80.0) and st["cond_out"]["phase"] == "saturated_liquid"
    assert st["subcool_out"]["T"] == pytest.approx(75.0) and st["subcool_out"]["phase"] == "liquid"
    assert st["subcool_out"]["h"] - st["ihx_out"]["h"] == pytest.approx(r["q_ihx_kJkg"])   # IHX after the subcooler
    assert st["evap_in"]["h"] == pytest.approx(st["ihx_out"]["h"])            # isenthalpic valve
    assert st["evap_in"]["phase"] == "two_phase" and st["evap_in"]["T"] == pytest.approx(10.0)
    dh_is = st["isentropic_discharge"]["h"] - st["suction"]["h"]
    assert st["discharge"]["h"] - st["suction"]["h"] == pytest.approx(dh_is / 0.7)
    assert st["discharge"]["T"] > 80.0
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-9)
    assert r["q_heat_kJkg"] == pytest.approx(r["q_cond_kJkg"] + r["q_subcool_kJkg"])
    assert r["q_subcool_kJkg"] == pytest.approx(toy.cpl * 5.0)
    assert res.COP == pytest.approx(r["q_heat_kJkg"] / (dh_is / 0.7 / 0.9))
    assert res.COP - res.COP_cooling == pytest.approx(r["w_to_gas_kJkg"] / r["w_electric_kJkg"])
    assert 1.0 < res.COP < r["COP_carnot"]
    assert res.VHC == pytest.approx(r["q_heat_kJkg"] * st["suction"]["D"] / 1000.0)
    assert "m_dot_kgs" not in r                                            # no size given


def test_subcooling_raises_heating_cop_and_can_be_excluded(toy):
    cops = [model(toy, dT_subcool=dT).solve().COP for dT in (0.0, 5.0, 15.0, 30.0)]
    assert all(b > a for a, b in zip(cops, cops[1:]))
    no_credit = model(toy, dT_subcool=15.0, subcooler_heat_useful=False).solve()
    assert no_credit.report["q_heat_kJkg"] == pytest.approx(no_credit.report["q_cond_kJkg"])
    assert no_credit.COP < cops[2]


def test_zero_superheat_and_subcooling_use_saturated_states(toy):
    res = model(toy, dT_superheat=0.0, dT_subcool=0.0).solve()
    assert res.success, res.message
    assert res.states["suction"]["phase"] == "saturated_vapour"
    assert res.report["q_subcool_kJkg"] == pytest.approx(0.0)


def test_pure_fluid_basis_does_not_matter(toy):
    base = model(toy).solve()
    for basis in ("bubble", "mean"):
        other = model(toy, evap_basis=basis, cond_basis=basis).solve()
        assert other.COP == pytest.approx(base.COP, rel=1e-9)


def test_wet_compression_is_rejected(toy):
    # the toy fluid's saturated-vapour entropy rises with temperature (an overhanging dome, like
    # butane), so a big lift with an efficient compressor and no superheat ends two-phase
    res = model(toy, T_evap=0.0, T_cond=120.0, dT_superheat=0.0, compressor=EfficiencyCurve(eta_is=1.0)).solve()
    assert not res.success and "wet compression" in res.message


def test_defaults():
    cfg = ConventionalHPConfig()
    assert (cfg.refrigerant, cfg.T_source, cfg.dTmin_source, cfg.T_steam, cfg.dTmin_steam) == ("Butane", 55.0, 5.0, 120.0, 5.0)
    assert (cfg.T_water_in, cfg.T_water_out, cfg.dTmin_water) == (80.0, 95.0, 5.0)
    assert (cfg.T_evap, cfg.T_cond, cfg.dT_subcool) == (None, None, None)   # set from the source, steam and hot water
    assert (cfg.dT_superheat, cfg.dT_discharge_superheat) == (15.0, None)
    assert not hasattr(cfg, "loss_calibration")                  # UA is always given, never calibrated
    assert cfg.compressor_Hz is None and cfg.compressor.frequency_Hz == 50.0 and cfg.compressor.rpm == 1450.0
    assert cfg.dTmin_ihx == 5.0
    assert not hasattr(cfg, "ihx") and not hasattr(cfg, "dT_evap_superheat")   # IHX on whenever superheat > 0
    assert not hasattr(cfg, "ihx_position")                       # the IHX always sits after the subcooler
    assert (cfg.loss_UA_kWK, cfg.T_ambient, cfg.loss_dT_basis) == (0.015, 15.0, "mean")
    assert cfg.compressor.name == "Recip-M"
    assert cfg.compressor.displacement_m3h == pytest.approx(56.25)
    assert recip_r600a("Recip-M", displacement_m3h=None).displacement_m3h is None
    assert recip_r600a("Recip-L").displacement_m3h == pytest.approx(151.6)    # every model now has its swept volume


def test_discharge_superheat_sets_compressor_heat_loss(toy):
    free = model(toy, T_evap=20.0, T_cond=90.0, dT_superheat=15.0).solve()
    res = model(toy, T_evap=20.0, T_cond=90.0, dT_superheat=15.0, dT_discharge_superheat=10.0).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    assert st["discharge"]["T"] == pytest.approx(100.0) and r["dT_discharge_superheat"] == pytest.approx(10.0)
    assert free.states["discharge"]["T"] > 100.0                  # adiabatic gas would be hotter
    assert r["q_compressor_loss_kJkg"] == pytest.approx(r["w_electric_kJkg"] - (st["discharge"]["h"] - st["suction"]["h"]))
    assert r["q_compressor_loss_kJkg"] > free.report["q_compressor_loss_kJkg"] >= 0.0
    assert r["w_electric_kJkg"] == pytest.approx(free.report["w_electric_kJkg"])   # same compressor input
    assert res.COP < free.COP                                    # heat lost from the shell
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-9)
    hot = model(toy, T_evap=20.0, T_cond=90.0, dT_superheat=15.0, dT_discharge_superheat=80.0).solve()
    assert not hot.success and "more than the compressor" in hot.message


SIZED = dict(T_evap=20.0, T_cond=90.0, dT_superheat=15.0)


def ua_model(toy, **kw):
    # cycle mechanics of the shell loss: read the map at the toy fluid's own (to, tc) (the toy's vapour
    # densities are not n-butane's, so the density basis would move it far across the map)
    return model(toy, **(SIZED | dict(compressor=recip_r600a("Recip-M", displacement_m3h=60.0, basis="temperature")) | kw))


def test_ua_loss_balance(toy):
    res = ua_model(toy, loss_UA_kWK=0.01, loss_dT_basis="discharge").solve()
    assert res.success, res.message
    r, st = res.report, res.states
    assert r["loss_mode"] == "UA" and r["UA_loss_kWK"] == pytest.approx(0.01)
    assert r["loss_dT_K"] == pytest.approx(st["discharge"]["T"] - 15.0)
    assert r["Q_compressor_loss_kW"] == pytest.approx(0.01 * (st["discharge"]["T"] - 15.0), rel=1e-6)
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-9)
    adiabatic = ua_model(toy, loss_UA_kWK=0.0, loss_dT_basis="discharge").solve()
    assert adiabatic.report["q_compressor_loss_kJkg"] == pytest.approx(0.0, abs=1e-6)
    lossy = ua_model(toy, loss_UA_kWK=0.02, loss_dT_basis="discharge").solve()
    assert lossy.states["discharge"]["T"] < res.states["discharge"]["T"] < adiabatic.states["discharge"]["T"]
    mean = ua_model(toy, loss_UA_kWK=0.01).solve()                       # default basis: mean
    T_wall = 0.5 * (mean.states["suction"]["T"] + mean.states["discharge"]["T"])
    assert mean.report["Q_compressor_loss_kW"] == pytest.approx(0.01 * (T_wall - 15.0), rel=1e-6)


def test_fixed_ua_carries_off_design(toy):
    hp = ConventionalHeatPump(ua_model(toy, loss_UA_kWK=0.01).c, props=toy)
    assert hp.loss_mode == "UA" and hp.loss_UA() == 0.01
    assert not hasattr(hp, "calibration") and not hasattr(hp, "calibration_note")
    off = hp.rate(T_evap=10.0, T_cond=80.0)
    assert off.success and off.report["UA_loss_kWK"] == pytest.approx(0.01)
    T_wall = 0.5 * (off.states["suction"]["T"] + off.states["discharge"]["T"])
    assert off.report["Q_compressor_loss_kW"] == pytest.approx(0.01 * (T_wall - 15.0), rel=1e-6)
    records = hp.rating_map([10.0, 20.0], [80.0, 90.0], key="T_cond", low_key="T_evap")
    assert all(rec["success"] and rec["UA_loss_kWK"] == pytest.approx(0.01) for rec in records)


def test_no_ua_uses_the_curve_discharge(toy):
    hp = ConventionalHeatPump(ua_model(toy, loss_UA_kWK=None).c, props=toy)
    assert hp.loss_mode == "curve"
    with pytest.raises(ValueError, match="loss_UA_kWK"):
        hp.loss_UA()
    with pytest.raises(TypeError):
        ConventionalHPConfig(loss_calibration=None)              # the calibration option is gone


def test_compressor_speed_scales_swept_volume(toy):
    base = dict(T_evap=20.0, T_cond=90.0, dT_superheat=15.0, compressor=recip_r600a("Recip-M", displacement_m3h=60.0))
    r50 = model(toy, **base).solve().report
    r35 = model(toy, **base, compressor_Hz=35.0).solve().report
    assert (r50["compressor_Hz"], r50["compressor_rpm"]) == (50.0, 1450.0)
    assert (r35["compressor_Hz"], r35["compressor_rpm"]) == (35.0, pytest.approx(1015.0))
    assert r35["m_dot_kgs"] == pytest.approx(r50["m_dot_kgs"] * 0.7)
    assert r35["V_suction_m3h"] == pytest.approx(r50["V_suction_m3h"] * 0.7)
    assert r35["COP"] == pytest.approx(r50["COP"])                 # no heat loss: per-kg cycle unchanged


def test_speed_with_ua_loss_changes_the_cycle(toy):
    hp = ConventionalHeatPump(ua_model(toy, loss_UA_kWK=0.03, T_ambient=25.0).c, props=toy)   # warm room: dry at 30 Hz
    slow, fast = hp.rate(compressor_Hz=30.0), hp.rate(compressor_Hz=70.0)
    assert slow.success and fast.success, (slow.message, fast.message)
    assert slow.report["UA_loss_kWK"] == fast.report["UA_loss_kWK"] == pytest.approx(0.03)
    # a fixed UA loses a larger share at low flow: lower discharge temperature and COP
    assert slow.report["compressor_loss_fraction"] > fast.report["compressor_loss_fraction"]
    assert slow.states["discharge"]["T"] < fast.states["discharge"]["T"]
    assert slow.COP < fast.COP


def test_speed_needs_a_speed_aware_curve(toy):
    curve = MapCurve(mass_flow=lambda to, tc: 3600.0, power=lambda to, tc: 100_000.0)
    res = model(toy, compressor=curve, compressor_Hz=40.0).solve()
    assert not res.success and "speed" in res.message


def test_ua_needs_a_sized_compressor(toy):
    res = model(toy, loss_UA_kWK=0.01).solve()
    assert not res.success and "sized compressor" in res.message


def test_ihx_after_subcooler_moves_superheat_off_the_evaporator(toy):
    kw = dict(T_evap=10.0, T_cond=80.0, dT_superheat=15.0, dT_subcool=20.0)
    ihx = model(toy, **kw).solve()
    assert ihx.success, ihx.message
    r, st = ihx.report, ihx.states
    assert st["evap_out"]["phase"] == "saturated_vapour" and st["evap_out"]["T"] == pytest.approx(10.0)
    assert st["suction"]["T"] == pytest.approx(25.0)
    q_ihx = st["suction"]["h"] - st["evap_out"]["h"]
    assert r["q_ihx_kJkg"] == pytest.approx(q_ihx)
    assert st["subcool_out"]["h"] - st["ihx_out"]["h"] == pytest.approx(q_ihx)      # liquid gives it up
    assert st["evap_in"]["h"] == pytest.approx(st["ihx_out"]["h"])                  # valve after the IHX
    assert st["ihx_out"]["T"] < st["subcool_out"]["T"] == pytest.approx(60.0)
    assert r["ihx_min_approach_K"] >= 5.0
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-9)
    assert r["ihx"] is True


def test_no_superheat_means_no_ihx(toy):
    res = model(toy, T_evap=10.0, T_cond=80.0, dT_superheat=0.0, dT_subcool=20.0).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    assert r["ihx"] is False and r["q_ihx_kJkg"] == 0.0 and r["ihx_min_approach_K"] is None
    assert st["suction"]["h"] == pytest.approx(st["evap_out"]["h"])               # saturated suction
    assert st["evap_in"]["h"] == pytest.approx(st["subcool_out"]["h"])            # subcooler straight to the valve


def test_ihx_sits_after_the_subcooler(toy):
    res = model(toy, T_evap=10.0, T_cond=80.0, dT_superheat=15.0, dT_subcool=20.0).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    assert st["subcool_out"]["h"] - st["ihx_out"]["h"] == pytest.approx(r["q_ihx_kJkg"])   # 3 -> 4 -> 4i
    assert r["q_subcool_kJkg"] == pytest.approx(st["cond_out"]["h"] - st["subcool_out"]["h"])
    assert st["evap_in"]["h"] == pytest.approx(st["ihx_out"]["h"]) and r["T_valve_in"] == pytest.approx(st["ihx_out"]["T"])
    assert "ihx_position" not in r
    with pytest.raises(TypeError):
        ConventionalHPConfig(ihx_position="before_subcooler")     # no longer an option


def test_ihx_approach_limits_the_superheat(toy):
    # liquid leaves the subcooler at 80 - 60 = 20 C: too cold to superheat 10 C vapour to 25 C
    res = model(toy, T_evap=10.0, T_cond=80.0, dT_superheat=15.0, dT_subcool=60.0).solve()
    assert not res.success and "IHX approach" in res.message


def test_invalid_inputs_fail_cleanly(toy):
    assert not model(toy, T_cond=5.0).solve().success
    assert "critical" in model(toy, T_cond=120.0, dP_crit_min=2.9).solve().message
    assert not model(toy, dT_subcool=-1.0).solve().success


def test_heat_duty_scaling(toy):
    res = model(toy, Q_heat_kW=100.0).solve()
    r = res.report
    assert r["sized_by"] == "Q_heat_kW"
    assert r["Q_heat_kW"] == pytest.approx(100.0)
    assert r["W_electric_kW"] == pytest.approx(100.0 / res.COP)
    assert r["Q_evap_kW"] + r["m_dot_kgs"] * r["w_to_gas_kJkg"] == pytest.approx(r["Q_cond_kW"] + r["Q_subcool_kW"])


def test_manufacturer_map_sizes_the_cycle(toy):
    comp = recip_r600a("Recip-M", displacement_m3h=60.0)
    res = model(toy, T_evap=20.0, T_cond=90.0, compressor=comp).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    to_m, tc_m = r["map_to"], r["map_tc"]                    # R600a point with the same vapour densities
    assert r["eta_vol"] == pytest.approx(comp.eta_vol(to_m, tc_m))
    assert r["eta_is_overall"] == pytest.approx(comp.eta_is(to_m, tc_m))      # electrical basis
    assert r["m_dot_kgs"] == pytest.approx(r["eta_vol"] * 60.0 / 3600.0 * st["suction"]["D"])
    assert r["V_suction_m3h"] == pytest.approx(r["eta_vol"] * 60.0)
    dh_is = st["isentropic_discharge"]["h"] - st["suction"]["h"]
    assert r["W_electric_kW"] == pytest.approx(r["m_dot_kgs"] * dh_is / r["eta_is_overall"])
    # suction-gas-cooled semi-hermetic: all electrical input ends up in the gas
    assert r["w_to_gas_kJkg"] == pytest.approx(r["w_electric_kJkg"])
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-9)
    assert r["in_envelope"] and not res.warnings


def test_manufacturer_map_warns_outside_envelope(toy):
    res = model(toy, T_evap=40.0, T_cond=45.0, dT_superheat=0.0,       # no IHX: 5 K subcooling leaves the liquid too cold
                compressor=recip_r600a("Recip-M", displacement_m3h=None)).solve()
    assert res.success
    assert not res.report["in_envelope"]
    assert any("envelope" in w for w in res.warnings)


def test_rating_map(toy):
    hp = model(toy, compressor=recip_r600a("Recip-M", displacement_m3h=60.0))
    records = hp.rating_map([0.0, 10.0, 20.0], [60.0, 80.0], key="T_cond", low_key="T_evap")
    assert len(records) == 6 and all(rec["success"] for rec in records)
    cop = {(rec["T_evap"], rec["T_cond"]): rec["COP"] for rec in records}
    assert cop[(20.0, 60.0)] > cop[(0.0, 80.0)]


# ---------------------------------------------------------------- compressor curves


def test_bicubic_matches_workbook_check_cell():
    # workbook cells C22/C23 for Recip-M at to = -10 C, tc = 100 C
    c = R600A_COEFFICIENTS["Recip-M"]
    assert bicubic_polynomial(c["volumetric"])(-10.0, 100.0) == pytest.approx(0.4889602100683594, rel=1e-12)
    assert bicubic_polynomial(c["isentropic"])(-10.0, 100.0) == pytest.approx(0.6567477345580446, rel=1e-12)


def test_envelope_polygon():
    assert in_envelope(20.0, 80.0, R600A_ENVELOPE)
    assert in_envelope(-10.0, 50.0, R600A_ENVELOPE)          # on the left edge
    assert not in_envelope(40.0, 40.0, R600A_ENVELOPE)       # below the lower boundary
    assert not in_envelope(0.0, 115.0, R600A_ENVELOPE)       # above the upper-left boundary
    assert in_envelope(5.0, 50.0, ((0.0, 10.0), (40.0, 60.0)))   # rectangle form


def test_en12900_polynomial_and_fit_roundtrip():
    coeffs = [1000.0, -20.0, 30.0, 0.1, -0.2, 0.3, 1e-3, -2e-3, 3e-3, -1e-4]
    f = en12900_polynomial(coeffs)
    assert f(0.0, 0.0) == pytest.approx(1000.0)
    to = [t for t in (-10, 0, 10, 20) for _ in range(4)]
    tc = [c for _ in range(4) for c in (40, 55, 70, 85)]
    assert fit_en12900(to, tc, [f(a, b) for a, b in zip(to, tc)]) == pytest.approx(coeffs, rel=1e-6, abs=1e-9)


def test_map_curve_superheat_correction(toy):
    # constant map rated at 10 K superheat; run at 5 K -> mass flow scales with suction density
    curve = MapCurve(mass_flow=lambda to, tc: 3600.0, power=lambda to, tc: 100_000.0, rated_superheat=10.0)
    res = model(toy, compressor=curve).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    assert 0.5 < r["eta_is_overall"] < 1.0
    rated = toy.ref_PT(res.P_L, st["evap_dew"]["T"] + 10.0)
    assert r["m_dot_kgs"] == pytest.approx(1.0 * st["suction"]["D"] / rated["D"])
    assert r["W_electric_kW"] == pytest.approx(100.0)
    at_rating = model(toy, compressor=curve, dT_superheat=10.0).solve()
    assert at_rating.report["m_dot_kgs"] == pytest.approx(1.0)
    en = EN12900Curve([100_000.0] + [0.0] * 9, [3600.0] + [0.0] * 9, rated_superheat=10.0)
    assert model(toy, compressor=en).solve().COP == pytest.approx(res.COP)


def test_efficiency_curve_forms(toy):
    base = model(toy, compressor=EfficiencyCurve(eta_is=0.7, motor_eff=0.9)).solve()
    poly = model(toy, compressor=EfficiencyCurve(eta_is=[0.7, 0.0], motor_eff=0.9)).solve()
    func = model(toy, compressor=EfficiencyCurve(eta_is=lambda pr: 0.7, motor_eff=0.9)).solve()
    assert base.COP == pytest.approx(poly.COP) == pytest.approx(func.COP)
    pareto = model(toy, compressor=EfficiencyCurve()).solve()
    assert pareto.success and 0.6 < pareto.report["eta_is_shaft"] < 0.9


# ---------------------------------------------------------------- REFPROP-backed


@pytest.fixture
def refprop(hpl):
    return hpl


def test_default_case(refprop):
    res = ConventionalHeatPump(ConventionalHPConfig()).solve()
    assert res.success, res.message
    r, st = res.report, res.states
    assert r["loss_mode"] == "UA" and r["UA_loss_kWK"] == pytest.approx(0.015)          # fixed 15 W/K
    T_wall = 0.5 * (st["suction"]["T"] + st["discharge"]["T"])
    assert r["Q_compressor_loss_kW"] == pytest.approx(0.015 * (T_wall - 15.0), rel=1e-5)
    # 55 C source, 120 C steam, hot water 80 -> 95 C: evaporating 50 C, condenser outlet 125 C, liquid out 85 C
    assert (r["T_source"], r["T_steam"], r["T_water_in"], r["T_water_out"]) == (55.0, 120.0, 80.0, 95.0)
    assert r["T_evap_dew"] == pytest.approx(50.0, abs=1e-6) and r["T_cond_bubble"] == pytest.approx(125.0, abs=1e-6)
    assert r["dT_subcool"] == pytest.approx(40.0, abs=1e-6) and r["subcooler_min_approach_K"] == pytest.approx(5.0, abs=1e-6)
    assert r["m_water_kgs"] == pytest.approx(r["Q_subcool_kW"] / (4.19 * 15.0)) and r["m_water_kgs"] == pytest.approx(0.284, rel=0.02)
    assert st["suction"]["T"] == pytest.approx(65.0, abs=1e-6)
    assert st["discharge"]["T"] > 130.0
    assert st["subcool_out"]["T"] == pytest.approx(85.0, abs=1e-6)
    assert st["evap_out"]["T"] == pytest.approx(50.0, abs=1e-6) and st["evap_out"]["phase"] == "saturated_vapour"
    assert r["q_ihx_kJkg"] > 0 and st["ihx_out"]["T"] < 85.0 and r["ihx_min_approach_K"] >= 5.0
    assert 0.0 < r["compressor_loss_fraction"] < 0.2 and r["q_compressor_loss_kJkg"] > 0.0
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-6)
    # density basis: n-butane at 50/125 C, 15 K suction superheat, reads the R600a map at about
    # 35.9/108.6 C, inside the envelope
    assert r["in_envelope"] and r["map_to"] == pytest.approx(35.87, abs=0.1) and r["map_tc"] == pytest.approx(108.59, abs=0.1)
    assert r["V_suction_m3h"] == pytest.approx(r["eta_vol"] * 56.25)
    assert r["Q_heat_kW"] == pytest.approx(r["m_dot_kgs"] * r["q_heat_kJkg"])


def test_mixture_uses_the_same_fixed_ua(refprop):
    # 10 K subcooling: with the default 40 K the liquid is too cold for the IHX to give 15 K suction superheat
    mix = ConventionalHeatPump(ConventionalHPConfig(refrigerant="Propane;Isopentane|0.2;0.8|mass", T_evap=40.0,
                                                    T_cond=110.0, cond_basis="mean", dT_subcool=10.0))
    assert mix.loss_UA() == pytest.approx(0.015)
    res = mix.solve()
    assert res.success, res.message
    assert res.report["UA_loss_kWK"] == pytest.approx(0.015)


def test_default_ua_off_design(refprop):
    hp = ConventionalHeatPump(ConventionalHPConfig())
    res = hp.rate(T_evap=40.0, T_cond=110.0)
    assert res.success, res.message
    r = res.report
    assert r["UA_loss_kWK"] == pytest.approx(hp.loss_UA())
    T_wall = 0.5 * (r["T_suction"] + r["T_discharge"])
    assert r["Q_compressor_loss_kW"] == pytest.approx(r["UA_loss_kWK"] * (T_wall - 15.0), rel=1e-6, abs=1e-9)
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-6)


def test_butane_with_recip_m(refprop):
    cfg = ConventionalHPConfig(refrigerant="Butane", T_evap=20.0, T_cond=90.0, dT_subcool=10.0,
                               compressor=recip_r600a("Recip-M", displacement_m3h=60.0))
    res = ConventionalHeatPump(cfg).solve()
    assert res.success, res.message
    r = res.report
    assert r["T_evap_dew"] == pytest.approx(20.0, abs=1e-6)
    assert r["T_cond_dew"] == pytest.approx(90.0, abs=1e-6)
    assert r["evap_glide_K"] == pytest.approx(0.0, abs=1e-6)
    assert r["energy_balance_residual_kJkg"] == pytest.approx(0.0, abs=1e-6)
    assert 2.0 < res.COP < r["COP_carnot"]
    assert r["Q_heat_kW"] == pytest.approx(r["m_dot_kgs"] * r["q_heat_kJkg"])
    assert r["in_envelope"]


def test_zeotropic_mixture_saturation_bases(refprop):
    fluid = "Propane;Isopentane|0.2;0.8|mass"
    base = dict(refrigerant=fluid, T_evap=30.0, T_cond=110.0, dT_subcool=5.0, dT_superheat=20.0,
                loss_UA_kWK=None, compressor=EfficiencyCurve(eta_is=0.7))
    dew = ConventionalHeatPump(ConventionalHPConfig(**base)).solve()
    mean = ConventionalHeatPump(ConventionalHPConfig(**base, evap_basis="mean", cond_basis="mean")).solve()
    assert dew.success and mean.success, (dew.message, mean.message)
    assert dew.report["T_evap_dew"] == pytest.approx(30.0, abs=1e-6)
    assert dew.report["evap_glide_K"] > 5.0
    m = mean.report
    assert 0.5 * (m["T_evap_dew"] + m["T_evap_bubble"]) == pytest.approx(30.0, abs=1e-5)
    assert 0.5 * (m["T_cond_dew"] + m["T_cond_bubble"]) == pytest.approx(110.0, abs=1e-5)
    assert mean.P_L > dew.P_L                                     # mean is above the dew pressure
