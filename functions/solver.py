__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from systems import TRHP
from components import GasCooler
from streams import MassStream
from scipy.optimize import root_scalar, minimize
from functions.readWriteOutput import print_result_to_csv
import numpy as np
from functions.decorator_functions import time_func

# @time_func
def solve_heat_pump_model(hp: TRHP) -> TRHP:
    #Solve heat pump cycle
    # try:
    hp = _find_hp_with_min_discharge_pressure(hp)
    if (abs(hp.evaporator.sink_in.temperature - hp.Tref_evap_inlet) > 0.1 and _check_hp_model_convergance(hp)):
        r = root_scalar(
                        f=_find_evap_inlet_temperature_offset,
                        method='brentq',
                        x0=hp.Tref_evap_outlet,
                        bracket=(hp.Tref_evap_outlet, hp.T_source_in - hp.dTmin_phx),
                        xtol=hp.tol,
                        args=(hp),
                        )
        hp.Tref_evap_outlet = r['root']
        hp.evaporator.sink_out.props_in_ab('TQ', hp.Tref_evap_outlet, 1.0)
        hp = _find_hp_with_min_discharge_pressure(hp)

    hp.converged = _check_hp_model_convergance(hp)
    hp.COP = _get_COP(hp)
    hp.cycle_type = _get_hp_cycle_type(hp)
    # except:
    #     _print_failed_to_solve_report(hp)
    #     print('Heat pump solver failed with an exception. The failed case has been recorded.')
    #     hp.converged = False       
    return hp

def _check_hp_model_convergance(hp: TRHP) -> bool:
    return True if ( sum([c.power for c in hp.compressors]) > 0.0 and hp.converged == True and hp.feasible == True ) else False

def _get_COP(hp: TRHP) -> float:
    return sum([g.duty for g in hp.gas_coolers]) / sum([c.power for c in hp.compressors]) if hp.converged else -np.inf    

def _get_hp_cycle_type(hp: TRHP) -> str:
    return 'Transcritical' if hp.comp_discharge.pressure > hp.comp_discharge.P_crit else 'Subcritical'

def _find_evap_inlet_temperature_offset(Tref_evap_outlet: float, hp: TRHP) -> float:
    hp.evaporator.sink_out.props_in_ab('TQ', Tref_evap_outlet, 1.0)
    hp = _find_hp_with_min_discharge_pressure(hp)
    return hp.evaporator.sink_in.temperature - hp.Tref_evap_inlet

def _find_hp_with_min_discharge_pressure(hp: TRHP) -> TRHP:
    bnds, x0, hp.feasible, bnd_val = _get_discharge_pressure_bounds(hp)
    if hp.feasible == True:
        r = root_scalar(f=_find_current_gc_dtmin_offset, #find the gas discharge pressure resulting in the target dtmin
                        method='brentq', 
                        bracket=bnds, 
                        x0=x0, 
                        xtol=hp.tol,
                        args=(hp),
                        )
        hp = _do_sequential_iterative_cycle_analysis(r['root'], hp)
        hp = _update_pinch_props_for_gcs(P_comp_out=r['root'], hp=hp)    
        hp.converged = r['converged']
    else:
        _print_failed_to_solve_report(hp)
        raise Exception("Infeasible pressure bounds. The failed case has been recorded.")
    return hp

def _get_discharge_pressure_bounds(hp: TRHP) -> list[tuple[float, float], float, bool, tuple[float, float]]:
    feasible: bool = False
    P_lb: float
    P_ub: float
    dt_neg: float
    dt_plus: float
    P_i: float

    # Positive dtmin is P_ub and negative dtmin is P_lb    
    P_lb = P_ub = hp.comp_discharge.pressure 
    dt_plus = dt_neg = _find_current_gc_dtmin_offset(P_lb, hp)
    if -1.0 < dt_plus < 1.0:
        P_lb = P_ub = hp.comp_discharge.pressure * 0.95
        dt_plus = dt_neg = _find_current_gc_dtmin_offset(P_lb, hp)

    def _recursive_find_bnds(P_suction: float, P_0: float, dt_0: float, P_1: float, dt_1: float, step: float, sign: int, i: int):
        if dt_0 * sign > 1.0:
            dt_1 = dt_0
            P_1 = P_0
        P_0 = max(P_suction + (P_0 - P_suction) * step, P_suction)
        dt_0 = _find_current_gc_dtmin_offset(P_0, hp)
        if (dt_0 * sign < -1.0 or i >= 100):
            return P_0, dt_0, P_1, dt_1  
        else: 
            return _recursive_find_bnds(P_suction, P_0, dt_0, P_1, dt_1, step, sign, i+1)

    if dt_plus > 0.0:
        P_lb, dt_neg, P_ub, dt_plus = _recursive_find_bnds(hp.comp1.flow_in.pressure, P_lb, dt_neg, P_ub, dt_plus, 0.5, 1, 0)
    else:
        P_ub, dt_plus, P_lb, dt_neg = _recursive_find_bnds(hp.comp1.flow_in.pressure, P_ub, dt_plus, P_lb, dt_neg, 2.0, -1, 0)

    if abs(dt_plus * dt_neg) == -dt_plus * dt_neg and abs(dt_plus * dt_neg) > 0.0:
        feasible = True
        P_i = P_lb - dt_neg / (dt_plus - dt_neg) * (P_ub - P_lb)

    return (P_lb, P_ub), P_i, feasible, (dt_plus, dt_neg)

def _find_current_gc_dtmin_offset(P_comp_out: float, hp: TRHP) -> float:
    dtmin_ls: list[list[float,float]]

    dtmin_ls = _get_all_possible_gc_dtmin_pos(P_comp_out, hp)

    dtmin_ls.pop(-1)            # Ignore the cold-side dtmin of the inlet air to the gas cooler
    if dtmin_ls[-1][1] < 1e-3:  # Ignore the mid exchanger dtmin of the inlet air to the gas cooler if it is the cold-side
        dtmin_ls.pop(-1)

    return min(dtmin_ls)[0] - hp.dT_GC_min_tar   

def _get_all_possible_gc_dtmin_pos(P_comp_out: float, hp: TRHP) -> list[list[float,float]]:
    dtmin_ls: list[float]
    gc: GasCooler

    hp = _do_sequential_iterative_cycle_analysis(P_comp_out, hp)
    dtmin_ls = []
    for gc in hp.gas_coolers:
        dtmin_ls.append([gc.source_in.temperature - gc.sink_out.temperature, 1.0])  # Pinch on the hot side of the gc
        dtmin_ls.append(_find_gc_dt_approaches(gc))                       # Pinch in the middle of the gc
        dtmin_ls.append([gc.source_out.temperature - gc.sink_in.temperature, 0.0])  # Pinch on the cold side of the gc

    return dtmin_ls

def _do_sequential_iterative_cycle_analysis(P_comp_out: float, hp: TRHP) -> TRHP:
    hp.comp_discharge.pressure = P_comp_out
    i: int = 0
    h: float 
    P: float 
    W: float 
    Q: float 
    most_neg_residual: float = 0.0

    while ((i == 0 or most_neg_residual < 0.0) and i < 10): #Check for model convergance
        # Saves pressure and temperature properties to check for solution convergance
        h = hp.exp_valve.flow_in.enthalpy
        P = hp.comp_discharge.pressure - hp.gas_cooler1.source_out.pressure
        W = sum([c.power for c in hp.compressors])
        Q = sum([g.duty for g in hp.gas_coolers])
        i += 1

        #Solve the performance of each component and update each special relationship
        [obj.solve() for obj in hp.operations]

        most_neg_residual = min(
            hp.tol * h - abs(hp.exp_valve.flow_in.enthalpy - h),
            hp.tol * P - abs(hp.comp_discharge.pressure - hp.gas_cooler1.source_out.pressure - P),
            hp.tol * W - abs(sum([c.power for c in hp.compressors]) - W),
            hp.tol * Q - abs(sum([g.duty for g in hp.gas_coolers]) - Q),    
        )
    hp = _set_gas_coolers_sink_temperatures(hp)
    return hp

def _set_gas_coolers_sink_temperatures(hp: TRHP) -> TRHP:
    gc: GasCooler

    hp.gas_coolers[0].sink_out.temperature = hp.T_air_tar
    hp.gas_coolers[0].sink_out.props_in_mPT()

    m_air = sum([gc.duty for gc in hp.gas_coolers]) / (hp.gas_coolers[0].sink_out.enthalpy - hp.gas_coolers[-1].sink_in.enthalpy)
    hp.gas_coolers[0].sink_out.mass_flow = m_air
    hp.gas_coolers[0].sink_out.heat_flow = m_air * hp.gas_coolers[0].sink_out.enthalpy

    for gc in hp.gas_coolers:
        gc.sink_pinch.mass_flow = m_air
        gc.sink_in.mass_flow = m_air
        gc.sink_in.enthalpy = gc.sink_out.enthalpy - gc.duty / m_air
        gc.sink_in.props_in_mPH()

    return hp

def _find_gc_dt_approaches(gc: GasCooler) -> list[float, float]:
    pinch_point: list[float, float]
    r: dict

    gc.source_pinch.pressure = gc.source_out.pressure
    gc.sink_pinch.pressure = gc.sink_in.pressure
    
    if (_calc_jac_at_given_point(gc.source_in, gc.sink_out) > 0.0 # Implies converging temperature profiles from hot-side
        or 
        _calc_jac_at_given_point(gc.source_out, gc.sink_in) < 0.0 # Implies converging temperature profiles from cold-side
        ):
        r = minimize(fun=_calc_pinch_dtmin_and_jac, args=(gc), method='TNC', x0=[gc.duty / 2], bounds=[(0.0, gc.duty)], jac=True) # l-bfgs-b args=(gc), 
        pinch_point = [r['fun'], r['x'] / gc.duty]
    else:
        pinch_point = min(
            [_calc_dtmin_at_given_point(gc.source_in, gc.sink_out), 1.0], 
            [_calc_dtmin_at_given_point(gc.source_out, gc.sink_in), 0.0]
            )
    
    return pinch_point

def _calc_pinch_dtmin_and_jac(Q: float, gc: GasCooler):
    gc = _calc_pinch_props_at_Q(Q, gc)
    return (
            _calc_dtmin_at_given_point(gc.source_pinch, gc.sink_pinch), #objective
            _calc_jac_at_given_point(gc.source_pinch, gc.sink_pinch) #Jacobian
            )

def _update_pinch_props_for_gcs(P_comp_out: float, hp: TRHP) -> TRHP:
    dtmin_ls = _get_all_possible_gc_dtmin_pos(P_comp_out, hp)
    i = 1
    for gc in hp.gas_coolers:
        gc = _calc_pinch_props_at_Q(Q=dtmin_ls[i][1], gc=gc)
        i += 3
    return hp

def _calc_pinch_props_at_Q(Q: float, gc: GasCooler) -> GasCooler:
    gc.source_pinch.enthalpy = gc.source_out.enthalpy + Q / gc.source_out.mass_flow
    gc.source_pinch.props_in_mPH()
    gc.sink_pinch.enthalpy = gc.sink_in.enthalpy + Q / gc.sink_in.mass_flow
    gc.sink_pinch.props_in_mPH()
    return gc

def _calc_dtmin_at_given_point(source: MassStream, sink: MassStream) -> float:
    return source.temperature - sink.temperature

def _calc_jac_at_given_point(source: MassStream, sink: MassStream) -> float:
    return 1 / (source.heat_capacity * source.mass_flow) - 1 / (sink.heat_capacity * sink.mass_flow)

def _print_failed_to_solve_report(hp: TRHP, filename: str ='results/Cases-Unsolved') -> None:
    print('Heat pump model failed to solve.')
    hp.get_summary()
    print_result_to_csv(hp.output, file_name=filename)