__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

from systems import HeatPump
from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator, SetRelationship, CustomRelationship
from streams import MassStream, PropertyPackage
from scipy.optimize import root_scalar, minimize
from functions.readWriteOutput import print_result_to_csv
import numpy as np

def solve_heat_pump_model(hp: HeatPump) -> HeatPump:
    #Solve heat pump cycle
    if hp.feasible:
        hp, pressure_bounds = _get_solve_evap_out_temperature(hp.Tsat_evap_tar, hp)

        # if (abs(r) > 1e-6 and abs(hp.evaporator.sink_in.temperature - hp.evaporator.sink_out.temperature) > 0.1):
        #     r = root_scalar(
        #                     f=_get_solve_evap_out_temperature,
        #                     method='secant',
        #                     x0=hp.Tsat_evap_tar,
        #                     bracket=(hp.Tsat_evap_tar, hp.T_source_in),
        #                     xtol=hp.tol,
        #                     args=(hp, pressure_bounds)
        #                 )

    if _check_hp_model_convergance(hp):
        hp.converged = True
        hp.COP = sum([g.duty for g in hp.gas_coolers]) / sum([c.power for c in hp.compressors])
    else:
        hp.converged = False
        hp.COP = -np.inf

    hp.cycle_type = _get_hp_cycle_type(hp)

    # if hp.disp or 1:
    #     hp.get_cycle_plot()

    return hp

def _get_solve_evap_out_temperature(T_evap_sat_out: float, hp: HeatPump, prev_pressure_bounds: tuple[float, float] = None) -> HeatPump:
    #Solve heat pump cycle
    hp.converged = False
    hp.COP = -np.inf
    hp = _set_evap_out_temperature(hp, T_evap_sat_out) #Is this the best place????

    try:
        pressure_bounds, P_i, hp.feasible = _get_discharge_pressure_bounds(hp, prev_pressure_bounds)
        if hp.feasible == True:
            hp = _find_hp_with_min_discharge_pressure(hp=hp, x0=P_i, bnds=pressure_bounds)
    except:
        _print_failed_to_solve_report(hp)
        print('Function _get_solve_gas_pressure exited with an exception.')
        hp.converged = False
        pressure_bounds = (0.0, 0.0)

    return hp, pressure_bounds

def _check_hp_model_convergance(hp: HeatPump) -> bool:
    return True if ( sum([c.power for c in hp.compressors]) > 0.0 and hp.converged == True ) else False

def _get_hp_cycle_type(hp: HeatPump) -> str:
    return 'Transcritical' if hp.comp_discharge.pressure > hp.RP_ref.P_crit else 'Subcritical'

def _set_evap_out_temperature(hp: HeatPump, T_evap_sat_out) -> HeatPump:
    hp.evaporator.sink_out.temperature = T_evap_sat_out
    hp.evaporator.sink_out.props_in_mTX()
    return hp

def _get_discharge_pressure_bounds(hp: HeatPump, prev_pressure_bounds: tuple[float, float] =None) -> list[tuple[float, float], float, bool]:
    feasible: bool
    P_lb: float
    P_ub: float
    P_i: float
    dt_plus: float
    dt_neg: float
    P_hist: np.array
    dt_hist: np.array
    i: int
    
    feasible = False
    if prev_pressure_bounds == None:
        P_lb = hp.comp_discharge.pressure
        P_ub = P_lb
        P_i = P_lb
        dt_plus = _find_current_gc_dtmin(P_lb, hp)
        dt_neg = dt_plus
    else:
        P_lb = prev_pressure_bounds[0]
        P_ub = prev_pressure_bounds[1]
        dt_plus = _find_current_gc_dtmin(P_lb, hp)
        dt_neg = _find_current_gc_dtmin(P_ub, hp)

    P_hist = np.array(P_lb)
    dt_hist = np.array(dt_plus)

    i = 0
    m_step = 0.75
    if dt_neg > 0.0:
        #Find a lower bound pressure resulting in an air temperature less than the target
        while dt_neg > 0.0 and i < 100:
            dt_plus = dt_neg
            P_ub = P_lb
            P_lb = hp.comp1.flow_in.pressure + (P_lb - hp.comp1.flow_in.pressure) * m_step
            dt_neg = _find_current_gc_dtmin(P_lb, hp)
            P_hist = np.append(P_hist, P_lb)
            dt_hist = np.append(dt_hist, dt_neg)
            i += 1

    else:
        #Find a upper bound pressure resulting in an air temperature greater than the target
        while dt_plus < 0.0 and i < 100:
            dt_neg = dt_plus
            P_lb = P_ub
            P_ub = hp.comp1.flow_in.pressure + (P_ub - hp.comp1.flow_in.pressure) / m_step
            dt_plus = _find_current_gc_dtmin(P_ub, hp)
            P_hist = np.append(P_hist, P_ub)
            dt_hist = np.append(dt_hist, dt_plus)
            i += 1

    if P_lb < hp.comp1.flow_in.pressure:
        #Check the low pressure bound remains above the evaporator pressure
        P_lb = hp.comp1.flow_in.pressure
        dt_neg = _find_current_gc_dtmin(P_lb, hp)
        if P_ub < P_lb:
            #Check the upper pressure bound stays above the low pressure bound
            dt_neg = dt_plus
    
    if abs(dt_plus * dt_neg) == -dt_plus * dt_neg and abs(dt_plus * dt_neg) > 0.0:
        feasible = True
        P_i = P_lb - dt_neg / (dt_plus - dt_neg) * (P_ub - P_lb)
        if P_i * 1.01 < P_ub:
            dt = _find_current_gc_dtmin(P_i * 1.01, hp)
            if dt > 0.0:
                P_ub = P_i * 1.01
                dt_plus = dt
        if P_i * 0.99 > P_lb:
            dt = _find_current_gc_dtmin(P_i * 0.99, hp)
            if dt < 0.0:
                P_lb = max(P_i * 0.99, P_lb)
                dt_neg = dt
        
        P_i = P_lb - dt_neg / (dt_plus - dt_neg) * (P_ub - P_lb)

    return (P_lb, P_ub), P_i, feasible

def _find_hp_with_min_discharge_pressure(hp: HeatPump, x0: float, bnds: tuple[float, float]) -> list[float, float]:
    r = root_scalar(f=_find_current_gc_dtmin, #find the gas discharge pressure resulting in the target dtmin
                    method='brentq', 
                    bracket=bnds, 
                    x0=x0, 
                    xtol=hp.tol,
                    args=(hp),
                    )
    hp = _do_sequential_iterative_cycle_analysis(r['root'], hp)
    hp = _update_pinch_props_for_gcs(P_comp_out=r['root'], hp=hp)    
    hp.converged = r['converged']
    return hp

def _find_current_gc_dtmin(P_comp_out: float, hp: HeatPump) -> float:
    dtmin_ls: list[list[float,float]]

    dtmin_ls = _get_all_possible_gc_dtmin_pos(P_comp_out, hp)

    dtmin_ls.pop(-1)            # Ignore the cold-side dtmin of the inlet air to the gas cooler
    if dtmin_ls[-1][1] < 1e-3:  # Ignore the mid exchanger dtmin of the inlet air to the gas cooler if it is the cold-side
        dtmin_ls.pop(-1)

    return min(dtmin_ls)[0] - hp.dT_GC_min_tar    

def _get_all_possible_gc_dtmin_pos(P_comp_out: float, hp: HeatPump) -> list[list[float,float]]:
    dtmin_ls: list[float]
    gc: GasCooler

    hp = _do_sequential_iterative_cycle_analysis(P_comp_out, hp)
    dtmin_ls = []
    for gc in hp.gas_coolers:
        dtmin_ls.append([gc.source_in.temperature - gc.sink_out.temperature, 1.0])  # Pinch on the hot side of the gc
        dtmin_ls.append(_find_gc_dt_approaches(gc))                       # Pinch in the middle of the gc
        dtmin_ls.append([gc.source_out.temperature - gc.sink_in.temperature, 0.0])  # Pinch on the cold side of the gc

    return dtmin_ls

def _do_sequential_iterative_cycle_analysis(P_comp_out: float, hp: HeatPump) -> HeatPump:
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

def _set_gas_coolers_sink_temperatures(hp: HeatPump) -> HeatPump:
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

def _update_pinch_props_for_gcs(P_comp_out: float, hp: HeatPump) -> HeatPump:
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
    try: # if gc.source_pinch.heat_capacity is float:
        return 1 / (source.heat_capacity * source.mass_flow) - 1 / (sink.heat_capacity * sink.mass_flow)
    except: # else:
        return - 1 / (sink.heat_capacity * sink.mass_flow)

def _print_failed_to_solve_report(hp: HeatPump, filename: str ='results/Cases-Unsolved') -> None:
    print('Heat pump model failed to solve.')
    hp.get_summary()
    print_result_to_csv(hp.output, file_name=filename)