__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

from systems.heat_pump import HeatPump
from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator, SetRelationship, CustomRelationship
from streams import MassStream, PropertyPackage

def set_trhp_model(hp: HeatPump, version: int = 1) -> HeatPump:
    hp = _get_trhp_structure(hp)
    if _check_input_param_feasibility(hp) == True:
        hp.feasible = True
        hp = _get_trhp_components(hp)
        hp = _get_trhp_streams(hp)
        if version == 2:
            hp = _get_trhp_logic_blocks(hp)
            hp = _get_trhp_solver_pathway(hp)
        else:
            hp = _get_trhp_logic_blocks_v1(hp)
            hp = _get_trhp_solver_pathway_v1(hp)        
    else:
        hp.feasible = False
    return hp

def _check_input_param_feasibility(hp: HeatPump):
    return True if ((hp.T_air_tar - hp.T_air_inlet) * 2/3 >= hp.dT_phx1_cold + hp.dT_phx2_cold + hp.dT_phx3_cold) else False

def _get_trhp_structure(hp: HeatPump) -> HeatPump:
    #Check the PHE combo
    hp.dT_phx1_cold = hp.dT_phx1_cold if hp.dT_phx1_cold > 0.0 else 0.0
    hp.dT_phx2_cold = hp.dT_phx2_cold if hp.n_comp > 1 and hp.dT_phx2_cold > 0.0 else 0.0
    hp.dT_phx3_cold = hp.dT_phx3_cold if hp.n_comp > 2 and hp.dT_phx3_cold > 0.0 else 0.0
    
    #Check the compressor combo
    if hp.comp_PR_uniform or hp.comp2_power_frac + hp.comp3_power_frac >= 1.0:
        hp.comp2_power_frac = 1 / hp.n_comp if hp.n_comp > 1 else 0.0
        hp.comp3_power_frac = 1 / hp.n_comp if hp.n_comp > 2 else 0.0
    else:
        hp.comp2_power_frac = hp.comp2_power_frac if hp.n_comp > 1 else 0.0
        hp.comp3_power_frac = hp.comp3_power_frac if hp.n_comp > 2 else 0.0
    
    if hp.comp2_power_frac == 0.0:
        hp.comp3_power_frac = 0.0
    
    if hp.n_comp == 3 and hp.comp3_power_frac == 0.0:
        hp.n_comp = 2
    
    if hp.n_comp == 2 and hp.comp2_power_frac == 0.0:
        hp.n_comp = 1

    hp.comp1_power_frac = 1 - (hp.comp2_power_frac + hp.comp3_power_frac)

    #Binary variables for the heat pump structure
    hp.comp2_on = 1  if hp.n_comp > 1 else 0
    hp.comp3_on = 1  if hp.n_comp > 2 else 0
    hp.ihx_on = 1    if hp.T_source_in < hp.T_gc_out else 0
            
    hp.phx1_on = 1   if hp.dT_phx1_cold > 0.0 else 0
    hp.phx2_on = 1   if hp.dT_phx2_cold > 0.0 and hp.comp2_on == 1 else 0
    hp.phx3_on = 1   if hp.dT_phx3_cold > 0.0 and hp.comp3_on == 1 else 0

    return hp

def _get_trhp_components(hp: HeatPump) -> HeatPump:
    #Initialisation of unit operations

    hp.comp2_on = 1  if hp.n_comp > 1 else 0
    hp.comp3_on = 1  if hp.n_comp > 2 else 0
    hp.ihx_on = 1    if hp.T_source_in < hp.T_gc_out else 0
            
    hp.phx1_on = 1   if hp.dT_phx1_cold > 0.0 else 0
    hp.phx2_on = 1   if hp.dT_phx2_cold > 0.0 else 0
    hp.phx3_on = 1   if hp.dT_phx3_cold > 0.0 else 0

    hp.comp1 = Compressor(name='Compressor', motor_eff=hp.eff_motor, model_index=1)
    if hp.comp2_on == 1:
        hp.comp2 = Compressor(name='Compressor', motor_eff=hp.eff_motor, model_index=1)
        hp.gas_cooler3 = GasCooler(name='Gas cooler', dTmin=hp.dT_GC_min_tar)
        if hp.comp3_on == 1:
            hp.comp3 = Compressor(name='Compressor', motor_eff=hp.eff_motor, model_index=1) 
            hp.gas_cooler4 = GasCooler(name='Gas cooler', dTmin=hp.dT_GC_min_tar)
    
    hp.gas_cooler1 = GasCooler(name='Gas cooler', dTmin=hp.dT_GC_min_tar)
    hp.gas_cooler2 = GasCooler(name='Gas cooler', dTmin=hp.dT_GC_min_tar)
    
    hp.ihx = PlateHeatExchanger(name='Internal heat exchanger', dTmin=hp.dTmin, model='dTmin_terminal')
    if hp.phx1_on == 1:
        hp.phx1 = PlateHeatExchanger(name='Recuperator', dT_c=hp.dT_phx1_cold, model='fixed cold delta temperature')
    if hp.phx2_on == 1:
        hp.phx2 = PlateHeatExchanger(name='Recuperator', dT_c=hp.dT_phx2_cold, model='fixed cold delta temperature') 
    if hp.phx3_on == 1:       
        hp.phx3 = PlateHeatExchanger(name='Recuperator', dT_c=hp.dT_phx3_cold, model='fixed cold delta temperature')
    
    hp.exp_valve = ExpansionValve(name='Expansion valve')
    hp.evaporator = Evaporator(name='Evaporator', model='saturated gas outlet evaporator')
    hp.evap_sup = PlateHeatExchanger(name='Evaporator', model='fixed cold temperature outlet')
    
    return hp

def _get_trhp_streams(hp: HeatPump) -> HeatPump:
    temp_ref_out: MassStream

    #Evaporator superheat section
    hp.evap_sup.source_in = MassStream( name='source In', fluid=hp.RP_src, pressure=0.101, temperature=hp.T_source_in, mass_flow=10^6)
    hp.evap_sup.source_out = MassStream( name='source Out', fluid=hp.RP_src, pressure=0.101, temperature=hp.T_source_in, mass_flow=10^6)
    hp.evap_sup.sink_in = MassStream( name='LP Ref: Evap to Super', fluid=hp.RP_ref, temperature=hp.Tsat_evap, quality=1.0)
    P_L = hp.evap_sup.sink_in.pressure
    hp.evap_sup.sink_out = MassStream( name='LP Ref: Super to IHX', fluid=hp.RP_ref, pressure=P_L, temperature=(hp.T_source_in - hp.dTmin))
    temp_ref_out = hp.evap_sup.sink_out

    #Internal heat exchanger, if present
    if hp.ihx_on == 1:
        hp.ihx.sink_in = temp_ref_out
        hp.ihx.sink_out = MassStream( name='LP Ref: IHX to PHX1', fluid=hp.RP_ref, pressure=P_L, temperature=(hp.T_gc_out - hp.dTmin))
        temp_ref_out = hp.ihx.sink_out

    #Plate heat exchanger 1 - cold-side, if present
    if hp.phx1_on == 1:
        hp.phx1.sink_in = temp_ref_out
        hp.phx1.sink_out = MassStream( name='LP Ref: PHX1 to Comp1', fluid=hp.RP_ref, pressure=P_L, temperature=(hp.phx1.sink_in.temperature + hp.dT_phx1_cold))
        temp_ref_out = hp.phx1.sink_out

    def _get_discharge_pressure_ratio_estimate(P_H_est: float, P_L: float, n_comp: int, T_tar: float, dT_min: float, s_out: float, RP: PropertyPackage) -> float:
        #Final discharge pressure
        if P_H_est > 0.0:
            P_H = P_H_est
        else:
            P_H = RP.get_properties('P', 'TS', T_tar + dT_min, s_out)[0]
            if P_H < P_L:
                P_H = P_L * 16
            
        return (P_H / P_L) ** (1 / n_comp)

    PR = _get_discharge_pressure_ratio_estimate(P_H_est=hp.P_H_est, P_L=P_L, n_comp=hp.n_comp, T_tar=hp.T_air_tar, dT_min=hp.dT_GC_min_tar, s_out=temp_ref_out.entropy, RP=hp.RP_ref)
    
    #Compressor 1, always present
    hp.comp1.flow_in = temp_ref_out
    hp.comp1.flow_out = MassStream( name='Comp1_out', fluid=hp.RP_ref, pressure=(P_L * PR), entropy=hp.phx1.sink_out.entropy)
    temp_ref_out = hp.comp1.flow_out
    hp.compressors = [hp.comp1]

    #Plate heat exchanger 2 - cold-side
    if hp.phx2_on == 1:
        hp.phx2.sink_in = temp_ref_out
        hp.phx2.sink_out = MassStream( name='Comp2_in', fluid=hp.RP_ref, pressure=hp.phx2.sink_in.pressure, temperature=hp.phx2.sink_in.temperature + hp.dT_phx2_cold)
        temp_ref_out = hp.phx2.sink_out
    
    #Compressor 2
    if hp.comp2_on == 1:
        hp.comp2.flow_in = temp_ref_out
        hp.comp2.flow_out = MassStream( name='Comp2_out', fluid=hp.RP_ref, pressure=(P_L * PR * PR), entropy=hp.phx1.sink_out.entropy)  
        temp_ref_out = hp.comp2.flow_out
        hp.compressors.append(hp.comp2)

    #Plate heat exchanger 3 - cold-side
    if hp.phx3_on == 1:
        hp.phx3.sink_in = temp_ref_out
        hp.phx3.sink_out = MassStream( name='Comp3_in', fluid=hp.RP_ref, pressure=hp.phx3.sink_in.pressure, temperature=hp.phx3.sink_in.temperature + hp.dT_phx3_cold)
        temp_ref_out = hp.phx3.sink_out

    #Compressor 3
    if hp.comp3_on == 1:
        hp.comp3.flow_in = temp_ref_out
        hp.comp3.flow_out = MassStream( name='Ref: Comp3 to GC4', fluid=hp.RP_ref, pressure=(P_L * PR * PR * PR), entropy=hp.phx1.sink_out.entropy)
        temp_ref_out = hp.comp3.flow_out
        hp.compressors.append(hp.comp3)

    hp.comp_discharge = temp_ref_out
    P_H = temp_ref_out.pressure
    hp.gas_coolers = []

    #Gas cooler 4 and Plate heat exchanger 3 - hot-sides
    hp.air_out = MassStream( name='Air: Outlet', fluid=hp.RP_sink, pressure=0.101, temperature=hp.T_air_tar )
    if hp.phx3_on == 1:
        hp.gas_cooler4.sink_out = hp.air_out
        hp.gas_cooler4.sink_pinch = MassStream( name='Air: GC4 pinch point', fluid=hp.RP_sink, pressure=0.101, temperature=hp.T_air_tar )
        hp.gas_cooler4.sink_in = MassStream( name='Air: GC3 to GC4', fluid=hp.RP_sink, pressure=0.101, temperature=hp.T_air_tar )
        temp_air_in = hp.gas_cooler4.sink_in

        hp.gas_coolers.append(hp.gas_cooler4)

        hp.gas_cooler4.source_in = temp_ref_out
        hp.gas_cooler4.source_pinch = MassStream( name='Ref: GC4 pinch point', fluid=hp.RP_ref )
        hp.gas_cooler4.source_out = MassStream( name='Ref: GC4 to PHX3', fluid=hp.RP_ref )
        hp.phx3.source_in = hp.gas_cooler4.source_out
        hp.phx3.source_out = MassStream( name='Ref: PHX3 to GC3', fluid=hp.RP_ref )     
        temp_ref_out = hp.phx3.source_out
    else:
        temp_air_in = hp.air_out
    
    #Gas cooler 3 and Plate heat exchanger 2 - hot-sides
    if hp.phx2_on == 1:
        hp.gas_cooler3.sink_out = temp_air_in
        hp.gas_cooler3.sink_pinch = MassStream( name='Air: GC3 pinch point', fluid=hp.RP_sink )
        hp.gas_cooler3.sink_in = MassStream( name='Air: GC2 to GC3', fluid=hp.RP_sink, pressure=0.101, temperature=hp.T_air_tar )
        temp_air_in = hp.gas_cooler3.sink_in

        hp.gas_coolers.append(hp.gas_cooler3)

        hp.gas_cooler3.source_in = temp_ref_out
        hp.gas_cooler3.source_pinch = MassStream( name='Ref: GC3 pinch point', fluid=hp.RP_ref )
        hp.gas_cooler3.source_out = MassStream( name='Ref: GC3 to PHX2', fluid=hp.RP_ref )
        hp.phx2.source_in = hp.gas_cooler3.source_out
        hp.phx2.source_out = MassStream( name='Ref: PHX2 to GC2', fluid=hp.RP_ref )
        temp_ref_out = hp.phx2.source_out

    #Gas cooler 2 and Plate heat exchanger 1 - hot-sides
    if hp.phx1_on == 1:
        hp.gas_cooler2.sink_out = temp_air_in
        hp.gas_cooler2.sink_pinch = MassStream( name='Air: GC2 pinch point', fluid=hp.RP_sink )
        hp.gas_cooler2.sink_in = MassStream( name='Air: GC1 to GC2', fluid=hp.RP_sink, pressure=0.101, temperature=hp.T_air_tar )
        temp_air_in = hp.gas_cooler2.sink_in

        hp.gas_coolers.append(hp.gas_cooler2)

        hp.gas_cooler2.source_in = temp_ref_out
        hp.gas_cooler2.source_pinch = MassStream( name='Ref: GC2 pinch point', fluid=hp.RP_ref )
        hp.gas_cooler2.source_out = MassStream( name='Ref: GC2 to PHX1', fluid=hp.RP_ref )
        hp.phx1.source_in = hp.gas_cooler2.source_out
        hp.phx1.source_out = MassStream( name='Ref: PHX1 to GC1', fluid=hp.RP_ref )
        temp_ref_out = hp.phx1.source_out

    #Gas cooler 1
    hp.gas_cooler1.sink_out = temp_air_in
    hp.gas_cooler1.sink_pinch = MassStream( name='Air: GC1 pinch point', fluid=hp.RP_sink )
    hp.gas_cooler1.sink_in = MassStream( name='Air: Inlet to GC1', fluid=hp.RP_sink, pressure=0.101, temperature=hp.T_air_inlet )
    
    hp.gas_coolers.append(hp.gas_cooler1)

    hp.gas_cooler1.source_in = temp_ref_out
    hp.gas_cooler1.source_pinch = MassStream( name='Ref: GC1 pinch point', fluid=hp.RP_ref )
    hp.gas_cooler1.source_out = MassStream( name='Ref: GC1 to IHX', fluid=hp.RP_ref, pressure=P_H, temperature=hp.T_gc_out )
    temp_ref_out = hp.gas_cooler1.source_out

    #Internal heat exchanger
    if hp.ihx_on == 1:
        hp.ihx.source_in = temp_ref_out
        hp.ihx.source_out = MassStream( name='Ref: IHX to Exp V', fluid=hp.RP_ref, pressure=P_H, temperature=hp.T_gc_out )
        temp_ref_out = hp.ihx.source_out
    
    #Expansion valve
    hp.exp_valve.flow_in = temp_ref_out
    hp.exp_valve.flow_out = MassStream( name='Evap_in', fluid=hp.RP_ref, pressure=P_L, enthalpy=hp.exp_valve.flow_in.enthalpy )

    #Evaporator
    hp.evaporator.sink_in = hp.exp_valve.flow_out
    hp.evaporator.sink_out = hp.evap_sup.sink_in
    hp.evaporator.source_in = MassStream( name='Exhaust air in', fluid=hp.RP_src, pressure=0.101, temperature=hp.T_source_in, mass_flow=1000)
    hp.evaporator.source_out = MassStream( name='Exhaust air out', fluid=hp.RP_src, pressure=0.101, temperature=min(hp.evaporator.sink_in.temperature + hp.dTmin, hp.T_source_in), mass_flow=1000)

    # hp.streams = [hp.evaporator.sink_out, 
    #                 hp.comp1.flow_out, 
    #                 hp.comp2.flow_out, 
    #                 hp.comp3.flow_out, 
    #                 hp.gas_cooler1.sink_out, 
    #                 hp.gas_cooler2.sink_out, 
    #                 hp.gas_cooler3.sink_out, 
    #                 hp.gas_cooler4.sink_out, 
    #                 hp.ihx.source_out, 
    #                 hp.exp_valve.flow_out]
    return hp

def _get_trhp_logic_blocks(hp: HeatPump) -> HeatPump:
    #Set special relationships between points in the cycle 
    if hp.phx3_on == 1:
        hp.set_phx3_dT_cold = SetRelationship( source_obj=hp.phx3.sink_in, target_obj=[hp.phx3.sink_out], c=hp.dT_phx3_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc4_hot_out = SetRelationship( source_obj=hp.phx3.sink_out, target_obj=[hp.phx3.source_in], c=hp.dTmin, prop='temperature', prop2='pressure' ) 

    if hp.phx2_on == 1:    
        hp.set_phx2_dT_cold = SetRelationship( source_obj=hp.phx2.sink_in, target_obj=[hp.phx2.sink_out], c=hp.dT_phx2_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc3_hot_out = SetRelationship( source_obj=hp.phx2.sink_out, target_obj=[hp.phx2.source_in], c=hp.dTmin, prop='temperature', prop2='pressure' ) 
    
    if hp.phx1_on == 1:
        hp.set_phx1_dT_cold = SetRelationship( source_obj=hp.phx1.sink_in, target_obj=[hp.phx1.sink_out], c=hp.dT_phx1_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc2_hot_out = SetRelationship( source_obj=hp.phx1.sink_out, target_obj=[hp.phx1.source_in], c=hp.dTmin, prop='temperature', prop2='pressure' ) 

    if hp.n_comp == 2:
        hp.set_P_ints = CustomRelationship( source_obj=[hp.comp1.flow_in, hp.comp_discharge, hp.comp2_power_frac, hp.comp3_power_frac], target_obj=[hp.comp1.flow_out, hp.comp_discharge], relation=1, tol=hp.tol )
    elif hp.n_comp == 3:
        hp.set_P_ints = CustomRelationship( source_obj=[hp.comp1.flow_in, hp.comp_discharge, hp.comp2_power_frac, hp.comp3_power_frac], target_obj=[hp.comp1.flow_out, hp.comp2.flow_out], relation=1, tol=hp.tol )

    # hp.gas_coolers.reverse()
    # hp.set_gc_air_massflow = CustomRelationship( source_obj=hp.gas_coolers, target_obj=hp.gas_coolers, relation=2, tol=hp.tol )
    # hp.set_evap_air_massflow = CustomRelationship( source_obj=[hp.evaporator, hp.evap_sup], target_obj=[hp.evaporator, hp.evap_sup], relation=3)
    return hp

def _get_trhp_solver_pathway(hp: HeatPump) -> HeatPump:
    #Solution pathway set-up
    hp.operations = [hp.evap_sup]

    if hp.n_comp > 1:
        hp.operations.append(hp.set_P_ints)
    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)
    if hp.phx1_on == 1:
        hp.operations.append(hp.set_phx1_dT_cold)
        # hp.operations.append(hp.set_gc2_hot_out)

    hp.operations.append(hp.comp1)

    if hp.phx2_on == 1:
        hp.operations.append(hp.set_phx2_dT_cold)
        # hp.operations.append(hp.set_gc3_hot_out)
    if hp.comp2_on == 1:
        hp.operations.append(hp.comp2)
    if hp.phx3_on == 1:
        hp.operations.append(hp.set_phx3_dT_cold)
        # hp.operations.append(hp.set_gc4_hot_out)
    if hp.comp3_on == 1:
        hp.operations.append(hp.comp3)
    
    if hp.phx3_on == 1:
        hp.operations.append(hp.gas_cooler4)
        hp.operations.append(hp.phx3)
    if hp.phx2_on == 1:
        hp.operations.append(hp.gas_cooler3)
        hp.operations.append(hp.phx2)
    if hp.phx1_on == 1:
        hp.operations.append(hp.gas_cooler2)        
        hp.operations.append(hp.phx1)

    hp.operations.append(hp.gas_cooler1)

    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)

    hp.operations.append(hp.exp_valve)
    hp.operations.append(hp.evaporator)

    return hp

def _get_trhp_logic_blocks_v1(hp: HeatPump) -> HeatPump:
    #Set special relationships between points in the cycle 
    if hp.phx3_on == 1:
        hp.set_phx3_dT_cold = SetRelationship( source_obj=hp.phx3.sink_in, target_obj=[hp.phx3.sink_out], c=hp.dT_phx3_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc4_hot_out = SetRelationship( source_obj=hp.phx3.sink_out, target_obj=[hp.phx3.source_in], c=hp.dTmin, prop='temperature', prop2='pressure' ) 

    if hp.phx2_on == 1:    
        hp.set_phx2_dT_cold = SetRelationship( source_obj=hp.phx2.sink_in, target_obj=[hp.phx2.sink_out], c=hp.dT_phx2_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc3_hot_out = SetRelationship( source_obj=hp.phx2.sink_out, target_obj=[hp.phx2.source_in], c=hp.dTmin, prop='temperature', prop2='pressure' ) 
    
    if hp.phx1_on == 1:
        hp.set_phx1_dT_cold = SetRelationship( source_obj=hp.phx1.sink_in, target_obj=[hp.phx1.sink_out], c=hp.dT_phx1_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc2_hot_out = SetRelationship( source_obj=hp.phx1.sink_out, target_obj=[hp.phx1.source_in], c=hp.dTmin, prop='temperature', prop2='pressure' ) 

    if hp.n_comp == 2:
        hp.set_P_ints = CustomRelationship( source_obj=[hp.comp1.flow_in, hp.comp_discharge, hp.comp2_power_frac, hp.comp3_power_frac], target_obj=[hp.comp1.flow_out, hp.comp_discharge], relation=1, tol=hp.tol )
    elif hp.n_comp == 3:
        hp.set_P_ints = CustomRelationship( source_obj=[hp.comp1.flow_in, hp.comp_discharge, hp.comp2_power_frac, hp.comp3_power_frac], target_obj=[hp.comp1.flow_out, hp.comp2.flow_out], relation=1, tol=hp.tol )

    hp.gas_coolers.reverse()
    hp.set_gc_air_massflow = CustomRelationship( source_obj=hp.gas_coolers, target_obj=hp.gas_coolers, relation=2, tol=hp.tol )
    #hp.set_evap_air_massflow = CustomRelationship( source_obj=[hp.evaporator, hp.evap_sup], target_obj=[hp.evaporator, hp.evap_sup], relation=3)
    return hp

def _get_trhp_solver_pathway_v1(hp: HeatPump) -> HeatPump:
    #Solution pathway set-up
    hp.operations = [hp.evap_sup]

    if hp.n_comp > 1:
        hp.operations.append(hp.set_P_ints)
    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)
    if hp.phx1_on == 1:
        hp.operations.append(hp.set_phx1_dT_cold)
        hp.operations.append(hp.set_gc2_hot_out)

    hp.operations.append(hp.comp1)

    if hp.phx2_on == 1:
        hp.operations.append(hp.set_phx2_dT_cold)
        hp.operations.append(hp.set_gc3_hot_out)
    if hp.comp2_on == 1:
        hp.operations.append(hp.comp2)
    if hp.phx3_on == 1:
        hp.operations.append(hp.set_phx3_dT_cold)
        hp.operations.append(hp.set_gc4_hot_out)
    if hp.comp3_on == 1:
        hp.operations.append(hp.comp3)
    
    if hp.phx3_on == 1:
        hp.operations.append(hp.gas_cooler4)
        hp.operations.append(hp.phx3)
    if hp.phx2_on == 1:
        hp.operations.append(hp.gas_cooler3)
        hp.operations.append(hp.phx2)
    if hp.phx1_on == 1:
        hp.operations.append(hp.gas_cooler2)        
        hp.operations.append(hp.phx1)

    hp.operations.append(hp.gas_cooler1)

    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)

    hp.operations.append(hp.exp_valve)
    hp.operations.append(hp.evaporator)

    return hp
 
def set_quality_relationship_v1(hp: HeatPump) -> HeatPump:
    e = []
    gcs = [hp.gas_cooler1, hp.gas_cooler2, hp.gas_cooler3, hp.gas_cooler4]
    for g in gcs:
        T_sat = hp.RP_ref.get_properties('T', 'PQ', g.source_out.pressure, 0.0)[0]
        e.append(abs(T_sat - g.source_out.temperature))
    
    hp.q_gc_out = [1.0]

    n = e.index(min(e))
    new_rel = SetRelationship( source_obj=hp.q_gc_out, target_obj=[gcs[n].source_out], prop='quality', prop2='pressure' )        
    if n == 1:
        old_rel = hp.set_gc2_hot_out
        hp.set_gc2_hot_out = new_rel
        hp.phx_fixed_src_quality = hp.phx1
    elif n == 2:
        old_rel = hp.set_gc3_hot_out
        hp.set_gc3_hot_out = new_rel
        hp.phx_fixed_src_quality = hp.phx2
    elif n == 3:
        old_rel = hp.set_gc4_hot_out
        hp.set_gc4_hot_out = new_rel
        hp.phx_fixed_src_quality = hp.phx3
    else:
        ValueError('Error in the model.')

    for i in range(len(hp.operations)):
        if hp.operations[i] == old_rel:
            hp.operations[i] = new_rel

    gcs[n].model = 'fixed hot quality outlet'
    return gcs[n]