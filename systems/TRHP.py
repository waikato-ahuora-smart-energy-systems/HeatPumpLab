__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

from systems.heat_pump import HeatPump
from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator, SetRelationship, CustomRelationship
from streams import MassStream

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
    hp.comp1 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.eff_motor, model_index=1)
    if hp.comp2_on == 1:
        hp.comp2 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.eff_motor, model_index=1)
        hp.gas_cooler3 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
        if hp.comp3_on == 1:
            hp.comp3 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.eff_motor, model_index=1) 
            hp.gas_cooler4 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
    
    hp.gas_cooler1 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
    hp.gas_cooler2 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
    
    hp.ihx = PlateHeatExchanger(name='Internal heat exchanger', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dTmin=hp.dTmin, model='dTmin_terminal')
    if hp.phx1_on == 1:
        hp.phx1 = PlateHeatExchanger(name='Recuperator', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dT_c=hp.dT_phx1_cold, model='fixed cold delta temperature')
    if hp.phx2_on == 1:
        hp.phx2 = PlateHeatExchanger(name='Recuperator', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dT_c=hp.dT_phx2_cold, model='fixed cold delta temperature') 
    if hp.phx3_on == 1:       
        hp.phx3 = PlateHeatExchanger(name='Recuperator', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dT_c=hp.dT_phx3_cold, model='fixed cold delta temperature')
    
    hp.exp_valve = ExpansionValve(name='Expansion valve', fluid=hp.refrigerant)
    hp.evaporator = Evaporator(name='Evaporator', source_fluid=hp.source, sink_fluid=hp.refrigerant, model='saturated gas outlet evaporator')
    hp.evap_sup = PlateHeatExchanger(name='Evaporator', source_fluid=hp.source, sink_fluid=hp.refrigerant, model='fixed cold temperature outlet')
    
    return hp

def _get_trhp_streams(hp: HeatPump) -> HeatPump:
    temp_ref_out: MassStream
    P_H: float
    P_L: float

    #Evaporator superheat section
    hp.evap_sup.source_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_source_in, mass_flow=1e6, name='Source in')
    hp.evap_sup.source_out.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_source_in, mass_flow=1e6, name='Source out')
    hp.evap_sup.sink_in.props_in_ab(input_str='TQ', prop1=hp.Tsat_evap, prop2=1.0, name='Evaporator sink outlet')
    P_L = hp.evap_sup.sink_in.pressure
    hp.evap_sup.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.T_source_in - hp.dTmin), name='Superheater sink outlet')
    temp_ref_out = hp.evap_sup.sink_out

    #Internal heat exchanger, if present
    if hp.ihx_on == 1:
        hp.ihx.sink_in = temp_ref_out
        hp.ihx.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.T_gc_out - hp.dTmin), name='IHX sink outlet')
        temp_ref_out = hp.ihx.sink_out

    #Plate heat exchanger 1 - cold-side, if present
    if hp.phx1_on == 1:
        hp.phx1.sink_in = temp_ref_out
        hp.phx1.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.phx1.sink_in.temperature + hp.dT_phx1_cold), name='PHX1 sink outlet')
        temp_ref_out = hp.phx1.sink_out

    def _get_discharge_pressure_ratio_estimate(P_H_est: float, P_L: float, n_comp: int, T_tar: float, dT_min: float, s_out: float, fld: str) -> float:
        #Final discharge pressure
        if P_H_est > 0.0:
            P_H = P_H_est
        else:
            P_H = MassStream(fluid=fld, temperature=T_tar + dT_min, entropy=s_out).pressure
            if P_H < P_L:
                P_H = P_L * 16
            
        return (P_H / P_L) ** (1 / n_comp)

    PR = _get_discharge_pressure_ratio_estimate(P_H_est=hp.P_H_est, P_L=P_L, n_comp=hp.n_comp, T_tar=hp.T_air_tar, dT_min=hp.dT_GC_min_tar, s_out=temp_ref_out.entropy, fld=hp.refrigerant)
    
    #Compressor 1, always present
    hp.comp1.flow_in = temp_ref_out
    hp.comp1.flow_out.props_in_ab(input_str='PS', prop1=(P_L * PR), prop2=hp.comp1.flow_in.entropy, name='Compressor 1 outlet') 
    temp_ref_out = hp.comp1.flow_out
    hp.compressors = [hp.comp1]

    #Plate heat exchanger 2 - cold-side
    if hp.phx2_on == 1:
        hp.phx2.sink_in = temp_ref_out
        hp.phx2.sink_out.props_in_ab(input_str='PT', prop1=hp.phx2.sink_in.pressure, prop2=(hp.phx2.sink_in.temperature + hp.dT_phx2_cold), name='PHX2 sink outlet') 
        temp_ref_out = hp.phx2.sink_out
    
    #Compressor 2
    if hp.comp2_on == 1:
        hp.comp2.flow_in = temp_ref_out
        hp.comp2.flow_out.props_in_ab(input_str='PS', prop1=(P_L * PR * PR), prop2=hp.comp2.flow_in.entropy, name='Compressor 2 outlet') 
        temp_ref_out = hp.comp2.flow_out
        hp.compressors.append(hp.comp2)

    #Plate heat exchanger 3 - cold-side
    if hp.phx3_on == 1:
        hp.phx3.sink_in = temp_ref_out
        hp.phx3.sink_out.props_in_ab(input_str='PT', prop1=hp.phx3.sink_in.pressure, prop2=(hp.phx3.sink_in.temperature + hp.dT_phx3_cold), name='PHX3 sink outlet') 
        temp_ref_out = hp.phx3.sink_out

    #Compressor 3
    if hp.comp3_on == 1:
        hp.comp3.flow_in = temp_ref_out
        hp.comp3.flow_out.props_in_ab(input_str='PS', prop1=(P_L * PR * PR * PR), prop2=hp.comp3.flow_in.entropy, name='Compressor 3 outlet') 
        temp_ref_out = hp.comp3.flow_out
        hp.compressors.append(hp.comp3)

    hp.comp_discharge = temp_ref_out
    P_H = temp_ref_out.pressure
    hp.gas_coolers = []

    #Gas cooler 4 and Plate heat exchanger 3 - hot-sides
    hp.air_out = MassStream( name='Air sink outlet', fluid=hp.sink, pressure=0.101, temperature=hp.T_air_tar )
    if hp.phx3_on == 1:
        hp.gas_cooler4.sink_out = hp.air_out
        hp.gas_cooler4.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC4 sink pinch point') 
        hp.gas_cooler4.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC4 sink inlet') 
        temp_air_in = hp.gas_cooler4.sink_in

        hp.gas_coolers.append(hp.gas_cooler4)

        hp.gas_cooler4.source_in = temp_ref_out
        hp.gas_cooler4.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC4 source pinch point') 
        hp.gas_cooler4.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC4 source outlet')

        hp.phx3.source_in = hp.gas_cooler4.source_out
        hp.phx3.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='PHX3 source outlet')
        temp_ref_out = hp.phx3.source_out
    else:
        temp_air_in = hp.air_out
    
    #Gas cooler 3 and Plate heat exchanger 2 - hot-sides
    if hp.phx2_on == 1:
        hp.gas_cooler3.sink_out = temp_air_in
        hp.gas_cooler3.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC3 sink pinch point') 
        hp.gas_cooler3.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC3 sink inlet') 
        temp_air_in = hp.gas_cooler3.sink_in

        hp.gas_coolers.append(hp.gas_cooler3)

        hp.gas_cooler3.source_in = temp_ref_out
        hp.gas_cooler3.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC3 source pinch point') 
        hp.gas_cooler3.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC3 source outlet')

        hp.phx2.source_in = hp.gas_cooler3.source_out
        hp.phx2.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='PHX2 source outlet')
        temp_ref_out = hp.phx2.source_out

    #Gas cooler 2 and Plate heat exchanger 1 - hot-sides
    if hp.phx1_on == 1:
        hp.gas_cooler2.sink_out = temp_air_in
        hp.gas_cooler2.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC2 sink pinch point') 
        hp.gas_cooler2.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC2 sink inlet') 
        temp_air_in = hp.gas_cooler2.sink_in

        hp.gas_coolers.append(hp.gas_cooler2)

        hp.gas_cooler2.source_in = temp_ref_out
        hp.gas_cooler2.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC2 source pinch point') 
        hp.gas_cooler2.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC2 source outlet')

        hp.phx1.source_in = hp.gas_cooler2.source_out
        hp.phx1.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='PHX1 source outlet')
        temp_ref_out = hp.phx1.source_out

    #Gas cooler 1
    hp.gas_cooler1.sink_out = temp_air_in
    hp.gas_cooler1.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_inlet, name='GC1 sink pinch point') 
    hp.gas_cooler1.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_inlet, name='GC1 sink inlet')
    
    hp.gas_coolers.append(hp.gas_cooler1)

    hp.gas_cooler1.source_in = temp_ref_out
    hp.gas_cooler1.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.T_gc_out, name='GC1 source pinch point') 
    hp.gas_cooler1.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.T_gc_out, name='GC1 source outlet')

    temp_ref_out = hp.gas_cooler1.source_out

    #Internal heat exchanger
    if hp.ihx_on == 1:
        hp.ihx.source_in = temp_ref_out
        hp.ihx.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.T_gc_out, name='IHX source outlet')
        temp_ref_out = hp.ihx.source_out
    
    #Expansion valve
    hp.exp_valve.flow_in = temp_ref_out
    hp.exp_valve.flow_out.props_in_ab(input_str='PH', prop1=P_L, prop2=hp.exp_valve.flow_in.enthalpy, name='Expansion valve outlet')

    #Evaporator
    hp.evaporator.sink_in = hp.exp_valve.flow_out
    hp.evaporator.sink_out = hp.evap_sup.sink_in

    hp.evaporator.source_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_source_in, mass_flow=1e6, name='Evaporator source inlet')
    hp.evaporator.source_out.props_in_ab(input_str='PT', prop1=0.101, prop2=min(hp.evaporator.sink_in.temperature + hp.dTmin, hp.T_source_in), mass_flow=1e6, name='Evaporator source outlet')

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
        T_sat = MassStream(fluid=hp.refrigerant, pressure=g.source_out.pressure, quality=0.0).temperature
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