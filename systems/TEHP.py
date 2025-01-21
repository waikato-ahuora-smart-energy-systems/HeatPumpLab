__author__ = 'Tim Walmsley'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2025'

"""
High temperature transcritical-economised heat pump model
By Tim Walmsley, University of Waikato, 2025 (tim.walmsley@waikato.ac.nz)
"""

from systems.common.commonHeatPump import HeatPump
from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator, SetRelationship, CustomRelationship, Economiser, Mixer, Splitter
from streams import MassStream

class TEHP(HeatPump):
    dTmin_econ: float =5.0
    n_comp: float =2
    y_econ: float =0.0
    comp_PR_uniform: bool =False
    comp1_power_frac: float =0.5
    comp2_power_frac: float =0.5

    def __init__(self, **d_args):
        for key, value in d_args.items():
            setattr(self, key, value)
   
    def get_case_parameters(self) -> dict[str, float | bool]:
        return {
            'refrigerant' : self.refrigerant,
            'dT_GC_min_tar' : self.dT_GC_min_tar,
            'T_air_tar' : self.T_air_tar,
            'T_air_inlet' : self.T_air_inlet,
            'T_source_in' : self.T_source_in,
            'comp_PR_uniform' : self.comp_PR_uniform,            
            'comp2_power_frac' : self.comp2_power_frac,
            }
    
    def get_summary(self) -> dict[str, float | bool]:
        self.output = {
            #Headline results
            'cycle_type' : self.cycle_type,
            'COP' : self.COP,
            'SEC' : self.SEC,
            'VHC' : self.VHC,
            'Pmax' : self.comp_discharge.pressure * 1000,
            'converged' : self.converged,
            #Parameters
            'refrigerant' : self.refrigerant,
            'n_comp' : self.n_comp,
            'T_source_in' : self.T_source_in,
            'dT_GC_min_tar' : self.dT_GC_min_tar,
            'T_air_inlet' : self.T_air_inlet,
            'T_air_tar' : self.T_air_tar,
            'comp_PR_uniform' : self.comp_PR_uniform,
            'Tref_evap_inlet' : self.Tref_evap_inlet,
            #Optmisation variables
            'y_econ' : self.y_econ,
            'comp1_power_frac' : self.comp1_power_frac, 
            'comp2_power_frac' : self.comp2_power_frac,
            'dT_gc_out_min_offset' : self.dT_gc_out_min_offset,
            #Additional results
            'T_comp_discharge' : self.comp_discharge.temperature,
            'T_gc_out' : self.T_gc_out,
            'Tref_evap_outlet' : self.Tref_evap_outlet,
            'P_evap_sink_in' : self.evaporator.sink_in.pressure * 1000,
            'Comp_1_P_out' : self.comp1.flow_out.pressure * 1000,
            'Comp_2_P_out' : self.comp2.flow_out.pressure * 1000,
            'PR_1' : self.comp1.pressure_ratio,
            'PR_2' : self.comp2.pressure_ratio,
            'Overall_PR' : self.comp_discharge.pressure / self.evaporator.sink_in.pressure,
            }
        return self.output

    def get_components_sets_to_graph(self) -> list:
        cycle_ls: list[list[str, str]]
        sink_ls: list[list[str, str]]
        
        cycle_ls = [
            ['self.ihx', 'source'],
            ['self.exp_valve', 'siso'],
            ['self.evaporator', 'sink'],
            ['self.evap_sup', 'sink'],
            ['self.ihx', 'sink'],
            ['self.comp1', 'siso'],
            ['self.exp_valve_sm', 'siso'],
            ['self.econ', 'sink'],
            ['self.comp2', 'siso'],
            ['self.gas_cooler2', 'source'],
            ['self.econ', 'source'],
            ['self.gas_cooler1', 'source'],       
            ]
        
        sink_ls = [
            ['self.gas_cooler2', 'sink'],
            ['self.gas_cooler1', 'sink'],
            ]

        return [
            self._get_component_objs_for_graphing(ls) for ls in [cycle_ls, sink_ls]
            ] 
    
    def _get_component_objs_for_graphing(self, ls: list) -> list:
        obj_ls: list = []
        obj_id: list[str, str]
        obj_name: str
        att_name: str
        component: object
        
        for obj_id in ls:
            obj_name, att_name = obj_id[0].split('.')
            if hasattr(locals()[obj_name], att_name):
                component = getattr(locals()[obj_name], att_name)
                if component.converged:
                    obj_ls.append([component, obj_id[1]])     
        return obj_ls 

def set_tehp_model(hp: TEHP) -> TEHP:
    hp = _get_tehp_structure(hp)
    if _check_input_param_feasibility(hp) == True:
        hp.feasible = True
        hp = _get_tehp_components(hp)
        hp = _get_tehp_streams(hp)
        hp = _get_tehp_logic_blocks(hp)
        hp = _get_tehp_solver_pathway(hp)
    else:
        hp.feasible = False
    return hp

def _check_input_param_feasibility(hp: TEHP):
    return True if (
        0.0 <= hp.dT_gc_out_min_offset <= 50.0 and
        hp.comp2_power_frac <= 1.0
        ) else False

def _get_tehp_structure(hp: TEHP) -> TEHP:
    #Check the y_econ value is not negative
    hp.y_econ = hp.y_econ if hp.y_econ > 0.0 else 0.0
    
    #Check the compressor combo
    if hp.n_comp != 2:
        hp.n_comp = 2
    
    if hp.comp_PR_uniform:
        hp.comp2_power_frac = 0.5
    else:
        hp.comp2_power_frac = hp.comp2_power_frac if hp.comp2_power_frac > 0.0 else 0.0

    hp.comp1_power_frac = 1 - hp.comp2_power_frac

    #Binary variables for the heat pump structure
    hp.comp2_on = 1
    hp.T_gc_out = hp.T_air_inlet + hp.dT_GC_min_tar + hp.dT_gc_out_min_offset
    hp.ihx_on = 1    if hp.T_source_in < hp.T_gc_out and hp.ihx_on == 1 else 0

    return hp

def _get_tehp_components(hp: TEHP) -> TEHP:
    #Initialisation of unit operations
    hp.comp1 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.motor_eff, model_index=hp.comp_model_index)
    hp.comp2 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.motor_eff, model_index=hp.comp_model_index)
  
    hp.gas_cooler1 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
    hp.gas_cooler2 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, model='fixed hot duty', dTmin=hp.dT_GC_min_tar)

    hp.exp_valve_sm = ExpansionValve(name='Expansion valve', fluid=hp.refrigerant)
    hp.econ = Economiser(name='Economiser', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dTmin=hp.dTmin_econ, y_econ=hp.y_econ)
    hp.mixer = Mixer(name='Mixer', fluid=hp.refrigerant)
    hp.splitter = Splitter(name='Splitter', fluid=hp.refrigerant, r_ratio=hp.y_econ)

    hp.ihx = PlateHeatExchanger(name='Internal heat exchanger', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dTmin=hp.dTmin_econ, model='dTmin_terminal')

    hp.exp_valve = ExpansionValve(name='Expansion valve', fluid=hp.refrigerant)
    hp.evaporator = Evaporator(name='Evaporator', source_fluid=hp.source, sink_fluid=hp.refrigerant, model='saturated gas outlet evaporator')
    hp.evap_sup = PlateHeatExchanger(name='Evaporator', source_fluid=hp.source, sink_fluid=hp.refrigerant, model='fixed cold temperature outlet')
    
    return hp

def _get_tehp_streams(hp: TEHP) -> TEHP:
    temp_ref_out: MassStream
    P_H: float
    P_L: float

    #Evaporator superheat section
    hp.evap_sup.source_in.props_in_ab(input_str='PT', prop1=0.5, prop2=hp.T_source_in, mass_flow=1e6, name='Source in', pp_instance='src')
    hp.evap_sup.source_out.props_in_ab(input_str='PT', prop1=0.5, prop2=hp.T_source_in, mass_flow=1e6, name='Source out', pp_instance='src')
    hp.Tref_evap_outlet = hp.Tref_evap_inlet
    hp.evap_sup.sink_in.props_in_ab(input_str='TQ', prop1=hp.Tref_evap_outlet, prop2=1.0, name='Evaporator sink outlet')
    P_L = hp.evap_sup.sink_in.pressure
    hp.evap_sup.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.T_source_in - hp.dTmin_econ), name='Superheater sink outlet')
    temp_ref_out = hp.evap_sup.sink_out

    #Internal heat exchanger, if present
    if hp.ihx_on == 1:
        hp.ihx.sink_in = temp_ref_out
        hp.ihx.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.T_gc_out - hp.dTmin_econ), name='IHX sink outlet')
        temp_ref_out = hp.ihx.sink_out

    def _get_discharge_pressure_ratio_estimate(P_H_est: float, P_L: float, n_comp: int, T_tar: float, dT_min: float, s_out: float, fld: str) -> float:
        #Final discharge pressure
        if P_H_est > 0.0:
            P_H = P_H_est
        else:
            P_H = MassStream(fluid=fld, temperature=T_tar * 1.5, entropy=s_out).pressure
            if P_H < P_L:
                P_H = P_L * 16
            
        return (P_H / P_L) ** (1 / n_comp)

    PR = _get_discharge_pressure_ratio_estimate(P_H_est=hp.P_H_est, P_L=P_L, n_comp=hp.n_comp, T_tar=hp.T_air_tar, dT_min=hp.dT_GC_min_tar, s_out=temp_ref_out.entropy, fld=hp.refrigerant)
    
    #Compressor 1, always present
    hp.comp1.flow_in = temp_ref_out
    hp.comp1.flow_out.props_in_ab(input_str='PS', prop1=(P_L * PR), prop2=hp.comp1.flow_in.entropy, name='Compressor 1 outlet') 

    hp.mixer.flow_in_1 = hp.comp1.flow_out
    hp.mixer.flow_out.props_in_ab(input_str='PH', prop1=hp.comp1.flow_out.pressure, prop2=hp.comp1.flow_out.enthalpy, name='Mixer outlet') 
    
    hp.comp2.flow_in = hp.mixer.flow_out
    hp.comp2.flow_out.props_in_ab(input_str='PS', prop1=(P_L * PR * PR), prop2=hp.comp2.flow_in.entropy, name='Compressor 2 outlet') 
    hp.comp_discharge = hp.comp2.flow_out

    hp.compressors = [hp.comp1, hp.comp2]

    hp.comp_discharge.pressure_bnds = None
    P_H = hp.comp_discharge.pressure
    hp.gas_coolers = [hp.gas_cooler2, hp.gas_cooler1]

    #Gas cooler 2
        #Air side
    hp.air_out = MassStream( name='Air sink outlet', fluid=hp.sink, pp_instance='snk', pressure=0.101, temperature=hp.T_air_tar )
    hp.gas_cooler2.sink_out = hp.air_out
    hp.gas_cooler2.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC2 sink pinch point', pp_instance='snk') 
    hp.gas_cooler2.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC2 sink inlet', pp_instance='snk') 
        #Refrigerant side
    hp.gas_cooler2.source_in = hp.comp_discharge
    hp.gas_cooler2.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC2 source pinch point') 
    hp.gas_cooler2.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC2 source outlet')

    #Economiser and Splitter
    hp.econ.source_in = hp.gas_cooler2.source_out
    hp.econ.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='Econ source outlet')

    hp.splitter.flow_in = hp.econ.source_out
    hp.splitter.flow_out_1.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC1 source inlet')
    hp.splitter.flow_out_2.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='Econ sink inlet')

    hp.exp_valve_sm.flow_in = hp.splitter.flow_out_2
    hp.exp_valve_sm.flow_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='Econ sink outlet')

    hp.econ.sink_in = hp.exp_valve_sm.flow_out
    hp.econ.sink_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='Econ sink outlet')
    
    hp.mixer.flow_in_2 = hp.econ.sink_out

    #Gas cooler 1
        #Air side
    hp.gas_cooler1.sink_out = hp.gas_cooler2.sink_in
    hp.gas_cooler1.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_inlet, name='GC1 sink pinch point', pp_instance='snk') 
    hp.gas_cooler1.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_inlet, name='GC1 sink inlet', pp_instance='snk')
        #Refrigerant side
    hp.gas_cooler1.source_in = hp.splitter.flow_out_1
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

    hp.evaporator.source_in.props_in_ab(input_str='PT', prop1=0.5, prop2=hp.T_source_in, mass_flow=1e6, name='Evaporator source inlet', pp_instance='src')
    hp.evaporator.source_out.props_in_ab(input_str='PT', prop1=0.5, prop2=min(hp.evaporator.sink_in.temperature + hp.dTmin_econ, hp.T_source_in), mass_flow=1e6, name='Evaporator source outlet', pp_instance='src')

    return hp

def _get_tehp_logic_blocks(hp: TEHP) -> TEHP:
    #Set special relationships between points in the cycle 
    hp.set_P_ints = CustomRelationship( source_obj=[hp.comp1.flow_in, hp.comp_discharge, hp.comp2_power_frac, 0.0], target_obj=[hp.comp1.flow_out, hp.comp_discharge], relation=1, tol=hp.tol )
    hp.econ.set_add_obj(comp1_flow_out=hp.comp1.flow_out, comp2_flow_out=hp.comp2.flow_out, gc_2_dP_h=hp.gas_cooler2.dP_h, T_gc_out=hp.T_gc_out)
    hp.exp_valve_sm.dTmin = hp.dTmin_ev_sm
    hp.set_exp_sm_T_out = CustomRelationship( source_obj=[hp.comp1, hp.econ, hp.comp_discharge, hp.gas_cooler2.dP_h], target_obj=[hp.exp_valve_sm, hp.splitter], relation=4)
    return hp

def _get_tehp_solver_pathway(hp: TEHP) -> TEHP:
    #Solution pathway set-up
    hp.operations = [hp.evap_sup]
    hp.operations.append(hp.set_P_ints)

    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)

    hp.operations.append(hp.comp1)
    hp.operations.append(hp.set_exp_sm_T_out)
    hp.operations.append(hp.exp_valve_sm)
    hp.operations.append(hp.econ)
    hp.operations.append(hp.mixer)
    hp.operations.append(hp.comp2)
    hp.operations.append(hp.gas_cooler2)        
    hp.operations.append(hp.splitter)
    hp.operations.append(hp.gas_cooler1)

    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)

    hp.operations.append(hp.exp_valve)
    hp.operations.append(hp.evaporator)

    return hp
    