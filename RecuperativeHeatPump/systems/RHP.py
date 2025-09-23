
from .common.heatPump import BaseHeatPump
from ..components import *
from ..streams import MassStream

class RHP(BaseHeatPump):
    
    n_comp: float =3
    comp_PR_uniform: float =False 
    comp1_power_frac: float =1.0
    comp2_power_frac: float =0.0 
    comp3_power_frac: float =0.0

    dT_rhx1_cold: float =0.0
    dT_rhx2_cold: float =0.0
    dT_rhx3_cold: float =0.0

    comp1: Compressor
    comp2: Compressor
    comp3: Compressor

    rhx1: PlateHeatExchanger
    rhx2: PlateHeatExchanger
    rhx3: PlateHeatExchanger

    gas_cooler1: GasCooler
    gas_cooler2: GasCooler
    gas_cooler3: GasCooler
    gas_cooler4: GasCooler

    ihx: PlateHeatExchanger
    exp_valve: ExpansionValve
    evaporator: Evaporator
    evap_sup: PlateHeatExchanger

    set_P_ints: CustomRelationship
    set_rhx1_dT_cold: SetRelationship
    set_gc2_hot_out: SetRelationship
    set_rhx2_dT_cold: SetRelationship
    set_gc3_hot_out: SetRelationship
    set_rhx3_dT_cold: SetRelationship
    set_gc4_hot_out: SetRelationship

    comp_discharge: MassStream

    def __init__(self, **d_args):
        self.args = d_args
        for key, value in d_args.items():
            setattr(self, key, value)
   
    def get_case_parameters(self) -> dict[str, float | bool]:
        return {
            'refrigerant' : self.refrigerant,
            'dT_GC_min_tar' : self.dT_GC_min_tar,
            'T_air_tar' : self.T_air_tar,
            'T_air_inlet' : self.T_air_inlet,
            'T_source_in' : self.T_source_in,
            'n_comp' : self.n_comp,
            'dT_rhx1_cold' : self.dT_rhx1_cold,
            'dT_rhx2_cold' : self.dT_rhx2_cold,
            'dT_rhx3_cold' : self.dT_rhx3_cold,
            'comp2_power_frac' : self.comp2_power_frac,
            'comp3_power_frac' : self.comp3_power_frac,
            'comp_PR_uniform' : self.comp_PR_uniform
            }
    
    def get_summary(self, dt=0) -> dict[str, float | bool]:
        self.output = {
            #Headline results
            'cycle_type' : self.cycle_type,
            'COP' : self.COP,
            'SEC' : self.SEC,
            'VHC' : self.VHC,
            'eDTL' : self.eDTL,
            'SLE' : self.SLE,
            'dQ/dW': self.dQdW,
            'LMTD_gcs' : self.LMTD_gcs,
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
            'dTmin_rhx' : self.dTmin_rhx,
            'ihx_on': self.ihx_on,
            'comp_model_index': self.comp_model_index,
            'motor_eff': self.motor_eff,
            'sink': self.sink,
            'source': self.source,
            #Optmisation variables
            'dT_rhx1_cold' : self.dT_rhx1_cold,
            'dT_rhx2_cold' : self.dT_rhx2_cold,
            'dT_rhx3_cold' : self.dT_rhx3_cold,
            'dT_rhx_cold_sum' : self.dT_rhx1_cold + self.dT_rhx2_cold + self.dT_rhx3_cold,
            'comp1_power_frac' : self.comp1_power_frac, 
            'comp2_power_frac' : self.comp2_power_frac, 
            'comp3_power_frac' : self.comp3_power_frac,  
            'dT_gc_out_min_offset' : self.dT_gc_out_min_offset,
            #Additional results
            'T_comp_discharge' : self.comp_discharge.temperature,
            'T_gc_out' : self.T_gc_out,
            'Tref_evap_outlet' : self.Tref_evap_outlet,
            'P_evap_sink_in' : self.evaporator.sink_in.pressure * 1000,
            'Comp_1_P_out' : self.comp1.flow_out.pressure * 1000,
            'Comp_2_P_out' : self.comp2.flow_out.pressure * 1000 if self.n_comp > 1 else self.comp_discharge.pressure * 1000,
            'Comp_3_P_out' : self.comp_discharge.pressure * 1000,
            'PR_1' : self.comp1.pressure_ratio,
            'PR_2' : self.comp2.pressure_ratio if self.n_comp > 1 else 1.0,
            'PR_3' : self.comp3.pressure_ratio if self.n_comp > 2 else 1.0,
            'Overall_PR' : self.comp_discharge.pressure / self.evaporator.sink_in.pressure,
            'Comp_1_T_out' : self.comp1.flow_out.temperature,
            'Comp_2_T_out' : self.comp2.flow_out.temperature if self.n_comp > 1 else self.comp_discharge.temperature,
            'Comp_3_T_out' : self.comp_discharge.temperature,
            'n_is_1' : self.comp1.isentropic_eff,
            'n_is_2' : self.comp2.isentropic_eff if self.n_comp > 1 else 1.0,
            'n_is_3' : self.comp3.isentropic_eff if self.n_comp > 2 else 1.0,
            'n_vol_1': self.comp1.isentropic_eff,
            'n_vol_2': self.comp2.isentropic_eff if self.n_comp > 1 else 1.0,
            'n_vol_3': self.comp3.isentropic_eff if self.n_comp > 2 else 1.0,
            'Solve time' : dt,
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
            ['self.rhx1', 'sink'],
            ['self.comp1', 'siso'],
            ['self.rhx2', 'sink'],
            ['self.comp2', 'siso'],
            ['self.rhx3', 'sink'],
            ['self.comp3', 'siso'],
            ['self.gas_cooler4', 'source'],
            ['self.rhx3', 'source'],
            ['self.gas_cooler3', 'source'],
            ['self.rhx2', 'source'],
            ['self.gas_cooler2', 'source'],
            ['self.rhx1', 'source'],
            ['self.gas_cooler1', 'source'], 
            ]
        
        sink_ls = [
            ['self.gas_cooler4', 'sink'],
            ['self.gas_cooler3', 'sink'],
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


def set_RHP_model(hp: RHP) -> RHP:
    hp = _get_RHP_structure(hp)
    if _check_input_param_feasibility(hp) == True:
        hp.feasible = True
        hp = _get_RHP_components(hp)
        hp = _get_RHP_streams(hp)
        hp = _get_RHP_logic_blocks(hp)
        hp = _get_RHP_solver_pathway(hp)
    else:
        hp.feasible = False
    return hp

def _check_input_param_feasibility(hp: RHP):
    return True if (
        (hp.T_air_tar - hp.T_air_inlet) * 2/3 >= hp.dT_rhx1_cold + hp.dT_rhx2_cold + hp.dT_rhx3_cold and
        0.0 <= hp.dT_gc_out_min_offset <= 50.0 and
        hp.comp2_power_frac + hp.comp3_power_frac <= 1.0
        ) else False

def _get_RHP_structure(hp: RHP) -> RHP:
    #Check the PHE combo
    hp.dT_rhx1_cold = hp.dT_rhx1_cold if hp.dT_rhx1_cold > 0.0 else 0.0
    hp.dT_rhx2_cold = hp.dT_rhx2_cold if hp.n_comp > 1 and hp.dT_rhx2_cold > 0.0 else 0.0
    hp.dT_rhx3_cold = hp.dT_rhx3_cold if hp.n_comp > 2 and hp.dT_rhx3_cold > 0.0 else 0.0
    
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
    hp.T_gc_out = hp.T_air_inlet + hp.dT_GC_min_tar + hp.dT_gc_out_min_offset
    hp.ihx_on = 1    if hp.T_source_in < hp.T_gc_out and hp.ihx_on == 1 else 0
            
    hp.rhx1_on = 1   if hp.dT_rhx1_cold > 0.0 else 0
    hp.rhx2_on = 1   if hp.dT_rhx2_cold > 0.0 and hp.comp2_on == 1 else 0
    hp.rhx3_on = 1   if hp.dT_rhx3_cold > 0.0 and hp.comp3_on == 1 else 0

    return hp

def _get_RHP_components(hp: RHP) -> RHP:
    #Initialisation of unit operations
    hp.comp1 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.motor_eff, model_index=hp.comp_model_index)
    if hp.comp2_on == 1:
        hp.comp2 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.motor_eff, model_index=hp.comp_model_index)
        hp.gas_cooler3 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
        if hp.comp3_on == 1:
            hp.comp3 = Compressor(name='Compressor', fluid=hp.refrigerant, motor_eff=hp.motor_eff, model_index=hp.comp_model_index) 
            hp.gas_cooler4 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
    
    hp.gas_cooler1 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
    hp.gas_cooler2 = GasCooler(name='Gas cooler', source_fluid=hp.refrigerant, sink_fluid=hp.sink, dTmin=hp.dT_GC_min_tar)
    
    hp.ihx = PlateHeatExchanger(name='Internal heat exchanger', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dTmin=hp.dTmin_rhx, model='dTmin_terminal')
    if hp.rhx1_on == 1:
        hp.rhx1 = PlateHeatExchanger(name='Recuperator', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dT_c=hp.dT_rhx1_cold, model='fixed cold delta temperature')
    if hp.rhx2_on == 1:
        hp.rhx2 = PlateHeatExchanger(name='Recuperator', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dT_c=hp.dT_rhx2_cold, model='fixed cold delta temperature') 
    if hp.rhx3_on == 1:       
        hp.rhx3 = PlateHeatExchanger(name='Recuperator', source_fluid=hp.refrigerant, sink_fluid=hp.refrigerant, dT_c=hp.dT_rhx3_cold, model='fixed cold delta temperature')
    
    hp.exp_valve = ExpansionValve(name='Expansion valve', fluid=hp.refrigerant)
    hp.evaporator = Evaporator(name='Evaporator', source_fluid=hp.source, sink_fluid=hp.refrigerant, model='saturated gas outlet evaporator')
    hp.evap_sup = PlateHeatExchanger(name='Evaporator', source_fluid=hp.source, sink_fluid=hp.refrigerant, model='fixed cold temperature outlet')
    
    return hp

def _get_RHP_streams(hp: RHP) -> RHP:
    temp_ref_out: MassStream
    P_H: float
    P_L: float
    hp.gas_coolers = []
    hp.evaporators = []

    #Evaporator superheat section
    hp.evap_sup.source_in.props_in_ab(input_str='PT', prop1=0.5, prop2=hp.T_source_in, mass_flow=1e6, name='Source in', pp_instance='src')
    hp.evap_sup.source_out.props_in_ab(input_str='PT', prop1=0.5, prop2=hp.T_source_in, mass_flow=1e6, name='Source out', pp_instance='src')
    hp.Tref_evap_outlet = hp.Tref_evap_inlet
    hp.evap_sup.sink_in.props_in_ab(input_str='TQ', prop1=hp.Tref_evap_outlet, prop2=1.0, name='Evaporator sink outlet')
    P_L = hp.evap_sup.sink_in.pressure
    hp.evap_sup.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.T_source_in - hp.dTmin_rhx), name='Superheater sink outlet')
    temp_ref_out = hp.evap_sup.sink_out

    hp.evaporators.append(hp.evap_sup)

    #Internal heat exchanger, if present
    if hp.ihx_on == 1:
        hp.ihx.sink_in = temp_ref_out
        hp.ihx.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.T_gc_out - hp.dTmin_rhx), name='IHX sink outlet')
        temp_ref_out = hp.ihx.sink_out

    #Plate heat exchanger 1 - cold-side, if present
    if hp.rhx1_on == 1:
        hp.rhx1.sink_in = temp_ref_out
        hp.rhx1.sink_out.props_in_ab(input_str='PT', prop1=P_L, prop2=(hp.rhx1.sink_in.temperature + hp.dT_rhx1_cold), name='rhx1 sink outlet')
        temp_ref_out = hp.rhx1.sink_out

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
    if hp.rhx2_on == 1:
        hp.rhx2.sink_in = temp_ref_out
        hp.rhx2.sink_out.props_in_ab(input_str='PT', prop1=hp.rhx2.sink_in.pressure, prop2=(hp.rhx2.sink_in.temperature + hp.dT_rhx2_cold), name='rhx2 sink outlet') 
        temp_ref_out = hp.rhx2.sink_out
    
    #Compressor 2
    if hp.comp2_on == 1:
        hp.comp2.flow_in = temp_ref_out
        hp.comp2.flow_out.props_in_ab(input_str='PS', prop1=(P_L * PR * PR), prop2=hp.comp2.flow_in.entropy, name='Compressor 2 outlet') 
        temp_ref_out = hp.comp2.flow_out
        hp.compressors.append(hp.comp2)

    #Plate heat exchanger 3 - cold-side
    if hp.rhx3_on == 1:
        hp.rhx3.sink_in = temp_ref_out
        hp.rhx3.sink_out.props_in_ab(input_str='PT', prop1=hp.rhx3.sink_in.pressure, prop2=(hp.rhx3.sink_in.temperature + hp.dT_rhx3_cold), name='rhx3 sink outlet') 
        temp_ref_out = hp.rhx3.sink_out

    #Compressor 3
    if hp.comp3_on == 1:
        hp.comp3.flow_in = temp_ref_out
        hp.comp3.flow_out.props_in_ab(input_str='PS', prop1=(P_L * PR * PR * PR), prop2=hp.comp3.flow_in.entropy, name='Compressor 3 outlet') 
        temp_ref_out = hp.comp3.flow_out
        hp.compressors.append(hp.comp3)

    hp.comp_discharge = temp_ref_out
    hp.comp_discharge.pressure_bnds = None
    P_H = temp_ref_out.pressure

    #Gas cooler 4 and Plate heat exchanger 3 - hot-sides
    hp.air_out = MassStream( name='Air sink outlet', fluid=hp.sink, pp_instance='snk', pressure=0.101, temperature=hp.T_air_tar )
    if hp.rhx3_on == 1:
        hp.gas_cooler4.sink_out = hp.air_out
        hp.gas_cooler4.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC4 sink pinch point', pp_instance='snk') 
        hp.gas_cooler4.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC4 sink inlet', pp_instance='snk') 
        temp_air_in = hp.gas_cooler4.sink_in

        hp.gas_coolers.append(hp.gas_cooler4)

        hp.gas_cooler4.source_in = temp_ref_out
        hp.gas_cooler4.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC4 source pinch point') 
        hp.gas_cooler4.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC4 source outlet')

        hp.rhx3.source_in = hp.gas_cooler4.source_out
        hp.rhx3.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='rhx3 source outlet')
        temp_ref_out = hp.rhx3.source_out
    else:
        temp_air_in = hp.air_out
    
    #Gas cooler 3 and Plate heat exchanger 2 - hot-sides
    if hp.rhx2_on == 1:
        hp.gas_cooler3.sink_out = temp_air_in
        hp.gas_cooler3.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC3 sink pinch point', pp_instance='snk') 
        hp.gas_cooler3.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC3 sink inlet', pp_instance='snk') 
        temp_air_in = hp.gas_cooler3.sink_in

        hp.gas_coolers.append(hp.gas_cooler3)

        hp.gas_cooler3.source_in = temp_ref_out
        hp.gas_cooler3.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC3 source pinch point') 
        hp.gas_cooler3.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC3 source outlet')

        hp.rhx2.source_in = hp.gas_cooler3.source_out
        hp.rhx2.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='rhx2 source outlet')
        temp_ref_out = hp.rhx2.source_out

    #Gas cooler 2 and Plate heat exchanger 1 - hot-sides
    if hp.rhx1_on == 1:
        hp.gas_cooler2.sink_out = temp_air_in
        hp.gas_cooler2.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC2 sink pinch point', pp_instance='snk') 
        hp.gas_cooler2.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_tar, name='GC2 sink inlet', pp_instance='snk') 
        temp_air_in = hp.gas_cooler2.sink_in

        hp.gas_coolers.append(hp.gas_cooler2)

        hp.gas_cooler2.source_in = temp_ref_out
        hp.gas_cooler2.source_pinch.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC2 source pinch point') 
        hp.gas_cooler2.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='GC2 source outlet')

        hp.rhx1.source_in = hp.gas_cooler2.source_out
        hp.rhx1.source_out.props_in_ab(input_str='PT', prop1=P_H, prop2=hp.comp_discharge.temperature, name='rhx1 source outlet')
        temp_ref_out = hp.rhx1.source_out

    #Gas cooler 1
    hp.gas_cooler1.sink_out = temp_air_in
    hp.gas_cooler1.sink_pinch.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_inlet, name='GC1 sink pinch point', pp_instance='snk') 
    hp.gas_cooler1.sink_in.props_in_ab(input_str='PT', prop1=0.101, prop2=hp.T_air_inlet, name='GC1 sink inlet', pp_instance='snk')
    
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

    hp.evaporator.source_in.props_in_ab(input_str='PT', prop1=0.5, prop2=hp.T_source_in, mass_flow=1e6, name='Evaporator source inlet', pp_instance='src')
    hp.evaporator.source_out.props_in_ab(input_str='PT', prop1=0.5, prop2=min(hp.evaporator.sink_in.temperature + hp.dTmin_rhx, hp.T_source_in), mass_flow=1e6, name='Evaporator source outlet', pp_instance='src')
    
    hp.evaporators.append(hp.evaporator)

    return hp

def _get_RHP_logic_blocks(hp: RHP) -> RHP:
    #Set special relationships between points in the cycle 
    if hp.rhx3_on == 1:
        hp.set_rhx3_dT_cold = SetRelationship( source_obj=hp.rhx3.sink_in, target_obj=[hp.rhx3.sink_out], c=hp.dT_rhx3_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc4_hot_out = SetRelationship( source_obj=hp.rhx3.sink_out, target_obj=[hp.rhx3.source_in], c=hp.dTmin_rhx, prop='temperature', prop2='pressure' ) 

    if hp.rhx2_on == 1:    
        hp.set_rhx2_dT_cold = SetRelationship( source_obj=hp.rhx2.sink_in, target_obj=[hp.rhx2.sink_out], c=hp.dT_rhx2_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc3_hot_out = SetRelationship( source_obj=hp.rhx2.sink_out, target_obj=[hp.rhx2.source_in], c=hp.dTmin_rhx, prop='temperature', prop2='pressure' ) 
    
    if hp.rhx1_on == 1:
        hp.set_rhx1_dT_cold = SetRelationship( source_obj=hp.rhx1.sink_in, target_obj=[hp.rhx1.sink_out], c=hp.dT_rhx1_cold, prop='temperature', prop2='pressure' ) 
        hp.set_gc2_hot_out = SetRelationship( source_obj=hp.rhx1.sink_out, target_obj=[hp.rhx1.source_in], c=hp.dTmin_rhx, prop='temperature', prop2='pressure' ) 

    if hp.n_comp == 2:
        hp.set_P_ints = CustomRelationship( source_obj=[hp.comp1.flow_in, hp.comp_discharge, hp.comp2_power_frac, hp.comp3_power_frac], target_obj=[hp.comp1.flow_out, hp.comp_discharge], relation=1, tol=hp.tol )
    elif hp.n_comp == 3:
        hp.set_P_ints = CustomRelationship( source_obj=[hp.comp1.flow_in, hp.comp_discharge, hp.comp2_power_frac, hp.comp3_power_frac], target_obj=[hp.comp1.flow_out, hp.comp2.flow_out], relation=1, tol=hp.tol )

    return hp

def _get_RHP_solver_pathway(hp: RHP) -> RHP:
    #Solution pathway set-up
    hp.operations = [hp.evap_sup]

    if hp.n_comp > 1:
        hp.operations.append(hp.set_P_ints)
    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)
    if hp.rhx1_on == 1:
        hp.operations.append(hp.set_rhx1_dT_cold)
        hp.operations.append(hp.set_gc2_hot_out)

    hp.operations.append(hp.comp1)

    if hp.rhx2_on == 1:
        hp.operations.append(hp.set_rhx2_dT_cold)
        hp.operations.append(hp.set_gc3_hot_out)
    if hp.comp2_on == 1:
        hp.operations.append(hp.comp2)
    if hp.rhx3_on == 1:
        hp.operations.append(hp.set_rhx3_dT_cold)
        hp.operations.append(hp.set_gc4_hot_out)
    if hp.comp3_on == 1:
        hp.operations.append(hp.comp3)
    
    if hp.rhx3_on == 1:
        hp.operations.append(hp.gas_cooler4)
        hp.operations.append(hp.rhx3)
    if hp.rhx2_on == 1:
        hp.operations.append(hp.gas_cooler3)
        hp.operations.append(hp.rhx2)
    if hp.rhx1_on == 1:
        hp.operations.append(hp.gas_cooler2)        
        hp.operations.append(hp.rhx1)

    hp.operations.append(hp.gas_cooler1)

    if hp.ihx_on == 1:
        hp.operations.append(hp.ihx)

    hp.operations.append(hp.exp_valve)
    hp.operations.append(hp.evaporator)

    return hp
    