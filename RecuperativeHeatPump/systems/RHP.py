from ..components import *
from ..streams import MassStream
from .common.heatPump import BaseHeatPump
from .RHP_constructor import *
from .solver.model_solver import solve_heat_pump_model


class RHP(BaseHeatPump):
    
    n_comp: float =2
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

        set_RHP_model(self)

    def solve(self) -> "RHP":
        if self.feasible:
            solve_heat_pump_model(self)
        return self
   
    def get_model_parameters(self) -> dict[str, float | bool]:
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
    
    def get_report(self, dt=0) -> dict[str, float | bool]:
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
