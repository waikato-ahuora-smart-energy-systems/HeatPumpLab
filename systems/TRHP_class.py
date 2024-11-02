__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

from systems.common.commonHeatPump import HeatPump

class TRHP(HeatPump):
    T_air_tar: float =200.0 
    T_air_inlet: float =15.0 
    dT_GC_min_tar: float =20.0 
    T_gc_out_min: float =0.0
    T_gc_out: float 
    T_source_in: float =40.0 
    Tref_evap_inlet: float =10.0
    Tref_evap_outlet: float =10.0
    dTmin_phx: float =5.0
    
    min_dT_ihx : float = 5.0
    n_comp: float =3 
    comp_PR_uniform: float =False 
    eff_motor: float =0.96 
    dT_phx1_cold: float =0.0 
    dT_phx2_cold: float =0.0 
    dT_phx3_cold: float =0.0 
    comp1_power_frac: float =1.0
    comp2_power_frac: float =0.0 
    comp3_power_frac: float =0.0
    P_H_est: float =0.0
    result_parameters: float =None    

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
            'n_comp' : self.n_comp,
            'dT_phx1_cold' : self.dT_phx1_cold,
            'dT_phx2_cold' : self.dT_phx2_cold,
            'dT_phx3_cold' : self.dT_phx3_cold,
            'comp2_power_frac' : self.comp2_power_frac,
            'comp3_power_frac' : self.comp3_power_frac,
            'comp_PR_uniform' : self.comp_PR_uniform
            }

    def get_components_sets_to_graph(self) -> list:
        cycle_ls: list[list[str, str]]
        sink_ls: list[list[str, str]]
        
        cycle_ls = [
            ['self.ihx', 'source'],
            ['self.exp_valve', 'siso'],
            ['self.evaporator', 'sink'],
            ['self.evap_sup', 'sink'],
            ['self.ihx', 'sink'],
            ['self.phx1', 'sink'],
            ['self.comp1', 'siso'],
            ['self.phx2', 'sink'],
            ['self.comp2', 'siso'],
            ['self.phx3', 'sink'],
            ['self.comp3', 'siso'],
            ['self.gas_cooler4', 'source'],
            ['self.phx3', 'source'],
            ['self.gas_cooler3', 'source'],
            ['self.phx2', 'source'],
            ['self.gas_cooler2', 'source'],
            ['self.phx1', 'source'],
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
    
    def get_summary(self) -> dict[str, float | bool]:
        self.output = {
            #Headline results
            'Type:' : self.cycle_type,
            'COP:' : self.COP,
            'Converged:' : self.converged,
            #Parameters
            'Refrigerant:' : self.refrigerant,
            'n comp:' : self.n_comp,
            'T source in:' : self.evaporator.source_in.temperature,
            'dT GC min tar:' : self.dT_GC_min_tar,
            'T air inlet:' : self.gas_cooler1.sink_in.temperature,
            'T air tar:' : self.T_air_tar,
            'Comp PR uniform:' : self.comp_PR_uniform,
            'T evap inlet:' : self.Tref_evap_inlet,
            #Optmisation variables
            'dT phx max 1:' : self.dT_phx1_cold,
            'dT phx max 2:' : self.dT_phx2_cold,
            'dT phx max 3:' : self.dT_phx3_cold,
            'Comp 1 power frac:' : self.comp1_power_frac, 
            'Comp 2 power frac:' : self.comp2_power_frac, 
            'Comp 3 power frac:' : self.comp3_power_frac,  
            'GC refrigerant outlet T:' : self.gas_cooler1.source_out.temperature,
            #Additional results
            'T_evap_outlet:' : self.Tref_evap_outlet,
            'Evap P:' : self.evaporator.sink_in.pressure * 1000,
            'Comp 1 P out:' : self.comp1.flow_out.pressure * 1000,
            'Comp 2 P out:' : self.comp2.flow_out.pressure * 1000 if self.n_comp > 1 else self.comp_discharge.pressure * 1000,
            'Comp 3 P out:' : self.comp_discharge.pressure * 1000,
            'PR 1:' : self.comp1.pressure_ratio,
            'PR 2:' : self.comp2.pressure_ratio if self.n_comp > 1 else 1.0,
            'PR 3:' : self.comp3.pressure_ratio if self.n_comp > 2 else 1.0,
            'Overall P ratio:' : self.comp_discharge.pressure / self.evaporator.sink_in.pressure,
            'Comp discharge T:' : self.comp_discharge.temperature,
            }
        return self.output

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