__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

# from streams import PropertyPackage, MassStream
from components import Compressor, GasCooler, PlateHeatExchanger, ExpansionValve, Evaporator, SetRelationship, CustomRelationship
from functions.graph import Graph
from functions.readWriteOutput import print_result_to_csv
from scipy.optimize import root_scalar 
import numpy as np


class HeatPump():
    gas_coolers: list[GasCooler]
    compressors: list[Compressor]
    feasible: bool

    def __init__(
            self, 
            refrigerant='R1336mzz(Z)', 
            source='Water',
            sink='Air',
            T_air_tar=200.0, 
            T_air_inlet=15.0, 
            dT_GC_min_tar=20.0, 
            T_gc_out_min=0.0, 
            T_source_in=40.0, 
            Tsat_evap_inlet=10.0,
            dTmin_phx=5.0,
            min_dT_ihx = 5.0,
            n_comp=3, 
            comp_PR_uniform=False, 
            eff_motor=0.96, 
            dT_phx1_cold=0.0, 
            dT_phx2_cold=0.0, 
            dT_phx3_cold=0.0, 
            comp2_power_frac=0.0, 
            comp3_power_frac=0.0,
            P_H_est=0.0,
            tol=1e-12,
            display=False,
            result_parameters=None
        ):
        
        self.refrigerant = refrigerant
        self.source = source
        self.sink = sink
        self.cycle_type = None
        self.T_air_tar = float(T_air_tar)
        self.dT_GC_min_tar = float(dT_GC_min_tar)
        self.T_air_inlet = float(T_air_inlet)
        self.T_source_in = float(T_source_in)
        self.T_gc_out = max(self.T_air_inlet + self.dT_GC_min_tar, T_gc_out_min)
        self.Tsat_evap_outlet = float(Tsat_evap_inlet)
        self.Tsat_evap_inlet = float(Tsat_evap_inlet)
        self.n_comp = int(n_comp)
        self.comp2_power_frac = float(comp2_power_frac)
        self.comp3_power_frac = float(comp3_power_frac)
        self.comp_PR_uniform = bool(comp_PR_uniform)
        self.dTmin = float(dTmin_phx)
        self.eff_motor = float(eff_motor)
        self.min_dT_ihx = float(min_dT_ihx)
        self.dT_phx1_cold = float(dT_phx1_cold) #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.dT_phx2_cold = float(dT_phx2_cold) #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.dT_phx3_cold = float(dT_phx3_cold) #if dT_phx1_cold >= self.min_dT_ihx else 0.0
        self.P_H_est = float(P_H_est)
        self.tol = float(tol)
        self.converged = False
        self.COP = -np.inf
        self.disp = display
   
    def print_failed_to_solve_report(self, filename='results/Cases-Unsolved'):
        print('Heat pump model failed to solve.')
        self.get_summary()
        print_result_to_csv(self.output, file_name=filename)        

    def print_currrent_case(self):
        print({
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
            })

    def get_summary(self):
        self.output = {
            #Parameters
            'Refrigerant:' : self.refrigerant,
            'Type:' : self.cycle_type,
            'n_comp:' : self.n_comp,
            'T_source_in:' : self.evaporator.source_in.temperature,
            'dT_min Gas Cooler:' : self.dT_GC_min_tar,
            'T_air_inlet:' : self.gas_cooler1.sink_in.temperature,
            'T_air_tar:' : self.T_air_tar,
            'Comp_PR_uniform:' : self.comp_PR_uniform,
            #Key results
            'COP:' : self.COP,
            'Tol:' : abs(self.T_air_tar - self.air_out.temperature) / self.T_air_tar,
            'Converged:' : self.converged,
            #Additional results
            'Actual air T:' : self.air_out.temperature,            
            'Evap T:' : self.Tsat_evap_inlet,
            'Evap P:' : self.evaporator.sink_in.pressure * 1000,
            'Comp 1 P out:' : self.comp1.flow_out.pressure * 1000,
            'Comp 2 P out:' : self.comp2.flow_out.pressure * 1000 if self.n_comp > 1 else self.comp1.flow_out.pressure * 1000,
            'Comp 3 P out:' : self.comp_discharge.pressure * 1000,
            'Comp 1 power frac:' : self.comp1_power_frac, 
            'Comp 2 power frac:' : self.comp2_power_frac, 
            'Comp 3 power frac:' : self.comp3_power_frac,            
            'dT phx max 1:' : self.dT_phx1_cold,
            'dT phx max 2:' : self.dT_phx2_cold,
            'dT phx max 3:' : self.dT_phx3_cold,
            'PR 1:' : self.comp1.pressure_ratio,
            'PR 2:' : self.comp2.pressure_ratio if self.n_comp > 1 else 1.0,
            'PR 3:' : self.comp3.pressure_ratio if self.n_comp > 2 else 1.0,
            'Overall P ratio:' : self.comp_discharge.pressure / self.evaporator.sink_in.pressure,
            'Gas cooler(s) inlet T:' : self.comp_discharge.temperature,
            'Gas cooler(s) outlet T:' : self.gas_cooler1.source_out.temperature
        }

    def get_cycle_plot(self, x_axis='h', y_axis='T'):
        cycle_ls: list[list[str, str]]
        cycle_obj: list[list[object, str]]
        sink_ls: list[list[str, str]]
        sink_obj: list[list[object, str]]

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
        
        sink_ls=[
                ['self.gas_cooler4', 'sink'],
                ['self.gas_cooler3', 'sink'],
                ['self.gas_cooler2', 'sink'],
                ['self.gas_cooler1', 'sink'],
            ]
        
        cycle_obj = []
        for cyc_component in cycle_ls:
            obj, att = cyc_component[0].split('.')
            if hasattr(locals()[obj], att):
                component = getattr(locals()[obj], att)
                if component.converged:
                    cycle_obj.append([ component, cyc_component[1] ])

        sink_obj = []
        for cyc_component in sink_ls:
            obj, att = cyc_component[0].split('.')
            if hasattr(locals()[obj], att):
                component = getattr(locals()[obj], att)
                if component.converged:
                    sink_obj.append([ component, cyc_component[1] ])

        g = Graph()
        g.plot_cycle(cycle_obj=cycle_obj, 
                sink_obj=sink_obj,
                x_axis=x_axis, y_axis=y_axis)

    def get_driving_force_plot(self, x_axis='h'):
        g = Graph()
        g.plot_driving_force(hx_obj=self.gas_coolers,
                        x_axis=x_axis,
                        dtmin=self.dT_GC_min_tar)
