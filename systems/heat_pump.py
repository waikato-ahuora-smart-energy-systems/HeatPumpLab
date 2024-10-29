__author__ = 'Tim Walmsley and Keri Adamson'
__email__ = 'tim.walmsley@waikato.ac.nz'
__date__ = '2021-2024'

"""
High temperature transcritical-recuperative heat pump model
By Tim Walmsley and Keri Adamson, University of Waikato, 2021-2024 (tim.walmsley@waikato.ac.nz)
"""

from streams import PropertyPackage, MassStream
from components import Compressor, PlateHeatExchanger, GasCooler, ExpansionValve, Evaporator, SetRelationship, CustomRelationship
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
            Tsat_evap=10.0,
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
        self.Tsat_evap_tar = float(Tsat_evap)
        self.Tsat_evap = float(Tsat_evap)
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
            'Evap T:' : self.Tsat_evap,
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
        for component in cycle_ls:
            obj, att = component[0].split('.')
            if hasattr(locals()[obj], att):
                cycle_obj.append([ getattr(locals()[obj], att), component[1] ])

        sink_obj = []
        for component in sink_ls:
            obj, att = component[0].split('.')
            if hasattr(locals()[obj], att):
                sink_obj.append([ getattr(locals()[obj], att), component[1] ])

        g = Graph()
        g.plot_cycle(cycle_obj=cycle_obj, 
                sink_obj=sink_obj,
                x_axis=x_axis, y_axis=y_axis)

    def get_driving_force_plot(self, x_axis='h'):
        g = Graph()
        g.plot_driving_force(hx_obj=self.gas_coolers,
                        x_axis=x_axis,
                        dtmin=self.dT_GC_min_tar)






    def get_pressure_bounds(self):
        feasible = 0
        P_lb = self.comp_discharge.pressure
        P_ub = P_lb
        P_i = P_lb

        dt_plus = self.get_solve_state(P_lb)
        dt_neg = dt_plus
        P = np.array(P_lb)
        DT = np.array(dt_plus)

        i = 0
        m_step = 0.75
        if dt_neg > 0.0:
            #Find a lower bound pressure resulting in an air temperature less than the target
            while dt_neg > 0.0 and i < 100:
                dt_plus = dt_neg
                P_ub = P_lb
                P_lb = self.comp1.flow_in.pressure + (P_lb - self.comp1.flow_in.pressure) * m_step
                dt_neg = self.get_solve_state(P_lb)
                P = np.append(P, P_lb)
                DT = np.append(DT, dt_neg)
                i += 1

        else:
            #Find a upper bound pressure resulting in an air temperature greater than the target
            while dt_plus < 0.0 and i < 100:
                dt_neg = dt_plus
                P_lb = P_ub
                P_ub = self.comp1.flow_in.pressure + (P_ub - self.comp1.flow_in.pressure) / m_step
                dt_plus = self.get_solve_state(P_ub)
                P = np.append(P, P_ub)
                DT = np.append(DT, dt_plus)
                i += 1

        if P_lb < self.comp1.flow_in.pressure:
            #Check the low pressure bound remains above the evaporator pressure
            P_lb = self.comp1.flow_in.pressure
            dt_neg = self.get_solve_state(P_lb)
            if P_ub < P_lb:
                #Check the upper pressure bound stays above the low pressure bound
                dt_neg = dt_plus
        
        if abs(dt_plus * dt_neg) == -dt_plus * dt_neg and abs(dt_plus * dt_neg) > 0.0:
            feasible = 1
            P_i = P_lb - dt_neg / (dt_plus - dt_neg) * (P_ub - P_lb)
            if P_i * 1.01 < P_ub:
                dt = self.get_solve_state(P_i * 1.01)
                if dt > 0.0:
                    P_ub = P_i * 1.01
                    dt_plus = dt
            if P_i * 0.99 > P_lb:
                dt = self.get_solve_state(P_i * 0.99)
                if dt < 0.0:
                    P_lb = max(P_i * 0.99, P_lb)
                    dt_neg = dt
            
            P_i = P_lb - dt_neg / (dt_plus - dt_neg) * (P_ub - P_lb)

        return P_lb, P_ub, P_i, feasible

    def get_p_test(self):
        P = np.linspace(start=self.P_lb / 1.5, stop=self.P_ub * 3, num=1000)
        dT = np.array([self.T_air_tar + self.get_solve_state(P_i) for P_i in P])
        Graph().plot_simple_xy(x=P, y=dT)
        ddT = np.diff(dT)
        i = ddT.argmax()
        P0, P1 = P[i], P[i+1]
        None
 
    def get_solve_state(self, P_comp_out):
        self.comp_discharge.pressure = P_comp_out
        i = 0; i_max = 10
        h = 0.0; m = 0.0; P = 0.0; W = 0.0
        dh = 0; dm = 0; dP = 0; dW = 0

        while ((i == 0 or min(dP, dW, dQ, dh, dm) < 0.0) and i < i_max): #Check for model convergance: dh, dm, 
            # Saves pressure and temperature properties to check for solution convergance
            h = self.exp_valve.flow_in.enthalpy
            P = self.comp_discharge.pressure - self.gas_cooler1.source_out.pressure
            W = sum([c.power for c in self.compressors])
            Q = sum([g.duty for g in self.gas_coolers])
            i += 1

            #Solve the performance of each component and update each special relationship
            for obj in self.operations:
                obj.solve()

            dh = self.tol * h - abs(self.exp_valve.flow_in.enthalpy - h)
            dP = self.tol * P - abs(self.comp_discharge.pressure - self.gas_cooler1.source_out.pressure - P)
            dW = self.tol * W - abs(sum([c.power for c in self.compressors]) - W)
            dQ = self.tol * Q - abs(sum([g.duty for g in self.gas_coolers]) - Q)

        self.set_gc_air_massflow.solve()

        if 0:
            self.get_cycle_plot(x_axis='h', y_axis='T')
            None

        if i < i_max:
            self.converged = True
            e = self.air_out.temperature - self.T_air_tar
        else:
            self.converged = False
            e = -self.T_air_tar

        return e

    def get_solve_gas_pressure(self, method, bnds, x0):
        r = root_scalar(f=self.get_solve_state, 
                        method=method, 
                        bracket=bnds, 
                        x0=x0, 
                        xtol=self.tol,
                        )
        r['error'] = 1e-4 - abs(self.T_air_tar - self.air_out.temperature) / self.T_air_tar
        if self.disp:
            # self.get_cycle_plot()
            print(r['error'])
            
        return r

    def get_phx_min_dT_for_solve_quality_method(self, QMASS_gc_out):
        self.q_gc_out[0] = QMASS_gc_out
        r = self.get_solve_gas_pressure(method='secant', 
                                        x0=self.P_i, 
                                        bnds=(self.P_lb, self.P_ub))
        hx = self.phx_fixed_src_quality
        e = hx.source_in.temperature - hx.sink_out.temperature - self.dTmin
        return e
    
    def get_solve_quality(self):
        e_one = self.get_phx_min_dT_for_solve_quality_method(1.0)
        P_one = self.comp3.flow_out.pressure
        e_zero = self.get_phx_min_dT_for_solve_quality_method(0.0)
        P_zero = self.comp3.flow_out.pressure
        if e_zero * e_one < 0:
            r = root_scalar(f=self.get_phx_min_dT_for_solve_quality_method,
                            method='brentq',
                            bracket=(0.0, 1.0),
                            x0=0.5,
                            xtol=1e-3
                            )
            r['error'] = 1e-4 - abs(self.T_air_tar - self.gas_cooler4.sink_out.temperature)
        else:
            if abs(P_zero - P_one) > 1:
                None
            r = {'error' : 1.0}
        return r

    def get_solve_T_evap(self, T_evap_sat_out):
        #Solve heat pump cycle
        self.evaporator.sink_out.temperature = T_evap_sat_out
        self.evaporator.sink_out.props_in_mTQ()
        self.converged = False
        self.COP = -np.inf
        
        if (self.T_air_tar - self.T_air_inlet) * 2/3 + 2 >= self.dT_phx1_cold + self.dT_phx2_cold + self.dT_phx3_cold:
            self.P_lb, self.P_ub, self.P_i, self.feasible = self.get_pressure_bounds()

            if self.feasible == 1:
                # self.get_p_test()
                try:
                    r = self.get_solve_gas_pressure(method='secant', x0=self.P_i, bnds=(self.P_lb, self.P_ub))
                except:
                    # self.print_failed_to_solve_report()
                    r = {'error' : -1,
                         'root' : self.P_i}

                if r['error'] < 0.0:
                    # Possibly unconverged due to a fixed hot outlet temperature of a gas cooler 
                    # being near the saturation temperature, causing instability
                    if r['root'] < self.evaporator.sink_out.P_crit:
                        self.set_quality_relationship()
                        self.P_i = r['root']
                        r = self.get_solve_quality()
                        if self.disp:
                            self.get_cycle_plot()
                    else:
                        dt_plus = self.get_solve_state(self.evaporator.sink_out.P_crit - 10e-4)
                        dt_neg = self.get_solve_state(self.evaporator.sink_out.P_crit + 10e-4)                    
                        if dt_neg < 0.0:
                            self.P_lb = self.evaporator.sink_out.P_crit + 10e-4
                            r = self.get_solve_gas_pressure(method='bisect', x0=max(r['root'], self.P_lb), bnds=(self.P_lb, self.P_ub))
                            if self.T_air_tar < self.gas_cooler4.sink_out.temperature:
                                if self.get_solve_state(r['root'] + 10e-6) > 0.0:
                                    r['root'] += 10e-6
                        
                        if dt_plus > 0.0:
                            self.P_ub = self.evaporator.sink_out.P_crit - 10e-4
                            r = self.get_solve_gas_pressure(method='bisect', x0=min(r['root'], self.P_ub), bnds=(self.P_lb, self.P_ub))
                            if self.T_air_tar < self.gas_cooler4.sink_out.temperature:
                                if self.get_solve_state(r['root'] + 10e-6) > 0.0:
                                    r['root'] += 10e-6
                                
                        if self.disp:
                            self.get_cycle_plot()

                    if r['error'] < 0.0 and self.T_air_tar > self.gas_cooler4.sink_out.temperature:
                        self.print_failed_to_solve_report()
                             
                if (r['error'] > 0.0 or self.T_air_tar > self.gas_cooler4.sink_out.temperature) and self.converged:
                    W = sum([c.power for c in self.compressors])
                    Q = sum([g.duty for g in self.gas_coolers])
                    self.COP = Q / W
                else:
                    self.converged = False
        
                self.cycle_type = 'Transcritical' if self.comp_discharge.pressure > self.evaporator.sink_out.P_crit else 'Subcritical'

        return self.evaporator.sink_in.temperature - self.Tsat_evap_tar

    def solve(self):
        #Solve heat pump cycle
        self.converged = False
        self.COP = -np.inf
        # try:
        r = self.get_solve_T_evap(self.Tsat_evap_tar)
        if abs(r) > 1e-6 and abs(self.evaporator.sink_in.temperature - self.evaporator.sink_out.temperature) > 0.1:
            r = root_scalar(
                            f=self.get_solve_T_evap,
                            method='secant',
                            x0=self.Tsat_evap_tar,
                            bracket=(self.Tsat_evap_tar, self.T_source_in),
                            xtol=self.tol
                        )
        # except:
        #     self.print_failed_to_solve_report()
                            
        return self.COP, self.comp_discharge.pressure
