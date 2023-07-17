__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

class Splitter: 
    def __init__(self, split_ratio_list):
        self.FlowIn = None
        self.FlowOut = [None]
        self.split_ratio_list = split_ratio_list

    def add_flow_out(self, flow_j):
        self.FlowOut.append(flow_j)

    def solve(self):
        i = 0
        for f in self.FlowOut:
            f.Mass_Flow = self.split_ratio_list[i]
            f.Pressure = self.FlowIn.Pressure
            f.Temperature = self.FlowIn.Temperature
            f.Enthalpy = self.FlowIn.Enthalpy
            f.Entropy = self.FlowIn.Entropy
            f.Quality = self.FlowIn.Quality
            f.Heat_Flow = self.FlowIn.Enthalpy * f.Mass_Flow
            i += 1        