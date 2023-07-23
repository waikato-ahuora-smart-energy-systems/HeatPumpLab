__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams import MassStream

class Splitter: 
    def __init__(self, Name='Splitter', split_ratio_list = []):
        self.Name = Name
        self.FlowIn = MassStream()
        self.FlowOut = [MassStream()]
        self.split_ratio_list = split_ratio_list

    def add_flow_out(self, flow_j):
        self.FlowOut.append(flow_j)

    def solve(self):
        i = 1
        for f in self.FlowOut:
            f.SetEqual = self.FlowIn
            f.Mass_Flow = self.split_ratio_list[i]
            f.Heat_Flow = self.FlowIn.Enthalpy * f.Mass_Flow
            i += 1
            