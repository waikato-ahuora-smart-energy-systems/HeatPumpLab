__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from streams import MassStream

class Splitter: 
    def __init__(self, name='Splitter', split_ratio_list = []):
        self.name = name
        self.flow_in = MassStream()
        self.flow_out = [MassStream()]
        self.split_ratio_list = split_ratio_list

    def add_flow_out(self, flow_j):
        self.flow_out.append(flow_j)

    def solve(self):
        i = 1
        for f in self.flow_out:
            f.set_equal = self.flow_in
            f.mass_flow = self.split_ratio_list[i]
            f.heat_flow = self.flow_in.enthalpy * f.mass_flow
            i += 1
            