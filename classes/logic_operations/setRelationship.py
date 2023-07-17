__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

class Set_Relationship:
    def __init__(self, source_obj, target_obj, m = 1, c = 0, prop = 'Temperature', prop2 = 0):
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.m = m
        self.c = c
        self.prop = prop
        self.prop2 = prop2
        self.solve()

    def solve(self):
        if self.prop == 'Temperature':
            for obj in self.target_obj: 
                obj.Temperature = self.m * self.source_obj.Temperature + self.c
                if self.prop2 == 'Pressure':
                    obj.PropsIn_mPT()
                elif self.prop2 == 'Quality':
                    obj.PropsIn_mTX()

        elif self.prop == 'Mass Flow':
            for obj in self.target_obj: 
                obj.Mass_Flow = self.m * self.source_obj.Mass_Flow + self.c
                obj.Heat_Flow = obj.Mass_Flow * obj.Enthalpy