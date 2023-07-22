__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']


class SetRelationship:
    
    def __init__(self, Name='', source_obj='', target_obj='', m = 1, c = 0, prop = 'Temperature', prop2 = 0):
        self.Name = Name
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.m = m
        self.c = c
        self.prop = prop
        self.prop2 = prop2
        self.solve() #Change this to being explicit


    def solve(self):
        if self.prop == 'Temperature':
            self.set_temperature()

        elif self.prop == 'Mass Flow':
            self.set_mass_flow()


    def set_temperature(self):
        for obj in self.target_obj: 
            obj.Temperature = self.m * self.source_obj.Temperature + self.c
            if self.prop2 == 'Pressure':
                obj.PropsIn_mPT()
            elif self.prop2 == 'Quality':
                obj.PropsIn_mTX()
            else:
                raise ValueError('Set relationship failed: ' + self.Name)   


    def set_mass_flow(self):
        for obj in self.target_obj: 
            obj.Mass_Flow = self.m * self.source_obj.Mass_Flow + self.c
            obj.Heat_Flow = obj.Mass_Flow * obj.Enthalpy