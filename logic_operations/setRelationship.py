__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

class SetRelationship:
    def __init__(self, name='', source_obj='', target_obj='', m = 1, c = 0, prop = 'temperature', prop2 = 0):
        self.name = name
        self.source_obj = source_obj
        self.target_obj = target_obj
        self.m = m
        self.c = c
        self.prop = prop
        self.prop2 = prop2
        self.solve()

    def solve(self):
        if self.prop == 'temperature':
            self.set_temperature()

        elif self.prop == 'quality':
            self.set_quality()

        elif self.prop == 'Mass Flow':
            self.set_mass_flow()

    def set_temperature(self):
        for obj in self.target_obj: 
            obj.temperature = self.m * self.source_obj.temperature + self.c
            if self.prop2 == 'pressure':
                obj.props_in_mPT()
            elif self.prop2 == 'quality':
                obj.props_in_mTX()
            else:
                raise ValueError('Set relationship failed: ' + self.name)   

    def set_quality(self):
        for obj in self.target_obj: 
            obj.quality = self.source_obj[0]
            if self.prop2 == 'pressure':
                obj.props_in_mPX()
            elif self.prop2 == 'temperature':
                obj.props_in_mTX()
            else:
                raise ValueError('Set relationship failed: ' + self.name)   

    def set_mass_flow(self):
        for obj in self.target_obj: 
            obj.mass_flow = self.m * self.source_obj.mass_flow + self.c
            obj.heat_flow = obj.mass_flow * obj.enthalpy