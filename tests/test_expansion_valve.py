import sys
import os
import unittest

# Add the parent directory to the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from components import ExpansionValve
from streams import MassStream

class TestCompressor(unittest.TestCase):
    def setUp(self):
        flow_in = MassStream(fluid='Air', pressure=100, temperature=20)
        flow_out = MassStream(fluid='Air', pressure=100, temperature=20)
        self.c = ExpansionValve(flow_in=flow_in, flow_out=flow_out)
        self.c.heat_loss = 0.0

    def tearDown(self):
        del self.c
    
    def test_should_return_outlet_enthalpy(self):
        result = self.c.get_outlet_enthalpy()
        self.assertAlmostEqual(result, self.c.flow_in.enthalpy)


if __name__ == '__main__':
    unittest.main()
