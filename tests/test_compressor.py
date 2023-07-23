import sys
import os
import unittest

# Add the parent directory to the Python path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from components import Compressor
from streams import MassStream

class TestCompressor(unittest.TestCase):
    def setUp(self):
        flow_in = MassStream(fluid='Air', pressure=100, temperature=20)
        flow_out = MassStream(fluid='Air', pressure=100, temperature=20)
        self.c = Compressor(flow_in=flow_in, flow_out=flow_out)
        self.c.isentropic_eff = 1.0
        self.c.motor_eff = 1.0

    def tearDown(self):
        del self.c
    
    def test_should_solve_compressor(self):
        self.c.solve()
        self.assertAlmostEqual(self.c.flow_out.enthalpy, self.c.flow_in.enthalpy)

    def test_should_return_compressor_power(self):
        result = self.c.get_compressor_power()
        self.assertAlmostEqual(result, 0.0)

    def test_should_return_isen_efficiency(self):
        result = self.c.get_isentropic_efficiency()
        self.assertEqual(result, 1.0)

    def test_should_return_model_coeff(self):
        result = self.c.set_model_coeff()
        self.assertEqual(result[0], 1.0)
        self.assertEqual(result[1], 0.0)
        self.assertEqual(result[2], 0.0)


if __name__ == '__main__':
    unittest.main()
