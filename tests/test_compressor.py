from ..components import Compressor
import unittest


class TestCompressor(unittest.TestCase):

    def setUp(self):
        """
        If there is any setup needed for your tests, it goes here.
        This method is called before each individual test method runs.
        """
        pass


    def tearDown(self):
        """
        Any cleanup after tests goes here.
        This method is called after each individual test method runs.
        """
        pass


    def test_should_return_number(self):
        """
        Here's where your actual test methods go. 
        Test methods must start with the word 'test'.
        """
        # Call the function you're testing
        result = Compressor.get_turbine_work()
        self.assertTrue(result.isnumeric())

    def template(self):
        # Use an assertion method to check the result.
        # There are many assertion methods you can use, like:
        # - assertEqual(a, b)
        # - assertTrue(x)
        # - assertFalse(x)
        # - assertIn(a, b)
        # - assertIsInstance(a, b)
        #         
        pass


if __name__ == '__main__':
    unittest.main()
