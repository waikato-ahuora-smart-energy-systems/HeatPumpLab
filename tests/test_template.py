import unittest

class TestMyFunction(unittest.TestCase):
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

    def test_case_1(self):
        """
        Here's where your actual test methods go. 
        Test methods must start with the word 'test'.
        """
        # Call the function you're testing
        result = my_function(input1)
        
        # Use an assertion method to check the result.
        # There are many assertion methods you can use, like:
        # - assertEqual(a, b)
        # - assertTrue(x)
        # - assertFalse(x)
        # - assertIn(a, b)
        # - assertIsInstance(a, b)
        self.assertEqual(result, expected_output_1)

    def test_case_2(self):
        """
        It's often a good idea to have multiple tests for each function, to test different cases.
        """
        result = my_function(input2)
        self.assertEqual(result, expected_output_2)


if __name__ == '__main__':
    unittest.main()
