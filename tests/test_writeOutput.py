import unittest
from unittest.mock import patch, MagicMock
from outputs.writeOutput import print_result_to_csv


class TestAppendCSV(unittest.TestCase):

    @patch('os.path.isfile', return_value=True)
    @patch('csv.DictWriter')
    @patch('builtins.open', new_callable=MagicMock)
    def test_append_csv_existing_file(self, mock_open, mock_csv, mock_isfile):
        data = {"name": "John Doe", "age": 32, "city": "New York"}
        print_result_to_csv(data, 'results\sample.csv')
        mock_open.assert_called_once_with('results\sample.csv', 'a', newline='')
        mock_csv().writerow.assert_called_once_with(data)
        

    @patch('os.path.isfile', return_value=False)
    @patch('csv.DictWriter')
    @patch('builtins.open', new_callable=MagicMock)
    def test_append_csv_new_file(self, mock_open, mock_csv, mock_isfile):
        data = {"name": "John Doe", "age": 32, "city": "New York"}
        print_result_to_csv(data, 'results\sample.csv')
        mock_open.assert_called_once_with('results\sample.csv', 'a', newline='')
        mock_csv().writeheader.assert_called_once()
        mock_csv().writerow.assert_called_once_with(data)


    def test_append_csv_wrong_data(self):
        with self.assertRaises(TypeError):
            print_result_to_csv('not a dict', 'results\sample.csv')


if __name__ == '__main__':
    unittest.main()
