__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
import csv

class Write_Output:
    def __init__(self, output='') -> None:
        self.output = output


    def terminal_print(self):
        for key, value in self.output.items():
            print(f'{key} {value}')


    def append_csv(self, file_name='TRHP_results_v1.csv'):
        file_exists = os.path.isfile(file_name)
        with open(file_name, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=self.output.keys())

            if not file_exists:
                writer.writeheader()  # file doesn't exist yet, write a header

            writer.writerow(self.output)

    def append_output_dict(self, result):
        return self.output.update({
            'Obj function:' : result.fun,
            'Iterations:' : result.nfev,
            'Message:' : result.message,
            'Success:' : result.success
        })
