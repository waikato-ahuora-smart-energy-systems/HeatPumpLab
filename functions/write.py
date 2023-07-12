__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
import csv


def terminal_print_result(output):
    for key, value in output.items():
        print(f'{key} {value}')


def append_result_to_csv(output, file_name='TRHP_results_v1.csv'):
    file_exists = os.path.isfile(file_name)
    with open(file_name, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=output.keys())

        if not file_exists:
            writer.writeheader()  # file doesn't exist yet, write a header

        writer.writerow(output)

def append_result_to_output(output, result):
    return output.update({
        'Obj function:' : result.fun,
        'Iterations:' : result.nfev,
        'Message:' : result.message,
        'Success:' : result.success
    })


