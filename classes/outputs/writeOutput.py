__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
import csv


def terminal_print(output):
    for key, value in output.items():
        print(f'{key} {value}')


def append_csv(output, file_name='TRHP_results_v2.csv'):
    file_exists = os.path.isfile(file_name)
    with open(file_name, 'a', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=output.keys())

        if not file_exists:
            writer.writeheader()  # file doesn't exist yet, write a header

        writer.writerow(output)
