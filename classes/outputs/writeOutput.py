__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
import csv

def terminal_print(output):
    if type(output) is dict:
        for key, value in output.items():
            print(f'{key} {value}')

    else:
        print('Failed to print due to wrong output format.')


def append_csv(output='', file_name='untitled.csv'):
    if type(output) is dict:
        file_exists = os.path.isfile(file_name)
        with open(file_name, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=output.keys())

            if not file_exists:
                writer.writeheader()  # file doesn't exist yet, write a header

            writer.writerow(output)
            
    else:
        print('Failed to appended CSV due to wrong output format.')
