__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os, csv

def print_result_to_terminal(output):
    if type(output) is dict:
        for key, value in output.items():
            print(f'{key} {value}')

    else:
        print('Failed to print due to wrong output format.')

def print_result_to_csv(output='', file_name='results/untitled.csv'):
    if type(output) is dict:
        file_exists = os.path.isfile(file_name)
        with open(file_name, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=output.keys())

            if not file_exists:
                writer.writeheader()  # file doesn't exist yet, write a header

            writer.writerow(output)
            
    else:
        print('Failed to appended CSV due to wrong output format.')
