__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os, csv

def print_result_to_terminal(output):
    if type(output) is dict:
        for key, value in output.items():
            print(f'{key} {value}')

    else:
        print('Failed to print due to wrong output format.')

def print_result_to_csv(output=[], file_name='results/untitled.csv'):
    try:
        if output != []:
            file_exists = os.path.isfile(file_name)
            with open(file_name, 'a', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=output.keys())

                if not file_exists:
                    writer.writeheader()  # file doesn't exist yet, write a header

                writer.writerow(output)
                
    except:
        print('Failed to appended result to CSV file.')

def read_csv_to_results(file_name='results/untitled.csv'):
    try:
        with open(file_name, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            return list(reader)

    except:
        print('Failed to open CSV file.')