__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os, csv

def print_result_to_terminal(output):
    if type(output) is dict:
        for key, value in output.items():
            print(f'{key} {value}')

    else:
        print('Failed to print due to wrong output format.')

def print_result_to_csv(output=[], file_name='results/temp/R', file_ext='.csv', suffix=''):
    try:
        if output != []:
            file_name = file_name + suffix + file_ext
            file_exists = os.path.isfile(file_name)
            if file_exists:
                with open(file_name, 'a', newline='') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=output.keys())
                    writer.writerow(output)

            else:
                with open(file_name, 'w', newline='') as csvfile:
                    writer = csv.DictWriter(csvfile, fieldnames=output.keys())
                    writer.writeheader()  
                    writer.writerow(output)  
                
    except:
        print('Failed to appended result to CSV file.')
  