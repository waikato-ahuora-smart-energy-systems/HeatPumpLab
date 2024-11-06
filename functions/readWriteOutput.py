__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os
from glob import glob
from csv import DictReader, DictWriter


def print_result_to_terminal(output: dict[str, any]):
    key: str
    value: any

    if type(output) is dict:
        for key, value in output.items():
            print(f'{key} {value}')
    else:
        print('Failed to print due to wrong output format.')

def print_result_to_csv(output: list[any] =[], file_name: str ='results/temp/R', file_ext: str ='.csv', suffix: str =''):
    try:
        if output != []:
            file_name = file_name + suffix + file_ext
            file_exists = os.path.isfile(file_name)
            if file_exists:
                with open(file_name, 'a', newline='') as csvfile:
                    writer = DictWriter(csvfile, fieldnames=output.keys())
                    writer.writerow(output)

            else:
                with open(file_name, 'w', newline='') as csvfile:
                    writer = DictWriter(csvfile, fieldnames=output.keys())
                    writer.writeheader()  
                    writer.writerow(output)  
                
    except:
        print('Failed to appended result to CSV file.')

def transfer_temp_files_to_results(file_name: str ='results/untitled') -> None:
    csv_file_list: glob
    f: str
    output: list[any]
    o: list[any]

    csv_file_list = glob(os.getcwd() + '/results/temp/*.csv')
    for f in csv_file_list:
        output = read_csv_to_results(file_name=f, file_ext='')
        for o in output:
            print_result_to_csv(output=o, file_name=file_name)
        os.remove(f)

def read_csv_to_results(file_name: str ='results/untitled', file_ext: str ='.csv') -> list[any]:
    reader: DictReader
    try:
        with open(file_name + file_ext, 'r') as csvfile:
            reader = DictReader(csvfile)
            return list(reader)
    except:
        print('Failed to open CSV file.')


