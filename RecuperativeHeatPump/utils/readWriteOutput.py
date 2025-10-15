"""Helpers for printing and persisting heat pump case results."""

import os
from glob import glob
import csv
# from csv import DictReader, DictWriter
from typing import Any

__all__ = ['print_result_to_terminal', 'print_result_to_csv', 'transfer_temp_files_to_results', 'read_csv_to_results']

#######################################################################################################
# Public API
#######################################################################################################

def print_result_to_terminal(output: dict[str, any]):
    """Pretty-print a result dictionary to the terminal."""
    key: str
    value: any

    if type(output) is dict:
        for key, value in output.items():
            print(f'{key} {value}')
    else:
        print('Failed to print due to wrong output format.')

def print_result_to_csv(output: list[any] =[], file_name: str ='R', file_ext: str ='.csv', suffix: str =''):
    """Append a result dictionary to a CSV file, creating it if needed."""
    try:
        if output != []:
            file_name = 'results/' + file_name + suffix + file_ext
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

def transfer_temp_files_to_results(file_name: str ='untitled') -> None:
    """Move temporary CSV files into the main results file."""
    csv_file_list: glob
    f: str
    fout: list[any]

    csv_file_list = glob(os.getcwd() + '/results/temp/*.csv')
    header_saved = False
    with open('results/' + file_name + '.csv','a') as fout:
        for f in reversed(csv_file_list):
            with open(f) as fin:
                header = next(fin)
                if not header_saved:
                    fout.write(header)
                    header_saved = True
                for line in fin:
                    fout.write(line)
            os.remove(f)

def read_csv_to_results(file_name: str = 'untitled', file_ext: str = '.csv') -> list[dict[str, Any]]:
    """Reads a CSV file from the results directory into a list of dictionaries."""
    file_path = os.path.join('results', file_name + file_ext)

    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return []

    try:
        with open(file_path, mode='r', newline='', encoding='utf-8') as csvfile:
            reader = csv.DictReader(csvfile)
            return list(reader)
    except (OSError, csv.Error) as e:
        print(f"Failed to read CSV file '{file_path}': {e}")
        return []
