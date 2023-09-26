__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import os, csv, glob

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
                    writer.writeheader()  # file doesn't exist yet, write a header
                    writer.writerow(output)  
                
    except:
        print('Failed to appended result to CSV file.')

def read_csv_to_results(file_name='results/untitled', file_ext='.csv'):
    try:
        with open(file_name + file_ext, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            return list(reader)

    except:
        print('Failed to open CSV file.')

def transfer_results(file_name='results/untitled'):
        csv_file_list = glob.glob(os.getcwd() + '/results/temp/*.csv')
        for f in csv_file_list:
            output = read_csv_to_results(file_name=f, file_ext='')
            print_result_to_csv(output=output[0], file_name=file_name)
            os.remove(f)
   