__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import csv

def write_to_file(filename='default.csv', output=''):
    with open(filename, 'a', newline= '') as f:
        writer = csv.writer(f)
        writer.writerow(output)


def write_headings(filename='default.csv'):
    output = ['Obj function:',
              'Iterations:',
              'Refrigerant:',
              'dT Gas Cooler:',
              'T_air_tar:',
              'T_air_inlet:',
              'T_air_exh:',
              'n_comp:',
              'dT phx max 1:',
              'dT phx max 2:',
              'PR 2:',
              'dT phx max 3:',
              'PR 3:',
              'Message:',
              'Success:']
              
    write_to_file(filename, output, 'w')


def write_output(filename, result, indep_var, n_comp, print_to_screen=True, print_to_file=False):
    #Print to file
    if print_to_file:
        output = []
        output.append(result.fun)
        output.append(result.nfev)
        output.append(indep_var[0])
        output.append(indep_var[1])
        output.append(indep_var[2])
        output.append(indep_var[3])
        output.append(indep_var[4])
        output.append(indep_var[5])
        output.append(result.x[0])
        if n_comp > 1:
            output.append(result.x[1])
            output.append(result.x[2])
        else:
            output.append('')
            output.append('')
        if n_comp > 2:
            output.append(result.x[3])
            output.append(result.x[4])
        else:
            output.append('')
            output.append('')
        output.append(result.message)
        output.append(result.success)

        write_to_file(filename=filename,output=output)

    #Print to screen
    if print_to_screen:
        output2 = []
        output2.append(['Obj function:', result.fun])
        output2.append(['Iterations:', result.nfev])
        output2.append(['Refrigerant:', indep_var[0]])
        output2.append(['dT Gas Cooler:', indep_var[1]])
        output2.append(['T_air_tar:', indep_var[2]])
        output2.append(['T_air_inlet:', indep_var[3]])
        output2.append(['T_air_exh:', indep_var[4]])
        output2.append(['n_comp:', indep_var[5]])
        output2.append(['dT phx max 1:', result.x[0]])
        if n_comp == 2:
            output2.append(['dT phx max 2:', result.x[1]])
            output2.append(['PR 2:', result.x[2]])
        if n_comp == 3:
            output2.append(['dT phx max 3:', result.x[3]])
            output2.append(['PR 3:', result.x[4]])
        output2.append(['Message:', result.message])
        output2.append(['Success:', result.success])

        print(output2)


def write_screen(heat_pump):
        print('COP:', heat_pump.COP)
        print('Refrigerant:', heat_pump.Refrigerant)
        print('Comp 3 P out:', heat_pump.comp3.FlowOut.Pressure)
        print('Comp 2 P out:', heat_pump.comp2.FlowOut.Pressure)
        print('Comp 1 P out:', heat_pump.comp1.FlowOut.Pressure)
        print('Evap P:', heat_pump.evaporator.ColdIn.Pressure)
        print('PHX1:', heat_pump.dT_phx1_max)
        print('PHX2:', heat_pump.dT_phx2_max)
        print('PHX3:', heat_pump.dT_phx3_max)
        print('Gas cooler(s) inlet T:', heat_pump.gas_cooler4.HotIn.Temperature)
        print('Gas cooler(s) outlet T:', heat_pump.gas_cooler1.HotOut.Temperature)
        print('Evap T:', heat_pump.T_evap_init)
        print('Target air T:', heat_pump.T_air_tar)
        print('Actual air T:', heat_pump.gas_cooler4.ColdOut.Temperature)
        print('dTmin:', heat_pump.dT_GC_min_tar)
        print('Overall P ratio:', heat_pump.gas_cooler4.HotIn.Pressure / heat_pump.evaporator.ColdIn.Pressure)
        print('Dryer exh air T:', heat_pump.evaporator.HotIn.Temperature)
        print('Tol:', abs(heat_pump.T_air_tar - heat_pump.gas_cooler4.ColdOut.Temperature) / heat_pump.T_air_tar )
