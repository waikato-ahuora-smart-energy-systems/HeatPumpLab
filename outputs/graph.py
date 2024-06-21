__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import matplotlib.pyplot as plt

class Graph:
    def __init__(self) -> None:
        pass
    
    def rescale(self, array_to_scale, array_new_scale):
        if abs(array_to_scale[0] - array_to_scale[-1]) > 0.0:
            m = (array_new_scale[-1] - array_new_scale[0]) / (array_to_scale[0] - array_to_scale[-1])
        else:
            m = 0.0
        
        c = array_new_scale[0]
        
        return [(x - array_to_scale[-1]) * m + c for x in array_to_scale]

    def cycle(self, cycle_obj, sink_obj, title='', x_axis='h', y_axis='T', res=100):
        cycle_obj = cycle_obj + [cycle_obj[0]]

        plt
        if x_axis == 'h' or x_axis == 'enthalpy':
            plt.xlabel('Enthalpy [kJ/kg]')
        elif x_axis == 's' or x_axis == 'S' or x_axis == 'entropy':
            plt.xlabel('Entropy [kJ/kgK]')
        elif x_axis == 'p' or x_axis == 'P' or x_axis == 'pressure':
            plt.xscale('log')
            plt.xlabel('Pressure [kPa]')
        elif x_axis == 't' or x_axis == 'T' or x_axis == 'temperature':
            plt.xlabel('Temperature [°C]')
        else:
            plt.xlabel('Specific Volume [m^3/kg]')

        if y_axis == 't' or y_axis == 'T' or y_axis == 'temperature':
            plt.ylabel('Temperature [°C]')
        else:
            plt.yscale('log')
            plt.ylabel('Pressure [kPa]') 

        for obj in cycle_obj:
            data = obj[0].get_graph_data(divisions=res)
            if x_axis == 'h' or x_axis == 'enthalpy':
                x = data[obj[1]]['h'] 
            elif x_axis == 's' or x_axis == 'S' or x_axis == 'entropy':
                x = data[obj[1]]['s']
            elif x_axis == 'p' or x_axis == 'P' or x_axis == 'pressure':
                x = data[obj[1]]['P']
            elif x_axis == 't' or x_axis == 'T' or x_axis == 'temperature':
                x = data[obj[1]]['T'] 
            else:
                x = data[obj[1]]['v']

            if y_axis == 't' or y_axis == 'T' or y_axis == 'temperature':
                y = data[obj[1]]['T']
            else:
                y = data[obj[1]]['P']

            if obj[0].name == 'Gas cooler':
                plt.plot(x, y, color='red')
            elif obj[0].name == 'Plate heat exchanger':
                plt.plot(x, y, color='green')
            elif obj[0].name == 'Evaporator':
                plt.plot(x, y, color='orange')
            else:
                plt.plot(x, y, color='black')
        
        if x_axis == 's' or x_axis == 'S' or x_axis == 'entropy' or x_axis == 'h' or x_axis == 'enthalpy':
            for obj in reversed(sink_obj):
                if obj[1] == 'sink':
                    data = obj[0].get_graph_data(divisions=res)
                    data[obj[1]]['h'] = self.rescale(data[obj[1]]['h'], data['source']['h'])
                    data[obj[1]]['s'] = self.rescale(data[obj[1]]['s'], data['source']['s'])

                    if x_axis == 'h' or x_axis == 'enthalpy':
                        x = data[obj[1]]['h'] 
                    elif x_axis == 's' or x_axis == 'S' or x_axis == 'entropy':
                        x = data[obj[1]]['s']
                    elif x_axis == 'p' or x_axis == 'P' or x_axis == 'pressure':
                        x = data[obj[1]]['P']
                    elif x_axis == 't' or x_axis == 'T' or x_axis == 'temperature':
                        x = data[obj[1]]['T'] 
                    else:
                        x = data[obj[1]]['v']

                    if y_axis == 't' or y_axis == 'T' or y_axis == 'temperature':
                        y = data[obj[1]]['T']
                    else:
                        y = data[obj[1]]['P']

                    plt.plot(x, y, color='blue')            

        if x_axis == 'h' or x_axis == 'enthalpy':
            plt.xlabel('Enthalpy [kJ/kg]')
        elif x_axis == 's' or x_axis == 'S' or x_axis == 'entropy':
            plt.xlabel('Entropy [kJ/kgK]')
        elif x_axis == 'p' or x_axis == 'P' or x_axis == 'pressure':
            plt.xscale('log')
            plt.xlabel('Pressure [kPa]')
        elif x_axis == 't' or x_axis == 'T' or x_axis == 'temperature':
            plt.xlabel('Temperature [°C]')
        else:
            plt.xlabel('Specific Volume [m^3/kg]')

        if y_axis == 't' or y_axis == 'T' or y_axis == 'temperature':
            plt.ylabel('Temperature [°C]')
        else:
            plt.yscale('log')
            plt.ylabel('Pressure [kPa]') 
        
        if title != '':
            plt.title(title)

        plt.show()           