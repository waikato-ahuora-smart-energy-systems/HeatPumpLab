__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import matplotlib.pyplot as plt

class Graph:
    def __init__(self) -> None:
        pass

    def cycle(self, graph_obj, title='', x_axis='h', y_axis='T', res=50):
        h = []; T = []; P = []; s = []; v = []
        cycle_obj = graph_obj + [graph_obj[0]]
        for obj in cycle_obj:
            data = obj[0].get_graph_data()
            P += data[obj[1]]['P']
            h += data[obj[1]]['h']  
            T += data[obj[1]]['T']  
            s += data[obj[1]]['s']  
            v += data[obj[1]]['v']
        
        plt
        if x_axis == 'h' or x_axis == 'enthalpy':
            x = h
            plt.xlabel('Enthalpy [kJ/kg]')
            
        elif x_axis == 's' or x_axis == 'S' or x_axis == 'entropy':
            x = s
            plt.xlabel('Entropy [kJ/kgK]')
        elif x_axis == 'p' or x_axis == 'P' or x_axis == 'pressure':
            x = P
            plt.xscale('log')
            plt.xlabel('Pressure [kPa]')
        elif x_axis == 't' or x_axis == 'T' or x_axis == 'temperature':
            x = T
            plt.xlabel('Temperature [°C]')
        else:
            x = v
            plt.xlabel('Specific Volume [m^3/kg]')

        if y_axis == 't' or y_axis == 'T' or y_axis == 'temperature':
            y = T
            plt.ylabel('Temperature [°C]')
        else:
            y = P
            plt.yscale('log')
            plt.ylabel('Pressure [kPa]')

        plt.plot(x, y, color='black')    

        if title != '':
            plt.title(title)

        plt.show() 

    # def heat_exchangers(self, graph_obj, title=''):
    #     H = []; T_hot = []; T_cold = []
    #     for obj in graph_obj:
    #         H_i, T_hot_i, T_cold_i = obj.get_graph_data()
    #         if len(H) > 0:
    #             j = len(H) - 1
    #             H_i = [H_i[i] + H[j] for i in range(len(H_i))]
    #             H += H_i 
    #         else:
    #             H += H_i

    #         T_hot += T_hot_i
    #         T_cold += T_cold_i

    #     plt.plot(H, T_hot, color='red', label = 'Hot stream')
    #     plt.plot(H, T_cold, color='blue', label = 'Cold stream')

    #     plt.xlabel('Heat Flow [kJ/kg(ref.)]')
    #     plt.ylabel('Temperature [°C]')
    #     if title != '':
    #         plt.title(title)
    #     plt.show()             