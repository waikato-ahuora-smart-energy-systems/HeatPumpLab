__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

import matplotlib.pyplot as plt

def graph_heat_exchangers(graph_obj):

        H = []; T_hot = []; T_cold = []
        for obj in graph_obj:
            H_i, T_hot_i, T_cold_i = obj.graph_data()
            if len(H) > 0:
                j = len(H) - 1
                H_i = [H_i[i] + H[j] for i in range(len(H_i))]
                H += H_i 
            else:
                H += H_i

            T_hot += T_hot_i
            T_cold += T_cold_i

        plt.plot(H, T_hot, color='red', label = 'Hot stream')
        plt.plot(H, T_cold, color='blue', label = 'Cold stream')

        plt.xlabel('Heat Flow [kJ/kg(ref.)]')
        plt.ylabel('Temperature [°C]')
        #plt.title('Gas Cooler')
        plt.show()   