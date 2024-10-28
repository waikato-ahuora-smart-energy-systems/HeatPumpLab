__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from matplotlib import pyplot
import numpy as np

class Graph:
    
    def __init__(self) -> None:
        pass
    
    def rescale(self, array_to_scale, array_new_scale) -> list[float]:
        m: float
        c: float
        if abs(array_to_scale[0] - array_to_scale[-1]) > 0.0:
            m = (array_new_scale[-1] - array_new_scale[0]) / (array_to_scale[0] - array_to_scale[-1])
        else:
            m = 0.0
        c = array_new_scale[0]
        return [(x - array_to_scale[-1]) * m + c for x in array_to_scale]

    def get_raw_data(self, ls, prop):
        x = []
        for s in ls:
            x += [getattr(s, prop)]
        return x

    def cycle(self, cycle_obj, sink_obj, title='', x_axis='h', y_axis='T', res=100):
        cycle_obj = cycle_obj + [cycle_obj[0]]

        plt = pyplot
        if title != '':
            plt.title(title)
        
        plt, x_axis = self._add_xaxis_labels(plt, x_axis)
        plt, y_axis = self._add_yaxis_labels(plt, y_axis)

        for obj in cycle_obj:
            data = obj[0].get_graph_data(divisions=res)
            x = self.get_raw_data(data[obj[1]], x_axis)
            y = self.get_raw_data(data[obj[1]], y_axis)

            if obj[0].name == 'Gas cooler':
                plt.plot(x, y, color='red')
            elif obj[0].name == 'Recuperator':
                plt.plot(x, y, color='green')
            elif obj[0].name == 'Internal heat exchanger':
                plt.plot(x, y, color='orange')
            elif obj[0].name == 'Evaporator':
                plt.plot(x, y, color='aqua')
            else:
                plt.plot(x, y, color='black')
        
        if x_axis == 'entropy' or x_axis == 'enthalpy':
            for obj in reversed(sink_obj):
                if obj[1] == 'sink':
                    data = obj[0].get_graph_data(divisions=res)
                    x = self.rescale(self.get_raw_data(data['sink'],   x_axis), 
                                     self.get_raw_data(data['source'], x_axis))
                    y = self.get_raw_data(data['sink'], y_axis)

                    if obj[0].name == 'Gas cooler':
                        plt.plot(x, y, color='blue')
                    else:
                        plt.plot(x, y, color='black')    

        plt.show()  
        del plt         

    def driving_force(self, hx_obj, title='', x_axis='h', y_axis='dT*', res=100, dtmin=0):
        plt = pyplot
        if title != '':
            plt.title(title)
        
        plt, x_axis = self._add_xaxis_labels(plt, x_axis)
        plt, y_axis = self._add_yaxis_labels(plt, y_axis)

        for obj in hx_obj:
            data = obj.get_graph_data(divisions=res)
            x = np.array(self.get_raw_data(data['source'], x_axis))
            y1 = np.array(self.get_raw_data(data['source'], y_axis)) - np.array(self.get_raw_data(data['sink'], y_axis))[::-1]
            y2 = np.repeat(dtmin, len(x))

            plt.plot(x, y1, color='red')
            plt.plot(x, y2, color='black')
        
        plt.show()
        del plt           

    def simple(self, x: list, y: list, title: str ='', x_axis: str ='', y_axis: str =''):
        plt = pyplot
        if title != '':
            plt.title(title)
        if x_axis != '':
            plt.xlabel(x_axis)
        if y_axis != '':
            plt.ylabel(y_axis)
        plt.plot(x, y, color='black')
        plt.show()
        del plt

    def _add_xaxis_labels(self, plt: pyplot, x_axis: str) -> list[pyplot, str]: # type: ignore
        if x_axis == 'h' or x_axis == 'enthalpy':
            plt.xlabel('Enthalpy [kJ/kg]')
            x_axis = 'enthalpy'
        elif x_axis == 's' or x_axis == 'S' or x_axis == 'entropy':
            plt.xlabel('Entropy [kJ/kgK]')
            x_axis = 'entropy'
        elif x_axis == 'p' or x_axis == 'P' or x_axis == 'pressure':
            plt.xscale('log')
            plt.xlabel('Pressure [kPa]')
            x_axis = 'pressure'
        elif x_axis == 't' or x_axis == 'T' or x_axis == 'temperature':
            plt.xlabel('Temperature [°C]')
            x_axis == 'temperature'
        else:
            plt.xlabel('Specific Volume [m^3/kg]')
            x_axis == 'specific_volume'
        return plt, x_axis

    def _add_yaxis_labels(self, plt: pyplot, y_axis: str) -> list[pyplot, str]: # type: ignore
        if y_axis == 't' or y_axis == 'T' or y_axis == 'temperature' or y_axis == 'Temperature':
            plt.ylabel('Temperature [°C]')
            y_axis = 'temperature'
        elif y_axis == 'p' or y_axis == 'P' or y_axis == 'pressure' or y_axis == 'Pressure':
            plt.yscale('log')
            plt.ylabel('Pressure [kPa]') 
            y_axis = 'pressure'
        elif y_axis == 'dT*':
            plt.ylabel(y_axis)
            y_axis = 'temperature'
        else:
            plt.ylabel('Temperature [°C]')
            y_axis = 'temperature'            
        return plt, y_axis
              