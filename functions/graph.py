__author__ = 'Tim Walmsley and Keri Adamson'
__credits__ = ['Alex Geary']

from matplotlib import pyplot
import numpy as np

class Graph:
    
    def __init__(self) -> None:
        pass
    
    def plot_cycle(self, cycle_obj, sink_obj, title='', x_axis='h', y_axis='T', res=100):
        cycle_obj = cycle_obj + [cycle_obj[0]]

        plt = pyplot
        if title != '':
            plt.title(title)
        
        xlabel, x_axis, plt = self._get_axis_labels(plt, x_axis)
        plt.xlabel(xlabel)
        ylabel, y_axis, plt = self._get_axis_labels(plt, y_axis)
        plt.ylabel(ylabel)

        for obj in cycle_obj:
            data = obj[0].get_graph_data(divisions=res)
            x = self._get_raw_data(data[obj[1]], x_axis)
            y = self._get_raw_data(data[obj[1]], y_axis)

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
                    x = self._rescale(self._get_raw_data(data['sink'],   x_axis), 
                                     self._get_raw_data(data['source'], x_axis))
                    y = self._get_raw_data(data['sink'], y_axis)

                    if obj[0].name == 'Gas cooler':
                        plt.plot(x, y, color='blue')
                    else:
                        plt.plot(x, y, color='black')    

        plt.show()
        filename = 'Heat pump - ' + y_axis + ' vs ' + x_axis + '.png'
        plt.savefig('results/figures/' + filename) 
        del plt         

    def plot_driving_force(self, hx_obj, title='', x_axis='h', y_axis='dT*', res=100, dtmin=0):
        plt = pyplot
        if title != '':
            plt.title(title)
        
        xlabel, x_axis, plt = self._get_axis_labels(plt, x_axis)
        plt.xlabel(xlabel)
        ylabel, y_axis, plt = self._get_axis_labels(plt, y_axis)
        plt.ylabel(ylabel)
        
        for obj in hx_obj:
            data = obj.get_graph_data(divisions=res)
            x = np.array(self._get_raw_data(data['source'], x_axis))
            y1 = np.array(self._get_raw_data(data['source'], y_axis)) - np.array(self._get_raw_data(data['sink'], y_axis))[::-1]
            y2 = np.repeat(dtmin, len(x))

            plt.plot(x, y1, color='red')
            plt.plot(x, y2, color='black')
        
        plt.show()
        del plt           

    def plot_simple_xy(self, x: list, y: list, title: str ='', x_axis: str ='', y_axis: str =''):
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
 
    def _rescale(self, array_to_scale, array_new_scale) -> list[float]:
        m: float
        c: float
        if abs(array_to_scale[0] - array_to_scale[-1]) > 0.0:
            m = (array_new_scale[-1] - array_new_scale[0]) / (array_to_scale[0] - array_to_scale[-1])
        else:
            m = 0.0
        c = array_new_scale[0]
        return [(x - array_to_scale[-1]) * m + c for x in array_to_scale]

    def _get_raw_data(self, ls, prop):
        x = []
        for s in ls:
            x += [getattr(s, prop)]
        return x

    def _get_axis_labels(self, plt: pyplot, axis: str) -> list[pyplot, str]: # type: ignore
        axis = axis.lower()
        if axis == 'h' or axis == 'enthalpy':
            label = 'Enthalpy [kJ/kg]'
            axis = 'enthalpy'
        elif axis == 's' or axis == 'entropy':
            label = 'Entropy [kJ/kgK]'
            axis = 'entropy'
        elif axis == 'p' or axis == 'pressure':
            plt.xscale('log')
            label = 'Pressure [kPa]'
            axis = 'pressure'
        elif axis == 't' or axis == 'temperature':
            label = 'Temperature [°C]'
            axis = 'temperature'
        elif axis == 'v' or axis == 'specific volume':
            label = 'Specific Volume [m^3/kg]'
            axis = 'specific_volume'
        elif axis == 'dT*':
            label = axis + ' [°C]'
            axis = 'temperature'
        else:
            label = 'Temperature [°C]'
            axis = 'temperature'            
        return label, axis, plt
              