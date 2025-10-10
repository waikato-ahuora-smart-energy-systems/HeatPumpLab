import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import pandas as pd
from pathlib import Path

__all__ = ['plot_cycle', 'plot_driving_force', 'plot_simple_xy', 'plot_heat_exchanger_profiles']

#######################################################################################################
# Public API
#######################################################################################################
    
def plot_cycle(cycle_obj, sink_obj, title='', x_axis='h', y_axis='T', res=50, refrigerant='default', include_sink=True):
    cycle_obj = cycle_obj + [cycle_obj[0]]
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots()

    if title != '':
        ax.set_title(title)
    
    xlabel, x_axis = _get_axis_labels(x_axis)
    ax.set_xlabel(xlabel)
    ylabel, y_axis = _get_axis_labels(y_axis)
    ax.set_ylabel(ylabel)

    if y_axis == 'p' or y_axis == 'pressure':
        ax.set_yscale('log')

    for obj in cycle_obj:
        data = obj[0].get_graph_data(divisions=res)
        x = _get_raw_data(data[obj[1]], x_axis)
        y = _get_raw_data(data[obj[1]], y_axis)
        if y_axis == 'p' or y_axis == 'pressure':
            y = [yi * 1000 for yi in y]

        if obj[0].name == 'Gas cooler':
            ax.plot(x, y, color='red')
        elif obj[0].name == 'Recuperator':
            ax.plot(x, y, color='green')
        elif obj[0].name == 'Internal heat exchanger':
            ax.plot(x, y, color='orange')
        elif obj[0].name == 'Evaporator':
            ax.plot(x, y, color='aqua')
        else:
            ax.plot(x, y, color='black')
    
    if (x_axis == 'entropy' or x_axis == 'enthalpy') and include_sink:

        for obj in reversed(sink_obj):
            if obj[1] == 'sink':
                data = obj[0].get_graph_data(divisions=res)
                x = _rescale(_get_raw_data(data['sink'],   x_axis), 
                                    _get_raw_data(data['source'], x_axis))
                y = _get_raw_data(data['sink'], y_axis)

                if obj[0].name == 'Gas cooler':
                    ax.plot(x, y, color='blue')
                else:
                    ax.plot(x, y, color='black')    

    fig.show()
    filename = 'Cycle--' + y_axis + '_vs_' + x_axis + '.png'
    output_dir = Path(__file__).resolve().parents[2] / "results" / "cycles" / refrigerant
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename)
    return fig, ax        

def plot_driving_force(hx_obj, title='', x_axis='h', y_axis='dT*', res=50, dtmin=0, refrigerant='default'):
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots()

    if title != '':
        ax.set_title(title)
    
    xlabel, x_axis = _get_axis_labels(x_axis)
    ax.set_xlabel(xlabel)
    ylabel, y_axis = _get_axis_labels(y_axis)
    ax.set_ylabel(ylabel)
    
    for obj in hx_obj:
        data = obj.get_graph_data(divisions=res)
        x = np.array(_get_raw_data(data['source'], x_axis))
        y1 = np.array(_get_raw_data(data['source'], y_axis)) - np.array(_get_raw_data(data['sink'], y_axis))[::-1]
        y2 = np.repeat(dtmin, len(x))

        ax.plot(x, y1, color='red')
        ax.plot(x, y2, color='black')
    
    fig.show()
    filename = 'Driving_force_plot.png'
    output_dir = Path(__file__).resolve().parents[2] / "results" / "cycles" / refrigerant
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename)
    return fig, ax           

def plot_heat_exchanger_profiles(hx_obj=[], title='', x_axis='h', y_axis='T', res=50, refrigerant='default'):
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots()

    if title != '':
        ax.set_title(title)
    
    xlabel, x_axis = _get_axis_labels(x_axis)
    ax.set_xlabel(xlabel)
    ylabel, y_axis = _get_axis_labels(y_axis)
    ax.set_ylabel(ylabel)
    
    for obj in hx_obj:
        data = obj.get_graph_data(divisions=res)
        x1 = np.array(_get_raw_data(data['source'], x_axis))
        x2 = np.array(_get_raw_data(data['sink'], x_axis)) 
        y1 = np.array(_get_raw_data(data['source'], y_axis))
        y2 = np.array(_get_raw_data(data['sink'], y_axis))

        # x1 = (x1[:] - x1[0]) * data['source'].mass_flow
        # x2 = (x2[:] - x2[0]) * data['source'].mass_flow
        
        ax.plot(x1, y1, color='red')
        ax.plot(x2, y2, color='blue')
    
    fig.show()
    filename = 'Heat_exchanger_profiles.png'
    output_dir = Path(__file__).resolve().parents[2] / "results" / "cycles" / refrigerant
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename)   
    return fig, ax           

def plot_simple_xy(x: list, y: list, title: str ='', x_axis: str ='', y_axis: str ='', refrigerant='default'):
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots()

    if title != '':
        ax.set_title(title)
    if x_axis != '':
        ax.set_xlabel(x_axis)
    if y_axis != '':
        ax.set_ylabel(y_axis)
    ax.plot(x, y, color='black')

    fig.show()
    filename = 'Cycle--' + y_axis + '_vs_' + x_axis + '.png'
    output_dir = Path(__file__).resolve().parents[2] / "results" / "cycles" / refrigerant
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename)
    return fig, ax


def plot_heat_map(title, x_val, x_label, y_val, y_label, z_val, z_label):
    df = pd.DataFrame({"x": x_val, "y": y_val, "z": z_val})
    fig = px.density_heatmap(df, x="x", y="y", z="z", nbinsx=20, nbinsy=20, color_continuous_scale="Viridis")
    fig.update_layout(
        title=title,
        xaxis_title="X-axis",
        yaxis_title="Y-axis"
    )
    fig.show()


#######################################################################################################
# Helper functions
#######################################################################################################

def _rescale(array_to_scale, array_new_scale) -> list[float]:
    m: float
    c: float
    if abs(array_to_scale[0] - array_to_scale[-1]) > 0.0:
        m = (array_new_scale[-1] - array_new_scale[0]) / (array_to_scale[0] - array_to_scale[-1])
    else:
        m = 0.0
    c = array_new_scale[0]
    return [(x - array_to_scale[-1]) * m + c for x in array_to_scale]

def _get_raw_data(ls, prop):
    x = []
    for s in ls:
        x += [getattr(s, prop)]
    return x

def _get_axis_labels(axis_var: str) -> list[str]: # type: ignore
    axis = axis_var.lower()
    if axis == 'h' or axis == 'enthalpy':
        label = 'Enthalpy (kJ/kg)'
        axis = 'enthalpy'
    elif axis == 's' or axis == 'entropy':
        label = 'Entropy (kJ/kgK)'
        axis = 'entropy'
    elif axis == 'p' or axis == 'pressure':
        label = 'Pressure (kPa)'
        axis = 'pressure'
    elif axis == 't' or axis == 'temperature':
        label = 'Temperature (°C)'
        axis = 'temperature'
    elif axis == 'v' or axis == 'specific volume':
        label = r'Specific Volume (m$^3$/kg)'
        axis = 'specific_volume'
    elif axis == 'dT*':
        label = axis + ' (°C)'
        axis = 'temperature'
    elif axis == 'heat flow':
        label = 'Heat Flow' + ' (kW)'
        axis = 'heat_flow'        
    else:
        label = 'Temperature (°C)'
        axis = 'temperature'            
    return label, axis
            