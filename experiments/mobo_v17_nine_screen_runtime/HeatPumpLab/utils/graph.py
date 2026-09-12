"""Plotting helpers for heat pump cycle outputs."""

import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
import pandas as pd
from pathlib import Path

__all__ = [
    "plot_cycle",
    "plot_driving_force",
    "plot_simple_xy",
    "plot_heat_exchanger_profiles",
]

DEFAULT_FIGSIZE: tuple[float, float] = (4.5, 3.75)
DEFAULT_DPI: int = 900

# Allows customization of legend names on a per-component basis.
LEGEND_NAME_MAP: dict[str, str] = {
    "Gas cooler": "GC",
    "Recuperator": "RHX",
    "Internal heat exchanger": "IHX",
    "Evaporator": "EVAP",
    "Compressor": "COMP",
    "Expansion valve": "VAL",
    "Economiser": "Economiser",
    "Plate heat exchanger": "Plate heat exchanger",
    "Heat exchanger": "Heat exchanger",
    "Mixer": "Mixer",
    "Splitter": "Splitter",
    "Turbine": "Turbine",
    "Sink - GC": "Sink",
}

#######################################################################################################
# Public API
#######################################################################################################


def plot_cycle(
    cycle_obj,
    sink_obj,
    title="",
    x_axis="h",
    y_axis="T",
    res=50,
    refrigerant="default",
    include_sink=True,
    show: bool = True,
    output_dir: str | Path | None = None,
):
    """Plot a cycle diagram with optional sink curves and save the figure."""
    cycle_obj = cycle_obj + [cycle_obj[0]]
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    fig.set_dpi(DEFAULT_DPI)
    plotted_labels: set[str] = set()
    sink_labels: set[str] = set()

    if title != "":
        ax.set_title(title)

    xlabel, x_axis = _get_axis_labels(x_axis)
    ax.set_xlabel(xlabel)
    ylabel, y_axis = _get_axis_labels(y_axis)
    ax.set_ylabel(ylabel)

    if y_axis == "p" or y_axis == "pressure":
        ax.set_yscale("log")

    ax.grid(True, which="major", linestyle="--", linewidth=0.5, color="0.85")

    for obj in cycle_obj:
        data = obj[0].get_graph_data(divisions=res)
        x = _get_raw_data(data[obj[1]], x_axis)
        y = _get_raw_data(data[obj[1]], y_axis)
        if y_axis == "p" or y_axis == "pressure":
            y = [yi * 1000 for yi in y]

        component_name = getattr(obj[0], "name", "Component")
        if component_name is None:
            component_name = "Component"
        component_name = str(component_name)
        color = "black"
        if component_name == "Gas cooler":
            color = "red"
        elif component_name == "Recuperator":
            color = "green"
        elif component_name == "Internal heat exchanger":
            color = "orange"
        elif component_name == "Evaporator":
            color = "aqua"
        elif component_name == "Expansion valve":
            color = "grey"

        display_name = _map_legend_name(component_name)
        label = display_name if display_name not in plotted_labels else None
        if label is not None:
            plotted_labels.add(display_name)
        ax.plot(x, y, color=color, label=label)

    if (x_axis == "entropy" or x_axis == "enthalpy") and include_sink:
        for obj in reversed(sink_obj):
            if obj[1] == "sink":
                data = obj[0].get_graph_data(divisions=res)
                x = _rescale(
                    _get_raw_data(data["sink"], x_axis),
                    _get_raw_data(data["source"], x_axis),
                )
                y = _get_raw_data(data["sink"], y_axis)

                if obj[0].name == "Gas cooler":
                    color = "blue"
                else:
                    color = "black"

                sink_component = getattr(obj[0], "name", "Component")
                if sink_component is None:
                    sink_component = "Component"
                sink_component = str(sink_component)
                sink_label = "AIR*"
                if sink_label in sink_labels:
                    sink_label = None
                else:
                    sink_labels.add(sink_label)

                ax.plot(x, y, color=color, label=sink_label)

    _add_legend(ax)
    _apply_layout(fig)
    if show:
        fig.show()
    filename = "Cycle--" + y_axis + "_vs_" + x_axis + ".png"
    output_dir = _get_output_dir(refrigerant, output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename, bbox_inches="tight", dpi=DEFAULT_DPI)
    return fig, ax


def plot_driving_force(
    hx_obj,
    title="",
    x_axis="h",
    y_axis="dT*",
    res=50,
    dtmin=0,
    refrigerant="default",
    show: bool = True,
    output_dir: str | Path | None = None,
):
    """Plot driving-force curves for a list of heat exchangers."""
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    fig.set_dpi(DEFAULT_DPI)
    plotted_labels: set[str] = set()
    dtmin_label_added = False

    if title != "":
        ax.set_title(title)

    xlabel, x_axis = _get_axis_labels(x_axis)
    ax.set_xlabel(xlabel)
    ylabel, y_axis = _get_axis_labels(y_axis)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", linestyle="--", linewidth=0.5, color="0.85")

    for obj in hx_obj:
        data = obj.get_graph_data(divisions=res)
        x = np.array(_get_raw_data(data["source"], x_axis))
        y1 = (
            np.array(_get_raw_data(data["source"], y_axis))
            - np.array(_get_raw_data(data["sink"], y_axis))[::-1]
        )
        y2 = np.repeat(dtmin, len(x))

        component_label = getattr(obj, "name", None)
        if component_label:
            component_label = str(component_label)
        mapped_label = _map_legend_name(component_label) if component_label else None
        label = (
            mapped_label
            if mapped_label and mapped_label not in plotted_labels
            else None
        )
        if label is not None:
            plotted_labels.add(mapped_label)
        ax.plot(x, y1, color="red", label=label)

        dtmin_label = r"$ΔT_{min}$"
        if dtmin_label_added:
            dtmin_label = None
        else:
            dtmin_label_added = True
        ax.plot(x, y2, color="black", label=dtmin_label)

    _add_legend(ax)
    _apply_layout(fig)
    if show:
        fig.show()
    filename = "Driving_force_plot.png"
    output_dir = _get_output_dir(refrigerant, output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename, bbox_inches="tight", dpi=DEFAULT_DPI)
    return fig, ax


def plot_heat_exchanger_profiles(
    hx_obj=None,
    title="",
    x_axis="h",
    y_axis="T",
    res=50,
    refrigerant="default",
    show: bool = True,
    output_dir: str | Path | None = None,
):
    """Plot source and sink temperature profiles for heat exchangers."""
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    fig.set_dpi(DEFAULT_DPI)
    source_labels: set[str] = set()
    sink_labels: set[str] = set()
    if hx_obj is None:
        hx_obj = []

    if title != "":
        ax.set_title(title)

    xlabel, x_axis = _get_axis_labels(x_axis)
    ax.set_xlabel(xlabel)
    ylabel, y_axis = _get_axis_labels(y_axis)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", linestyle="--", linewidth=0.5, color="0.85")

    for obj in hx_obj:
        data = obj.get_graph_data(divisions=res)
        x1 = np.array(_get_raw_data(data["source"], x_axis))
        x2 = np.array(_get_raw_data(data["sink"], x_axis))
        y1 = np.array(_get_raw_data(data["source"], y_axis))
        y2 = np.array(_get_raw_data(data["sink"], y_axis))

        # x1 = (x1[:] - x1[0]) * data['source'].mass_flow
        # x2 = (x2[:] - x2[0]) * data['source'].mass_flow

        component_name = getattr(obj, "name", None)
        if component_name:
            component_name = str(component_name)
        mapped_component_name = (
            _map_legend_name(component_name) if component_name else None
        )
        base_source_label = (
            f"Source - {mapped_component_name}" if mapped_component_name else "Source"
        )
        source_label = (
            base_source_label if base_source_label not in source_labels else None
        )
        if source_label is not None:
            source_labels.add(base_source_label)

        base_sink_label = (
            f"Sink - {mapped_component_name}" if mapped_component_name else "Sink"
        )
        sink_label = base_sink_label if base_sink_label not in sink_labels else None
        if sink_label is not None:
            sink_labels.add(base_sink_label)

        ax.plot(x1, y1, color="red", label=source_label)
        ax.plot(x2, y2, color="blue", label=sink_label)

    _add_legend(ax)
    _apply_layout(fig)
    if show:
        fig.show()
    filename = "Heat_exchanger_profiles.png"
    output_dir = _get_output_dir(refrigerant, output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename, bbox_inches="tight", dpi=DEFAULT_DPI)
    return fig, ax


def plot_simple_xy(
    x: list,
    y: list,
    title: str = "",
    x_axis: str = "",
    y_axis: str = "",
    refrigerant="default",
    show: bool = True,
    output_dir: str | Path | None = None,
):
    """Plot a simple 2D curve and write it to the figures directory."""
    fig: plt.Figure
    ax: plt.Axes

    fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)
    fig.set_dpi(DEFAULT_DPI)

    if title != "":
        ax.set_title(title)
    if x_axis != "":
        ax.set_xlabel(x_axis)
    if y_axis != "":
        ax.set_ylabel(y_axis)
    ax.grid(True, which="major", linestyle="--", linewidth=0.5, color="0.85")
    label = title if title else "Series"
    ax.plot(x, y, color="black", label=label)

    _add_legend(ax)
    _apply_layout(fig)
    if show:
        fig.show()
    filename = "Cycle--" + y_axis + "_vs_" + x_axis + ".png"
    output_dir = _get_output_dir(refrigerant, output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / filename, bbox_inches="tight", dpi=DEFAULT_DPI)
    return fig, ax


def plot_heat_map(title, x_val, x_label, y_val, y_label, z_val, z_label, show=True):
    """Render a heat map using Plotly and return the figure."""
    df = pd.DataFrame({"x": x_val, "y": y_val, "z": z_val})
    fig = px.density_heatmap(
        df, x="x", y="y", z="z", nbinsx=20, nbinsy=20, color_continuous_scale="Viridis"
    )
    fig.update_layout(title=title, xaxis_title=x_label, yaxis_title=y_label)
    fig.update_coloraxes(colorbar_title=z_label)
    if show:
        fig.show()
    return fig


#######################################################################################################
# Helper functions
#######################################################################################################


def _add_legend(ax: plt.Axes) -> None:
    """Add a legend if the axes currently have labelled artists."""
    _, labels = ax.get_legend_handles_labels()
    if any(label for label in labels):
        ax.legend()


def _map_legend_name(name: str | None) -> str:
    """Convert component names to shorter legend-friendly labels."""
    if not name:
        return "Component"
    return LEGEND_NAME_MAP.get(name, name)


def _apply_layout(fig: plt.Figure) -> None:
    """Apply tight layout for consistent sizing."""
    try:
        fig.tight_layout(pad=1.4)
    except Exception:
        pass


def _get_output_dir(refrigerant: str, output_dir: str | Path | None) -> Path:
    """Return the directory used for figure output."""
    if output_dir is not None:
        return Path(output_dir)
    return (
        Path(__file__).resolve().parents[2]
        / "results"
        / "figures"
        / "cycles"
        / refrigerant
    )


def _rescale(array_to_scale, array_new_scale) -> list[float]:
    """Rescale an array to match the span of another array."""
    m: float
    c: float
    if abs(array_to_scale[0] - array_to_scale[-1]) > 0.0:
        m = (array_new_scale[-1] - array_new_scale[0]) / (
            array_to_scale[0] - array_to_scale[-1]
        )
    else:
        m = 0.0
    c = array_new_scale[0]
    return [(x - array_to_scale[-1]) * m + c for x in array_to_scale]


def _get_raw_data(ls, prop):
    """Extract a list of property values from stream objects."""
    x = []
    for s in ls:
        x += [getattr(s, prop)]
    return x


def _get_axis_labels(axis_var: str) -> list[str]:  # type: ignore
    """Return axis label text and canonical property key."""
    axis = axis_var.lower()
    if axis == "h" or axis == "enthalpy":
        label = "Enthalpy [kJ/kg]"
        axis = "enthalpy"
    elif axis == "s" or axis == "entropy":
        label = "Entropy [kJ/kgK]"
        axis = "entropy"
    elif axis == "p" or axis == "pressure":
        label = "Pressure [kPa]"
        axis = "pressure"
    elif axis == "t" or axis == "temperature":
        label = "Temperature [°C]"
        axis = "temperature"
    elif axis == "v" or axis == "specific volume":
        label = r"Specific Volume [m$^3$/kg]"
        axis = "specific_volume"
    elif axis == "dt*":
        label = r"$\Delta$T [K]"
        axis = "temperature"
    elif axis == "heat flow":
        label = "Heat Flow" + " [kW]"
        axis = "heat_flow"
    else:
        label = "Temperature [°C]"
        axis = "temperature"
    return label, axis
