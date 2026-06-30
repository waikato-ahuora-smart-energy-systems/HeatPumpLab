from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.figure  # noqa: E402

from HeatPumpLab.utils.graph import (  # noqa: E402
    _apply_layout,
    _get_axis_labels,
    _get_output_dir,
    _map_legend_name,
    _rescale,
    plot_cycle,
    plot_driving_force,
    plot_heat_exchanger_profiles,
    plot_heat_map,
    plot_simple_xy,
)


def state(temperature, enthalpy, pressure=1.0, entropy=1.0, heat_flow=1.0):
    return SimpleNamespace(
        temperature=temperature,
        enthalpy=enthalpy,
        pressure=pressure,
        entropy=entropy,
        specific_volume=1.0,
        heat_flow=heat_flow,
    )


class FakeHeatExchanger:
    name = "Gas cooler"

    def get_graph_data(self, divisions=50):
        return {
            "source": [
                state(100.0, 500.0, pressure=2.0, entropy=2.0, heat_flow=10.0),
                state(70.0, 400.0, pressure=1.8, entropy=1.8, heat_flow=8.0),
            ],
            "sink": [
                state(20.0, 100.0, pressure=1.0, entropy=1.0, heat_flow=2.0),
                state(60.0, 250.0, pressure=1.0, entropy=1.2, heat_flow=5.0),
            ],
        }


class NamedFakeHeatExchanger(FakeHeatExchanger):
    def __init__(self, name):
        self.name = name


def test_simple_xy_plot_returns_axes_and_writes_to_requested_output_dir(tmp_path):
    fig, ax = plot_simple_xy(
        [1, 2, 3],
        [4, 5, 6],
        title="Performance",
        x_axis="Load",
        y_axis="COP",
        output_dir=tmp_path,
        show=False,
    )

    assert ax.get_title() == "Performance"
    assert ax.get_xlabel() == "Load"
    assert ax.get_ylabel() == "COP"
    assert len(ax.lines) == 1
    assert (tmp_path / "Cycle--COP_vs_Load.png").exists()
    fig.clf()


def test_simple_xy_plot_allows_empty_title_and_axis_labels(tmp_path):
    fig, ax = plot_simple_xy(
        [1, 2],
        [3, 4],
        output_dir=tmp_path,
        show=False,
    )

    assert ax.get_title() == ""
    assert ax.get_xlabel() == ""
    assert ax.get_ylabel() == ""
    assert ax.get_legend_handles_labels()[1] == ["Series"]
    assert (tmp_path / "Cycle--_vs_.png").exists()
    fig.clf()


def test_heat_exchanger_profile_plot_uses_source_and_sink_labels(tmp_path):
    fig, ax = plot_heat_exchanger_profiles(
        [FakeHeatExchanger()],
        title="HX",
        output_dir=tmp_path,
        show=False,
    )

    labels = ax.get_legend_handles_labels()[1]
    assert "Source - GC" in labels
    assert "Sink - GC" in labels
    assert len(ax.lines) == 2
    assert (tmp_path / "Heat_exchanger_profiles.png").exists()
    fig.clf()


def test_heat_exchanger_profile_plot_accepts_no_exchangers(tmp_path):
    fig, ax = plot_heat_exchanger_profiles(
        hx_obj=None,
        title="Empty HX",
        output_dir=tmp_path,
        show=False,
    )

    assert ax.get_title() == "Empty HX"
    assert len(ax.lines) == 0
    assert (tmp_path / "Heat_exchanger_profiles.png").exists()
    fig.clf()


def test_driving_force_plot_includes_component_and_dtmin_lines(tmp_path):
    fig, ax = plot_driving_force(
        [FakeHeatExchanger()],
        title="Driving force",
        dtmin=5.0,
        output_dir=tmp_path,
        show=False,
    )

    labels = ax.get_legend_handles_labels()[1]
    assert "GC" in labels
    assert r"$ΔT_{min}$" in labels
    assert len(ax.lines) == 2
    assert (tmp_path / "Driving_force_plot.png").exists()
    fig.clf()


def test_driving_force_plot_labels_repeated_components_once(tmp_path):
    fig, ax = plot_driving_force(
        [FakeHeatExchanger(), FakeHeatExchanger()],
        dtmin=5.0,
        output_dir=tmp_path,
        show=False,
    )

    labels = ax.get_legend_handles_labels()[1]
    assert labels.count("GC") == 1
    assert labels.count(r"$ΔT_{min}$") == 1
    assert len(ax.lines) == 4
    fig.clf()


def test_cycle_plot_draws_cycle_and_sink_curves(tmp_path):
    hx = FakeHeatExchanger()

    fig, ax = plot_cycle(
        cycle_obj=[[hx, "source"]],
        sink_obj=[[hx, "sink"]],
        title="Cycle",
        output_dir=tmp_path,
        show=False,
    )

    labels = ax.get_legend_handles_labels()[1]
    assert "GC" in labels
    assert "AIR*" in labels
    assert len(ax.lines) == 3
    assert (tmp_path / "Cycle--temperature_vs_enthalpy.png").exists()
    fig.clf()


def test_cycle_plot_handles_pressure_axis_names_and_repeated_sink_labels(tmp_path):
    gas_cooler = NamedFakeHeatExchanger("Gas cooler")
    recuperator = NamedFakeHeatExchanger("Recuperator")
    ihx = NamedFakeHeatExchanger("Internal heat exchanger")
    evaporator = NamedFakeHeatExchanger("Evaporator")
    valve = NamedFakeHeatExchanger("Expansion valve")
    unnamed = NamedFakeHeatExchanger(None)
    generic_sink = NamedFakeHeatExchanger("Heat exchanger")

    fig, ax = plot_cycle(
        cycle_obj=[
            [gas_cooler, "source"],
            [recuperator, "source"],
            [ihx, "source"],
            [evaporator, "source"],
            [valve, "source"],
            [unnamed, "source"],
        ],
        sink_obj=[
            [gas_cooler, "sink"],
            [generic_sink, "sink"],
        ],
        x_axis="h",
        y_axis="P",
        output_dir=tmp_path,
        show=False,
    )

    labels = ax.get_legend_handles_labels()[1]
    assert ax.get_yscale() == "log"
    assert "GC" in labels
    assert "RHX" in labels
    assert "IHX" in labels
    assert "EVAP" in labels
    assert "VAL" in labels
    assert "Component" in labels
    assert labels.count("AIR*") == 1
    assert (tmp_path / "Cycle--pressure_vs_enthalpy.png").exists()
    fig.clf()


def test_cycle_plot_handles_sink_without_component_name_and_show_branch(
    monkeypatch, tmp_path
):
    shown = []
    monkeypatch.setattr(
        matplotlib.figure.Figure,
        "show",
        lambda self: shown.append(type(self).__name__),
    )
    gas_cooler = NamedFakeHeatExchanger("Gas cooler")
    unnamed_sink = NamedFakeHeatExchanger(None)

    fig, ax = plot_cycle(
        cycle_obj=[[gas_cooler, "source"]],
        sink_obj=[[unnamed_sink, "sink"]],
        x_axis="h",
        y_axis="T",
        output_dir=tmp_path,
        show=True,
    )

    assert shown == ["Figure"]
    assert ax.get_legend_handles_labels()[1].count("AIR*") == 1
    fig.clf()


def test_plot_show_branches_call_figure_show(monkeypatch, tmp_path):
    shown = []
    monkeypatch.setattr(
        matplotlib.figure.Figure,
        "show",
        lambda self: shown.append(type(self).__name__),
    )

    fig1, _ = plot_driving_force(
        [FakeHeatExchanger()],
        output_dir=tmp_path,
        show=True,
    )
    fig2, _ = plot_heat_exchanger_profiles(
        [FakeHeatExchanger()],
        output_dir=tmp_path,
        show=True,
    )
    fig3, _ = plot_simple_xy(
        [1, 2],
        [3, 4],
        output_dir=tmp_path,
        show=True,
    )

    assert shown == ["Figure", "Figure", "Figure"]
    fig1.clf()
    fig2.clf()
    fig3.clf()


def test_heat_map_returns_plotly_figure_with_domain_labels():
    fig = plot_heat_map(
        title="SEC map",
        x_val=[1, 2, 3],
        x_label="Lift",
        y_val=[4, 5, 6],
        y_label="Pressure",
        z_val=[0.4, 0.5, 0.6],
        z_label="SEC",
        show=False,
    )

    assert fig.layout.title.text == "SEC map"
    assert fig.layout.xaxis.title.text == "Lift"
    assert fig.layout.yaxis.title.text == "Pressure"
    assert fig.layout.coloraxis.colorbar.title.text == "SEC"


def test_heat_map_show_branch_calls_plotly_show(monkeypatch):
    shown = []
    monkeypatch.setattr(
        "plotly.graph_objects.Figure.show",
        lambda self: shown.append(self.layout.title.text),
    )

    fig = plot_heat_map(
        title="Visible map",
        x_val=[1],
        x_label="x",
        y_val=[2],
        y_label="y",
        z_val=[3],
        z_label="z",
        show=True,
    )

    assert fig.layout.title.text == "Visible map"
    assert shown == ["Visible map"]


def test_graph_axis_helpers_cover_supported_aliases_and_fallback():
    assert _get_axis_labels("s") == ("Entropy [kJ/kgK]", "entropy")
    assert _get_axis_labels("v") == (
        r"Specific Volume [m$^3$/kg]",
        "specific_volume",
    )
    assert _get_axis_labels("dt*") == (r"$\Delta$T [K]", "temperature")
    assert _get_axis_labels("heat flow") == ("Heat Flow [kW]", "heat_flow")
    assert _get_axis_labels("unknown") == ("Temperature [°C]", "temperature")


def test_legend_and_output_directory_helpers_cover_defaults(tmp_path):
    assert _map_legend_name(None) == "Component"
    assert _get_output_dir("Butane", tmp_path) == tmp_path
    assert _get_output_dir("Butane", None).parts[-4:] == (
        "results",
        "figures",
        "cycles",
        "Butane",
    )


def test_rescale_handles_flat_source_axis():
    assert _rescale([2.0, 2.0], [10.0, 20.0]) == [10.0, 10.0]


def test_layout_errors_are_ignored_for_plotting_fallbacks():
    class BadFigure:
        def tight_layout(self, pad):
            raise RuntimeError("layout unavailable")

    assert _apply_layout(BadFigure()) is None
