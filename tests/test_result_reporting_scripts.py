from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from HeatPumpLab.cycles.RecuperativeHeatPump.utils import (  # noqa: E402
    plot_comp2_power_frac_analysis as comp2_plot,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.utils import (  # noqa: E402
    plot_heatmap_metrics as heatmap,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.utils import (  # noqa: E402
    summarize_best_results as summary,
)


def sample_results() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "refrigerant": "Butane",
                "Success": True,
                "converged": True,
                "COP": 2.0,
                "SEC": 0.5,
                "VHC": 10.0,
                "eDTL": 50.0,
                "SLE": 0.01,
                "LMTD_gcs": 18.0,
                "Pmax": 1000.0,
                "dT_rhx1_cold": 0.0,
                "dT_rhx2_cold": 0.0,
                "dT_rhx3_cold": 0.0,
                "comp2_power_frac": 0.5,
                "dT_gc_out_min_offset": 0.0,
                "cycle_type": "Subcritical",
            },
            {
                "refrigerant": "Butane",
                "Success": True,
                "converged": True,
                "COP": 2.2,
                "SEC": 0.45,
                "VHC": 11.0,
                "eDTL": 48.0,
                "SLE": 0.009,
                "LMTD_gcs": 16.0,
                "Pmax": 900.0,
                "dT_rhx1_cold": 5.0,
                "dT_rhx2_cold": 5.0,
                "dT_rhx3_cold": 0.0,
                "comp2_power_frac": 0.5,
                "dT_gc_out_min_offset": 2.0,
                "cycle_type": "Subcritical",
            },
            {
                "refrigerant": "Propane",
                "Success": True,
                "converged": True,
                "COP": 3.0,
                "SEC": 0.33,
                "VHC": 12.0,
                "eDTL": 45.0,
                "SLE": 0.008,
                "LMTD_gcs": 15.0,
                "Pmax": 1200.0,
                "dT_rhx1_cold": 0.0,
                "dT_rhx2_cold": 0.0,
                "dT_rhx3_cold": 0.0,
                "comp2_power_frac": 0.5,
                "dT_gc_out_min_offset": 1.0,
                "cycle_type": "Subcritical",
            },
        ]
    )


def test_summary_load_results_filters_success_and_converged_rows(tmp_path):
    df = sample_results()
    df.loc[1, "Success"] = False
    df.loc[2, "converged"] = False
    path = tmp_path / "results.csv"
    df.to_csv(path, index=False)

    loaded = summary.load_results(path, filter_success=True, filter_converged=True)

    assert loaded["refrigerant"].tolist() == ["Butane"]


def test_summary_parse_args_exposes_cli_defaults(monkeypatch):
    monkeypatch.setattr("sys.argv", ["summary"])

    args = summary.parse_args()

    assert "COP:max" in args.metrics
    assert "cycle_type" in args.include_columns
    assert args.float_precision == 4
    assert args.refrigerants is None
    assert args.comp2_values is None


def test_summary_metric_specs_parse_modes_and_reject_invalid_values():
    assert summary.parse_metric_specs(["COP:min", "VHC"]) == [
        ("COP", "min"),
        ("VHC", "max"),
    ]
    with pytest.raises(ValueError):
        summary.parse_metric_specs([":max"])
    with pytest.raises(ValueError):
        summary.parse_metric_specs(["COP:median"])


def test_summary_unique_sorted_ignores_non_numeric_values():
    assert summary.unique_sorted([0.5, "bad", float("nan"), 0.4]) == [0.4, 0.5]


def test_summary_basecase_metrics_are_relative_to_each_refrigerant():
    df = summary.add_basecase_metrics(sample_results())
    improved = df[(df["refrigerant"] == "Butane") & (df["COP"] == 2.2)].iloc[0]

    assert improved["COP_pct_increase_vs_base"] == pytest.approx(10.0)
    assert improved["Pmax_reduction_vs_base"] == pytest.approx(10.0)


def test_summary_basecase_metrics_raise_when_required_columns_are_missing():
    with pytest.raises(ValueError):
        summary.add_basecase_metrics(pd.DataFrame({"refrigerant": ["Butane"]}))


def test_summary_basecase_metrics_warn_when_no_basecase_exists(capsys):
    df = sample_results()
    df["comp2_power_frac"] = 0.4

    enriched = summary.add_basecase_metrics(df)

    assert enriched["COP_pct_increase_vs_base"].isna().all()
    assert "No base-case rows" in capsys.readouterr().out


def test_summary_filters_refrigerants_and_comp2_values():
    df = sample_results()

    filtered = summary.filter_dataframe(df, refrigerants=["Butane"], comp2_values=[0.5])

    assert filtered["refrigerant"].unique().tolist() == ["Butane"]
    with pytest.raises(ValueError):
        summary.filter_dataframe(df, refrigerants=["Missing"], comp2_values=None)
    with pytest.raises(ValueError):
        summary.filter_dataframe(df, refrigerants=None, comp2_values=[0.1])


def test_summary_filter_reports_missing_comp2_column():
    with pytest.raises(ValueError, match="comp2_power_frac"):
        summary.filter_dataframe(
            pd.DataFrame({"refrigerant": ["Butane"]}),
            refrigerants=None,
            comp2_values=[0.5],
        )


def test_summary_best_cases_selects_best_metric_per_refrigerant():
    df = summary.add_basecase_metrics(sample_results())

    best = summary.summarise_best_cases(
        df,
        metrics=[("COP", "max"), ("Pmax", "min")],
        include_columns=["Pmax", "COP_pct_increase_vs_base"],
    )

    butane_cop = best[(best["metric"] == "COP") & (best["refrigerant"] == "Butane")]
    butane_pmax = best[(best["metric"] == "Pmax") & (best["refrigerant"] == "Butane")]
    assert butane_cop.iloc[0]["metric_value"] == pytest.approx(2.2)
    assert butane_pmax.iloc[0]["metric_value"] == pytest.approx(900.0)


def test_summary_best_cases_skips_missing_or_invalid_metrics(capsys):
    best = summary.summarise_best_cases(
        sample_results(),
        metrics=[("Missing", "max")],
        include_columns=["COP"],
    )

    assert best.empty
    assert "Skipping metric 'Missing'" in capsys.readouterr().out


def test_summary_best_cases_handles_non_numeric_metrics_and_default_columns(capsys):
    df = sample_results()
    df["bad_metric"] = "not numeric"

    best = summary.summarise_best_cases(
        df,
        metrics=[("bad_metric", "max"), ("COP", "max")],
        include_columns=None,
    )

    assert best["metric"].tolist() == ["COP", "COP"]
    assert "Skipping metric 'bad_metric'" in capsys.readouterr().out


def test_best_cases_summary_writes_filtered_summary(monkeypatch, tmp_path):
    input_csv = tmp_path / "results.csv"
    output_csv = tmp_path / "summary.csv"
    sample_results().to_csv(input_csv, index=False)
    monkeypatch.setattr(
        summary,
        "parse_args",
        lambda: SimpleNamespace(
            metrics=["COP:max"],
            include_columns=["Pmax"],
            refrigerants=["Butane"],
            comp2_values=None,
            skip_success_filter=False,
            skip_converged_filter=False,
            float_precision=2,
        ),
    )

    summary.get_best_cases_summary(filepath=input_csv, output=output_csv)

    written = pd.read_csv(output_csv)
    assert written["metric"].tolist() == ["COP"]
    assert written["refrigerant"].tolist() == ["Butane"]
    assert written["metric_value"].tolist() == [2.2]


def test_best_cases_summary_reports_empty_summary(monkeypatch, tmp_path, capsys):
    input_csv = tmp_path / "results.csv"
    sample_results().to_csv(input_csv, index=False)
    monkeypatch.setattr(
        summary,
        "parse_args",
        lambda: SimpleNamespace(
            metrics=["Missing:max"],
            include_columns=["Pmax"],
            refrigerants=["Butane"],
            comp2_values=None,
            skip_success_filter=False,
            skip_converged_filter=False,
            float_precision=2,
        ),
    )

    summary.get_best_cases_summary(filepath=input_csv, output=None)

    assert "No summary results generated" in capsys.readouterr().out


def test_comp2_parse_args_exposes_cli_defaults(monkeypatch):
    monkeypatch.setattr("sys.argv", ["comp2"])

    args = comp2_plot.parse_args()

    assert "COP_pct_increase_vs_base:max" in args.metrics
    assert "Butane" in args.refrigerants
    assert args.comp2_values is None
    assert args.dpi == 300


def test_comp2_load_results_filters_success_and_converged_rows(tmp_path):
    df = sample_results()
    df.loc[1, "Success"] = False
    df.loc[2, "converged"] = False
    path = tmp_path / "results.csv"
    df.to_csv(path, index=False)

    loaded = comp2_plot.load_results(
        path,
        filter_success=True,
        filter_converged=True,
    )

    assert loaded["refrigerant"].tolist() == ["Butane"]


def test_comp2_helpers_parse_sanitise_and_filter_data():
    assert comp2_plot.unique_sorted([0.5, 0.500000001, np.nan, 0.4]) == [0.4, 0.5]
    assert comp2_plot.sanitize_metric_name("dQ/dW value") == "dQ_per_dW_value"
    assert comp2_plot.parse_metric_specs(["Pmax:min", "COP"]) == [
        ("Pmax", "min"),
        ("COP", "max"),
    ]
    with pytest.raises(ValueError):
        comp2_plot.parse_metric_specs(["Pmax:middle"])
    with pytest.raises(ValueError):
        comp2_plot.parse_metric_specs([":max"])

    filtered = comp2_plot.filter_dataframe(
        sample_results(),
        refrigerants=["Butane"],
        comp2_values=[0.5],
    )
    assert filtered["refrigerant"].unique().tolist() == ["Butane"]


def test_comp2_basecase_metrics_reports_missing_columns_and_no_basecase(capsys):
    with pytest.raises(ValueError):
        comp2_plot.add_basecase_metrics(pd.DataFrame({"refrigerant": ["Butane"]}))

    df = sample_results()
    df["comp2_power_frac"] = 0.4
    enriched = comp2_plot.add_basecase_metrics(df)

    assert enriched["COP_pct_increase_vs_base"].isna().all()
    assert "No base-case rows" in capsys.readouterr().out


def test_comp2_basecase_metrics_and_plot_file(tmp_path):
    df = comp2_plot.add_basecase_metrics(sample_results())

    comp2_plot.plot_metric_for_refrigerants(
        df=df,
        metric="COP_pct_increase_vs_base",
        mode="max",
        include_basecase_line=True,
        fig_width=3.0,
        fig_height=2.5,
        dpi=80,
        output_dir=tmp_path,
    )

    assert df.loc[1, "COP_pct_increase_vs_base"] == pytest.approx(10.0)
    assert (tmp_path / "COP_pct_increase_vs_base_best_vs_comp2_power_frac.png").exists()


def test_comp2_plot_supports_minimisation_mode(tmp_path):
    comp2_plot.plot_metric_for_refrigerants(
        df=sample_results(),
        metric="Pmax",
        mode="min",
        include_basecase_line=False,
        fig_width=3.0,
        fig_height=2.5,
        dpi=80,
        output_dir=tmp_path,
    )

    assert (tmp_path / "Pmax_best_vs_comp2_power_frac.png").exists()


def test_comp2_plot_skips_missing_metric(tmp_path, capsys):
    comp2_plot.plot_metric_for_refrigerants(
        df=sample_results(),
        metric="Missing",
        mode="max",
        include_basecase_line=False,
        fig_width=3.0,
        fig_height=2.5,
        dpi=80,
        output_dir=tmp_path,
    )

    assert "column not found" in capsys.readouterr().out


def test_comp2_plot_reports_when_no_refrigerant_or_metric_data_remains(tmp_path, capsys):
    empty_refrigerant_df = sample_results().iloc[0:0]
    comp2_plot.plot_metric_for_refrigerants(
        df=empty_refrigerant_df,
        metric="COP",
        mode="max",
        include_basecase_line=False,
        fig_width=3.0,
        fig_height=2.5,
        dpi=80,
        output_dir=tmp_path,
    )
    assert "No refrigerant data available" in capsys.readouterr().out

    no_metric_df = sample_results()
    no_metric_df["COP"] = np.nan
    comp2_plot.plot_metric_for_refrigerants(
        df=no_metric_df,
        metric="COP",
        mode="min",
        include_basecase_line=False,
        fig_width=3.0,
        fig_height=2.5,
        dpi=80,
        output_dir=tmp_path,
    )
    assert "No data available to plot" in capsys.readouterr().out


def test_comp2_comparison_analysis_loads_filters_and_writes_plots(monkeypatch, tmp_path):
    input_csv = tmp_path / "results.csv"
    output_dir = tmp_path / "figures"
    sample_results().to_csv(input_csv, index=False)
    monkeypatch.setattr(
        comp2_plot,
        "parse_args",
        lambda: SimpleNamespace(
            csv=input_csv,
            output_dir=output_dir,
            metrics=["COP_pct_increase_vs_base:max"],
            refrigerants=["Butane"],
            comp2_values=[0.5],
            skip_success_filter=False,
            skip_converged_filter=False,
            include_basecase_line=True,
            dpi=80,
            fig_width=3.0,
            fig_height=2.5,
        ),
    )

    comp2_plot.get_comp2_power_frac_comparison_analysis()

    assert (
        output_dir / "COP_pct_increase_vs_base_best_vs_comp2_power_frac.png"
    ).exists()


def test_comp2_comparison_analysis_reports_missing_columns_and_empty_filters(
    monkeypatch, tmp_path
):
    missing_csv = tmp_path / "missing.csv"
    sample_results().drop(columns=["comp2_power_frac"]).to_csv(missing_csv, index=False)
    monkeypatch.setattr(
        comp2_plot,
        "parse_args",
        lambda: SimpleNamespace(
            csv=missing_csv,
            output_dir=tmp_path,
            metrics=["COP:max"],
            refrigerants=["Butane"],
            comp2_values=None,
            skip_success_filter=False,
            skip_converged_filter=False,
            include_basecase_line=False,
            dpi=80,
            fig_width=3.0,
            fig_height=2.5,
        ),
    )
    with pytest.raises(ValueError, match="missing columns"):
        comp2_plot.get_comp2_power_frac_comparison_analysis()

    filtered_csv = tmp_path / "filtered.csv"
    sample_results().to_csv(filtered_csv, index=False)
    monkeypatch.setattr(
        comp2_plot,
        "parse_args",
        lambda: SimpleNamespace(
            csv=filtered_csv,
            output_dir=tmp_path,
            metrics=["COP:max"],
            refrigerants=["Missing"],
            comp2_values=None,
            skip_success_filter=False,
            skip_converged_filter=False,
            include_basecase_line=False,
            dpi=80,
            fig_width=3.0,
            fig_height=2.5,
        ),
    )
    with pytest.raises(ValueError, match="No data remains"):
        comp2_plot.get_comp2_power_frac_comparison_analysis()

    already_enriched_csv = tmp_path / "already_enriched.csv"
    enriched = sample_results().drop(columns=["comp2_power_frac"])
    enriched["COP_pct_increase_vs_base"] = 0.0
    enriched["Pmax_reduction_vs_base"] = 0.0
    enriched.to_csv(already_enriched_csv, index=False)
    monkeypatch.setattr(
        comp2_plot,
        "parse_args",
        lambda: SimpleNamespace(
            csv=already_enriched_csv,
            output_dir=tmp_path,
            metrics=["COP:max"],
            refrigerants=["Butane"],
            comp2_values=None,
            skip_success_filter=False,
            skip_converged_filter=False,
            include_basecase_line=False,
            dpi=80,
            fig_width=3.0,
            fig_height=2.5,
        ),
    )
    with pytest.raises(ValueError, match="required columns"):
        comp2_plot.get_comp2_power_frac_comparison_analysis()


def test_heatmap_load_results_removes_header_rows_and_filters_success(tmp_path):
    path = tmp_path / "results.csv"
    df = sample_results()
    bad_row = {column: column for column in df.columns}
    failed_row = df.iloc[0].copy()
    failed_row["Success"] = False
    pd.concat([df.iloc[[0]], pd.DataFrame([bad_row, failed_row])]).to_csv(
        path, index=False
    )

    loaded = heatmap._load_results(path)

    assert len(loaded) == 1
    assert loaded.iloc[0]["refrigerant"] == "Butane"


def test_heatmap_parse_args_exposes_cli_defaults(monkeypatch):
    monkeypatch.setattr("sys.argv", ["heatmap"])

    args = heatmap._parse_args()

    assert "COP:max" in args.metrics
    assert args.comp2_values == [0.5]
    assert "Butane" in args.refrigerants
    assert args.dpi == 900


def test_heatmap_helpers_parse_sanitise_and_build_grid():
    assert heatmap._unique_sorted([0.5, "0.4", None, float("nan")]) == [0.4, 0.5]
    assert heatmap._sanitize_float(-1.25) == "neg1p25"
    assert heatmap._sanitize_metric_name("dQ/dW value") == "dQ_per_dW_value"
    assert heatmap._format_axis_label("COP") == "COP"
    assert heatmap._parse_metric_specs(["COP:min"]) == [("COP", "min")]
    with pytest.raises(ValueError):
        heatmap._parse_metric_specs(["COP:median"])

    x_vals, y_vals, grid = heatmap._build_grid(
        sample_results(),
        metric="COP",
        x_col="dT_rhx1_cold",
        y_col="dT_rhx2_cold",
    )

    assert x_vals.tolist() == [0.0, 5.0]
    assert y_vals.tolist() == [0.0, 5.0]
    assert grid.shape == (2, 2)


def test_heatmap_unique_sorted_ignores_unusable_values():
    assert heatmap._unique_sorted(["bad", None, float("nan"), "0.5"]) == [0.5]


def test_heatmap_helpers_reject_empty_metric_name_and_empty_grid():
    with pytest.raises(ValueError, match="Metric name"):
        heatmap._parse_metric_specs([":max"])

    x_vals, y_vals, grid = heatmap._build_grid(
        pd.DataFrame({"dT_rhx1_cold": [np.nan], "dT_rhx2_cold": [np.nan], "COP": [1.0]}),
        metric="COP",
        x_col="dT_rhx1_cold",
        y_col="dT_rhx2_cold",
    )
    assert x_vals.size == 0
    assert y_vals.size == 0
    assert grid.size == 0

    _, _, all_nan_grid = heatmap._build_grid(
        pd.DataFrame(
            {
                "dT_rhx1_cold": [0.0],
                "dT_rhx2_cold": [0.0],
                "COP": [np.nan],
            }
        ),
        metric="COP",
        x_col="dT_rhx1_cold",
        y_col="dT_rhx2_cold",
    )
    assert all_nan_grid.size == 0


def test_heatmap_basecase_metrics_and_missing_columns():
    df = heatmap._add_basecase_metrics(sample_results())

    assert df.loc[1, "Pmax_reduction_vs_base"] == pytest.approx(10.0)
    with pytest.raises(ValueError):
        heatmap._add_basecase_metrics(pd.DataFrame({"refrigerant": ["Butane"]}))


def test_heatmap_basecase_metrics_reports_missing_basecase(capsys):
    df = sample_results()
    df["comp2_power_frac"] = 0.4

    enriched = heatmap._add_basecase_metrics(df)

    assert enriched["COP_pct_increase_vs_base"].isna().all()
    assert "No base-case rows" in capsys.readouterr().out


def test_heatmap_creation_writes_expected_nested_file(tmp_path):
    df = heatmap._add_basecase_metrics(sample_results())

    heatmap._create_heatmaps(
        df=df,
        metrics=[("COP", "max"), ("Pmax", "min")],
        output_dir=tmp_path,
        comp2_values=[0.5],
        refrigerants_filter=["Butane"],
        cmap="viridis",
        dpi=80,
        use_log_pmax=True,
    )

    assert (tmp_path / "COP" / "Butane" / "comp2_0p5" / "heatmap.png").exists()
    assert (tmp_path / "Pmax" / "Butane" / "comp2_0p5" / "heatmap.png").exists()


def test_heatmap_creation_reports_invalid_filters(tmp_path):
    df = heatmap._add_basecase_metrics(sample_results())

    with pytest.raises(ValueError):
        heatmap._create_heatmaps(
            df=df,
            metrics=[("COP", "max")],
            output_dir=tmp_path,
            comp2_values=[0.5],
            refrigerants_filter=["Missing"],
            cmap="viridis",
            dpi=80,
            use_log_pmax=False,
        )


def test_heatmap_creation_reports_missing_comp2_filter(tmp_path):
    df = heatmap._add_basecase_metrics(sample_results())

    with pytest.raises(ValueError, match="comp2_power_frac"):
        heatmap._create_heatmaps(
            df=df,
            metrics=[("COP", "max")],
            output_dir=tmp_path,
            comp2_values=[0.1],
            refrigerants_filter=["Butane"],
            cmap="viridis",
            dpi=80,
            use_log_pmax=False,
        )


def test_heatmap_creation_skips_sparse_refrigerant_slices(tmp_path):
    df = pd.DataFrame(
        [
            {
                "refrigerant": "Butane",
                "Success": True,
                "COP": 2.0,
                "Pmax": 1000.0,
                "dT_rhx1_cold": 0.0,
                "dT_rhx2_cold": 0.0,
                "comp2_power_frac": 0.5,
            },
            {
                "refrigerant": "Propane",
                "Success": True,
                "COP": 3.0,
                "Pmax": 1200.0,
                "dT_rhx1_cold": 0.0,
                "dT_rhx2_cold": 0.0,
                "comp2_power_frac": 0.6,
            },
        ]
    )

    heatmap._create_heatmaps(
        df=df,
        metrics=[("COP", "max")],
        output_dir=tmp_path,
        comp2_values=[0.5, 0.6],
        refrigerants_filter=["Butane", "Propane"],
        cmap="viridis",
        dpi=80,
        use_log_pmax=False,
    )

    assert (tmp_path / "COP" / "Butane" / "comp2_0p5" / "heatmap.png").exists()
    assert (tmp_path / "COP" / "Propane" / "comp2_0p6" / "heatmap.png").exists()


def test_heatmap_creation_skips_empty_grids_and_falls_back_from_unusable_log_scale(
    tmp_path,
):
    df = pd.DataFrame(
        [
            {
                "refrigerant": "Butane",
                "Success": True,
                "Pmax": np.nan,
                "dT_rhx1_cold": 0.0,
                "dT_rhx2_cold": 0.0,
                "comp2_power_frac": 0.5,
            },
            {
                "refrigerant": "Propane",
                "Success": True,
                "Pmax": 10.0,
                "dT_rhx1_cold": 0.0,
                "dT_rhx2_cold": 0.0,
                "comp2_power_frac": 0.5,
            },
            {
                "refrigerant": "DME",
                "Success": True,
                "Pmax": 0.0,
                "dT_rhx1_cold": 0.0,
                "dT_rhx2_cold": 0.0,
                "comp2_power_frac": 0.5,
            },
        ]
    )

    heatmap._create_heatmaps(
        df=df,
        metrics=[("Pmax", "min")],
        output_dir=tmp_path,
        comp2_values=[0.5],
        refrigerants_filter=["Butane", "Propane", "DME"],
        cmap="viridis",
        dpi=80,
        use_log_pmax=True,
    )

    assert not (tmp_path / "Pmax" / "Butane" / "comp2_0p5" / "heatmap.png").exists()
    assert (tmp_path / "Pmax" / "Propane" / "comp2_0p5" / "heatmap.png").exists()
    assert (tmp_path / "Pmax" / "DME" / "comp2_0p5" / "heatmap.png").exists()


def test_heatmap_creation_reports_missing_data_and_skips_invalid_metrics(
    tmp_path, capsys
):
    with pytest.raises(ValueError, match="No refrigerant data"):
        heatmap._create_heatmaps(
            df=pd.DataFrame({"comp2_power_frac": [0.5]}),
            metrics=[("COP", "max")],
            output_dir=tmp_path,
            comp2_values=[0.5],
            refrigerants_filter=["Butane"],
            cmap="viridis",
            dpi=80,
            use_log_pmax=False,
        )

    with pytest.raises(ValueError, match="comp2_power_frac"):
        heatmap._create_heatmaps(
            df=pd.DataFrame({"refrigerant": ["Butane"]}),
            metrics=[("COP", "max")],
            output_dir=tmp_path,
            comp2_values=[0.5],
            refrigerants_filter=["Butane"],
            cmap="viridis",
            dpi=80,
            use_log_pmax=False,
        )

    df = sample_results()
    df["Invalid"] = np.nan
    heatmap._create_heatmaps(
        df=df,
        metrics=[("Missing", "max"), ("Invalid", "max")],
        output_dir=tmp_path,
        comp2_values=[0.5],
        refrigerants_filter=["Butane"],
        cmap="viridis",
        dpi=80,
        use_log_pmax=False,
    )

    out = capsys.readouterr().out
    assert "Skipping metric 'Missing'" in out
    assert "Skipping metric 'Invalid'" in out


def test_heatmap_results_comparison_uses_cli_defaults_and_writes_files(
    monkeypatch, tmp_path
):
    input_csv = tmp_path / "results.csv"
    output_dir = tmp_path / "heatmaps"
    sample_results().to_csv(input_csv, index=False)
    monkeypatch.setattr("sys.argv", ["heatmap"])

    heatmap.get_heatmap_results_comparison(
        filepath=input_csv,
        output_dir=output_dir,
        metrics_ls=[("COP", "max")],
        refrigerant_ls=["Butane"],
        comp2_load_frac_ls=[0.5],
    )

    assert (output_dir / "COP" / "Butane" / "comp2_0p5" / "heatmap.png").exists()


def test_heatmap_and_summary_agree_on_metric_directions():
    heat = dict(heatmap.DEFAULT_METRICS_AND_OBJ_PAIRS)
    best = dict(summary.DEFAULT_METRICS)
    assert {m: heat[m] for m in heat.keys() & best.keys()} == {m: best[m] for m in heat.keys() & best.keys()}
    assert best["LMTD_gcs"] == "min"


def test_vhc_label_matches_the_report_unit():
    # the report's VHC is kW per m3/s of suction, i.e. kJ/m3
    assert "kJ/m" in heatmap.AXIS_LABELS["VHC"]
