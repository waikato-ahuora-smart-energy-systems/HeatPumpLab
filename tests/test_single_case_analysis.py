import pytest

from HeatPumpLab.cycles.RecuperativeHeatPump.utils import single_case_analysis


def test_single_case_analysis_requires_arguments():
    with pytest.raises(ValueError):
        single_case_analysis.analyse_single_case(object, args=None)


def test_single_case_analysis_reports_and_plots_converged_case(monkeypatch):
    terminal_outputs = []
    csv_outputs = []

    class FakeHP:
        def __init__(self, **args):
            self.args = args
            self.converged = False
            self.feasible = True
            self.output = {}
            self.gas_coolers = ["gc"]
            self.dT_GC_min_tar = 20.0
            self.plots = []

        def solve(self):
            self.converged = True

        def get_report(self):
            self.output = {"COP": 3.0}
            return self.output

        def get_driving_force_plot(self, *args, **kwargs):
            self.plots.append(("driving", args, kwargs))

        def get_cycle_plot(self, *args, **kwargs):
            self.plots.append(("cycle", args, kwargs))

    monkeypatch.setattr(
        single_case_analysis,
        "print_result_to_terminal",
        lambda output: terminal_outputs.append(output),
    )
    monkeypatch.setattr(
        single_case_analysis,
        "print_result_to_csv",
        lambda output, file_name: csv_outputs.append((output, file_name)),
    )

    hp = single_case_analysis.analyse_single_case(
        FakeHP, args={"refrigerant": "Butane"}
    )

    assert hp is None
    assert terminal_outputs == [{"COP": 3.0}]
    assert csv_outputs == [({"COP": 3.0}, "single-result")]


def test_single_case_analysis_reports_infeasible_case(capsys):
    class FakeHP:
        def __init__(self, **args):
            self.converged = False
            self.feasible = False

        def solve(self):
            return None

    single_case_analysis.analyse_single_case(FakeHP, args={"refrigerant": "Butane"})

    assert "infeasible" in capsys.readouterr().out


def test_single_case_analysis_reports_failed_solve(capsys):
    class FakeHP:
        def __init__(self, **args):
            self.converged = False
            self.feasible = True

        def solve(self):
            return None

    single_case_analysis.analyse_single_case(FakeHP, args={"refrigerant": "Butane"})

    assert "failed to solve" in capsys.readouterr().out
