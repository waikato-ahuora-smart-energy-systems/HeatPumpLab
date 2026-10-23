import json

import matplotlib

matplotlib.use("Agg")

import pytest  # noqa: E402

from HeatPumpLab.cycles.RecuperativeHeatPump.utils import (  # noqa: E402
    optimiser_benchmark_report as optimiser_report,
)
from HeatPumpLab.cycles.RecuperativeHeatPump.utils import (  # noqa: E402
    origin_benchmark_report as origin_report,
)

REF = "Propane;Pentane|0.3;0.7|mass"
BASE_CASE = {"T_air_tar": 200.0, "T_air_inlet": 15.0, "dT_rhx3_cold": 0.0}
COMPOSITION = [0.3, 0.7]


def _table_rows(markdown, header_start):
    """Return the data rows of the first markdown table whose header starts with header_start."""
    lines = markdown.splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(header_start))
    rows = []
    for line in lines[start + 2:]:
        if not line.startswith("|"):
            break
        rows.append([cell.strip() for cell in line.strip("|").split("|")])
    return rows


# --- case_label -----------------------------------------------------------

@pytest.mark.parametrize("refrigerant, label", [
    ("Propane;Pentane|0.3;0.7|mass", "R290/R601 30/70"),
    ("Propane;Isobutane|0.5;0.5|mass", "R290/R600A 50/50"),
    ("Propane;Isopentane|0.25;0.75|mass", "R290/R601A 25/75"),
    ("Propane;Butane|0.9;0.1|mass", "R290/R600 90/10"),
    ("Propane;Hexane|0.4;0.6|mass", "R290/Hexane 40/60"),
])
def test_case_label_assumes_propane_is_first_component(refrigerant, label):
    # Current behaviour: the first fluid is always labelled R290 and only the
    # second fluid's name is mapped; percentages keep the given order.
    assert optimiser_report.case_label(refrigerant) == label


# --- origin_benchmark_report ----------------------------------------------

def _origin_trial(method, seed, cops, *, failed=(), budget=None):
    """cops[0] is the origin COP; failed indexes are unsuccessful evaluations."""
    budget = budget or len(cops)
    records = []
    for i, cop in enumerate(cops):
        success = i not in failed
        records.append({
            "x": [float(i), 0.0, 0.5] + COMPOSITION,
            "success": success,
            "COP": cop if success else None,
            "SEC": 1 / cop if success else None,
            "evaluation_wall_seconds": 10.0 * seed,
            "proposal_seconds": 1.0,
        })
    best = max((r for r in records if r["success"]), key=lambda r: r["COP"])
    return {
        "metadata": {
            "refrigerant": REF, "method": method, "seed": seed, "budget": budget,
            "origin": [0.0, 0.0, 0.5], "base_case": BASE_CASE, "baseline_COP": 2.4,
            "search_source_sha256": "abc", "versions": {"numpy": "1"},
        },
        "records": records,
        "replay": {"success": True, "x": list(best["x"]), "COP": best["COP"],
                   "VHC_MJ_per_m3": 3.25, "Pmax_MPa": 4.5},
    }


def _origin_trials():
    return [
        _origin_trial("rbf", 1, [2.0, 2.9985, 2.1]),
        _origin_trial("gp_ei", 2, [2.0, 2.1, 3.0]),
        _origin_trial("gp_ei", 1, [2.0, 2.5, 2.4]),
        _origin_trial("gp_ei", 3, [2.0, 2.2, 9.9], failed={2}),
    ]


def _write_trials(tmp_path, trials):
    for i, trial in enumerate(trials):
        (tmp_path / f"trial_{i}.json").write_text(json.dumps(trial))
    return [(t["metadata"]["refrigerant"], t["metadata"]["method"], t["metadata"]["seed"],
             t["metadata"]["budget"], "extra") for t in trials]


def test_origin_report_summary_uses_replayed_cops_and_best_observed_reference(tmp_path):
    expected = _write_trials(tmp_path, _origin_trials())
    path = origin_report.write_origin_report(tmp_path, expected)

    summary = {row["method"]: row for row in json.loads((tmp_path / "summary.json").read_text())}
    gp = summary["gp_ei"]
    assert gp["runs"] == 3
    assert gp["median_COP"] == pytest.approx(2.5)
    assert (gp["min_COP"], gp["max_COP"]) == (pytest.approx(2.2), pytest.approx(3.0))
    assert gp["median_COP_at_64"] == pytest.approx(2.5)  # last running best for short budgets
    assert gp["best_observed_COP"] == pytest.approx(3.0)
    assert gp["within_0_1_percent_best_observed"] == 1
    assert gp["beats_v11"] == 2  # 2.5 and 3.0 exceed baseline 2.4
    assert gp["failed_evaluations"] == 1
    assert gp["origin_COP"] == 2.0
    assert gp["median_wall_seconds"] == pytest.approx(3 * (20.0 + 1.0))  # seeds 1,2,3 -> median seed 2
    # 2.9985 is within 0.1% of the best observed 3.0 found by a different method.
    assert summary["rbf"]["within_0_1_percent_best_observed"] == 1
    assert summary["rbf"]["best_observed_COP"] == pytest.approx(3.0)

    markdown = path.read_text()
    consistency = _table_rows(markdown, "| Composition by mass | Method | Seeds")
    assert consistency[0] == ["R290/R601 30/70", "BO (GP-EI)", "3", "2.500000", "2.500000",
                              "2.200000–3.000000", "1/3", "63.0", "1"]
    assert consistency[1][:3] == ["R290/R601 30/70", "RBF", "1"]

    trials_table = _table_rows(markdown, "| Composition | Method | Seed")
    assert [(row[1], row[2]) for row in trials_table] == [
        ("BO (GP-EI)", "1"), ("BO (GP-EI)", "2"), ("BO (GP-EI)", "3"), ("RBF", "1")]
    assert trials_table[1][4:7] == ["3.000000", "3.2500", "4.5000"]
    assert trials_table[1][7:10] == ["2.00000", "0.00000", "0.500000"]

    validation = json.loads((tmp_path / "validation.json").read_text())
    assert validation["trials"] == 4
    assert validation["candidate_evaluations"] == 12
    assert validation["failed_evaluations"] == 1
    assert (tmp_path / "origin_convergence.png").stat().st_size > 0
    assert "Inner-solver sensitivity checks" not in markdown  # no diagnostic files present


def test_origin_report_ignores_unexpected_trials_and_reports_missing(tmp_path):
    trials = _origin_trials()
    expected = _write_trials(tmp_path, trials)
    with pytest.raises(ValueError, match="Missing 1 expected trials"):
        origin_report.write_origin_report(tmp_path, expected + [(REF, "cma_es", 1, 3)])


def test_origin_report_rejects_differing_origin_results(tmp_path):
    trials = _origin_trials()
    trials[0]["records"][0]["COP"] = 2.01
    trials[0]["records"][0]["SEC"] = 1 / 2.01
    with pytest.raises(ValueError, match="Origin results differ"):
        origin_report.write_origin_report(tmp_path, _write_trials(tmp_path, trials))


def test_validate_trial_requires_replay_of_lowest_sec_design():
    trial = _origin_trial("gp_ei", 1, [2.0, 2.5, 2.4])
    origin_report.validate_trial(trial)  # well-formed trial passes
    trial["replay"]["x"] = trial["records"][2]["x"]
    with pytest.raises(ValueError, match="Replay used a different design"):
        origin_report.validate_trial(trial)


@pytest.mark.parametrize("mutate, message", [
    (lambda t: t["records"].pop(), "has not completed its budget"),
    (lambda t: t["records"][0]["x"].__setitem__(0, 1.0), "not the declared origin"),
    (lambda t: t["records"][1]["x"].__setitem__(3, 0.31), "Composition varied"),
    (lambda t: t["records"][1]["x"].__setitem__(1, 150.0), "outside individual bounds"),
    (lambda t: t["records"][1]["x"].__setitem__(1, 99.0) or t["records"][1]["x"].__setitem__(0, 30.0),
     "combined lift limit"),
    (lambda t: t["replay"].__setitem__("COP", 2.6), "failed the fresh-process replay"),
    (lambda t: t["replay"].__setitem__("VHC_MJ_per_m3", 0.0), "Missing positive VHC"),
], ids=["budget", "origin", "composition", "bounds", "lift", "replay-cop", "vhc"])
def test_validate_trial_rejects_inconsistent_trials(mutate, message):
    trial = _origin_trial("gp_ei", 1, [2.0, 2.5, 2.4])
    mutate(trial)
    with pytest.raises(ValueError, match=message):
        origin_report.validate_trial(trial)


# --- optimiser_benchmark_report -------------------------------------------

def _pilot_trial(method, seed, cops, *, replay_cop, warm_start=False, failed=(), replay_success=True):
    records = []
    for i, cop in enumerate(cops):
        success = i not in failed
        records.append({"x": [float(i), 1.0, 0.5], "success": success,
                        "COP": cop if success else None, "SEC": 1 / cop if success else None,
                        "elapsed_seconds": 2.0, "proposal_seconds": 0.5})
    return {
        "metadata": {"refrigerant": REF, "method": method, "seed": seed, "budget": len(cops),
                     "baseline_COP": 2.0, "baseline_evaluations": 120, "warm_start": warm_start},
        "records": records,
        "replay": {"success": replay_success, "COP": replay_cop, "VHC_MJ_per_m3": 3.5,
                   "x": [1.25, 2.5, 0.45], "Pmax_MPa": 4.25},
    }


def test_benchmark_report_tables_use_replay_baseline_and_record_columns(tmp_path):
    trials = {
        "a_rbf.json": _pilot_trial("rbf", 7, [1.9, 2.05, 2.2], replay_cop=2.211),
        "b_gp.json": _pilot_trial("gp_ei", 7, [1.9, 2.1, 2.0], replay_cop=2.1, failed={2}),
        "c_warm.json": _pilot_trial("pattern", 0, [2.0, 2.3], replay_cop=2.3, warm_start=True),
        "d_unreplayed.json": {**_pilot_trial("space_filling", 7, [2.0], replay_cop=2.0), "replay": {}},
        "e_other.json": {"not": "a trial"},
    }
    for name, trial in trials.items():
        (tmp_path / name).write_text(json.dumps(trial))

    path = optimiser_report.write_benchmark_report(tmp_path)
    markdown = path.read_text()

    cold = _table_rows(markdown, "| Composition by mass | Method | Calls")
    assert cold == [
        ["R290/R601 30/70", "RBF", "3", "120", "2.000000", "2.211000", "+10.550%", "3.5000", "6.0 / 1.5", "0"],
        ["R290/R601 30/70", "GP expected improvement", "3", "120", "2.000000", "2.100000", "+5.000%",
         "3.5000", "6.0 / 1.5", "1"],
    ]
    warm = _table_rows(markdown, "| Composition by mass | Extra calls")
    assert warm == [["R290/R601 30/70", "2", "2.000000", "2.300000", "+15.000%"]]

    designs = _table_rows(markdown, "| Composition by mass | Method | RHX1")
    assert [row[1] for row in designs] == ["RBF", "GP expected improvement", "Warm joint polish"]
    assert designs[0][2:6] == ["1.250000", "2.500000", "0.450000", "4.2500"]
    # Drift compares replay COP with the lowest-SEC successful record.
    assert designs[0][6] == "+0.500000%"  # 2.211 vs best record 2.2
    assert designs[1][6] == "+0.000000%"
    assert (tmp_path / "optimiser_convergence.png").stat().st_size > 0
    assert not (tmp_path / "pentane30_landscape.png").exists()


def test_benchmark_report_shows_nan_for_failed_cold_replay(tmp_path):
    (tmp_path / "rbf.json").write_text(json.dumps(
        _pilot_trial("rbf", 1, [2.0, 2.1], replay_cop=None, replay_success=False)))
    markdown = optimiser_report.write_benchmark_report(tmp_path).read_text()
    cold = _table_rows(markdown, "| Composition by mass | Method | Calls")
    assert cold[0][5:7] == ["nan", "+nan%"]
    assert _table_rows(markdown, "| Composition by mass | Method | RHX1") == []


def test_benchmark_report_requires_cold_start_trial(tmp_path):
    (tmp_path / "warm.json").write_text(json.dumps(
        _pilot_trial("pattern", 0, [2.0], replay_cop=2.0, warm_start=True)))
    with pytest.raises(ValueError, match="No completed cold-start trials"):
        optimiser_report.write_benchmark_report(tmp_path)


def test_benchmark_report_plots_landscape_when_present(tmp_path):
    (tmp_path / "rbf.json").write_text(json.dumps(_pilot_trial("rbf", 1, [2.0, 2.1], replay_cop=2.1)))
    (tmp_path / "pentane30_landscape.json").write_text(json.dumps({"records": [
        {"interpolation_fraction": f, "fast_inner_objective": True, "success": True, "COP": 2 + f, "Pmax_MPa": 4}
        for f in (0.0, 0.5, 1.0)
    ]}))
    optimiser_report.write_benchmark_report(tmp_path)
    assert (tmp_path / "pentane30_landscape.png").stat().st_size > 0
