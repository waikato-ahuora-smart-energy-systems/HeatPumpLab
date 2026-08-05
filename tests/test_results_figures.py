from collections import Counter
from pathlib import Path

import numpy as np

from paper.figures.generate_results_figures import (
    allocation,
    build_base_ternary_comparisons,
    build_nesting_comparisons,
    build_target_envelopes,
    generate_all,
    load_results,
    primary_results,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "results" / "RHPopt" / "RHPopt_ncomponent_results_v11.csv"


def test_results_figure_data_matches_section_four_audit():
    df = load_results(SOURCE)
    direct, reverse = build_nesting_comparisons(df)
    ternary, global_leader = build_base_ternary_comparisons(df)
    envelope = build_target_envelopes(df)

    assert len(df) == 1650
    assert direct["order"].value_counts().to_dict() == {3: 875, 2: 525}
    assert len(reverse) == 526
    assert int(reverse["sec_change_pct"].lt(-0.1).sum()) == 55
    assert int(reverse["sec_change_pct"].lt(-1.0).sum()) == 13

    assert len(ternary) == 35
    assert int(ternary["retained_three"].sum()) == 18
    assert int((ternary["retained_three"] & ternary["cop_change_pct"].gt(0.0)).sum()) == 4
    assert global_leader["refrigerant_components"] == "Propane;Pentane"

    warning_targets = envelope.loc[
        envelope["role"].eq("pooled") & envelope["warning"], "target"
    ].tolist()
    assert warning_targets == [140, 155, 160]
    assert envelope.groupby("role").size().to_dict() == {
        "binary": 25,
        "pooled": 25,
        "pure": 25,
        "three": 25,
    }

    base = primary_results(df)
    base = base[np.isclose(base["T_air_tar"], 200.0)]
    counts = Counter(allocation(row) for _, row in base.iterrows())
    assert counts == Counter({"RHX2 edge": 31, "RHX1 edge": 17, "Zero": 13, "Interior": 2})


def test_results_figure_generator_writes_all_four_graphs(tmp_path):
    outputs = generate_all(SOURCE, tmp_path, formats=("png",), dpi=72)

    assert {path.name for path in outputs} == {
        "results-nesting-audit.png",
        "results-ternary-ranking.png",
        "results-target-temperature.png",
        "results-recuperator-allocation.png",
    }
    assert all(path.stat().st_size > 10_000 for path in outputs)
