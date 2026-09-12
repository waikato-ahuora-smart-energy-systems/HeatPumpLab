"""Compare the nine validated screening fronts with the completed propane runs."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from main_RHPmobo import OUTPUT, PAIRS, frontier_indices

ROOT = Path(__file__).resolve().parent
REFERENCES = ROOT.parents[1] / "results/RHPcomposition/mobo_v17_4000"


def hypervolume(records):
    front = sorted(
        (r for r in records if r["success"] and not r.get("validation_rejected")),
        key=lambda r: -r["COP"],
    )
    total, high = 0.0, 0.0
    for r in front:
        if r["VHC_MJ_per_m3"] > high:
            total += (r["VHC_MJ_per_m3"] - high) * r["COP"]
            high = r["VHC_MJ_per_m3"]
    return total


def build_report(output=OUTPUT):
    output = Path(output)
    references = []
    for partner in ("Butane", "Isopentane", "Pentane", "Isobutane"):
        data = json.loads((REFERENCES / f"Propane_{partner}.json").read_text())
        assert data["complete"]
        references.extend(data["records"][i] for i in data["frontier_indices"])
    reference_hv = hypervolume(references)
    datasets = {k: json.loads((output / f"{k}.json").read_text()) for k in PAIRS}
    summaries = []
    fronts = {}
    for key, data in datasets.items():
        assert data["complete"] and len(data["records"]) == 800
        assert frontier_indices(data["records"]) == data["frontier_indices"]
        selected = []
        for i in data["frontier_indices"]:
            original = data["records"][i]
            checked = data["frontier_replays"][str(i)]
            assert checked["success"] and np.allclose(
                checked["x"], original["x"], rtol=0, atol=1e-10
            )
            assert all(
                np.isclose(checked[k], original[k], rtol=1e-6, atol=1e-9)
                for k in (
                    "COP",
                    "VHC_MJ_per_m3",
                    "SEC",
                    "Pmax_MPa",
                    "max_reported_refrigerant_C",
                    "PR_1",
                    "PR_2",
                )
            )
            selected.append(dict(checked, evaluation=i + 1))
        fronts[key] = selected
        good = [
            r
            for r in data["records"]
            if r["success"] and not r.get("validation_rejected")
        ]
        summaries.append(
            dict(
                pair=key,
                components=list(PAIRS[key]),
                estimated_mixture_parameters=data["metadata"]["property_model_status"][
                    "estimated_mixture_parameters"
                ],
                evaluations=800,
                failures=sum(not r["success"] for r in data["records"]),
                validated_front_points=len(selected),
                replay_evaluations=len(data["frontier_replays"]),
                validation_rejections=sum(
                    bool(r.get("validation_rejected")) for r in data["records"]
                ),
                highest_COP=max(selected, key=lambda r: r["COP"]) if selected else None,
                highest_VHC=max(selected, key=lambda r: r["VHC_MJ_per_m3"])
                if selected
                else None,
                hypervolume=hypervolume(good),
                gain_over_reference_front_percent=100
                * (hypervolume(references + selected) / reference_hv - 1),
                median_proposal_seconds=float(
                    np.median([r["proposal_seconds"] for r in data["records"][32:]])
                ),
                maximum_reported_temperature_C=max(
                    (r["max_reported_refrigerant_C"] for r in good), default=None
                ),
            )
        )
    shortlist = sorted(
        summaries,
        key=lambda s: (s["gain_over_reference_front_percent"], s["hypervolume"]),
        reverse=True,
    )[:3]
    (output / "summary.json").write_text(
        json.dumps(
            dict(
                families=summaries,
                shortlist=[s["pair"] for s in shortlist],
                reference_hypervolume=reference_hv,
                validated_fronts=fronts,
            ),
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    colors = plt.get_cmap("tab10").colors
    fig, ax = plt.subplots(figsize=(11, 7), layout="constrained")
    old = [references[i] for i in frontier_indices(references)]
    old.sort(key=lambda r: r["VHC_MJ_per_m3"])
    ax.plot(
        [r["VHC_MJ_per_m3"] for r in old],
        [r["COP"] for r in old],
        "k--",
        lw=2,
        label="Previous four families (4,000 each)",
    )
    for (key, rows), color in zip(fronts.items(), colors):
        rows = sorted(rows, key=lambda r: r["VHC_MJ_per_m3"])
        ax.plot(
            [r["VHC_MJ_per_m3"] for r in rows],
            [r["COP"] for r in rows],
            ".-",
            ms=3,
            lw=1,
            color=color,
            label="/".join(PAIRS[key]),
        )
    ax.set(
        xlabel="Total volumetric heating capacity (MJ/m³)",
        ylabel="Heating COP",
        title="Nine additional binary families · 800 global BO evaluations each",
    )
    ax.grid(alpha=0.2)
    ax.legend(fontsize=8)
    fig.savefig(output / "pareto_fronts.png", dpi=180)
    fig.savefig(output / "pareto_fronts.svg")
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    for summary, color in zip(summaries, colors):
        r = summary["highest_COP"]
        if r:
            axes[0].scatter(
                r["max_reported_refrigerant_C"],
                r["COP"],
                color=color,
                label="/".join(summary["components"]),
            )
            axes[1].scatter(r["Pmax_MPa"], r["COP"], color=color)
    axes[0].set(
        xlabel="Maximum reported refrigerant temperature (°C)",
        ylabel="COP",
        title="Each family’s maximum-COP design",
    )
    axes[1].set(
        xlabel="Maximum pressure (MPa)",
        ylabel="COP",
        title="Pressure at the same designs",
    )
    for ax in axes:
        ax.grid(alpha=0.2)
    axes[0].legend(fontsize=7)
    fig.savefig(output / "best_case_conditions.png", dpi=180)
    plt.close(fig)
    lines = [
        "# Nine additional binary-fluid screens",
        "",
        "All nine families completed 800 global evaluations each (7,200 total). All retained COP–VHC_total frontier designs were independently replayed. "
        "The previous four propane-family fronts are references with 4,000 evaluations each; this is not an equal-budget algorithm comparison.",
        "",
        "Cyclopentane blends use alternative mixing rules because fitted binary interaction parameters are unavailable in the installed REFPROP database. Their predictions are exploratory, even when replays match.",
        "",
        "| Pair | Failed / 800 | Validated front | Best COP | VHC at best COP (MJ/m³) | Max reported T (°C) | Max P (MPa) | PR1 / PR2 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for s in summaries:
        r = s["highest_COP"]
        if r:
            lines.append(
                f"| {'/'.join(s['components'])} | {s['failures']} | {s['validated_front_points']} | {r['COP']:.6f} | {r['VHC_MJ_per_m3']:.4f} | {r['max_reported_refrigerant_C']:.1f} | {r['Pmax_MPa']:.3f} | {r['PR_1']:.2f} / {r['PR_2']:.2f} |"
            )
        else:
            lines.append(
                f"| {'/'.join(s['components'])} | {s['failures']} | 0 | — | — | — | — | — |"
            )
    lines += [
        "",
        "![Validated screening fronts](pareto_fronts.png)",
        "",
        "![Temperatures and pressures at maximum COP](best_case_conditions.png)",
        "",
        "## Local-refinement shortlist",
        "",
        "Cyclopentane candidates require property-model validation before treating small predicted advantages as established.",
        "",
        "Provisional ranking by incremental hypervolume beyond the completed reference front, with family hypervolume breaking ties. "
        "Temperatures and pressure ratios must also be considered when selecting a design; this ranking imposes no unprovided equipment or chemical-life limits.",
        "",
    ]
    for s in shortlist:
        lines.append(
            f"- {'/'.join(s['components'])}: additional reference-front hypervolume {s['gain_over_reference_front_percent']:.5f}%."
        )
    lines += [
        "",
        "The 128-point local method is available for the selected fronts. No additional 4,000-evaluation runs have been launched automatically.",
        "",
        "[Protocol](PROTOCOL.md) · [Full numerical results](summary.json) · [Frozen source](source_manifest.json)",
        "",
        "Maximum reported refrigerant temperature is the maximum of reported compressor outlets, final discharge, gas-cooler outlet and evaporator outlet. "
        "It is not a heat-exchanger wall temperature or a chemical-stability limit. Chemical degradation, lubricant compatibility and lifetime are not modelled. "
        "Pure endpoints are included; composition x[3] always refers to the first named component. Lines connect sampled designs.",
    ]
    (output / "REPORT.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    build_report()
