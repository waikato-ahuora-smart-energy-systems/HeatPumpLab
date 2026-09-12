"""Plot read-only observed COP/VHC or COP/pressure Pareto snapshots."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT/'experiments/mobo_v17_clean_runtime'))
from main_RHPmobo import PARTNERS, frontier_indices  # noqa: E402

COLORS = ('#087f8c', '#db7c26', '#6554c0', '#b83d77')
NAMES = {'Butane': 'Butane (R600)', 'Isopentane': 'Isopentane (R601a)',
         'Pentane': 'Pentane (R601)', 'Isobutane': 'Isobutane (R600a)'}


def pressure_front_indices(records):
    """Maximise COP and minimise pressure; preserve exact ties."""
    valid = [i for i,r in enumerate(records) if r['success'] and not r.get('validation_rejected')
             and r.get('Pmax_MPa') is not None and r.get('COP') is not None
             and np.isfinite([r['Pmax_MPa'],r['COP']]).all() and min(r['Pmax_MPa'],r['COP']) > 0]
    if not valid:
        return []
    pressure = np.asarray([records[i]['Pmax_MPa'] for i in valid])
    cop = np.asarray([records[i]['COP'] for i in valid])
    return [i for j,i in enumerate(valid) if not np.any(
        (pressure <= pressure[j]) & (cop >= cop[j]) & ((pressure < pressure[j]) | (cop > cop[j])))]


def main(source, objective='vhc'):
    pressure_mode = objective == 'pressure'
    key = 'Pmax_MPa' if pressure_mode else 'VHC_MJ_per_m3'
    select_front = pressure_front_indices if pressure_mode else frontier_indices
    filename = 'pareto_cop_pressure' if pressure_mode else 'pareto_fronts'
    stamp = datetime.datetime.now(datetime.timezone.utc)
    output = source/'snapshots'/stamp.strftime('%Y%m%dT%H%M%SZ')
    output.mkdir(parents=True, exist_ok=False)
    datasets, provenance, fronts = {}, {}, {}
    for partner in PARTNERS:
        path = source/f'Propane_{partner}.json'
        raw = path.read_bytes()  # Checkpoints are atomically replaced by the search.
        data = json.loads(raw)
        datasets[partner] = data
        (output/path.name).write_bytes(raw)
        provenance[partner] = dict(evaluations=len(data['records']), sha256=hashlib.sha256(raw).hexdigest())
        indices = select_front(data['records'])
        if not pressure_mode:
            assert indices == data['frontier_indices'], partner
        fronts[partner] = sorted((data['records'][i] for i in indices), key=lambda r:r[key])
    total = sum(len(d['records']) for d in datasets.values())
    all_front = [r for front in fronts.values() for r in front]
    combined = [all_front[i] for i in select_front(all_front)]
    combined = sorted({(r[key],r['COP']) for r in combined})
    best = max(all_front,key=lambda r:r['COP'])

    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11, 'axes.labelcolor':'#263446',
                         'text.color':'#263446','axes.spines.top':False,'axes.spines.right':False})
    fig, ax = plt.subplots(figsize=(11,7))
    fig.subplots_adjust(left=.09,right=.97,top=.83,bottom=.29)
    fig.text(.09,.945,'COP–maximum pressure Pareto fronts' if pressure_mode else 'COP–VHC_total Pareto fronts',
             fontsize=21,weight='bold')
    fig.text(.09,.9,f'Interim snapshot · {total:,} / 16,000 evaluations · 200°C sink',fontsize=12,color='#596779')
    for partner, color in zip(PARTNERS,COLORS):
        front = fronts[partner]
        ax.plot([r[key] for r in front],[r['COP'] for r in front],
                color=color,lw=1.8,marker='o',ms=3,alpha=.95,zorder=3,
                label=f"Propane / {NAMES[partner]} · n={len(datasets[partner]['records']):,}")
    ax.plot(*np.asarray(combined).T, color='#263446',lw=1,ls=(0,(4,4)),zorder=4,
            label='Combined observed front')
    ax.annotate(f"Highest observed COP: {best['COP']:.5f}",
                xy=(best[key],best['COP']),xytext=((-20 if pressure_mode else 20),15),
                ha='right' if pressure_mode else 'left',textcoords='offset points',
                fontsize=10,color='#4b3d99',arrowprops=dict(arrowstyle='-',color='#6554c0',lw=.8))
    ax.set(xlabel='Maximum cycle pressure (MPa) — lower is better' if pressure_mode else
           'Total volumetric heating capacity, VHC_total (MJ/m³)',ylabel='Heating COP — higher is better')
    ax.margins(x=.04,y=.13)
    ax.grid(color='#dce2e9',lw=.7,alpha=.75)
    ax.set_axisbelow(True)
    fig.legend(*ax.get_legend_handles_labels(),loc='upper left',bbox_to_anchor=(.085,.195),
               ncol=2,frameon=False,fontsize=9.5)
    fig.text(.09,.075,'Re-ranked all successful COP–VHC search samples; pressure was not an optimisation objective.'
             if pressure_mode else 'Observed nondominated candidates; final fresh-process validation is pending.',
             fontsize=10,color='#596779')
    fig.text(.09,.043,'Observed fronts; fresh validation pending. Lines connect sampled designs.' if pressure_mode else
             'Lines connect sampled designs. VHC_total = Q_H_total / (V_suc_1 + V_suc_2).',fontsize=9,color='#596779')
    fig.text(.97,.043,stamp.strftime('%d %b %Y · %H:%M UTC'),ha='right',fontsize=9,color='#596779')
    for extension in ('png','svg'):
        fig.savefig(output/f'{filename}.{extension}',dpi=200,facecolor='white')
    plt.close(fig)
    summary=dict(snapshot_utc=stamp.isoformat(),total_evaluations=total,validation='pending',
                 families=provenance,front_points={p:len(f) for p,f in fronts.items()},
                 highest_observed_COP=best['COP'],highest_COP_design=best,
                 objectives={'COP':'maximise',key:'minimise' if pressure_mode else 'maximise'},
                 family_fronts=fronts,
                 combined_front=[{key:v,'COP':c} for v,c in combined])
    (output/'snapshot.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(output),total_evaluations=total,families=provenance),indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,default=ROOT/'results/RHPcomposition/mobo_v17_4000')
    parser.add_argument('--objective',choices=('vhc','pressure'),default='vhc')
    args=parser.parse_args()
    main(args.source.resolve(),args.objective)
