"""Summarise the completed, independently replayed local GP subset trial."""
import hashlib
import json
from pathlib import Path
import statistics

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from rhp_local_pareto_search import hypervolume

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results/RHPcomposition/local_pareto_subset_trial'
PARTNERS=('Butane','Isopentane','Pentane','Isobutane')
METHODS=('128','256','full')
COLORS=('#087f8c','#db7c26','#6554c0')


def main():
    baseline_manifest=json.loads((OUT/'baseline_manifest.json').read_text())
    main_out=ROOT/'results/RHPcomposition/mobo_v17_4000'
    for name,digest in baseline_manifest.items():
        assert hashlib.sha256((OUT/'baseline'/name).read_bytes()).hexdigest()==digest
        assert hashlib.sha256((main_out/name).read_bytes()).hexdigest()==digest, 'Main search changed during trial'
    assert json.loads((main_out/'runner_state.json').read_text())['status']=='paused'
    rows=[]
    traces={}
    for partner in PARTNERS:
        baseline=json.loads((OUT/'baseline'/f'Propane_{partner}.json').read_text())
        values=[[r['COP'],r['VHC_MJ_per_m3']] if r['success'] and not r.get('validation_rejected') else [0.,0.] for r in baseline['records']]
        initial=hypervolume(values)
        for method in METHODS:
            data=json.loads((OUT/f'{partner}_{method}.json').read_text())
            assert data['complete'] and len(data['records'])==data['metadata']['budget']==12
            for i in data['new_front_indices']:
                assert data['replays'][str(i)]['matches']
            records=data['records']
            prediction_errors=[]
            for r in records:
                diag=r['search_diagnostic']
                if r['success'] and 'predicted_objectives' in diag:
                    prediction_errors.append(np.abs(np.asarray(diag['predicted_objectives'])/[r['COP'],r['VHC_MJ_per_m3']]-1)*100)
            t=[0.]
            history=list(values)
            for r in records:
                history.append([r['COP'],r['VHC_MJ_per_m3']] if r['success'] and not r.get('validation_rejected') else [0.,0.])
                t.append(100*(hypervolume(history)/initial-1))
            traces[(partner,method)]=t
            row=dict(partner=partner,method=method,evaluations=len(records),
                     successful=sum(r['success'] for r in records),validated_new_front_points=len(data['new_front_indices']),
                     replay_evaluations=len(data['replays']),
                     rejected_replays=sum(not r['matches'] for r in data['replays'].values()),
                     mean_proposal_seconds=statistics.mean(r['proposal_seconds'] for r in records),
                     median_proposal_seconds=statistics.median(r['proposal_seconds'] for r in records),
                     total_proposal_seconds=sum(r['proposal_seconds'] for r in records),
                     total_oracle_seconds=sum(r['evaluation_wall_seconds'] for r in records),
                     peak_parent_RSS_MiB=data['peak_rss_MiB'],hypervolume_gain_percent=t[-1],
                     fit_seconds=sum(r['search_diagnostic'].get('fit_seconds',0) for r in records),
                     pool_score_seconds=sum(r['search_diagnostic'].get('pool_score_seconds',0) for r in records),
                     acquisition_opt_seconds=sum(r['search_diagnostic'].get('acquisition_opt_seconds',0) for r in records),
                     fit_count=sum(r['search_diagnostic'].get('hyperparameters_fitted',False) for r in records),
                     highest_COP_before=max(v[0] for v in values),highest_COP_after=max(v[0] for v in history),
                     prediction_mean_absolute_percent_error=np.mean(prediction_errors,axis=0).tolist() if prediction_errors else None)
            rows.append(row)
    for row in rows:
        full=next(r for r in rows if r['partner']==row['partner'] and r['method']=='full')
        row['speedup_vs_full']=full['total_proposal_seconds']/row['total_proposal_seconds']
    summary=dict(rows=rows,main_search_unchanged=True,search_evaluations=sum(r['evaluations'] for r in rows),
                 replay_evaluations=sum(r['replay_evaluations'] for r in rows))
    (OUT/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
    for ax,key,label in zip(axes,('mean_proposal_seconds','peak_parent_RSS_MiB','hypervolume_gain_percent'),
                           ('Mean candidate-selection time (s)','Peak search-process memory (MiB)','Hypervolume gain (%)')):
        for j,(method,color) in enumerate(zip(METHODS,COLORS)):
            y=[next(r[key] for r in rows if r['partner']==p and r['method']==method) for p in PARTNERS]
            ax.bar(np.arange(4)+(j-1)*.25,y,width=.23,color=color,label=f'{method}-point GP' if method!='full' else 'Full-data GP')
        ax.set_xticks(np.arange(4),PARTNERS,rotation=25,ha='right')
        ax.set_ylabel(label)
        ax.grid(axis='y',alpha=.2)
        ax.set_axisbelow(True)
        if key=='mean_proposal_seconds':
            ax.set_yscale('log')
    axes[0].legend(fontsize=9)
    fig.suptitle('Local Pareto refinement · identical regions · 12 new candidates per arm/family')
    fig.savefig(OUT/'comparison.png',dpi=180)
    fig.savefig(OUT/'comparison.svg')
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for ax,partner in zip(axes.flat,PARTNERS):
        for method,color in zip(METHODS,COLORS):
            ax.plot(range(13),traces[(partner,method)],'o-',ms=3,color=color,label=method)
        ax.set(title=f'Propane/{partner}',xlabel='New candidate evaluations',ylabel='Hypervolume gain (%)')
        ax.grid(alpha=.2)
        ax.legend(title='GP training points',fontsize=9)
    fig.suptitle('Improvement over the full paused baseline; new retained points replay-validated')
    fig.savefig(OUT/'hypervolume_progress.png',dpi=180)
    plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    for ax,partner in zip(axes.flat,PARTNERS):
        baseline=json.loads((OUT/'baseline'/f'Propane_{partner}.json').read_text())
        data=json.loads((OUT/f'{partner}_128.json').read_text())
        front=[baseline['records'][i] for i in baseline['frontier_indices']]
        ax.scatter([100*r['x'][3] for r in front],[r['x'][1] for r in front],s=14,color='#aab2bd',label='Observed frontier')
        for j,r in enumerate(data['records'][:6]):
            lower,upper=np.asarray(r['search_diagnostic']['bounds'])
            centre=np.asarray(r['search_diagnostic']['centre'])
            ax.add_patch(Rectangle((100*lower[3],100*lower[1]),100*(upper[3]-lower[3]),100*(upper[1]-lower[1]),
                                   facecolor='#087f8c',alpha=.12,edgecolor='#087f8c'))
            ax.plot(100*centre[3],100*centre[1],'+',color='#087f8c')
            ax.annotate(str(j+1),(100*centre[3],100*centre[1]),xytext=(4,4),textcoords='offset points',fontsize=8)
        ax.set(title=f'Propane/{partner}',xlabel='Propane mass fraction (%)',ylabel='RHX2 cold-side lift (K)',xlim=(-2,102))
        ax.grid(alpha=.15)
    fig.suptitle('Six local search boxes per mixture (projection)\nRHX1 and compressor split are also bounded within each box')
    fig.savefig(OUT/'local_regions.png',dpi=180)
    plt.close(fig)
    text=['# Local Pareto GP subset trial','',
          f"Completed {summary['search_evaluations']} new candidate evaluations and {summary['replay_evaluations']} extra frontier replays. "
          'The main searches remain paused and their checkpoints are byte-for-byte unchanged.','',
          'Six small fixed neighbourhoods along each mixture front were searched, twice each. '
          'Subsets use the nearest 128 or 256 successful designs. The full-data comparator searches the same boxes. '
          'All observations remain in the hypervolume baseline and failure-distance heuristic.','',
          '| Propane partner | GP training set | Mean proposal s | Speedup vs full | Peak MiB | Successful / 12 | Validated new front points | HV gain % |',
          '|---|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:
        text.append(f"| {r['partner']} | {r['method']} | {r['mean_proposal_seconds']:.2f} | {r['speedup_vs_full']:.2f}× | "
                    f"{r['peak_parent_RSS_MiB']:.0f} | {r['successful']} | {r['validated_new_front_points']} | {r['hypervolume_gain_percent']:.5f} |")
    text+=['','![Runtime, memory and improvement](comparison.png)','','![Hypervolume progress](hypervolume_progress.png)','',
           '![Local search regions](local_regions.png)','',
           '## Interpretation and limits','',
           'This is a one-seed, 12-candidate local refinement trial per arm and family. Different arms generate different candidates; '
           'the common starting data, anchors, boxes, budgets and acquisition settings are controlled. It does not measure global discovery '
           'or guarantee that the same subset size will remain best over a longer run. Candidate selection includes model fitting and acquisition optimisation. '
           'The full-data model refits every five BO visits; local models have separate per-anchor fit caches, as detailed in the protocol.','',
           'Scoring batches are eight and restart batches are one in every arm. Consequently the full-data comparator also incorporates a memory-oriented '
           'batch-size change relative to the original global run. Peak RSS measures the search process, excluding its separate native evaluation worker. '
           'Fresh-worker startup is included in oracle wall time. New retained frontier points are replayed; the historical baseline has not undergone '
           'complete frontier validation in this trial.','',
           '[Detailed protocol](PROTOCOL.md) · [Machine-readable summary](summary.json)']
    (OUT/'REPORT.md').write_text('\n'.join(text)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    main()
