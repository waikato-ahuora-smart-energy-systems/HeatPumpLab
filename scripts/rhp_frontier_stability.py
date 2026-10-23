"""Progress checks for a fixed-temperature COP/total-VHC frontier."""
import numpy as np


def frontier(records):
    good = [dict(evaluation=i+1, **r) for i, r in enumerate(records)
            if r['success'] and not r.get('validation_rejected')
            and r.get('COP') is not None and r.get('VHC_MJ_per_m3') is not None
            and np.isfinite([r['COP'], r['VHC_MJ_per_m3']]).all()
            and min(r['COP'], r['VHC_MJ_per_m3']) > 0]
    ordered = sorted(good, key=lambda r: (-r['VHC_MJ_per_m3'], -r['COP']))
    selected, best = [], -np.inf
    for record in ordered:
        if record['COP'] > best:
            selected.append(record)
            best = record['COP']
    return selected[::-1]


def assess(records, previous=None, *, tolerance=.001, minimum_evaluations=2000, hv_only=False,
           required_blocks=3):
    """Use a fixed reference and VHC thresholds; never interpolate observations.

    Thresholds are fixed by the first check. Endpoint changes also guard against
    growth outside that original range. Three stable blocks are required.
    """
    points = frontier(records)
    vhc = np.array([p['VHC_MJ_per_m3'] for p in points])
    cop = np.array([p['COP'] for p in points])
    hv = float(np.dot(np.diff(np.r_[0., vhc]), cop)) if points else 0.
    thresholds = (previous['vhc_thresholds'] if previous else
                  np.linspace(0, float(vhc.max()) if points else 1., 101).tolist())
    envelope = [float(cop[vhc >= t].max()) if np.any(vhc >= t) else None for t in thresholds]
    gain, curve_change, extent_change = None, None, None
    stable = False
    best_cop = float(cop.max()) if points else None
    max_vhc = float(vhc.max()) if points else None
    if previous and previous['hypervolume'] > 0 and points:
        gain = hv/previous['hypervolume']-1
        changes = [abs(new/old-1) if new is not None and old is not None else
                   (0. if new is None and old is None else 1.)
                   for new, old in zip(envelope, previous['cop_at_vhc_thresholds'])]
        curve_change = max(changes)
        extent_change = abs(max_vhc/previous['max_vhc']-1)
        recent = records[previous['evaluations']:]
        # Sparse or failing searches do not establish a stable curve.
        healthy = len(points) >= 10 and recent and sum(r['success'] for r in recent)/len(recent) >= .8
        stable = bool(healthy and abs(gain) < tolerance and (
            hv_only or (curve_change < tolerance and extent_change < tolerance)))
    streak = previous['stable_blocks']+1 if previous and stable else 0
    return dict(evaluations=len(records), hypervolume=hv, reference_point=[0., 0.],
                hypervolume_relative_gain=gain, curve_max_relative_change=curve_change,
                vhc_extent_relative_change=extent_change, vhc_thresholds=thresholds,
                cop_at_vhc_thresholds=envelope, best_cop=best_cop, max_vhc=max_vhc,
                frontier_points=len(points), frontier=points,
                failures=sum(not r['success'] for r in records), stable_block=stable,
                stable_blocks=streak, stop_for_stability=bool(streak >= required_blocks and len(records) >= minimum_evaluations),
                tolerance=tolerance, label='Observed fixed-temperature frontier; final replay validation pending')


def render(snapshot, path, partner, first='Propane'):
    import plotly.graph_objects as go
    points = snapshot['frontier']
    figure = go.Figure(go.Scatter(
        x=[r['VHC_MJ_per_m3'] for r in points], y=[r['COP'] for r in points],
        mode='lines+markers', marker=dict(color=[r['x'][3] for r in points],
            colorscale='Viridis', cmin=0, cmax=1, colorbar=dict(title=f'{first} mass fraction')),
        customdata=[[r['evaluation'], r['x'][3], r.get('Pmax_MPa')] for r in points],
        hovertemplate='VHC: %{x:.5f} MJ/m³<br>COP: %{y:.5f}<br>'+first+': %{customdata[1]:.5f}'
                      '<br>Pressure: %{customdata[2]:.3f} MPa<br>Evaluation: %{customdata[0]}<extra></extra>'))
    figure.update_layout(title=f"{first} / {partner} at 200°C · {snapshot['evaluations']:,} evaluations",
                         xaxis_title='VHC_total (MJ/m³)', yaxis_title='COP', template='plotly_white',
                         annotations=[dict(text=snapshot['label']+'; connecting lines are visual guides',
                                           xref='paper', yref='paper', x=0, y=-.18, showarrow=False)],
                         margin=dict(b=120))
    figure.write_html(str(path), include_plotlyjs=True)
