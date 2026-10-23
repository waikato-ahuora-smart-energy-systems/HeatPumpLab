"""Interactive evaluated 3D fronts, with explicitly optional interpolation."""
from pathlib import Path

import numpy as np
import plotly.graph_objects as go
from scipy.spatial import Delaunay, QhullError


def render_front(snapshot, path):
    points = snapshot['frontier']
    figure = go.Figure()
    if points:
        temperature = np.array([r['T_air_tar'] for r in points])
        vhc = np.array([r['VHC_MJ_per_m3'] for r in points])
        cop = np.array([r['COP'] for r in points])
        if len(np.unique(temperature)) > 1 and len(points) >= 3:
            _, unique = np.unique(np.column_stack((temperature, vhc)), axis=0, return_index=True)
            try:
                triangles = unique[Delaunay(np.column_stack((temperature[unique], vhc[unique]))).simplices]
                # Never bridge a missing temperature slice. Faces are only a visual aid.
                triangles = triangles[np.ptp(temperature[triangles], axis=1) <= 5.000001]
                if len(triangles):
                    figure.add_trace(go.Mesh3d(x=temperature, y=vhc, z=cop,
                        i=triangles[:, 0], j=triangles[:, 1], k=triangles[:, 2],
                        opacity=.18, color='#386c91', name='Interpolated faces (not evaluated)',
                        hoverinfo='skip', visible='legendonly', showlegend=True))
            except QhullError:
                pass
        hover = [[r['evaluation'],r['x'][3],r.get('Pmax_MPa'),r['x'][0],r['x'][1],r['x'][2]] for r in points]
        figure.add_trace(go.Scatter3d(x=temperature, y=vhc, z=cop, mode='markers',
            marker=dict(size=3, color=[r['x'][3] for r in points], colorscale='Viridis',
                        cmin=0,cmax=1,colorbar=dict(title='Propane<br>mass fraction')),
            customdata=hover, name='Evaluated Pareto points',
            hovertemplate='Evaluation %{customdata[0]}<br>Air: %{x} °C<br>COP: %{z:.5f}<br>VHC: %{y:.5f} MJ/m³'
                          '<br>Propane: %{customdata[1]:.4f}<br>Pressure: %{customdata[2]:.3f} MPa'
                          '<br>RHX lifts: %{customdata[3]:.2f}, %{customdata[4]:.2f} K'
                          '<br>Stage allocation: %{customdata[5]:.4f}<extra></extra>'))
    figure.update_layout(title=f"{Path(path).parent.name.replace('_', ' / ')} · {snapshot['evaluations']:,} evaluations",
        scene=dict(xaxis_title='Target air temperature (°C)',yaxis_title='VHC_total (MJ/m³)',zaxis_title='COP'),
        margin=dict(l=20,r=20,b=20,t=75),height=750,
        annotations=[dict(text=snapshot['label'],xref='paper',yref='paper',x=0,y=1.07,showarrow=False)])
    figure.write_html(str(Path(path)),include_plotlyjs=True)
