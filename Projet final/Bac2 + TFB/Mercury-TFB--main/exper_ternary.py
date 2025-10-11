"""
This script generates an interactive ternary diagram showing experimental compositions
(Mer8 and Mer15) in relation to their partial melting rate (F). The plot features:
• Experimental points colored by F (%)
• Smoothed experimental trend curves
• Separate colorbars per pressure level
• Black axis borders and gridlines
All data are plotted using Plotly, with ternary coordinates normalized to 100%.
"""



import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from scipy.interpolate import UnivariateSpline
import plotly.graph_objects as go
import os

DATA_FOLDER = "data"

# ============================================================
# Utility functions
# ============================================================

# Inverse Distance Weighting (IDW) interpolation of scalar values
def idw_interpolation(query_points, points, values, power=2):
    interpolated = []
    for qp in query_points:
        d = np.linalg.norm(points - qp, axis=1)
        if np.any(d < 1e-8):  # Exact match found
            interpolated.append(values[d.argmin()])
        else:
            w = 1 / d ** power
            interpolated.append(np.sum(w * values) / np.sum(w))
    return np.array(interpolated)

# Fit a curve using PCA projection + spline interpolation, then reproject
def compute_curve(df_subset, n_points=200, power=2):
    X = df_subset[["Mg/Si", "Ca/Si", "Al/Si"]].values
    F_exp = df_subset["F"].values
    t = PCA(1).fit_transform(X).flatten()  # Reduce to 1D projection
    order = np.argsort(t)
    t_sorted = t[order]
    # Fit a spline to each ternary axis
    splines = [UnivariateSpline(t_sorted, X[order, i], s=len(t_sorted)) for i in range(3)]
    t_fine = np.linspace(t_sorted.min(), t_sorted.max(), n_points)
    curve = np.vstack([s(t_fine) for s in splines]).T
    F_curve = idw_interpolation(curve, X, F_exp, power=power)
    return curve, F_curve

# Convert compositions to ternary by normalizing so Mg+Ca+Al = 100%
def to_ternary(arr):
    total = arr.sum(axis=1, keepdims=True)
    total[total == 0] = np.nan
    return arr / total

# ============================================================
# Data loading
# ============================================================

datasets = {
    "8": pd.read_csv(os.path.join(DATA_FOLDER, "data_Mer8.csv")),
    "15": pd.read_csv(os.path.join(DATA_FOLDER, "data_Mer15.csv")),
}
pressures = [1.5, 3.5, 5.0]

# ============================================================
# Styling and colorbar placement
# ============================================================

# Define color map for each (dataset, pressure) pair
cmap = {
    ("8", 1.5): "YlGnBu",
    ("8", 3.5): "BuGn",
    ("8", 5.0): "GnBu",
    ("15", 1.5): "YlOrBr",
    ("15", 3.5): "OrRd",
    ("15", 5.0): "RdPu",
}
marker_shape = {"8": "circle", "15": "diamond"}

# Position of colorbars on the plot
bar_positions = {
    ("8", 1.5): 0.75,
    ("8", 3.5): 0.80,
    ("8", 5.0): 0.85,
    ("15", 1.5): 0.90,
    ("15", 3.5): 0.95,
    ("15", 5.0): 1.00,
}
bar_thickness = 12
bar_length = 0.50

# ============================================================
# Plot construction
# ============================================================

fig = go.Figure()
annotations = []

# Add general title above the colorbars
annotations.append(dict(
    x=1.0, y=0.82, xref="paper", yref="paper",
    text="Taux de fusion partielle (%)", showarrow=False,
    font=dict(size=20, color="black"),
))

# Loop over datasets and pressures
for mer, df in datasets.items():
    for p in pressures:
        sub = df[df["Pression"] == p]
        if sub.empty:
            continue

        xyz = sub[["Mg/Si", "Ca/Si", "Al/Si"]].values
        Fvals = sub["F"].values
        tern = to_ternary(xyz) * 100  # Normalize to % for display

        # Add experimental points
        fig.add_trace(
            go.Scatterternary(
                a=tern[:, 0], b=tern[:, 1], c=tern[:, 2],
                mode="markers",
                marker=dict(
                    symbol=marker_shape[mer],
                    size=7,
                    opacity=0.8,
                    line=dict(width=0.5, color="grey"),
                    color=Fvals,
                    colorscale=cmap[(mer, p)],
                    cmin=0, cmax=100,
                    showscale=True,
                    colorbar=dict(
                        len=bar_length,
                        thickness=bar_thickness,
                        x=bar_positions[(mer, p)],
                        y=0.5,
                        tickmode="array",
                        tickvals=[0, 100],
                        ticktext=["0", "100"],
                        outlinewidth=0,
                        tickfont=dict(size=12),
                    ),
                ),
                name=f"{p} GPa (Mer {mer})",
                legendgroup=f"M{mer}P{p}",
            )
        )

        # Add label below each colorbar
        annotations.append(dict(
            x=bar_positions[(mer, p)], y=0.60,
            xref="paper", yref="paper",
            text=f"{p}\u00A0GPa (Mer{mer})", textangle=-90,
            showarrow=False,
            font=dict(size=13, color="dimgray"),
        ))

        # Compute and plot the smoothed trend line
        curve, Fcurve = compute_curve(sub)
        tern_curve = to_ternary(curve) * 100
        fig.add_trace(go.Scatterternary(
            a=tern_curve[:, 0], b=tern_curve[:, 1], c=tern_curve[:, 2],
            mode="lines",
            line=dict(width=2, color="rgba(80,80,80,0.35)"),
            showlegend=False,
            legendgroup=f"M{mer}P{p}",
        ))

        # Overlay F values on curve with small markers (no legend, no hover)
        fig.add_trace(go.Scatterternary(
            a=tern_curve[:, 0], b=tern_curve[:, 1], c=tern_curve[:, 2],
            mode="markers",
            marker=dict(size=4, color=Fcurve,
                        colorscale=cmap[(mer, p)], cmin=0, cmax=100,
                        showscale=False),
            hoverinfo="skip",
            showlegend=False,
        ))

# ============================================================
# Layout with black triangle and grid
# ============================================================

fig.update_layout(
    title="Ternaire des tendances de compositions expérimentales (Mer8 et Mer15)<br>en corrélation avec le taux de fusion partielle",
    ternary=dict(
        sum=100,
        aaxis=dict(title="Mg/Si (%)", linecolor="black", gridcolor="black"),
        baxis=dict(title="Ca/Si (%)", linecolor="black", gridcolor="black"),
        caxis=dict(title="Al/Si (%)", linecolor="black", gridcolor="black"),
        bgcolor="rgb(250,250,250)",
    ),
    legend=dict(x=0.01, y=0.98, bgcolor="rgba(255,255,255,0.7)"),
    margin=dict(l=30, r=90, t=110, b=30),
    annotations=annotations + [
        dict(
            text="<b>📌 Astuce :</b><br>Cliquez sur un élément de la légende<br>(ci-dessus) pour le masquer.<br>Sélectionnez pour zoomer.<br>Double-clic pour dézoomer.",
            x=0, y=0.3, showarrow=False,
            align='left',
            bordercolor='black', borderwidth=1,
            bgcolor='white', font=dict(size=14),
            xref='paper', yref='paper'
        )
    ]
)

fig.show()

# ============================================================
# Export
# ============================================================

OUTPUT_FOLDER = "interactive_diagrams"
os.makedirs(OUTPUT_FOLDER, exist_ok=True)
fig.write_html(os.path.join(OUTPUT_FOLDER, "FP_ternary_diagram.html"), include_plotlyjs='cdn')