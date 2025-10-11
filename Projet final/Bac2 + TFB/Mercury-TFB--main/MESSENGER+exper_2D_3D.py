"""
This script generates three interactive diagrams (1 in 3D and 2 in 2D)
to compare experimental melt trajectories (Mer8 and Mer15) 
with regional surface compositions on Mercury.

The plots include:
- 3D plot with Mg/Si, Ca/Si, Al/Si axes
- 2D diagrams: Ca/Si vs Mg/Si and Al/Si vs Mg/Si
- Experimental curves colored by pressure (1.5, 3.5, 5 GPa)
- Mer8 (●) and Mer15 (◆) markers shown in bold black
- Customizable transparency of regional cloud points
- Full legend and annotation support

The script uses Plotly for interactive rendering, and 
outputs HTML files for visual exploration.
"""

import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from lasagne import regions_array  # Custom function that returns a (3, 6, N) array for 6 Mercury regions
NORTH_ONLY = True     
# ─────────────────────────────────────────────────────────────
# Input/output directories
# ─────────────────────────────────────────────────────────────
DATA_FOLDER = "data"
OUTPUT_FOLDER = "interactive_diagrams"
os.makedirs(OUTPUT_FOLDER, exist_ok=True)

# ─────────────────────────────────────────────────────────────
# Utility functions
# ─────────────────────────────────────────────────────────────

def get_pixels(data):
    """
    Converts the (3, 6, N) array into a dictionary {region_id: (N, 3)},
    where each value is a matrix of valid (Mg/Si, Ca/Si, Al/Si) points per region.
    """
    pixels_by_region = {}
    for region_index in range(6):
        mg = data[0, region_index]
        ca = data[1, region_index]
        al = data[2, region_index]
        valid = ~np.isnan(mg) & ~np.isnan(ca) & ~np.isnan(al)
        pixels_by_region[region_index] = np.vstack([mg[valid], ca[valid], al[valid]]).T
    return pixels_by_region

def load_exp(filename):
    """
    Loads experimental CSV data and groups points by the 'Groupe' column.
    Each group is associated with a specific pressure value.
    """
    df = pd.read_csv(os.path.join(DATA_FOLDER, filename))
    df = df.rename(columns={"Mg/Si": "MgSi", "Ca/Si": "CaSi", "Al/Si": "AlSi"})
    return {
        g: (gdf[["MgSi", "CaSi", "AlSi"]].values, float(gdf["Pression"].iloc[0]))
        for g, gdf in df.groupby("Groupe")
    }

def add_exp_2d(fig, symbol, filename, press_colors, dim_x, dim_y):
    """
    Adds both line and marker plots to a 2D figure.
    Lines are color-coded by pressure, and markers are large black symbols.
    """
    for coords, press in load_exp(filename).values():
        colour = press_colors.get(press, "grey")
        fig.add_trace(
            go.Scatter(x=coords[:, dim_x], y=coords[:, dim_y], mode="lines",
                       line=dict(color=colour, width=5), showlegend=False))
        fig.add_trace(
            go.Scatter(x=coords[:, dim_x], y=coords[:, dim_y], mode="markers",
                       marker=dict(symbol=symbol, size=11, color="black"), showlegend=False))

def add_exp_3d(fig, symbol, filename, press_colors):
    """
    Adds 3D line + marker plots for experimental melt paths.
    Lines and markers are color-coded by pressure.
    """
    for coords, press in load_exp(filename).values():
        colour = press_colors.get(press, "grey")
        fig.add_trace(
            go.Scatter3d(x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                         mode="lines+markers",
                         line=dict(color=colour, width=5),
                         marker=dict(symbol=symbol, size=6, color=colour),
                         showlegend=False))

def add_dummy_marker(fig, symbol, name, is3d=False):
    """
    Adds a fake black marker to the legend.
    Used for Mer8 and Mer15 symbols.
    """
    if is3d:
        fig.add_trace(go.Scatter3d(x=[None], y=[None], z=[None], mode="markers",
                                   marker=dict(symbol=symbol, size=8, color="black"),
                                   name=name))
    else:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="markers",
                                 marker=dict(symbol=symbol, size=11, color="black"),
                                 name=name))

def add_dummy_line(fig, colour, name, is3d=False):
    """
    Adds a fake line in a given color for the legend.
    Used to represent pressure values.
    """
    if is3d:
        fig.add_trace(go.Scatter3d(x=[None], y=[None], z=[None], mode="lines",
                                   line=dict(color=colour, width=3), name=name))
    else:
        fig.add_trace(go.Scatter(x=[None], y=[None], mode="lines",
                                 line=dict(color=colour, width=2), name=name))

# ─────────────────────────────────────────────────────────────
# Global plot parameters
# ─────────────────────────────────────────────────────────────
region_opacity = 0.6
pressure_colors = {
    1.5: "rgb(255,165,0)",  # Orange
    3.5: "rgb(128,0,128)",  # Violet
    5.0: "rgb(0,191,255)"   # Dark turquoise
}
dshapes = {"Mer8": ("circle", "sphère"), "Mer15": ("diamond", "losange")}
region_colors = [
    "rgb(255,0,0)",      # High-Mg
    "rgb(0,255,0)",      # Al-rich
    "rgb(0,0,255)",      # Caloris
    "rgb(255,255,0)",    # Rach
    "rgb(255,0,255)",    # High-al NVP
    "rgb(0,255,255)"     # Low-al NVP
]
region_names = ["High‑Mg", "Al‑rich", "Caloris", "Rach", "High‑Al NVP", "Low‑Al NVP"]
region_pixels = get_pixels(regions_array())

# ─────────────────────────────────────────────────────────────
# === 3D FIGURE ===
# ─────────────────────────────────────────────────────────────
fig3d = go.Figure()

for dset, csv in [("Mer8", "data_Mer8.csv"), ("Mer15", "data_Mer15.csv")]:
    add_exp_3d(fig3d, dshapes[dset][0], csv, pressure_colors)

for i, pts in region_pixels.items():
    fig3d.add_trace(go.Scatter3d(x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
                                 mode="markers", name=region_names[i],
                                 marker=dict(size=3, color=region_colors[i], opacity=region_opacity)))

add_dummy_marker(fig3d, "circle", "Mer8 (sphère)", is3d=True)
add_dummy_marker(fig3d, "diamond", "Mer15 (losange)", is3d=True)
for p, c in pressure_colors.items():
    add_dummy_line(fig3d, c, f"{p} GPa", is3d=True)

fig3d.update_layout(
    title="Évolution de la composition pour les différentes régions de Mercure, en superposition des chemins de fusion partielle expérimentale sous différentes pressions.",
    title_font=dict(size=22),
    scene=dict(
        xaxis=dict(title="Mg/Si", title_font=dict(size=18), tickfont=dict(size=14)),
        yaxis=dict(title="Ca/Si", title_font=dict(size=18), tickfont=dict(size=14)),
        zaxis=dict(title="Al/Si", title_font=dict(size=18), tickfont=dict(size=14)),
        bgcolor="rgb(250,250,250)"
    ),
    legend=dict(font=dict(size=20)),
    margin=dict(l=30, r=30, t=80, b=30),
    annotations=[
        dict(
            text="<b>📌 Astuce :</b><br>Cliquez sur un élément de la légende pour le masquer.<br>Sélectionnez pour zoomer.<br>Double-clic pour dézoomer.",
            x=1, y=1, showarrow=False, align='left', bordercolor='black',
            borderwidth=1, bgcolor='white', font=dict(size=14),
            xref='paper', yref='paper'
        )
    ]
)
fig3d.write_html(os.path.join(OUTPUT_FOLDER, "cloud_3D_diagram.html"), include_plotlyjs="cdn")

# ─────────────────────────────────────────────────────────────
# === 2D FIGURE: Ca/Si vs Mg/Si ===
# ─────────────────────────────────────────────────────────────
fig_ca = go.Figure()
for i, pts in region_pixels.items():
    fig_ca.add_trace(go.Scatter(x=pts[:, 0], y=pts[:, 1], mode="markers",
                                name=region_names[i],
                                marker=dict(size=5, color=region_colors[i], opacity=region_opacity)))
for dset, csv in [("Mer8", "data_Mer8.csv"), ("Mer15", "data_Mer15.csv")]:
    add_exp_2d(fig_ca, dshapes[dset][0], csv, pressure_colors, 0, 1)

add_dummy_marker(fig_ca, "circle", "Mer8 (sphère)")
add_dummy_marker(fig_ca, "diamond", "Mer15 (losange)")
for p, c in pressure_colors.items():
    add_dummy_line(fig_ca, c, f"{p} GPa")

fig_ca.update_layout(
    title="Évolution du rapport Ca/Si en fonction de Mg/Si pour les différentes régions de Mercure,<br>en superposition des chemins de fusion partielle expérimentale sous différentes pressions.",
    title_font=dict(size=22),
    xaxis=dict(title="Mg/Si", title_font=dict(size=18), tickfont=dict(size=14)),
    yaxis=dict(title="Ca/Si", title_font=dict(size=18), tickfont=dict(size=14)),
    legend=dict(font=dict(size=20)),
    annotations=[
        dict(
            text="<b>📌 Astuce :</b><br>Cliquez sur un élément de la légende pour le masquer.<br>Sélectionnez pour zoomer.<br>Double-clic pour dézoomer.",
            x=2, y=1, showarrow=False, align='left', bordercolor='black',
            borderwidth=1, bgcolor='white', font=dict(size=14),
            xref='paper', yref='paper'
        )
    ]
)
fig_ca.write_html(os.path.join(OUTPUT_FOLDER, "cloud_Ca_diagram.html"), include_plotlyjs="cdn")

# ─────────────────────────────────────────────────────────────
# === 2D FIGURE: Al/Si vs Mg/Si ===
# ─────────────────────────────────────────────────────────────
fig_al = go.Figure()
for i, pts in region_pixels.items():
    fig_al.add_trace(go.Scatter(x=pts[:, 0], y=pts[:, 2], mode="markers",
                                name=region_names[i],
                                marker=dict(size=5, color=region_colors[i], opacity=region_opacity)))
for dset, csv in [("Mer8", "data_Mer8.csv"), ("Mer15", "data_Mer15.csv")]:
    add_exp_2d(fig_al, dshapes[dset][0], csv, pressure_colors, 0, 2)

add_dummy_marker(fig_al, "circle", "Mer8 (sphère)")
add_dummy_marker(fig_al, "diamond", "Mer15 (losange)")
for p, c in pressure_colors.items():
    add_dummy_line(fig_al, c, f"{p} GPa")

fig_al.update_layout(
    title="Évolution du rapport Al/Si en fonction de Mg/Si pour les différentes régions de Mercure,<br>en superposition des chemins de fusion partielle expérimentale sous différentes pressions.",
    title_font=dict(size=22),
    xaxis=dict(title="Mg/Si", title_font=dict(size=18), tickfont=dict(size=14)),
    yaxis=dict(title="Al/Si", title_font=dict(size=18), tickfont=dict(size=14)),
    legend=dict(font=dict(size=20)),
    annotations=[
        dict(
            text="<b>📌 Astuce :</b><br>Cliquez sur un élément de la légende pour le masquer.<br>Sélectionnez pour zoomer.<br>Double-clic pour dézoomer.",
            x=2, y=1, showarrow=False, align='left', bordercolor='black',
            borderwidth=1, bgcolor='white', font=dict(size=14),
            xref='paper', yref='paper'
        )
    ]
)
fig_al.write_html(os.path.join(OUTPUT_FOLDER, "cloud_Al_diagram.html"), include_plotlyjs="cdn")

print("✔️  Figures interactives générées dans", OUTPUT_FOLDER)