
"""
This script performs cubic interpolation of temperature (Temp) as a function 
of degree of melting (F, in %) and pressure (in GPa) for two datasets: Mer8 and Mer15.

It generates a side-by-side (1×2) figure using a custom colormap that goes from black → red → yellow.
Each subplot shows a smooth interpolated temperature map for one mission.
"""

from pathlib import Path
from typing import Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from matplotlib.colors import LinearSegmentedColormap

# Set global font and layout styles for matplotlib
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 20,
})

# ─────────────────────────────────────────────────────────────
# PARAMETERS
# ─────────────────────────────────────────────────────────────

# CSV file paths for each mission
CSV_FILES: Dict[str, str] = {
    "8": "data/data_Mer8.csv",
    "15": "data/data_Mer15.csv",
}

INTERP_METHOD: str = "cubic"  # interpolation method
CMAP = LinearSegmentedColormap.from_list("black_red_yellow", ["black", "red", "yellow"])  # custom color map
FIGSIZE: tuple = (12, 5)  # figure size
GRID_RES_F: int = 250  # number of points for F axis
GRID_RES_P: int = 250  # number of points for Pressure axis

# ─────────────────────────────────────────────────────────────
# LOAD CSV DATA
# ─────────────────────────────────────────────────────────────

required_cols = {"Pression", "F", "Temp"}
data: Dict[str, pd.DataFrame] = {}

# Check that all required columns exist in each CSV
for mission, file in CSV_FILES.items():
    path = Path(file)
    if not path.is_file():
        raise FileNotFoundError(f"Fichier introuvable : {file}")
    df = pd.read_csv(path)
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Colonnes manquantes dans « {file} » : {', '.join(missing)}")
    data[mission] = df

# ─────────────────────────────────────────────────────────────
# CREATE COMMON INTERPOLATION GRID
# ─────────────────────────────────────────────────────────────

# Gather all F and Pressure values from both datasets
all_F = np.concatenate([d["F"].values for d in data.values()])
all_P = np.concatenate([d["Pression"].values for d in data.values()])

# Generate 1D grid vectors for F and Pressure axes
F_grid = np.linspace(all_F.min(), 100.0, GRID_RES_F)
P_grid = np.linspace(all_P.min(), all_P.max(), GRID_RES_P)

# Generate 2D meshgrid (F, P) for interpolation
F_mesh, P_mesh = np.meshgrid(F_grid, P_grid)

# ─────────────────────────────────────────────────────────────
# PLOT SETUP: 1 × 2 FIGURE
# ─────────────────────────────────────────────────────────────

# Create subplots for Mer8 and Mer15
fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, sharex=True, sharey=True, constrained_layout=True)

# These will track global colorbar range
vmin_global, vmax_global = None, None

# ─────────────────────────────────────────────────────────────
# INTERPOLATION AND RANGE TRACKING
# ─────────────────────────────────────────────────────────────

Ti_all = {}  # store interpolated results per mission

for mission in ["8", "15"]:
    df = data[mission]
    points = df[["F", "Pression"]].values
    values = df["Temp"].values

    # Perform cubic interpolation over the 2D grid
    Ti = griddata(points, values, (F_mesh, P_mesh), method=INTERP_METHOD)
    Ti_all[mission] = Ti

    # Update global min/max for color scale
    vmin = np.nanmin(Ti)
    vmax = np.nanmax(Ti)
    vmin_global = vmin if vmin_global is None else min(vmin_global, vmin)
    vmax_global = vmax if vmax_global is None else max(vmax_global, vmax)

# ─────────────────────────────────────────────────────────────
# DRAW EACH INTERPOLATED TEMPERATURE MAP
# ─────────────────────────────────────────────────────────────

for i, mission in enumerate(["8", "15"]):
    ax = axes[i]
    Ti = Ti_all[mission]
    im = ax.imshow(
        Ti,
        origin="lower",
        extent=[F_grid.min(), F_grid.max(), P_grid.min(), P_grid.max()],
        aspect="auto",
        cmap=CMAP,
        vmin=vmin_global,
        vmax=vmax_global,
    )

    # Invert y-axis to display increasing pressure downward
    ax.set_ylim(P_grid.max(), P_grid.min())
    ax.set_xlabel("F (%)")
    if i == 0:
        ax.set_ylabel("Pression (GPa)")
    label = "(a)" if mission == "8" else "(b)"
    ax.set_title(f"{label}")

# ─────────────────────────────────────────────────────────────
# SHARED COLORBAR BELOW THE FIGURE
# ─────────────────────────────────────────────────────────────

cbar_ax = fig.add_axes([0.15, -0.1, 0.7, 0.05])  # [left, bottom, width, height]
cbar = fig.colorbar(im, cax=cbar_ax, orientation="horizontal")
cbar.set_label("Température interpolée (°C)")


plt.show()