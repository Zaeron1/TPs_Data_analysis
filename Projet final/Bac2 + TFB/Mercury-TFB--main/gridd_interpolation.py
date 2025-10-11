"""
Comparison of experimental datasets Mer8 and Mer15:
– Performs a 2D interpolation (griddata) of geochemical ratios (Mg/Si, Ca/Si, Al/Si) 
  across a grid defined by pressure (P) and melt fraction (F)
– For each global pixel, finds the P–F combination with the smallest residual 
  (difference from observed planetary data)
– Generates interpolated maps and optimal P and F maps for both datasets
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from tqdm import tqdm
import lasagne  # Local module for accessing global planetary data
NORTH_ONLY = True     
from scipy.interpolate import griddata
import os

DATA_FOLDER = "data"

# —————————————————————————  Réglages esthétiques globaux  ————————————————————————— #
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 20,
})
deg_colors = LinearSegmentedColormap.from_list("black_red_yellow",
                                               ["black", "red", "yellow"])

# —————————————————————————————  Fonctions utilitaires  ————————————————————————————— #
def process_griddata(df, n=50):
    """
    Performs 2D interpolation (griddata) of Mg/Si, Ca/Si, and Al/Si on a P–F grid.
    """
    p, f = df["Pression"].values, df["F"].values
    comps = {c: df[c].values for c in ["Mg/Si", "Ca/Si", "Al/Si"]}

    Pm, Fm = np.meshgrid(np.linspace(p.min(), p.max(), n),
                         np.linspace(f.min(), 100.0, n))
    mesh = {c: griddata((p, f), v, (Pm, Fm), method="cubic") for c, v in comps.items()}

    df_int = pd.DataFrame(
        {"Pression": Pm.ravel(), "F": Fm.ravel(),
         **{c: mesh[c].ravel() for c in mesh}}
    ).dropna().reset_index(drop=True)
    return Pm, Fm, mesh, df_int


def compute_residuals(df_int):
    """
    For each pixel in the global dataset, finds the P–F combination from the 
    interpolated experimental dataset with the lowest relative residual.
    """
    g = lasagne.regions_array()
    mg, ca, al = g[0, 6].astype(float), g[1, 6].astype(float), g[2, 6].astype(float)
    h, w = mg.shape
    best_res = np.full((h, w), np.nan)
    best_idx = np.full((h, w), -1)
    vals = {c: df_int[c].values for c in ["Mg/Si", "Ca/Si", "Al/Si"]}

    for i in tqdm(range(h), desc="Calculating residuals", unit="ligne"):
        for j in range(w):
            obs = [mg[i, j], ca[i, j], al[i, j]]
            if any(np.isnan(obs)):
                continue
            with np.errstate(divide="ignore", invalid="ignore"):
                res = sum(np.abs(vals[c] - obs[k]) / obs[k]
                          for k, c in enumerate(["Mg/Si", "Ca/Si", "Al/Si"]))
            idx = np.nanargmin(res)
            best_res[i, j] = res[idx]
            best_idx[i, j] = idx

    Pmap = np.where(best_idx != -1, df_int["Pression"].values[best_idx], np.nan)
    Fmap = np.where(best_idx != -1, df_int["F"].values[best_idx], np.nan)
    return best_res, best_idx, Pmap, Fmap

# ————————————————————————————————  Pipeline principal  ——————————————————————————————— #
def compare():
    files = {"8": "data_Mer8.csv", "15": "data_Mer15.csv"}
    results = {}

    # 1) Pré-traitement & interpolation
    for k, f in files.items():
        print(f"\n── Processing Mer{k} ──")
        path = os.path.join(DATA_FOLDER, f)
        df = pd.read_csv(path)
        Pm, Fm, mesh, dfi = process_griddata(df)
        res, idx, Pmap, Fmap = compute_residuals(dfi)
        results[k] = {"P_mesh": Pm, "F_mesh": Fm, "mesh": mesh, "df_int": dfi,
                      "best_res": res, "best_idx": idx, "Pmap": Pmap, "Fmap": Fmap}

    # 2) Figure 1 : interpolation Mg/Si – Ca/Si – Al/Si
    comps = ["Mg/Si", "Ca/Si", "Al/Si"]
    fig1, axs1 = plt.subplots(2, 3, figsize=(15, 10), sharex="col", sharey="row")

    for r, (k, res) in enumerate(results.items()):
        for c, comp in enumerate(comps):
            ax = axs1[r, c]
            im = ax.imshow(
                res["mesh"][comp].T,
                extent=[res["F_mesh"].min(), res["F_mesh"].max(),
                        res["P_mesh"].max(), res["P_mesh"].min()],
                origin="upper", aspect="auto", cmap=deg_colors
            )
            ax.set_title(f"Mer{k} – {comp}")
            if r == 1:
                ax.set_xlabel("F (%)")
            if c == 0:
                ax.set_ylabel("Pression (GPa)")

    fig1.tight_layout(rect=[0, 0.10, 1, 0.95])
    cbar = fig1.colorbar(im, ax=axs1.flatten(), orientation="horizontal",
                         fraction=0.05, pad=0.12)
    cbar.set_label("Valeur interpolée")
    plt.show()

    # 3) Figure 2 : cartes de pression optimale
    fig2, axs2 = plt.subplots(2, 1, figsize=(9, 8))  # élargi verticalement
    vminP = np.nanmin([results[k]["Pmap"] for k in results])
    vmaxP = np.nanmax([results[k]["Pmap"] for k in results])

    lon_ticks = np.linspace(0, results["8"]["Pmap"].shape[1] - 1, 9)
    lat_ticks = np.linspace(0, results["8"]["Pmap"].shape[0] - 1, 3)
    lon_labels = [f"{x:.0f}°" for x in np.linspace(-180, 180, len(lon_ticks))]
    lat_labels = [f"{y:.0f}°" for y in np.linspace(90, 0, len(lat_ticks))]

    for ax, (k, res) in zip(axs2, results.items()):
        imP = ax.imshow(res["Pmap"], origin="upper", vmin=vminP, vmax=vmaxP,
                        cmap="viridis")
        ax.set_title(f"Mer{k} – Pression optimale")
        ax.set_xticks(lon_ticks); ax.set_xticklabels(lon_labels)
        ax.set_yticks(lat_ticks); ax.set_yticklabels(lat_labels)

    fig2.tight_layout(rect=[0, 0.08, 1, 0.97])
    cbar = fig2.colorbar(imP, ax=axs2, orientation="horizontal",
                         fraction=0.05, pad=0.08)
    cbar.set_label("Pression (GPa)")
    plt.show()

    # 4) Figure 3 : cartes de fraction fondue optimale
    fig3, axs3 = plt.subplots(2, 1, figsize=(9, 8))  # élargi verticalement
    vminF = np.nanmin([results[k]["Fmap"] for k in results])
    vmaxF = np.nanmax([results[k]["Fmap"] for k in results])

    for ax, (k, res) in zip(axs3, results.items()):
        imF = ax.imshow(res["Fmap"], origin="upper", vmin=vminF, vmax=vmaxF,
                        cmap="viridis")
        ax.set_title(f"Mer{k} – Fraction fondue optimale")
        ax.set_xticks(lon_ticks); ax.set_xticklabels(lon_labels)
        ax.set_yticks(lat_ticks); ax.set_yticklabels(lat_labels)

    fig3.tight_layout(rect=[0, 0.08, 1, 0.97])
    cbar = fig3.colorbar(imF, ax=axs3, orientation="horizontal",
                         fraction=0.05, pad=0.08)
    cbar.set_label("F (%)")
    plt.show()

if __name__ == "__main__":
    compare()