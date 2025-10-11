"""
Mercury Geochemical Analysis
======================================================
Author: Alexandre Michaux (University de Liège) – 2025

This single script reproduces the complete processing chain used in my TFB to
compare two experimental melt series (Mer8 and Mer15) with Mercury’s surface
composition maps.

Pipeline overview
-----------------
1. Load experimental CSV files (`data_Mer8.csv`, `data_Mer15.csv`).  Each row
   contains a pressure P (GPa), melt fraction F*(%), and the three elemental
   ratios Mg/Si, Ca/Si, Al/Si.
2. Polynomial interpolation (`numpy.polynomial.Polynomial`) of these three
   ratios as a function of *F* for every pressure level (default degree = 1).
3. Global pixel matching – for every pixel of Mercury’s global maps
   (delivered by my local module `lasagne.regions_array()`):
   - read observed Mg, Ca, Al values;
   - compute the absolute relative residual to every interpolated experimental
     record;
   - pick the (P, F) pair with the lowest sum of residuals;
   - store the optimal pressure map, melt‑fraction map and total residual.
4. Visualisation  – all figure titles, labels and legends are in French
   as required by the thesis guidelines:
   - scatter & line plots to check 1‑D interpolations;
   - 2‑D maps of optimal pressure and melt fraction (colour‑bar shared);
   - histograms and bar charts summarising residuals and region counts;
   - regional zoom comparing observed vs. modelled ratios with R² tables.
   - interactive 3D maps of optimal melt fraction and pressure
"""


import os
from typing import Dict, List

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.polynomial import Polynomial
import plotly.graph_objects as go
from tqdm import tqdm

import lasagne  # Local module for accessing global planetary data
NORTH_ONLY = True          # Show only the northern hemisphere in maps

# ---------------------------------------------------------------------------
# Global Matplotlib style: slightly larger fonts for thesis figures
# ---------------------------------------------------------------------------
plt.rcParams.update({
    "font.size": 14,      # base font size
    "axes.titlesize": 18, # subplot titles
    "axes.labelsize": 16, # axis labels
    "legend.fontsize": 20,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 20,  # top‑level figure title
})

DATA_FOLDER = "data"  # folder containing experimental CSV files

# ===========================================================================
# 1. PRE‑PROCESSING – interpolate experimental datasets
# ===========================================================================

def process_dataset(df: pd.DataFrame, deg: int = 1, label: str = "") -> Dict[str, np.ndarray]:
    """Return interpolated experimental values and best‑fit maps for *one*
    experimental dataset (Mer8 or Mer15).

    Parameters
    ----------
    df  : Raw experimental dataframe (one of the two CSVs).
    deg : Degree of the Polynomial fit (default = 1 ⇒ linear).

    Returns
    -------
    dict with the following keys (all NumPy arrays):
    - ``df_int``: the interpolated dataframe (long‑form)
    - ``best_res``: 2‑D array of minimal residuals per pixel
    - ``best_idx``: index of the winning experimental record per pixel
    - ``pressure_map`` / ``F_map``: best‑fit pressure & melt fraction maps
    - ``pressures``: sorted array of unique pressure levels in the CSV
    """

    # ---- 1‑D polynomial interpolation over F for each pressure -------------
    pressures = np.sort(df["Pression"].unique())
    records: List[Dict[str, float]] = []

    for p in pressures:
        sub = df[df["Pression"] == p].sort_values("F")
        if sub.empty:
            continue  # safeguard – should not happen
        Fv = np.linspace(sub["F"].min(), max(100, sub["F"].max()), 50)
        polys = {
            r: Polynomial.fit(sub["F"], sub[r], deg)
            for r in ("Mg/Si", "Ca/Si", "Al/Si")
        }
        records.extend(
            {"Pression": p, "F": fv, **{r: polys[r](fv) for r in polys}} for fv in Fv
        )

    df_int = pd.DataFrame(records)

    # ---- Load regional maps from *lasagne* ---------------------------------
    g = lasagne.regions_array()  # shape = (3, 6, Ny, Nx)
    maps = {r: np.asarray(g[i, 6], float) for i, r in enumerate(("Mg", "Ca", "Al"))}
    Ny, Nx = maps["Mg"].shape

    # ---- Pixel‑wise residual search ----------------------------------------
    exp = {r: df_int[r].values for r in ("Mg/Si", "Ca/Si", "Al/Si")}
    best_res = np.full((Ny, Nx), np.nan)
    best_idx = np.full((Ny, Nx), -1, int)

    for i in tqdm(range(Ny), desc=f"Processing Mer{label}...", unit=" rows"):
        for j in range(Nx):
            obs = {r: maps[r][i, j] for r in maps}
            if any(np.isnan(v) for v in obs.values()):
                continue  # invalid pixel in global map
            with np.errstate(divide="ignore", invalid="ignore"):
                res = (
                    np.abs(exp["Mg/Si"] - obs["Mg"]) / obs["Mg"]
                    + np.abs(exp["Ca/Si"] - obs["Ca"]) / obs["Ca"]
                    + np.abs(exp["Al/Si"] - obs["Al"]) / obs["Al"]
                )
            idx = int(np.nanargmin(res))  # index in interpolated dataframe
            best_res[i, j] = res[idx]
            best_idx[i, j] = idx

    pressure_map = np.where(best_idx >= 0, df_int["Pression"].values[best_idx], np.nan)
    F_map = np.where(best_idx >= 0, df_int["F"].values[best_idx], np.nan)

    return {
        "df_int": df_int,
        "best_res": best_res,
        "best_idx": best_idx,
        "pressure_map": pressure_map,
        "F_map": F_map,
        "pressures": pressures,
    }

# ===========================================================================
# 2. VISUALISATIONS – helper routines 
# ===========================================================================

def plot_interp(results: Dict[str, dict], data: Dict[str, pd.DataFrame], ratios: List[str]):
    """Scatter (raw) + line (interpolated) plots for each ratio and mission."""
    for ratio in ratios:
        fig, ax = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(6, 8))
        for a, k in zip(ax, results):
            raw, df_i = data[k], results[k]["df_int"]
            for p in results[k]["pressures"]:
                sr = raw[raw["Pression"] == p]
                si = df_i[df_i["Pression"] == p]
                a.scatter(sr["F"], sr[ratio], s=10)
                a.plot(si["F"], si[ratio], lw=1)
            a.set_title(f"Mer{k}")
        ax[1].set_xlabel("F (%)")
        ax[0].set_ylabel(ratio)
        ax[1].set_ylabel(ratio)
        plt.tight_layout()
        plt.show()


def plot_map(results: Dict[str, dict], key: str, title: str, cmap: str = "viridis", unit: str = ""):
    """Display two vertically‑stacked maps (Mer8 & Mer15) for *key*."""
    fig, ax = plt.subplots(2, 1, figsize=(9, 8))  # slightly taller for clarity

    # Unified colour scale
    if key == "F_map":
        vmin, vmax = 0, 60
    else:
        vmin = min(np.nanmin(results["8"][key]), np.nanmin(results["15"][key]))
        vmax = max(np.nanmax(results["8"][key]), np.nanmax(results["15"][key]))

    for a, k in zip(ax, results):
        im = a.imshow(results[k][key], origin="upper", vmin=vmin, vmax=vmax, cmap=cmap)
        # Latitude/longitude ticks
        lon_ticks = np.linspace(0, results[k][key].shape[1] - 1, 9)
        lat_ticks = np.linspace(0, results[k][key].shape[0] - 1, 3)
        lon_labels = [f"{x:.0f}°" for x in np.linspace(-180, 180, len(lon_ticks))]
        lat_labels = [f"{y:.0f}°" for y in np.linspace(90, 0, len(lat_ticks))]
        a.set_xticks(lon_ticks); a.set_xticklabels(lon_labels)
        a.set_yticks(lat_ticks); a.set_yticklabels(lat_labels)
        a.set_title(f"Mer{k} – {title}")

    # Shared colour‑bar
    fig.subplots_adjust(top=0.92, bottom=0.07, hspace=0.15)
    cax = fig.add_axes([0.2, 0.02, 0.6, 0.02])
    plt.colorbar(im, cax=cax, orientation="horizontal").set_label(f"{title} {unit}")
    plt.show()


def plot_hist_res(results):
    plt.figure()
    for k, style in zip(["8", "15"], [{"alpha": 0.5}, {"histtype": "step", "linestyle": "--"}]):
        res = results[k]["best_res"].ravel()
        plt.hist(res[~np.isnan(res)], bins=30, label=f"Mer{k}", **style)
    plt.xlabel("Somme des résidus relatifs")
    plt.ylabel("Nombre de pixels")
    plt.title("")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_frac_table(results):
    thr = [0.20, 0.30, 0.50]
    fracs = {
        k: [
            100
            * np.sum(results[k]["best_res"][~np.isnan(results[k]["best_res"])] <= t)
            / np.sum(~np.isnan(results[k]["best_res"]))
            for t in thr
        ]
        for k in ["8", "15"]
    }
    df = pd.DataFrame({"Threshold": thr, "Mer8 (%)": fracs["8"], "Mer15 (%)": fracs["15"]})
    print("\nTable des fractions de surface sous seuils (%):\n", df.to_string(index=False))
    x = np.arange(len(thr))
    plt.bar(x - 0.2, fracs["8"], 0.4, label="Mer8")
    plt.bar(x + 0.2, fracs["15"], 0.4, label="Mer15", hatch="///")
    plt.xticks(x, thr)
    plt.xlabel("Seuil de résidu")
    plt.ylabel("% surface")
    plt.title("")
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_region_counts(results):
    """
    Histogramme RELATIF des pixels par région et par pression (Mer8/Mer15).

    – Mer8 : barres pleines semi-transparentes ;
    – Mer15 : contours pointillés.
    – Chaque barre est normalisée par le nombre total de pixels valides de la région.
    """

    names = ["High-Mg", "Al-rich", "Caloris", "Rach", "High-al NVP"]
    ps = np.unique(np.concatenate([results["8"]["pressures"], results["15"]["pressures"]]))
    ga = lasagne.regions_array()
    cmap = plt.cm.get_cmap("tab20", len(names))
    fig, ax = plt.subplots(figsize=(9, 6))
    bw = 0.8 / len(names)

    for i, reg in enumerate(names):
        mask = ~np.isnan(np.asarray(ga[0, i], float))
        total_pixels = np.sum(mask)
        if total_pixels == 0:
            continue

        for k in ["8", "15"]:
            idxs = [
                results[k]["df_int"].at[idx, "Pression"]
                for (ii, jj), idx in np.ndenumerate(results[k]["best_idx"])
                if idx >= 0 and mask[ii, jj]
            ]
            counts = np.bincount([np.where(ps == p)[0][0] for p in idxs], minlength=len(ps))
            # normalisation relative
            counts = counts / total_pixels

            x = np.arange(len(ps)) + (i - len(names) / 2 + 0.5) * bw
            ax.bar(
                x,
                counts,
                bw,
                label=f"{reg} (Mer{k})",
                facecolor=cmap(i) if k == "8" else "none",
                edgecolor=cmap(i),
                alpha=0.6 if k == "8" else 1.0,
                linestyle="--",
                linewidth=1.2,
            )

    ax.set_xticks(np.arange(len(ps)))
    ax.set_xticklabels([f"{p} GPa" for p in ps])
    ax.set_ylabel("Fraction relative des pixels")
    ax.set_title("")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)
    plt.tight_layout()
    plt.show()


def plot_Fopt(results):
    bins = np.linspace(0, 100, 21)
    F8 = results["8"]["F_map"].ravel()[~np.isnan(results["8"]["F_map"].ravel())]
    F15 = results["15"]["F_map"].ravel()[~np.isnan(results["15"]["F_map"].ravel())]
    df = pd.DataFrame(
        {
            "F_bin_start": bins[:-1],
            "Mer8_count": np.histogram(F8, bins)[0],
            "Mer15_count": np.histogram(F15, bins)[0],
        }
    )
    print("\nTable des nombres de pixels par bin de F optimal :\n", df.to_string(index=False))
    plt.hist(F8, bins, alpha=0.7, label="Mer8")
    plt.hist(F15, bins, histtype="step", linestyle="--", label="Mer15")
    plt.xlabel("F optimal (%)")
    plt.ylabel("Nombre de pixels")
    plt.title("")
    plt.legend()
    plt.tight_layout()
    plt.show()

def analyze_regions(results: Dict[str, dict]):
    """Compare les rapports élémentaires observés vs modélisés pour plusieurs régions.

    – Trois sous‑graphes côte à côte (Mg/Si, Ca/Si, Al/Si), chacun avec ses propres limites.
    – Axes recadrés sur l’intervalle pertinent (plus de détour par l’origine 0, 0).
    – Diagonale pointillée *x = y* traversant toute la zone, identifiée dans la légende.
    """

    # Préparation des données --------------------------------------------------
    ga = lasagne.regions_array()
    region_names = [
        "High-Mg",
        "Al-rich",
        "Caloris",
        "Rach",
        "High-al NVP",
        "Low-al NVP",
    ]
    ratios = {"Mg/Si": 0, "Ca/Si": 1, "Al/Si": 2}
    colors = plt.cm.get_cmap("tab10", len(region_names))
    markers = {"8": "*", "15": "D"}
    marker_labels = {"8": "Mer8", "15": "Mer15"}

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    rows = []  # pour la table des R²
    region_handles = {}
    mission_handles = {}

    # Parcourir chaque ratio et sous‑graphe ------------------------------------
    for ax, (ratio, idx_r) in zip(axes, ratios.items()):
        obs_all, mod_all = [], []  # pour déterminer les limites d’axe
        for mission in markers:
            obs_m, mod_m = [], []
            for i, reg in enumerate(region_names):
                mask = ~np.isnan(np.asarray(ga[idx_r, i], float))
                obs = np.asarray(ga[idx_r, i], float)[mask]
                idxs = results[mission]["best_idx"][mask].astype(int)
                valid = idxs >= 0
                obs, idxs = obs[valid], idxs[valid]
                if obs.size == 0:
                    obs_m.append(np.nan)
                    mod_m.append(np.nan)
                    continue

                mod = results[mission]["df_int"][ratio].values[idxs]
                obs_mean = float(np.nanmean(obs))
                mod_mean = float(np.nanmean(mod))
                obs_m.append(obs_mean)
                mod_m.append(mod_mean)

                sc = ax.scatter(
                    mod_mean,
                    obs_mean,
                    marker=markers[mission],
                    s=110,
                    color=colors(i),
                    label=f"{reg}" if mission == "8" else "_"
                )
                obs_all.append(obs_mean)
                mod_all.append(mod_mean)

                # Handle légende couleur une seule fois
                if reg not in region_handles:
                    patch = plt.Line2D([0], [0], marker='s', color='none',
                                       markerfacecolor=colors(i), markersize=12,
                                       label=reg)
                    region_handles[reg] = patch

            # Handle légende mission une seule fois
            if mission not in mission_handles:
                handle = plt.Line2D([0], [0], marker=markers[mission], color='black',
                                    linestyle='None', markersize=10, label=marker_labels[mission])
                mission_handles[mission] = handle

            # Coefficient de corrélation pour le ratio/mission -----------------
            mask_valid = ~np.isnan(obs_m) & ~np.isnan(mod_m)
            if np.sum(mask_valid) >= 2:
                r = np.corrcoef(np.array(obs_m)[mask_valid], np.array(mod_m)[mask_valid])[0, 1]
                rows.append({"Ratio": ratio, "Mission": f"Mer{mission}", "R²": r ** 2})
            else:
                rows.append({"Ratio": ratio, "Mission": f"Mer{mission}", "R²": np.nan})

        # Ajustements esthétiques du sous‑graphe ------------------------------
        all_vals = np.array(obs_all + mod_all)
        vmin, vmax = np.nanmin(all_vals), np.nanmax(all_vals)
        pad = 0.05 * (vmax - vmin) if vmax > vmin else 1.0
        ax_lim_min = vmin - pad
        ax_lim_max = vmax + pad
        ax.plot([ax_lim_min, ax_lim_max], [ax_lim_min, ax_lim_max], "k--", linewidth=0.8,
                label=r"$x=y$" if ratio == "Mg/Si" else "_")
        ax.set_xlim(ax_lim_min, ax_lim_max)
        ax.set_ylim(ax_lim_min, ax_lim_max)
        ax.set_xlabel("Modélisé")
        ax.set_ylabel("Observé")
        ax.set_title(ratio)

    # Légende unique ----------------------------------------------------------
    all_handles = list(region_handles.values()) + list(mission_handles.values())
    all_labels = [h.get_label() for h in all_handles]
    fig.legend(all_handles, all_labels, loc="lower center", ncol=6, frameon=False, bbox_to_anchor=(0.5, -0.1))

    fig.suptitle("", y=1.02)
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    plt.show()

    # Affichage du tableau récapitulatif des R² -------------------------------
    df_r2 = pd.DataFrame(rows).pivot(index="Ratio", columns="Mission", values="R²").sort_index()
    print("\nTable des R² par ratio :\n", df_r2.to_string())
    


def plot_globe_map(results: Dict[str, dict], key: str, title: str, unit: str = "", colormap: str = "Viridis"):
    os.makedirs("interactive_diagrams", exist_ok=True)

    # Pour F_map, on limite la colorbar
    is_fusion = (key == "F_map")
    cmin_val = 0 if is_fusion else None
    cmax_val = 60 if is_fusion else None

    for mission in ["8", "15"]:
        data = results[mission][key]
        Ny, Nx = data.shape

        # grilles lat/lon (nord + sud fictif)
        lat_n = np.linspace(90, 0, Ny)
        lat_s = np.linspace(0, -90, Ny)
        lat = np.concatenate([lat_n, lat_s])
        lon = np.linspace(-180, 180, Nx)
        lon_grid, lat_grid = np.meshgrid(lon, lat)

        # matrice pleine : nord = données, sud = NaN
        data_full = np.full((2*Ny, Nx), np.nan)
        data_full[:Ny, :] = data

        # conversion sphérique → cartésien
        theta = np.radians(90 - lat_grid)
        phi   = np.radians(lon_grid)
        R = 1.0
        x_s = R * np.sin(theta) * np.cos(phi)
        y_s = R * np.sin(theta) * np.sin(phi)
        z_s = R * np.cos(theta)

        valid = ~np.isnan(data_full)
        x_v = x_s[valid];  y_v = y_s[valid];  z_v = z_s[valid]
        val = data_full[valid]

        fig = go.Figure()

        # 1) base sphère noire
        fig.add_trace(go.Surface(
            x=x_s, y=y_s, z=z_s,
            surfacecolor=np.zeros_like(z_s),
            colorscale=[[0,"black"],[1,"black"]],
            showscale=False,
            hoverinfo="skip",
        ))

        # 2) points valides colorés
        marker_cfg = dict(
            size=2.5,
            color=val,
            colorscale=colormap,
            showscale=True,
            opacity=1.0,
            colorbar=dict(title=f"{title} {unit}")
        )
        if cmin_val is not None: marker_cfg["cmin"] = cmin_val
        if cmax_val is not None: marker_cfg["cmax"] = cmax_val

        fig.add_trace(go.Scatter3d(
            x=x_v, y=y_v, z=z_v,
            mode="markers",
            marker=marker_cfg,
            name=f"Mer{mission}",
            hoverinfo="skip",
        ))

        fig.update_layout(
    title=f"Globe – {title.upper()} pour Mer{mission}",
    scene=dict(
        xaxis=dict(showticklabels=False, visible=False),
        yaxis=dict(showticklabels=False, visible=False),
        zaxis=dict(showticklabels=False, visible=False),
        aspectmode="data",
        bgcolor="white"
    ),
    margin=dict(l=0, r=0, t=60, b=0),
    annotations=[
        dict(
            text="● Données inexistantes ou non valides",
            showarrow=False,
            xref="paper",
            yref="paper",
            x=0.01,
            y=0.01,
            font=dict(size=14, color="black"),
            bgcolor="white",
            bordercolor="black",
            borderwidth=1,
            borderpad=4
        )
    ]
)

        # sauvegarde
        suffix = "fusion" if is_fusion else "pression"
        fname = f"globe_{suffix}_Mer{mission}.html"
        path = os.path.join("interactive_diagrams", fname)
        fig.write_html(path, include_plotlyjs="cdn")
        print(f"✔️ Globe Mer{mission} enregistré : {path}")

# ===========================================================================
# 3. MAIN EXECUTION BLOCK – run everything & generate figures
# ===========================================================================
if __name__ == "__main__":
    files = {"8": "data_Mer8.csv", "15": "data_Mer15.csv"}
    data = {k: pd.read_csv(os.path.join(DATA_FOLDER, fname)) for k, fname in files.items()}
    results = {k: process_dataset(df, label=k) for k, df in data.items()}
    ratios = ["Mg/Si", "Ca/Si", "Al/Si"]
    plot_interp(results, data, ratios)
    plot_map(results, "pressure_map", "Pression optimale", unit="(GPa)")
    plot_map(results, "F_map", "F optimale", unit="(%)")
    plot_hist_res(results)
    plot_frac_table(results)
    plot_region_counts(results)
    plot_Fopt(results)
    analyze_regions(results)
    plot_globe_map(results, "pressure_map", "Pression optimale", unit="(GPa)")
    plot_globe_map(results, "F_map", "F optimale", unit="(%)")