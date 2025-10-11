#!/usr/bin/env python3
"""
mercury_maps.py – Load elemental abundance maps (Mg/Si, Ca/Si, Al/Si) of Mercury,
build derived ratios and plot either the whole planet or only the northern
hemisphere.

Changelog (May 2025) ---------------------------------------------------------
* **show_borders** – overlay region outlines in solid black.
* **border_thick** – outline thickness.
* **auto-titles & colour-bar** – `plot_ratio_map()` now sets the figure title and
  colour-bar label automatically from *ratio* and *north_only*.

Quick example ::
    >>> from mercury_maps import regions_array, plot_ratio_map
    >>> cube = regions_array(north_only=True)
    >>> mg_si = cube[0, 6]              # Mg/Si ratio, global slice
    >>> plot_ratio_map(mg_si, ratio="Mg/Si", north_only=True)
"""



from __future__ import annotations
import os
from typing import Dict, Tuple

import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


plt.rcParams.update({
    "font.size": 14,      # base font size
    "axes.titlesize": 18, # subplot titles
    "axes.labelsize": 12, # axis labels
    "legend.fontsize": 20,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 20,  # top‑level figure title
})

# ----------------------------------------------------------------------------
# CONSTANTS & SETTINGS
# ----------------------------------------------------------------------------
DATA_FOLDER: str = "data"
FILE_NAMES: Tuple[str, ...] = ("mgsi.bmp", "casi.bmp", "alsi.bmp", "regions.bmp")
VARIABLES: Tuple[str, ...] = ("mg_si", "ca_si", "al_si", "regions")
REGION_VALUES: Tuple[int, ...] = (1, 2, 3, 4, 5, 6)
REGION_NAMES: Tuple[str, ...] = (
    "high-Mg",
    "Al-rich",
    "Caloris",
    "Rach",
    "high-Al NVP",
    "low-Al NVP",
)
SCALE_FACTORS: Dict[str, float] = {
    "mg_si": 0.860023,
    "ca_si": 0.318000,
    "al_si": 0.402477,
}
CHANNELS_ORDER: Tuple[str, ...] = (
    "mg_si",
    "ca_si",
    "al_si",
    "ca_mg",
    "al_mg",
    "ca_al",
)

# Map channel index → pretty ratio label (for auto-title)
CHANNEL_LABELS: Tuple[str, ...] = (
    "Mg/Si", "Ca/Si", "Al/Si", "Ca/Mg", "Al/Mg", "Ca/Al"
)

# ----------------------------------------------------------------------------
# CORE UTILITIES
# ----------------------------------------------------------------------------

def _load_raw_images() -> Dict[str, Image.Image]:
    missing = [f for f in FILE_NAMES if not os.path.exists(os.path.join(DATA_FOLDER, f))]
    if missing:
        raise FileNotFoundError(f"Missing image files in '{DATA_FOLDER}': {', '.join(missing)}")
    return {var: Image.open(os.path.join(DATA_FOLDER, f)).convert("L") for var, f in zip(VARIABLES, FILE_NAMES)}


def _slice_rows(arr: np.ndarray, north_only: bool) -> np.ndarray:
    return arr[: arr.shape[0] // 2] if north_only else arr


def _degree_per_pixel(mat: np.ndarray, north_only: bool) -> Tuple[float, float]:
    h, w = mat.shape
    return 360.0 / w, (90.0 if north_only else 180.0) / h


def _region_boundaries(mask: np.ndarray) -> np.ndarray:
    b = np.zeros_like(mask, bool)
    b[:-1, :] |= mask[:-1, :] != mask[1:, :]
    b[:, :-1] |= mask[:, :-1] != mask[:, 1:]
    return b


def _dilate(mask: np.ndarray, iters: int = 1) -> np.ndarray:
    m = mask.copy()
    for _ in range(iters):
        up    = np.pad(m[1:, : ], ((0,1),(0,0)), constant_values=False)
        down  = np.pad(m[:-1, :], ((1,0),(0,0)), constant_values=False)
        left  = np.pad(m[:, 1: ], ((0,0),(0,1)), constant_values=False)
        right = np.pad(m[:, :-1], ((0,0),(1,0)), constant_values=False)
        m |= up | down | left | right
    return m

# ----------------------------------------------------------------------------
# PUBLIC API
# ----------------------------------------------------------------------------

def regions_array(north_only: bool = True) -> np.ndarray:
    imgs = _load_raw_images()
    reg_mask_full = np.asarray(imgs["regions"], float)
    reg_mask = _slice_rows(reg_mask_full, north_only)

    data: Dict[str, np.ndarray] = {}
    for var in ("mg_si", "ca_si", "al_si"):
        full = np.asarray(imgs[var], float)
        sub = _slice_rows(full, north_only) * (SCALE_FACTORS[var] / 255.0)
        sub[sub == 0] = np.nan
        data[var] = sub

    valid = ~(np.isnan(data["mg_si"]) | np.isnan(data["ca_si"]) | np.isnan(data["al_si"]))
    for v in data:
        data[v] = np.where(valid, data[v], np.nan)

    data["ca_mg"] = data["ca_si"] / data["mg_si"]
    data["al_mg"] = data["al_si"] / data["mg_si"]
    data["ca_al"] = data["ca_si"] / data["al_si"]

    cube = np.stack([data[ch] for ch in CHANNELS_ORDER], axis=-1)

    mats: Dict[str, np.ndarray] = {}
    for val, name in zip(REGION_VALUES, REGION_NAMES):
        mask = (reg_mask == val) & valid
        mats[name] = np.where(mask[..., None], cube, np.nan)
    for i, ch in enumerate(CHANNELS_ORDER):
        mats[f"full_{ch}"] = cube[..., i]

    out = np.empty((len(CHANNELS_ORDER), len(REGION_NAMES) + 1), dtype=object)
    for i, ch in enumerate(CHANNELS_ORDER):
        for j, r in enumerate(REGION_NAMES):
            out[i, j] = mats[r][..., i]
        out[i, -1] = mats[f"full_{ch}"]
    return out


def plot_ratio_map(
    matrix: np.ndarray,
    *,
    ratio: str | None = None,
    north_only: bool = True,
    title: str | None = None,
    cmap: str = "jet",
    show_borders: bool = False,
    border_thick: int = 1,
    cbar_pad: float = 0.6,
    show: bool = True,
    save_path: str | None = None,
) -> None:
    """Display *matrix* with automatic title and colour-bar label.

    Parameters
    ----------
    ratio : str | None
        Label of the ratio (e.g. "Ca/Mg"). If *None*, defaults to "Ratio".
    north_only : bool
        True → northern hemisphere; False → global map.
    title : str | None
        Custom title. If *None*, auto-generated from *ratio* and *north_only*.
    """

    if matrix.ndim != 2:
        raise ValueError("`matrix` must be 2-D")

    ratio = ratio or "Ratio"
    if title is None:
        title = f""
    # Figure size ---------------------------------------------------------
    dx, dy = _degree_per_pixel(matrix, north_only)
    fig_w = 10.0
    fig_h = fig_w * (dy / dx)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    extent = [-180, 180, 0, 90] if north_only else [-180, 180, -90, 90]
    ax.imshow(matrix, cmap=cmap, extent=extent, origin="upper", interpolation="nearest")

    # Region outlines -----------------------------------------------------
    if show_borders:
        reg_slice = _slice_rows(np.asarray(_load_raw_images()["regions"], float), north_only)
        borders = _region_boundaries(reg_slice)
        if border_thick > 1:
            borders = _dilate(borders, border_thick - 1)
        ax.imshow(np.ma.masked_where(~borders, borders),
                  cmap="gray_r", vmin=0, vmax=1, alpha=1.0,
                  extent=extent, origin="upper", interpolation="nearest")

    # Axis formatting -----------------------------------------------------
    ax.set_title(title)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_xticks(np.linspace(-180, 180, 9))
    ax.set_yticks(np.linspace(0, 90, 4) if north_only else np.linspace(-90, 90, 7))
    ax.grid(True, linestyle="--", alpha=0.5)

    # Colour-bar ----------------------------------------------------------
    div = make_axes_locatable(ax)
    cax = div.append_axes("bottom", size="5%", pad=cbar_pad)
    cb = fig.colorbar(ax.images[0], cax=cax, orientation="horizontal")
    cb.set_label(ratio)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300)
    if show:
        plt.show()
    else:
        plt.close(fig)


# ----------------------------------------------------------------------------
# Quick test -----------------------------------------------------------------
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    SHOW_NORTH_ONLY = False
    SHOW_BORDERS   = True
    BORDER_THICK   = 2
    CHANNEL_IDX    = 2     # 0=Mg/Si, 1=Ca/Si, 2=Al/Si, ...

    cube = regions_array(north_only=SHOW_NORTH_ONLY)
    mat = cube[CHANNEL_IDX, 6]

    plot_ratio_map(mat,
                   ratio=CHANNEL_LABELS[CHANNEL_IDX],
                   north_only=SHOW_NORTH_ONLY,
                   show_borders=SHOW_BORDERS,
                   border_thick=BORDER_THICK)
