"""
This script displays a grayscale Digital Elevation Model (DEM) of Mercury
and overlays a semi-transparent region mask (from a BMP file).
Each region is assigned a unique transparent color. The map includes:
• Latitude/longitude axes
• Colored overlay of six geologic regions
• A custom legend
"""


from PIL import Image, ImageFile
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
import os

# ————————————————————————————————————————————————————————————————
# Display settings: globally increase font sizes for clarity
# ————————————————————————————————————————————————————————————————
plt.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 18,
    "axes.labelsize": 16,
    "legend.fontsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.titlesize": 20,
})

# Allow loading of corrupted or incomplete image files
ImageFile.LOAD_TRUNCATED_IMAGES = True

# Disable pixel size limit for very large TIFF files
Image.MAX_IMAGE_PIXELS = None

# Load the DEM (grayscale TIFF file)
DATA_folder = "data"
DEM_path = os.path.join(DATA_folder, 'DEM.tif')
DEM = Image.open(DEM_path)

# Load the geological region mask (BMP image)
region_path = os.path.join(DATA_folder, 'regions.bmp')
mask_image = Image.open(region_path)
mask_array = np.array(mask_image)

# Make sure the DEM and mask are the same size; resize if necessary
if DEM.size != mask_image.size:
    print("Resizing of DEM to fit with the mask.")
    DEM = DEM.resize(mask_image.size, Image.BILINEAR)

# Convert DEM image to a NumPy array
DEM_array = np.array(DEM)

# Define color mapping for the region mask
# Value 0 = fully transparent; values 1–6 = colored overlays
mask_colors = [
    (0, 0, 0, 0),       # 0 : transparent
    (1, 0, 0, 0.5),     # 1 : semi-transparent red
    (0, 1, 0, 0.5),     # 2 : semi-transparent green
    (0, 0, 1, 0.5),     # 3 : semi-transparent blue
    (1, 1, 0, 0.5),     # 4 : semi-transparent yellow
    (1, 0, 1, 0.5),     # 5 : semi-transparent magenta
    (0, 1, 1, 0.5)      # 6 : semi-transparent cyan
]
mask_cmap = ListedColormap(mask_colors)

# Define boundaries so that integer values (0–6) map correctly to colors
norm = BoundaryNorm(boundaries=[-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5], ncolors=mask_cmap.N)

# Define geographical extent and axis ticks (longitude and latitude)
extent = [-180, 180, -90, 90]
longitude_ticks = np.arange(-180, 181, 60)
latitude_ticks = np.arange(-90, 91, 30)

# Create the figure and axis
fig, axs = plt.subplots(figsize=(10, 8))

# Display the DEM as grayscale
axs.imshow(DEM_array, extent=extent, origin='upper', interpolation='none', cmap='gray')
axs.set_title("")
axs.set_xlabel("Longitude (°)")
axs.set_ylabel("Latitude (°)")
axs.set_xticks(longitude_ticks)
axs.set_yticks(latitude_ticks)
axs.grid(True)

# Overlay the region mask with transparent colors
axs.imshow(mask_array, extent=extent, origin='upper', interpolation='none', cmap=mask_cmap, norm=norm)

# Add a legend with region names and corresponding patch colors
legend_elements = [
    Patch(facecolor=(1, 0, 0, 0.5), edgecolor='r', label='High-Mg'),
    Patch(facecolor=(0, 1, 0, 0.5), edgecolor='g', label='Al-rich'),
    Patch(facecolor=(0, 0, 1, 0.5), edgecolor='b', label='Caloris'),
    Patch(facecolor=(1, 1, 0, 0.5), edgecolor='y', label='Rach'),
    Patch(facecolor=(1, 0, 1, 0.5), edgecolor='m', label='High-al NVP'),
    Patch(facecolor=(0, 1, 1, 0.5), edgecolor='c', label='Low-al NVP')
]
axs.legend(handles=legend_elements, loc='lower right', fontsize='small', title='Régions', title_fontsize='medium')

plt.show()