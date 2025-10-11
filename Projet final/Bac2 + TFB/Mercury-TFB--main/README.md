# TFB – Geochemical Analysis of Mercury's Surface

This repository contains the code and data related to my Bachelor's Thesis (TFB), conducted at the University of Liège. The project focuses on comparing spectral data from NASA's MESSENGER spacecraft with partial melting experimental data.

 
## 🧪 Project Objective

The goal is to compare the elemental ratios Mg/Si, Ca/Si, and Al/Si observed at the surface of Mercury with those obtained from partial melting experiments at different pressures (1.5, 3.5, and 5 GPa).  
This comparison helps to constrain the depth of origin of Mercury’s volcanic rocks and their mantle source conditions.

## 🛠 Modules Used

Confer requirements.txt

## 📊 Key Features

- 📈 Interpolation of experimental data
- 🗺 Generation of geochemical maps
- 🔍 Automatic residual calculation between experimental and planetary values
- 🔺 Interactive ternary diagrams

## 📥 Data

To run the scripts, please create a folder named data in the root directory of the repository and include the following files:

data/

├── data_Mer8.csv    # Experimental data for Mer8 composition (unpublished)

├── data_Mer15.csv   # Experimental data for Mer15 composition (unpublished)

├── mgsi.bmp         # Mg/Si ratio map from MESSENGER

├── casi.bmp         # Ca/Si ratio map from MESSENGER

├── alsi.bmp         # Al/Si ratio map from MESSENGER

├── regions.bmp      # Region classification mask

├── DEM.tif          # Digital Elevation Model

The `.bmp` files contain elemental ratio maps (Mg/Si, Ca/Si, Al/Si) extracted from MESSENGER observations. **Available on the following article** -> https://doi.org/10.1016/j.icarus.2020.113716

The `.csv` files contain the results of partial melting experiments at different pressures (Mer8, Mer15 compositions) but **are not given here** (unpublished data).


## 👨‍🔬 Author

Alexandre Michaux  
Geology student – University of Liège  
