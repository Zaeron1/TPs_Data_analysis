rm(list=ls())
# ==================== PACKAGES ====================
packages <- c("raster","bmp","png","rstudioapi","abind","fields","viridis","plotly","dplyr","ggplot2","tidyr","grid","patchwork")
installed_packages <- packages %in% installed.packages()[,"Package"]
if (any(!installed_packages)) install.packages(packages[!installed_packages])
lapply(packages, library, character.only = TRUE)

# ==================== CHEMINS ====================
base_path <- dirname(rstudioapi::getSourceEditorContext()$path)
directory_path <- file.path(base_path, "data")
result_file <- file.path(base_path, "result_array_with_uncertainty.rds")

# ==================== UTILS ====================
from03602180 <- function(mat) {
  middle <- dim(mat)[2] / 2
  cbind(mat[, (middle + 1):(2 * middle)], mat[, 1:middle])
}

resize_to_target <- function(mat, target_dim = c(720,1440)) {
  if (is.null(dim(mat))) stop("Matrix has no dimensions.")
  if (all(dim(mat) == target_dim)) return(mat)
  r <- raster::raster(mat)
  r_target <- raster::raster(nrows = target_dim[1], ncols = target_dim[2])
  r_resized <- raster::resample(r, r_target, method = "bilinear")
  as.matrix(r_resized)
}

# ==================== FONCTION DE LECTURE DES COUCHES ====================
process_files_to_3d_array <- function(directory_path, correction_factors = NULL, target_dim = c(720, 1440)) {
  # -------------------------------------------------------
  # Charge BMP / PNG / DAT / RDS d'un dossier,
  # applique une éventuelle correction par fichier,
  # redimensionne à target_dim et empile en tableau 3D.
  # (CSV exclus, tous les RDS inclus)
  # -------------------------------------------------------
  
  # Liste des fichiers, en excluant explicitement CSV et TIFF
  file_list <- list.files(directory_path, full.names = TRUE)
  file_list <- file_list[!grepl("\\.csv$", file_list, ignore.case = TRUE)]
  file_list <- file_list[!grepl("\\.tif(f)?$", file_list, ignore.case = TRUE)]
  # Ne garder que les extensions gérées
  file_list <- file_list[grepl("\\.(bmp|png|dat|rds)$", file_list, ignore.case = TRUE)]
  
  n_files <- length(file_list)
  if (n_files == 0) stop("Aucun fichier pris en charge (BMP/PNG/DAT/RDS) dans : ", directory_path)
  
  layer_names <- basename(file_list)
  array_3d <- array(NA_real_, dim = c(target_dim[1], target_dim[2], n_files))
  
  apply_correction <- function(m, cf) {
    if (is.null(cf)) return(m)
    m * cf
  }
  
  for (i in seq_along(file_list)) {
    file_path <- file_list[i]
    file_name <- layer_names[i]
    ext <- tolower(tools::file_ext(file_name))
    message("→ Lecture de ", file_name)
    
    if (ext == "bmp") {
      # BMP : lecture brute, orientation inchangée
      m <- read.bmp(file_path)
      mode(m) <- "numeric"
      m <- apply_correction(m, correction_factors[[file_name]])
      m <- resize_to_target(m, target_dim)
      array_3d[,,i] <- m
      
    } else if (ext == "png") {
      # PNG : via raster -> matrice ; flip vertical pour redresser
      m <- as.matrix(raster::raster(file_path))
      mode(m) <- "numeric"
      m <- apply_correction(m, correction_factors[[file_name]])
      m <- resize_to_target(m, target_dim)
      m <- m[nrow(m):1, ]  # redressement vertical
      array_3d[,,i] <- m
      
    } else if (ext == "dat") {
      # DAT : structure spéciale -> recollage + (dans ta version) inversion L/C
      m <- from03602180(as.matrix(read.table(file_path, header = FALSE)))
      mode(m) <- "numeric"
      m <- apply_correction(m, correction_factors[[file_name]])
      m <- resize_to_target(m, target_dim)
      array_3d[,,i] <- m
      
    } else if (ext == "rds") {
      # RDS : on prend TOUT fichier .rds
      m <- readRDS(file_path)
      mode(m) <- "numeric"
      m <- apply_correction(m, correction_factors[[file_name]])
      m <- resize_to_target(m, target_dim)
      array_3d[,,i] <- m
      
    } else {
      warning("⚠️ Type non supporté/ignoré : ", file_name)
      next
    }
  }
  
  attr(array_3d, "layer_names") <- layer_names
  array_3d
}
# ==================== FACTEURS DE CORRECTION ====================
correction_factors <- list(
  "mgsi.bmp"    = 0.860023 / 255.0,
  "alsi.bmp"    = 0.402477 / 255.0,
  "ssi.bmp"     = 0.161680 / 255.0,
  "fesi.bmp"    = 0.117737 / 255.0,
  "casi.bmp"    = 0.318000 / 255.0,
  "mgsierr.png" = 0.223226 / 255.0,
  "alsierr.png" = 0.153596 / 255.0,
  "ssierr.png"  = 0.0398775 / 255.0,
  "fesierr.png" = 0.0283532 / 255.0,
  "casierr.png" = 0.0809775 / 255.0
)

# ==================== CONSTRUCTION DU CUBE ====================
result_array <- process_files_to_3d_array(directory_path, correction_factors)
layer_names  <- attr(result_array, "layer_names")
cat("✅ Cube initial chargé :", dim(result_array)[3], "couches\n")

# ==================== MASQUE COMMUN (5 RAPPORTS .bmp) ====================
idx_mgsi <- grep("^mgsi\\.bmp$", layer_names, ignore.case = TRUE)
idx_alsi <- grep("^alsi\\.bmp$", layer_names, ignore.case = TRUE)
idx_casi <- grep("^casi\\.bmp$", layer_names, ignore.case = TRUE)
idx_fesi <- grep("^fesi\\.bmp$", layer_names, ignore.case = TRUE)
idx_ssi  <- grep("^ssi\\.bmp$",  layer_names, ignore.case = TRUE)

if (any(lengths(list(idx_mgsi,idx_alsi,idx_casi,idx_fesi,idx_ssi)) == 0)) {
  stop("❌ Introuvable : au moins une des couches mgsi/alsi/casi/fesi/ssi (.bmp).")
}

# On considère qu’une valeur < 0.0001 = pixel invalide
mask_common <- (result_array[,,idx_mgsi]  > 0.0001) &
  (result_array[,,idx_alsi]  > 0.0001) &
  (result_array[,,idx_casi]  > 0.0001) &
  (result_array[,,idx_fesi]  > 0.0001) &
  (result_array[,,idx_ssi]   > 0.0001)

# ==================== CUBE MASQUÉ (NA sur pixels invalides) ====================
masked_cube <- array(NA_real_, dim = dim(result_array))
for (i in 1:dim(result_array)[3]) {
  L <- result_array[,,i]
  L[!mask_common] <- NA_real_
  masked_cube[,,i] <- L
}

# ==================== CONCATÉNER : ORIGINAL + MASQUÉ ====================
result_array_full <- abind::abind(result_array, masked_cube, along = 3)
attr(result_array_full, "layer_names") <- c(layer_names, paste0(layer_names, "_masked"))

cat("✅ Cube final :", dim(result_array_full)[3], "couches (doublé)\n")

# ==================== STATISTIQUES DES COUCHES ====================
na_counts <- sapply(1:dim(result_array_full)[3], function(i) sum(is.na(result_array_full[,,i])))
layer_info <- data.frame(
  Index = seq_along(attr(result_array_full, "layer_names")),
  Layer = attr(result_array_full, "layer_names"),
  NA_Count = na_counts
)
cat("\n=== STATISTIQUES DES COUCHES ===\n")
print(layer_info)

# ==================== SAUVEGARDE ====================
saveRDS(result_array_full, file = result_file)
cat("\n💾 Sauvegardé :", result_file, "\n")

result_array_full[result_array_full == 0 | abs(result_array_full) < 1e-10] <- NA_real_
#nombre de NA dans chaque couche
na_counts <- sapply(1:dim(result_array_full)[3], function(i) sum(is.na(result_array_full[,,i])))
layer_info <- data.frame(Index = seq_along(attr(result_array_full, "layer_names")), Layer = attr(result_array_full, "layer_names"), NA_Count = na_counts)
print(layer_info) # Afficher les informations des couches avec le nombre de



get_layer_as_matrix <- function(result_array_full, layer_index) {
  
  if (layer_index < 1 || layer_index > dim(result_array_full)[3]) { #check si l'index est inf. à 1 (index minimum) ou si il est sup. au nombre de couche de result_array
    stop("Layer index out of bounds.")
  }
  return(result_array_full[,,layer_index])   # Extrait et retourne la couche spécifiée sous forme de matrice
}


library(ggplot2)
library(viridis)

plot_element_map <- function(result_array_full, index, legend_title, map_title,
                             palette = "viridis",  # <--- nouvel argument
                             save = FALSE, file_name = NULL) {
  
  # --- Extraction et orientation correcte ---
  L <- get_layer_as_matrix(result_array_full, index)
  L <- t(L)[, nrow(L):1]
  
  nx <- nrow(L)
  ny <- ncol(L)
  lon <- seq(-180, 180, length.out = nx)
  lat <- seq(-90, 90,  length.out = ny)
  
  df <- expand.grid(Longitude = lon, Latitude = lat)
  df$Z <- as.vector(L)
  
  # --- Grilles principales et secondaires ---
  lon_lines_major <- seq(-180, 180, by = 60)
  lat_lines_major <- seq(-90, 90,  by = 30)
  lon_lines_minor <- seq(-180, 180, by = 30)
  lat_lines_minor <- seq(-90, 90,  by = 15)
  
  # --- Choix dynamique de la palette ---
  if (palette %in% c("viridis", "plasma", "magma", "inferno", "cividis")) {
    fill_scale <- scale_fill_viridis(name = legend_title, option = palette)
  } else if (palette == "terrain") {
    fill_scale <- scale_fill_gradientn(
      name = legend_title,
      colours = terrain.colors(256)
    )
  } else {
    warning("Palette inconnue, utilisation de 'viridis' par défaut.")
    fill_scale <- scale_fill_viridis(name = legend_title, option = "C")
  }
  
  # --- Construction du graphique ---
  p <- ggplot(df, aes(x = Longitude, y = Latitude, fill = Z)) +
    geom_tile() +
    
    # Grille visible (au-dessus)
    geom_vline(xintercept = lon_lines_minor, color = "grey80", linetype = "dotted", linewidth = 0.3, alpha = 0.6) +
    geom_hline(yintercept = lat_lines_minor, color = "grey80", linetype = "dotted", linewidth = 0.3, alpha = 0.6) +
    geom_vline(xintercept = lon_lines_major, color = "grey60", linetype = "dotted", linewidth = 0.5, alpha = 0.8) +
    geom_hline(yintercept = lat_lines_major, color = "grey60", linetype = "dotted", linewidth = 0.5, alpha = 0.8) +
    
    # Palette choisie
    fill_scale +
    
    # Axes et proportion
    scale_x_continuous(breaks = lon_lines_major, expand = c(0, 0)) +
    scale_y_continuous(breaks = lat_lines_major, expand = c(0, 0)) +
    coord_equal(expand = FALSE, xlim = c(-180, 180), ylim = c(-90, 90)) +
    
    # Titres
    labs(
      title = map_title,
      x = "Longitude (°)",
      y = "Latitude (°)"
    ) +
    
    # Thème professionnel
    theme_minimal(base_size = 15) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "bottom",
      legend.key.width = unit(3, "cm"),
      legend.title = element_text(size = 13, face = "bold"),
      legend.text = element_text(size = 11),
      plot.margin = margin(10, 10, 15, 10)
    )
  
  print(p)
  
  # --- Sauvegarde optionnelle ---
  if (save) {
    if (is.null(file_name)) {
      file_name <- paste0("Carte_", gsub("/", "_", legend_title), "_", palette, ".png")
    }
    ggsave(
      filename = file_name,
      plot = p,
      width = 10, height = 5.5, dpi = 400
    )
    message("✅ Figure enregistrée sous : ", file_name)
  }
}

names <- c("Al/Si", "Ca/Si", "Fe/Si", "Mg/Si", "S/Si")
indices <- c(1, 3, 7, 10, 12)

for (k in seq_along(indices)) {
  i <- indices[k]
  plot_element_map(
    result_array_full,
    index = i,
    legend_title = names[k],
    map_title = paste("Carte élémentaire du", names[k]),
    palette = "magma",
    save = FALSE
  )
}




# ==================== BOXPLOTS CLASSIQUES + INCERTITUDES SUPERPOSÉES ====================

# --- 1. Données à partir de ton objet maps ---
ratios <- c("MgSi", "AlSi", "CaSi", "FeSi", "SSi")

df <- do.call(rbind, lapply(ratios, function(r) {
  vals <- as.vector(maps[[r]])
  errs <- as.vector(maps[[paste0(r, "_err")]])
  data.frame(
    Ratio = r,
    Value = vals,
    Error = errs
  )
}))

df <- df[complete.cases(df), ]

# --- 2. Calcul de statistiques par rapport ---
library(dplyr)
err_summary <- df %>%
  group_by(Ratio) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    mean_err = mean(Error, na.rm = TRUE),
    sd_err   = sd(Error, na.rm = TRUE)
  )

# --- 3. Couleurs pour les boxplots et les symboles d’erreur ---
cols_box <- c("darkorange", "steelblue", "seagreen3", "indianred3", "goldenrod")
cols_err <- c("chocolate3", "navy", "forestgreen", "darkred", "darkgoldenrod4")

# --- 4. Tracé des boxplots normaux ---
boxplot(Value ~ Ratio, data = df,
        col = cols_box,
        border = "gray30",
        ylab = "Rapport élémentaire (X/Si)",
        main = "Distribution des rapports élémentaires avec incertitudes",
        outline = FALSE,
        las = 1,
        cex.axis = 0.9)

# --- 5. Superposition des incertitudes moyennes (autre couleur / symbole) ---
# Positions x des boxplots (1:5)
xpos <- 1:length(ratios)

# Points = moyennes, barres = ± écart-type des erreurs
points(xpos, err_summary$mean_val, 
       pch = 21, bg = "white", col = "black", cex = 1.3, lwd = 1.5)
arrows(xpos,
       err_summary$mean_val - err_summary$mean_err,
       xpos,
       err_summary$mean_val + err_summary$mean_err,
       angle = 90, code = 3, length = 0.08, lwd = 2, col = cols_err)

# --- 6. Légende claire ---
legend("topleft",
       legend = c("Distribution (boxplot)", "Incertitude moyenne ±1σ(err)"),
       pch = c(15, 21),
       pt.bg = c(cols_box[1], "white"),
       col = c("gray30", "black"),
       lwd = c(0, 2),
       pt.cex = c(1.8, 1.3),
       bty = "n")

# ==================== DISTRIBUTIONS RÉGIONALES/GLOBALES — FIGURES SÉPARÉES ====================

# --- Packages ---
if (!require(patchwork)) install.packages("patchwork")
library(ggplot2)
library(dplyr)
library(patchwork)

# --- 1. Données ---
maps <- list(
  MgSi = get_layer_as_matrix(result_array_full, 24),
  AlSi = get_layer_as_matrix(result_array_full, 15),
  CaSi = get_layer_as_matrix(result_array_full, 17),
  FeSi = get_layer_as_matrix(result_array_full, 21),
  SSi  = get_layer_as_matrix(result_array_full, 26)
)

region_map   <- get_layer_as_matrix(result_array_full, 28)
region_vals  <- c(1,2,3,4,5,6)
region_names <- c("High-Mg","Al-rich","Caloris","Rach","High-Al NVP","Low-Al NVP")
Region <- factor(as.vector(region_map), levels = region_vals, labels = region_names)

ratios <- names(maps)
df <- do.call(rbind, lapply(ratios, function(r){
  data.frame(Ratio = r, Value = as.vector(maps[[r]]), Region = Region)
})) |> na.omit()

# --- 2. Palette géologique ---
cols_base <- c(
  "High-Mg"="#E4572E", "Al-rich"="#17BEBB", "Caloris"="#FFC914",
  "Rach"="#4B4E6D", "High-Al NVP"="#76B041", "Low-Al NVP"="#A23B72"
)

# --- 3. Fonction pour créer un couple avec légende ---
make_pair <- function(r){
  d  <- df[df$Ratio==r, ]
  xr <- range(d$Value, na.rm = TRUE)
  
  # Distribution régionale
  p_reg <- ggplot(d, aes(Value, fill=Region, color=Region)) +
    geom_density(alpha=.35, linewidth=1, adjust=1.2) +
    scale_fill_manual(values=cols_base, name="Régions géologiques") +
    scale_color_manual(values=cols_base, guide="none") +
    scale_x_continuous(limits = xr) +
    labs(title=paste("Distribution régionale de", r), x=NULL, y="Densité") +
    theme_minimal(base_size=12) +
    theme(
      plot.title = element_text(face="bold", hjust=.5),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  # Distribution globale
  p_glb <- ggplot(d, aes(Value)) +
    geom_density(fill="#CCCCCC", color="#444444", alpha=.6, linewidth=1.1, adjust=1.2) +
    scale_x_continuous(limits = xr) +
    labs(title=paste("Distribution globale de", r), x=r, y="Densité") +
    theme_minimal(base_size=12) +
    theme(
      plot.title = element_text(face="bold", hjust=.5),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
  
  # Empilement vertical avec légende commune
  p_reg / p_glb + plot_layout(heights = c(1, 1), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.text = element_text(size = 10),
      legend.title = element_text(face = "bold")
    )
}

# --- 4. Génération et affichage séparé des figures ---
pairs <- lapply(ratios, make_pair)

for (i in seq_along(ratios)) {
  print(pairs[[i]])
  cat("\n✅ Figure affichée pour", ratios[i], "\n")
}

# --- Option : Export individuel ---
# ggsave(paste0("distribution_", r, ".png"), pairs[[i]], width = 8, height = 6, dpi = 300)


# ==================== CARTES D'ANOMALIES NORMALISÉES (ggplot2 + grille 60°/30°) ====================



library(ggplot2)
library(dplyr)
library(viridis)

# --- 1. Indices et noms des couches ---
indices <- c(1, 3, 7, 10, 12)
names   <- c("Al/Si", "Ca/Si", "Fe/Si", "Mg/Si", "S/Si")

# --- 2. Fonction : calcul des anomalies Z-score bornées à ±3 ---
get_anomaly_layer <- function(mat) {
  mu  <- mean(mat, na.rm = TRUE)
  sdv <- sd(mat, na.rm = TRUE)
  Z   <- (mat - mu) / sdv
  Z[Z >  3] <-  3
  Z[Z < -3] <- -3
  return(Z)
}

# --- 3. Calcul des anomalies pour chaque couche sélectionnée ---
result_anomalies <- list()

for (k in seq_along(indices)) {
  i <- indices[k]
  layer <- get_layer_as_matrix(result_array_full, i)
  result_anomalies[[names[k]]] <- get_anomaly_layer(layer)
}

# --- 4. Boucle d’affichage avec ta fonction pro ---
for (k in seq_along(result_anomalies)) {
  cat("→ Génération de la carte :", names[k], "\n")
  
  # Transformation en tableau 3D (nécessaire pour plot_element_map)
  anomaly_array <- array(
    result_anomalies[[k]],
    dim = c(nrow(result_anomalies[[k]]),
            ncol(result_anomalies[[k]]),
            1)
  )
  
  # Appel de la fonction
  plot_element_map(
    result_array_full = anomaly_array,
    index = 1,
    legend_title = paste0(names[k], " (Z-score)"),
    map_title = paste("Anomalies normalisées du", names[k]),
    palette = "terrain",   # 🔥 choisis : "viridis", "magma", "terrain", etc.
    save = FALSE
  )
}

# ==================== DISTRIBUTIONS DENSITÉ & DEM (régions + globales séparées) ====================

# --- Installation et chargement des packages nécessaires ---
if (!require(patchwork)) install.packages("patchwork")
library(patchwork)
library(ggplot2)
library(dplyr)

# --- 1. Extraction des couches ---
Density <- get_layer_as_matrix(result_array_full, 20)  # DensityGrid.dat_masked
DEM     <- get_layer_as_matrix(result_array_full, 19)  # DEM.RDS_masked
region_map <- get_layer_as_matrix(result_array_full, 28)

# --- 2. Définition des régions ---
region_values <- c(1, 2, 3, 4, 5, 6)
region_names  <- c("High-Mg", "Al-rich", "Caloris", "Rach", "High-Al NVP", "Low-Al NVP")

# --- 3. Construction du dataframe principal ---
Region <- factor(as.vector(region_map), levels = region_values, labels = region_names)
df <- data.frame(
  Density = as.vector(Density),
  DEM     = as.vector(DEM),
  Region  = Region
) |> na.omit()

# --- 4. Palette harmonieuse ---
cols_base <- c(
  "High-Mg"     = "#E4572E",  # orange
  "Al-rich"     = "#17BEBB",  # turquoise
  "Caloris"     = "#FFC914",  # jaune
  "Rach"        = "#4B4E6D",  # bleu-gris
  "High-Al NVP" = "#76B041",  # vert clair
  "Low-Al NVP"  = "#A23B72"   # magenta
)

# === FIGURE 1 : DENSITY ===

## A. DENSITY par région
p_density_region <- ggplot(df, aes(x = Density, fill = Region, color = Region)) +
  geom_density(alpha = 0.35, linewidth = 1.0, adjust = 1.2) +
  scale_fill_manual(values = cols_base, name = "Régions géologiques") +
  scale_color_manual(values = cols_base, guide = "none") +
  labs(title = "Distribution régionale de la densité",
       x = NULL, y = "Densité de probabilité") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold")
  )

## B. DENSITY globale
p_density_global <- ggplot(df, aes(x = Density)) +
  geom_density(fill = "#A23B72", color = "#7A2853", alpha = 0.4, linewidth = 1.1) +
  labs(title = "Distribution globale de la densité",
       x = "Densité (unités relatives)", y = "Densité de probabilité") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

## Combinaison verticale avec légende commune
plot_density <- p_density_region / p_density_global +
  plot_layout(heights = c(1, 0.6), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

plot_density  # <-- Affiche la figure DENSITY avec légende


# === FIGURE 2 : DEM ===

## C. DEM par région
p_dem_region <- ggplot(df, aes(x = DEM, fill = Region, color = Region)) +
  geom_density(alpha = 0.35, linewidth = 1.0, adjust = 1.2) +
  scale_fill_manual(values = cols_base, name = "Régions géologiques") +
  scale_color_manual(values = cols_base, guide = "none") +
  labs(title = "Distribution régionale de l'altitude (DEM)",
       x = NULL, y = "Densité de probabilité") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "bottom",
    legend.justification = "center",
    legend.text = element_text(size = 10),
    legend.title = element_text(face = "bold")
  )

## D. DEM global
p_dem_global <- ggplot(df, aes(x = DEM)) +
  geom_density(fill = "#999999", color = "#555555", alpha = 0.4, linewidth = 1.1) +
  labs(title = "Distribution globale de l'altitude",
       x = "Altitude (m)", y = "Densité de probabilité") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

## Combinaison verticale avec légende commune
plot_dem <- p_dem_region / p_dem_global +
  plot_layout(heights = c(1, 0.6), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

plot_dem  # <-- Affiche la figure DEM avec légende



