rm(list=ls())

# ==================== PACKAGES ====================
packages <- c("raster","bmp","png","rstudioapi","abind")
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





# ==================== CALCUL DES CARTES DE PRESSION ET DE FUSION ====================

library(plotly)

# --- 1. Chargement des données expérimentales Mer8 et Mer15 ---
mer8_path  <- "/Users/alexandremichaux/Documents/UCA/Cours/Analyse des données/TP/TPs/Projet final/data/data_Mer8.csv"
mer15_path <- "/Users/alexandremichaux/Documents/UCA/Cours/Analyse des données/TP/TPs/Projet final/data/data_Mer15.csv"

data_Mer8  <- read.csv(mer8_path,  sep = ",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
data_Mer15 <- read.csv(mer15_path, sep = ",", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

# Combiner les deux jeux expérimentaux
exp_data <- rbind(data_Mer8, data_Mer15)
exp_data <- exp_data[complete.cases(exp_data[, c("Mg/Si","Al/Si","Ca/Si","Fe/Si","S/Si")]), ]

# --- 2. Extraction des couches masquées depuis result_array_full ---
maps <- list(
  MgSi     = get_layer_as_matrix(result_array_full, 24),
  MgSi_err = get_layer_as_matrix(result_array_full, 25),
  AlSi     = get_layer_as_matrix(result_array_full, 15),
  AlSi_err = get_layer_as_matrix(result_array_full, 16),
  CaSi     = get_layer_as_matrix(result_array_full, 17),
  CaSi_err = get_layer_as_matrix(result_array_full, 18),
  FeSi     = get_layer_as_matrix(result_array_full, 21),
  FeSi_err = get_layer_as_matrix(result_array_full, 22),
  SSi      = get_layer_as_matrix(result_array_full, 26),
  SSi_err  = get_layer_as_matrix(result_array_full, 27)
)

nx <- nrow(maps$MgSi)
ny <- ncol(maps$MgSi)

# --- 3. Matrices vides pour Pression et F ---
pressure_map <- matrix(NA, nrow = nx, ncol = ny)
fusion_map   <- matrix(NA, nrow = nx, ncol = ny)

# --- 4. Fonction de calcul du résidu pondéré ---
compute_residual <- function(M, sigma, E) {
  valid <- !is.na(M) & !is.na(sigma)
  if (sum(valid) < 3) return(NA)
  sqrt(sum(((M[valid] - E[valid]) / sigma[valid])^2))
}

# --- 5. Boucle principale : minimisation du résidu pondéré ---
pb <- txtProgressBar(min = 0, max = nx, style = 3)

for (x in 1:nx) {
  for (y in 1:ny) {
    M <- c(maps$MgSi[x,y], maps$AlSi[x,y], maps$CaSi[x,y],
           maps$FeSi[x,y], maps$SSi[x,y])
    sigma <- c(maps$MgSi_err[x,y], maps$AlSi_err[x,y], maps$CaSi_err[x,y],
               maps$FeSi_err[x,y], maps$SSi_err[x,y])
    if (all(is.na(M))) next
    
    residuals <- apply(exp_data[, c("Mg/Si","Al/Si","Ca/Si","Fe/Si","S/Si")], 1, function(E){
      compute_residual(M, sigma, E)
    })
    
    best_index <- which.min(residuals)
    if (is.na(best_index)) next
    
    pressure_map[x,y] <- exp_data$Pression[best_index]
    fusion_map[x,y]   <- exp_data$F[best_index]
  }
  setTxtProgressBar(pb, x)
}
close(pb)

lst <- list(
  pressure_map = pressure_map,
  fusion_map = fusion_map
)

for (name in lst) {
  image(
    z = t(apply(name, 2, rev)),   # rotation correcte pour l'affichage
    col = terrain.colors(100),
    main = "",
    axes = T,
    xlab = "Longitude",
    ylab = "Latitude"
    
  )}

# --- 7. Sauvegarde des résultats ---
saveRDS(pressure_map, file.path(base_path, "pressure_map_masked_weighted.rds"))
saveRDS(fusion_map,   file.path(base_path, "fusion_map_masked_weighted.rds"))
cat("\n✅ Cartes de pression et de fusion pondérées sauvegardées.\n")


  
