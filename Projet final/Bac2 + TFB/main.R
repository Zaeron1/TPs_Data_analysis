#------------------------------------------------------------------------------#
# Projet de Modélisation numérique, analyse satellite de Mercure      13/08/24 #
# Ecrit par: Alexandre Michaux, Anatole Devolder et Beniamino Orsini           #
# Ce script se base principalement sur l'utlisation du package "mercure.Ulg"   #
# les fonctions qui le compose ont été écrite dans le dossier "Fonctions"      #
# partagé avec le dossier principale "Projet Mercure"                          #
#------------------------------------------------------------------------------#


rm(list=ls()) # Supprime toutes les variables de l'environnement pour commencer avec un espace de travail vide

packages <- c("plotly", "raster", "bmp", "devtools", "rstudioapi")
installed_packages <- packages %in% installed.packages()[,"Package"]
if (any(!installed_packages)) {
  install.packages(packages[!installed_packages])
}

library(plotly); library(raster);library(bmp); library(devtools);library(rstudioapi)
base_path <- dirname(rstudioapi::getSourceEditorContext()$path)
directory_path <- file.path(base_path, "mercury1440x720")
result_file <- file.path(base_path, "result_array.rds")
devtools::install_github("Zaeron1/mercure.ulg", force = TRUE)
library(mercure.ulg)


#pour des raisons inconnues, la fonctions définies ci-dessous n'a pas pu être chargée dans notre package :/
process_files_to_3d_array <- function(directory_path, correction_factors, target_dim = c(720, 1440)) {
  file_list <- list.files(directory_path, full.names = TRUE) # Obtenir la liste des fichiers dans le répertoire
  n_files <- length(file_list) # Calculer le nombre total de fichiers
  array_3d <- array(0, dim = c(target_dim[1], target_dim[2], n_files)) # Initialiser un tableau 3D vide avec les dimensions cibles
  layer_names <- basename(file_list) # Initialiser un vecteur pour stocker les noms des fichiers (couches)
  apply_correction <- function(matrix, correction_factor) { # Définir une fonction d'aide pour appliquer un facteur de correction
    if (is.null(correction_factor)) correction_factor <- 1 # Si aucun facteur de correction n'est spécifié, utiliser 1 par défaut
    matrix * correction_factor # Appliquer le facteur de correction
  }
  for (i in seq_along(file_list)) { # Boucle pour traiter chaque fichier
    file_path <- file_list[i] # Chemin du fichier actuel
    file_name <- layer_names[i] # Nom du fichier actuel
    file_extension <- tolower(tools::file_ext(file_name)) # Obtenir l'extension du fichier en minuscules
    if (file_extension == "bmp") { # Traiter les fichiers BMP
      bmp_matrix <- read.bmp(file_path) # Lire le fichier BMP
      corrected_matrix <- apply_correction(bmp_matrix, correction_factors[[file_name]]) # Appliquer la correction
      array_3d[,,i] <- corrected_matrix # Stocker le résultat dans le tableau 3D
    } else if (file_extension == "dat") { # Traiter les fichiers DAT
      data_matrix <- from03602180(as.matrix(read.table(file_path, header = FALSE))) # Lire le fichier DAT
      array_3d[,,i] <- data_matrix # Stocker le résultat dans le tableau 3D
    } else if (file_extension == "tif") { # Traiter les fichiers TIF
      tif_raster <- raster(file_path) # Lire le fichier TIF en tant que raster
      resize_factor <- c(max(1, ncol(tif_raster) / target_dim[2]), max(1, nrow(tif_raster) / target_dim[1])) # Calculer le facteur de redimensionnement
      resized_raster <- aggregate(tif_raster, fact = resize_factor) # Redimensionner le raster
      resized_matrix <- as.matrix(resized_raster, ncol = target_dim[2], nrow = target_dim[1]) # Convertir le raster redimensionné en matrice
      array_3d[,,i] <- resized_matrix * 0.5 # Stocker le résultat dans le tableau 3D après correction
    } else if (file_extension == "csv") { # Traiter les fichiers CSV
      lines <- gsub(",\\s+", ",", gsub("(\\d)\\s+(\\d)", "\\1,\\2", readLines(file_path))) # Nettoyer les lignes du fichier CSV
      temp_file <- tempfile(fileext = ".csv") # Créer un fichier temporaire pour le CSV nettoyé
      writeLines(lines, temp_file) # Écrire les lignes nettoyées dans le fichier temporaire
      donnees <- read.csv(temp_file, header = FALSE) # Lire le fichier CSV temporaire
      matrix_csv <- from03602180(as.matrix(donnees[,-1])) # Convertir les données en matrice
      array_3d[,,i] <- matrix_csv # Stocker le résultat dans le tableau 3D
      unlink(temp_file) # Supprimer le fichier temporaire
    } else {
      warning(paste("File type not supported:", file_extension, "for file", file_name)) # Avertir si le type de fichier n'est pas supporté
    }
  }
  layer_info <- data.frame(Index = seq_along(layer_names), Layer = layer_names) # Créer un tableau de données avec les indices et noms des couches
  print(layer_info) # Afficher les informations des couches
  return(array_3d) # Retourner le tableau 3D
}

correction_factors <- list("mgsi.bmp" = 0.860023 / 255.0, "alsi.bmp" = 0.402477 / 255.0, "ssi.bmp" = 0.161680 / 255.0, "fesi.bmp" = 0.117737 / 255.0, "casi.bmp" = (0.318000/255.0)) # Facteurs de correction pour les fichiers BMP
result_array <- process_files_to_3d_array(directory_path, correction_factors) # Appeler la fonction pour traiter les fichiers et obtenir le tableau 3D
# IMPORTANT: le tableau généré dans la console renseigne les index qu'il faut utiliser dans les prochaines fonctions poru extraire les bonnes matrices du 3D array
saveRDS(result_array, file = result_file) # Sauvegarder le tableau 3D dans un fichier RDS
#result_array <- readRDS(result_file) # Lire le tableau 3D depuis le fichier RDS
rm(list=ls()[!ls() %in% c("result_array")]) # Supprimer toutes les variables sauf le tableau 3D

layer_matrix <- get_layer_as_matrix(result_array, 2) # EXEMPLE: Extraire une couche spécifique (index 2 ici) du tableau 3D sous forme de matrice
plot_matrix(4, "DEM de Mercure", "Élévation (m)", color_palette = colorRampPalette(c("blue", "cyan", "yellowgreen", "yellow", "orange", "orangered", "darkred"))(9)) # Tracer une matrice avec une palette de couleurs spécifiqu, en premier argument de la fonction on peut spécifier l'index ou alors le nom d'une matrice déja présente dans la mémoire (importer grace à la ligne précédente)
Classy <- classify_matrix(8) # Classifier une matrice spécifique en décile
plot_matrix(Classy, "Carte des déciles du rapport Mg/Si", "Décile", colorRampPalette(c("blue", "cyan", "yellowgreen", "yellow", "orange", "orangered", "darkred"))(9)) # Tracer la carte des déciles
plot_max_elevation(4) # Tracer l'élévation maximale pour une couche spécifique
plot_3D_surface(4, "DEM 3D") # Tracer une surface 3D pour une couche spécifique
analyse_regressions(5, c(1, 2, 6, 8, 9)) # Analyser les régressions entre différentes couches
correl(3, 7) # Calculer la corrélation entre deux couches spécifiques
plot_matrix(10, "Région de Mercure", "N° de région", color_palette = colorRampPalette(c("blue", "cyan", "yellowgreen", "yellow", "orange", "orangered", "darkred"))(7)) # Tracer la matrice de la région de Mercure
region(10, 8) # Analyser du rapport moyen de chaque région
generate_histograms(8, 10) # Générer des histogrammes
calculate_and_plot_composition() # Calculer et tracer la composition des éléments présents à la surface de Mercure



