# =======================================================
#  Affichage brut du CrustalThickness_ModelV1.csv
#  sans inversion ni redimensionnement
# =======================================================

# ==================== PACKAGES ====================
packages <- c("raster")
installed_packages <- packages %in% installed.packages()[, "Package"]
if (any(!installed_packages)) install.packages(packages[!installed_packages])
lapply(packages, library, character.only = TRUE)

# ==================== FONCTION D’AFFICHAGE ====================
open_device <- function() {
  if (.Platform$OS.type == "windows") windows(width = 8, height = 5)
  else if (capabilities("aqua")) quartz(width = 8, height = 5)  # macOS
  else dev.new(width = 8, height = 5)
}

# ==================== CHEMINS ====================
csv_path <- "/Users/alexandremichaux/Documents/UCA/Cours/Analyse des données/TP/TPs/Projet final/data/CrustalThickness_ModelV1.csv"
output_rds <- sub("\\.csv$", ".RDS", csv_path)

cat("→ Lecture du fichier CSV :", csv_path, "\n")

# ==================== LECTURE ET NETTOYAGE ====================
lines <- readLines(csv_path, warn = FALSE)
lines <- gsub(",\\s+", ",", gsub("(\\d)\\s+(\\d)", "\\1,\\2", lines))
temp_file <- tempfile(fileext = ".csv")
writeLines(lines, temp_file)

donnees <- read.csv(temp_file, header = FALSE)
unlink(temp_file)

# ==================== CONVERSION EN MATRICE ====================
if (ncol(donnees) > 2) {
  mat <- as.matrix(donnees[, -1])  # Ignore une colonne d’index éventuelle
} else {
  mat <- as.matrix(donnees)
}

mat <- mat[nrow(mat):1, ]
from03602180 <- function(mat) {
  middle <- dim(mat)[2] / 2
  mat1 <- mat[, 1:middle]
  mat2 <- mat[, (middle + 1):(2 * middle)]
  return(cbind(mat2, mat1))
}
mat <- from03602180(mat)

mat <- t(mat)
mode(mat) <- "numeric"

# ==================== SAUVEGARDE ====================
saveRDS(mat, file = output_rds)
cat("✅ Fichier RDS sauvegardé :", output_rds, "\n")

# ==================== AFFICHAGE ====================
open_device()

zmin <- min(mat, na.rm = TRUE)
zmax <- max(mat, na.rm = TRUE)

image(
  mat,
  col = terrain.colors(256),
  axes = FALSE,
  main = sprintf("Crustal Thickness 5?", zmin, zmax)
)

# Échelle de couleur
colorbar_ticks <- pretty(seq(zmin, zmax, length.out = 10))
legend(
  "bottomleft",
  legend = round(colorbar_ticks, 1),
  fill = terrain.colors(length(colorbar_ticks)),
  title = "Épaisseur crustale (km)",
  cex = 0.8
)

