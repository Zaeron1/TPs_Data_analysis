# ==================== PACKAGES ====================
library(raster)

# ==================== CHEMINS ====================
tif_path <- "/Users/alexandremichaux/Documents/UCA/Cours/Analyse des données/TP/TPs/Projet final/data/DEM.tif"
output_rds <- "/Users/alexandremichaux/Documents/UCA/Cours/Analyse des données/TP/TPs/Projet final/data/DEM.RDS"
target_dim <- c(720, 1440)

# ==================== VÉRIFICATION ====================
if (!file.exists(tif_path)) stop(paste("❌ DEM introuvable :", tif_path))
cat("✅ Fichier DEM trouvé :", tif_path, "\n")

# ==================== FONCTION from03602180 ====================
from03602180 <- function(mat) {
  middle <- dim(mat)[2] / 2
  mat1 <- mat[, 1:middle]
  mat2 <- mat[, (middle + 1):(2 * middle)]
  return(cbind(mat2, mat1))
}

# ==================== LECTURE & REDIMENSIONNEMENT ====================
cat("→ Lecture du DEM...\n")
tif_raster <- raster(tif_path)

resize_factor <- c(
  max(1, ncol(tif_raster) / target_dim[2]),
  max(1, nrow(tif_raster) / target_dim[1])
)
cat("→ Facteur de redimensionnement :", resize_factor, "\n")

resized_raster <- aggregate(tif_raster, fact = resize_factor)
resized_matrix <- as.matrix(resized_raster,
                            ncol = target_dim[2],
                            nrow = target_dim[1])

# ==================== CORRECTION FINALE ====================
# 🔁 Seulement la correction Est/Ouest
resized_matrix <- from03602180(resized_matrix)

# ❌ Pas d’inversion verticale ici !

# ==================== AFFICHAGE ====================
cat("→ Affichage du DEM final (Nord en haut)\n")
image(
  t(apply(resized_matrix, 2, rev)),
  col = terrain.colors(256),
  axes = FALSE,
  main = "DEM corrigé et redimensionné (tête à l’endroit)"
)

# ==================== SAUVEGARDE ====================
saveRDS(resized_matrix, file = output_rds)
cat("✅ DEM sauvegardé sous :", output_rds, "\n")
