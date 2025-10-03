plot_max_elevation <- function(layer_index){
  # Vérifie si l'index de couche est dans les limites
  if (layer_index < 1 || layer_index > dim(result_array)[3]) {
    stop("Layer index out of bounds.")
  }

  # Extrait la couche spécifiée sous forme de matrice
  layer_matrix <- get_layer_as_matrix(result_array, layer_index)

  # Trouver la position du maximum dans la matrice
  position_max <- which(layer_matrix == max(layer_matrix), arr.ind = TRUE)
  ligne_max <- position_max[1, 1]
  colonne_max <- position_max[1, 2]

  # Crée un heatmap avec plotly
  plot_heatmap <- plot_ly(
    z = layer_matrix,
    type = "heatmap",
    colors = colorRamp(c("blue", "cyan", "yellowgreen", "yellow", "orange", "orangered", "darkred")),
    showscale = TRUE,
    colorbar = list(title = "Élévation (m)")
  ) %>%
    layout(
      title = "Position des coupes sur le DEM",
      xaxis = list(
        title = "Longitude",
        tickvals = seq(0, ncol(layer_matrix), by = ncol(layer_matrix) / 12),
        ticktext = paste0(seq(-180, 180, length.out = 13), "°")
      ),
      yaxis = list(
        title = "Latitude",
        tickvals = seq(0, nrow(layer_matrix), by = nrow(layer_matrix) / 12),
        ticktext = paste0(seq(90, -90, length.out = 13), "°"), autorange = "reversed"
      )
    )

  # Ajoute des lignes pour les coupes
  plot_heatmap <- plot_heatmap %>%
    add_trace(
      x = c(0, ncol(layer_matrix)), y = ligne_max,
      type = 'scatter', mode = 'lines', line = list(color = 'blue', width = 2), name = 'Coupe 2'
    ) %>%
    add_trace(
      x = colonne_max, y = c(0, nrow(layer_matrix)),
      type = 'scatter', mode = 'lines', line = list(color = 'red', width = 2), name = 'Coupe 1'
    )

  # Afficher la position et l'élévation maximales
  cat("La valeur maximale est à la ligne", ligne_max, "et la colonne", colonne_max, ".\n")
  cat("L'élévation du point le plus haut est de", max(layer_matrix), "m.\n")

  # Extraire la ligne et la colonne maximales
  ligne_selectionnee2 <- layer_matrix[ligne_max, ]
  colonne_selectionnee2 <- layer_matrix[ , colonne_max]

  # Créer un plot pour la ligne (longitude)
  plot_colonne <- plot_ly(
    x = seq_along(colonne_selectionnee2),
    y = colonne_selectionnee2,
    type = 'scatter',
    mode = 'lines+markers',
    name = paste("Ligne", colonne_max)
  ) %>%
    layout(
      title = "Élévation le long de la coupe 1",
      xaxis = list(
        title = "Latitude",
        tickvals = seq(0, length(colonne_selectionnee2), by = length(colonne_selectionnee2) / 12), # Définit les ticks
        ticktext = paste0(seq(-90, 90, length.out = 13), "°")
      ),
      yaxis = list(title = "Élévation (m)")
    )

  # Créer un plot pour la colonne (latitude)
  plot_ligne <- plot_ly(
    x = seq_along(ligne_selectionnee2),
    y = ligne_selectionnee2,
    type = 'scatter',
    mode = 'lines+markers',
    name = paste("Colonne", ligne_max)
  ) %>%
    layout(
      title = "Élévation le long de la coupe 2",
      xaxis = list(
        title = "Longitude",
        tickvals = seq(0, length(ligne_selectionnee2), by = length(ligne_selectionnee2) / 12), # Définit les ticks
        ticktext = paste0(seq(-180, 180, length.out = 13), "°")
      ),
      yaxis = list(title = "Élévation (m)")
    )

  # Retourner les plots
  list(plot_heatmap = plot_heatmap, plot_ligne = plot_ligne, plot_colonne = plot_colonne)
}
