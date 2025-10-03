region <- function(index1, index2) {
  regions_matrix <- get_layer_as_matrix(result_array, index1)
  mg_si_matrix <- get_layer_as_matrix(result_array, index2)
  
  if (!all(dim(regions_matrix) == dim(mg_si_matrix))) {
    stop("Les matrices doivent avoir les mêmes dimensions.")
  }
  region_labels <- 0:6  # Initialiser un vecteur pour stocker la moyenne de Mg/Si par région
  mg_si_means <- numeric(length(region_labels))

  # Calculer la moyenne de Mg/Si pour chaque région
  for (region in region_labels) {
    region_indices <- which(regions_matrix == region, arr.ind = TRUE)# Trouver les indices des éléments appartenant à la région
    mg_si_values <- mg_si_matrix[region_indices]# Extraire les valeurs correspondantes de la matrice Mg/Si
    mg_si_means[region + 1] <- mean(mg_si_values, na.rm = TRUE)# Calculer la moyenne pour la région, en évitant les NA
  }

  # Créer un data frame pour le tracé
  df <- data.frame(
    Region = as.factor(region_labels),
    Mg_Si_Mean = mg_si_means
  )
  # Créer l'histogramme avec Plotly
  plot <- plot_ly(
    data = df,
    x = ~Region,
    y = ~Mg_Si_Mean,
    type = 'bar',
    name = 'Moyenne Mg/Si',
    marker = list(color = 'rgba(50, 171, 96, 0.6)')
  ) %>%
    layout(
      title = 'Moyenne de Mg/Si par région',
      xaxis = list(title = 'Région'),
      yaxis = list(title = 'Moyenne Mg/Si')
    )
  plot # Return the plot
}
