calculate_and_plot_composition <- function() {
  matrice_mgsi <- get_layer_as_matrix(result_array, layer_index = 8)
  matrice_casi <- get_layer_as_matrix(result_array, layer_index = 2)
  matrice_fesi <- get_layer_as_matrix(result_array, layer_index = 6)
  matrice_alsi <- get_layer_as_matrix(result_array, layer_index = 1)
  matrice_ssi <- get_layer_as_matrix(result_array, layer_index = 9)

  # Liste des matrices
  matrices <- list(
    mgsi = matrice_mgsi,
    casi = matrice_casi,
    fesi = matrice_fesi,
    alsi = matrice_alsi,
    ssi = matrice_ssi
  )

 
  mean_values <- sapply(matrices, function(matrix) mean(matrix, na.rm = TRUE))  # Calcul des moyennes

  # Fonction pour calculer la composition
  calculate_composition_from_means <- function(means) {
    total_moles <- 1 + sum(means)
    percentages <- c(1, means) / total_moles * 100
    names(percentages) <- c("Si", "Mg", "Ca", "Al", "Fe", "S")
    return(percentages)
  }

  composition_means <- calculate_composition_from_means(mean_values)# Composition
  print(composition_means)

  # Création du diagramme circulaire avec Plotly
  plot_ly(
    labels = names(composition_means),
    values = composition_means,
    type = 'pie',
    textinfo = 'label+percent',
    insidetextorientation = 'radial'
  ) %>%
    layout(
      title = "Composition de la surface de Mercure"
    )
}
