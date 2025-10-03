analyse_regressions <- function(density_index, compare_indices) {
  resultats_reg <- data.frame(Couche = character(), R_squared = numeric(), stringsAsFactors = FALSE) # Initialisation du tableau de résultat
  matrice_density <- get_layer_as_matrix(result_array, density_index)

  # Define custom titles for specific indices
  custom_titles <- list(
    "1" = "Al/Si",
    "2" = "Ca/Si",
    "6" = "Fe/Si",
    "8" = "Mg/Si",
    "9" =  "SSi"
  )

  for (index in compare_indices) {
    matrice_chim <- get_layer_as_matrix(result_array, index)
    matrice_density_filtered <- matrice_density[matrice_density > 0]  # ERRORE 02.09.24 - matrice_chim al posto di matrice_density
    matrice_chim_filtered <- matrice_chim[matrice_chim > 0]

    # Vérifier si les matrices sont non vides
    if (length(matrice_density_filtered) == 0 || length(matrice_chim_filtered) == 0) {
      warning(paste("Pas de valeurs positives pour l'indice", index))
      r_squared <- NA
    } else {
      # Déterminer la taille de l'échantillon
      sample_size <- min(10000, length(matrice_density_filtered))
      if (sample_size < 10000) {
        warning(paste("Sample size réduit à", sample_size, "pour l'indice", index))
      }

      # Échantillonnage
      positions <- sample(1:length(matrice_density_filtered), sample_size)
      valeurs_matrice_dens <- matrice_density_filtered[positions]
      valeurs_matrice_chim <- matrice_chim_filtered[positions]

      model <- lm(valeurs_matrice_chim ~ valeurs_matrice_dens)# Calcul du modèle de régression

      r_squared <- summary(model)$r.squared # Calcul du coefficient de détermination (R^2)

       #écriture du titre du graphe en fonction de l'index utilisé
      if (index %in% names(custom_titles)) {
        custom_text <- custom_titles[[as.character(index)]]
        plot_title <- paste("Régression linéaire pour le rapport", custom_text, "et la densité")
        yaxis_title <- paste("Rapport", custom_text)
      } else {
        plot_title <- paste("Régression linéaire pour l'indice", index)
        yaxis_title <- paste("Rapport Indice", index)
      }

      # Création du graphique 
      plot <- plot_ly() %>%
        add_trace(x = ~valeurs_matrice_dens, y = ~valeurs_matrice_chim, type = 'scatter', mode = 'markers',
                  marker = list(color = 'blue'), name = 'Data Points') %>%
        add_trace(x = ~valeurs_matrice_dens, y = fitted(model), mode = 'lines', line = list(color = 'red', width = 2), name = 'Regression Line') %>%
        layout(title = plot_title,
               xaxis = list(title = "Density"),
               yaxis = list(title = yaxis_title),
               annotations = list(
                 list(
                   x = max(valeurs_matrice_dens), y = max(valeurs_matrice_chim),
                   xref = 'x', yref = 'y',
                   text = paste("R^2 =", round(r_squared, 4)),
                   showarrow = TRUE,
                   arrowhead = 2,
                   ax = 20, ay = -30,
                   font = list(color = "black", size = 20)
                 )
               ))

      print(plot)
    }
    resultats_reg <- rbind(resultats_reg, data.frame(Couche = paste("Indice", index), R_squared = r_squared))     # Stockage des résultats
  }

  return(resultats_reg)
}
