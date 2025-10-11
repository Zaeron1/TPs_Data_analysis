generate_histograms <- function(corrected_mgsi, regions) {
  # Liste des noms des régions
  region_names <- c(
    "1: High-Mg region",
    "2: High-Al region",
    "3: Caloris Planitia",
    "4: Rachmaninoff",
    "5: High-Mg NVP",
    "6: Low-Mg NVP",
    "0: Unclassified"
  )
  corrected_mgsi <- get_layer_as_matrix(result_array, corrected_mgsi)
  regions <- get_layer_as_matrix(result_array, regions)
  
  for (i in seq_along(region_names)) {
    hist(
      corrected_mgsi[regions == (i %% 7)], # i %% 7 est utilisé pour gérer le cas du "Unclassified" avec index 0
      freq = FALSE,
      xlim = c(0, 1),
      col = "lightblue",
      xlab = "Rapport Mg/Si",
      ylab = "Densité",
      main = region_names[i]
    )
  }
}
