classify_matrix <- function(layer_index){
  if (layer_index < 1 || layer_index > dim(result_array)[3]) {
    stop("Layer index out of bounds.")
  }
  layer_matrix <- get_layer_as_matrix(result_array, layer_index)
  quantile <- quantile(layer_matrix, probs = seq(0, 1, by = 0.1))  # Calcul des quantile
  labels <- c("Lowest","Very Very Low", "Very Low", "Low", "Medium Low", "Medium", "Medium High", "High", "Very High", "Very Very High", "Highest")  # Étiquettes de classification
  classes <- cut(as.vector(layer_matrix), breaks = c(quantile, Inf), labels = labels, include.lowest = TRUE)  # Classification des valeurs
  class_matrix <- matrix(as.numeric(classes), nrow = nrow(layer_matrix), byrow = FALSE)  # Conversion en matrice numérique
  return(class_matrix)
}
