get_layer_as_matrix <- function(result_array, layer_index) {

  if (layer_index < 1 || layer_index > dim(result_array)[3]) { #check si l'index est inf. à 1 (index minimum) ou si il est sup. au nombre de couche de result_array
    stop("Layer index out of bounds.")
  }
  return(result_array[,,layer_index])   # Extrait et retourne la couche spécifiée sous forme de matrice
}
