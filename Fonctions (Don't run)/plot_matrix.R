plot_matrix <- function(layer_index_or_name, main_title, legend_name, color_palette) {
  # check si l'argument est une matrice dans l'environment ou un index à extraire de result_array
  if (is.matrix(layer_index_or_name)) {
    layer_matrix <- layer_index_or_name
  } else if (is.numeric(layer_index_or_name)) {
    if (layer_index_or_name < 1 || layer_index_or_name > dim(result_array)[3]) {
      stop("Layer index out of bounds.")
    }
    layer_matrix <- get_layer_as_matrix(result_array, layer_index_or_name)
  } else {
    stop("The input must be either a matrix or a numeric index of a layer.")
  }

  # Create a plotly heatmap
  plot_ly(
    z = layer_matrix,
    type = "heatmap",
    colors = colorRamp(c("blue", "cyan", "yellowgreen", "yellow", "orange", "orangered", "darkred")),
    showscale = TRUE,
    colorbar = list(title = legend_name)
  ) %>%
    layout(
      title = main_title,
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
}
