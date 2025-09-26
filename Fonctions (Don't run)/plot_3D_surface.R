plot_3D_surface <- function(layer_index, plot_title) {
  if (layer_index < 1 || layer_index > dim(result_array)[3]) {
    stop("Layer index out of bounds.")
  }
  layer_matrix <- get_layer_as_matrix(result_array, layer_index)
  reversed_matrix <- t(apply(layer_matrix, 1, rev))
  axx <- list(title = "Longitude")
  axy <- list(title = "Latitude")
  axz <- list(title = "Elevation")

  fig <- plot_ly(z = ~reversed_matrix) %>%
    add_surface(colorscale = list(
      list(0, 'rgb(0, 0, 255)'),       
      list(0.5, 'rgb(255, 255, 255)'), 
      list(1, 'rgb(255, 0, 0)')        
    )) %>%
    layout(
      title = plot_title,
      scene = list(
        xaxis = axx,
        yaxis = axy,
        zaxis = axz,
        aspectmode = 'manual',
        aspectratio = list(x = 2, y = 1, z = 0.07)
      )
    )
  fig
}
