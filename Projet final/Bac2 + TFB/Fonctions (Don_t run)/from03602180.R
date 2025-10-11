from03602180 <- function(mat) {
  middle <- dim(mat)[2] / 2
  mat1 <- mat[, 1:middle]
  mat2 <- mat[, (middle + 1):(2 * middle)]
  return(cbind(mat2, mat1))
}
