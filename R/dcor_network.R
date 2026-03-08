#' Compute pairwise distance correlation matrix
#'
#' Computes the distance correlation for every pair of variables using
#' [energy::dcor()]. The result is a symmetric matrix with zeros on the
#' diagonal.
#'
#' @param data A numeric matrix or data.frame with observations (samples)
#'   as rows and variables (MAGs) as columns.
#'
#' @return A symmetric numeric matrix of dCor values. Diagonal is 0.
#'
#' @importFrom energy dcor
#' @keywords internal
compute_dcor_matrix <- function(data) {
  data <- as.matrix(data)
  n_vars <- ncol(data)
  dcor_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)

  nms <- colnames(data)
  if (is.null(nms)) nms <- as.character(seq_len(n_vars))
  rownames(dcor_matrix) <- nms
  colnames(dcor_matrix) <- nms

  for (i in seq_len(n_vars - 1L)) {
    for (j in (i + 1L):n_vars) {
      val <- energy::dcor(data[, i], data[, j])
      dcor_matrix[i, j] <- val
      dcor_matrix[j, i] <- val
    }
  }

  dcor_matrix
}
