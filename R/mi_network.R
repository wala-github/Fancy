#' Compute pairwise mutual information matrix
#'
#' Estimates pairwise mutual information for all variable pairs using the
#' k-nearest-neighbour estimator from \pkg{knnmi}.
#'
#' @param data A numeric matrix or data.frame with observations (samples) as
#'   rows and variables (MAGs) as columns.
#' @param k Integer; number of neighbours for the kNN MI estimator.
#'
#' @return A symmetric numeric matrix of MI values with column names set.
#'
#' @importFrom knnmi mutual_inf_cc
#' @keywords internal
compute_mi_matrix <- function(data, k = 5L) {
  data <- as.matrix(data)
  n_vars <- ncol(data)
  mi_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)

  for (i in seq_len(n_vars - 1L)) {
    for (j in (i + 1L):n_vars) {
      val <- knnmi::mutual_inf_cc(data[, i], data[, j], k = k)
      mi_matrix[i, j] <- val
      mi_matrix[j, i] <- val
    }
  }

  nms <- colnames(data)
  if (is.null(nms)) nms <- as.character(seq_len(n_vars))
  rownames(mi_matrix) <- nms
  colnames(mi_matrix) <- nms

  mi_matrix
}

# -------------------------------------------------------------------

#' Apply MRNET to a mutual information matrix
#'
#' Runs the MRNET algorithm via [minet::mrnet()] and extracts non-zero
#' edges as a data.frame.
#'
#' @param mi_matrix A symmetric MI matrix (from [compute_mi_matrix()]).
#'
#' @return A data.frame with columns `source`, `target`, `weight`.
#'
#' @importFrom minet mrnet
#' @keywords internal
compute_mrnet_network <- function(mi_matrix) {
  network <- minet::mrnet(mi_matrix)
  edges <- which(network != 0, arr.ind = TRUE)
  data.frame(
    source = rownames(network)[edges[, 1]],
    target = colnames(network)[edges[, 2]],
    weight = network[edges],
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------

#' Extract all pairwise MI edges from the upper triangle
#'
#' Returns every non-zero MI pair from the upper triangle of the MI matrix
#' as an edge list. Used when `use_mrnet = FALSE` in [bootstrap_networks()]
#' as a non-sparse alternative to MRNET.
#'
#' @param mi_matrix A symmetric MI matrix (from [compute_mi_matrix()]).
#'
#' @return A data.frame with columns `source`, `target`, `weight`.
#'
#' @keywords internal
.extract_mi_edges <- function(mi_matrix) {
  idx <- which(upper.tri(mi_matrix) & mi_matrix != 0, arr.ind = TRUE)
  data.frame(
    source = rownames(mi_matrix)[idx[, 1]],
    target = colnames(mi_matrix)[idx[, 2]],
    weight = mi_matrix[idx],
    stringsAsFactors = FALSE
  )
}

# -------------------------------------------------------------------

#' Find optimal k for the kNN MI estimator
#'
#' Tests multiple k values by computing the full MI matrix for each k and
#' recording the average MI across all variable pairs. The resulting
#' elbow curve helps choose k: the point where MI stops improving steeply
#' (largest negative second derivative) is selected automatically when
#' `auto_select = TRUE`.
#'
#' Computation is parallelised across k values using \pkg{snowfall}.
#'
#' @param data A numeric matrix or data.frame with observations (samples)
#'   as rows and variables (MAGs) as columns.
#' @param k_values Integer vector of k values to evaluate.
#' @param cpus Integer; number of parallel workers.
#' @param auto_select Logical; if `TRUE` (default), returns the elbow
#'   point as `$best_k`.
#'
#' @return A list with:
#'   \describe{
#'     \item{`results`}{A `data.frame(k, avg_mi)` for plotting.}
#'     \item{`best_k`}{(when `auto_select = TRUE`) the k at the elbow
#'       point, determined by the most negative second derivative.}
#'   }
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' res <- find_optimal_k(t(fancy_tiny_clr), k_values = c(3L, 5L), cpus = 1)
#' plot(res$results$k, res$results$avg_mi, type = "b")
#' res$best_k
#' }
#'
#' @importFrom snowfall sfInit sfStop sfLibrary sfExport sfClusterApplyLB sfCpus
#' @importFrom knnmi mutual_inf_cc
find_optimal_k <- function(
    data,
    k_values = c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 10L, 15L, 18L, 20L),
    cpus = 4L,
    auto_select = TRUE
) {
  data <- as.matrix(data)
  k_values <- as.integer(k_values)

  snowfall::sfStop()
  snowfall::sfInit(parallel = TRUE, cpus = cpus)
  on.exit(snowfall::sfStop(), add = TRUE)
  snowfall::sfLibrary("knnmi", character.only = TRUE)

  .fok_data <- data
  .fok_k_values <- k_values
  snowfall::sfExport(".fok_data", ".fok_k_values", local = TRUE)

  worker_fn <- function(iter) {
    k <- .fok_k_values[iter]
    n_vars <- ncol(.fok_data)
    mi_matrix <- matrix(0, nrow = n_vars, ncol = n_vars)

    for (i in seq_len(n_vars - 1L)) {
      for (j in (i + 1L):n_vars) {
        val <- knnmi::mutual_inf_cc(.fok_data[, i], .fok_data[, j], k = k)
        mi_matrix[i, j] <- val
        mi_matrix[j, i] <- val
      }
    }

    avg_mi <- mean(mi_matrix[upper.tri(mi_matrix)])
    data.frame(k = k, avg_mi = avg_mi)
  }

  raw <- snowfall::sfClusterApplyLB(seq_along(k_values), worker_fn)
  results <- do.call(rbind, raw)
  results <- results[order(results$k), ]
  rownames(results) <- NULL

  out <- list(results = results)

  if (auto_select && nrow(results) >= 3L) {
    delta2 <- diff(diff(results$avg_mi))
    # The elbow is the k where the second derivative is most negative
    elbow_idx <- which.min(delta2) + 1L
    out$best_k <- results$k[elbow_idx]
  }

  out
}
