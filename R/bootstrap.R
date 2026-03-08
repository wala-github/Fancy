#' Initialise snowfall cluster
#'
#' Sets up a \pkg{snowfall} parallel cluster and loads required packages
#' on every worker.
#'
#' @param cpus Integer; number of workers.
#' @param slave_outfile Character path for worker log output, or `NULL`
#'   for the platform null device.
#'
#' @return Called for its side-effect (cluster initialisation).
#'
#' @importFrom snowfall sfInit sfStop sfLibrary
#' @keywords internal
.init_snowfall <- function(cpus, slave_outfile = NULL) {
  if (is.null(slave_outfile)) {
    slave_outfile <- if (.Platform$OS.type == "windows") "nul" else "/dev/null"
  }
  snowfall::sfStop()
  snowfall::sfInit(parallel = TRUE, cpus = cpus,
                   slaveOutfile = slave_outfile)
  snowfall::sfLibrary("minet",  character.only = TRUE)
  snowfall::sfLibrary("energy", character.only = TRUE)
  snowfall::sfLibrary("knnmi",  character.only = TRUE)
}

# -------------------------------------------------------------------

#' Single bootstrap iteration
#'
#' Resamples rows of `data` with replacement, computes the MI matrix,
#' optionally applies MRNET or extracts all MI edges, and computes dCor
#' edges.
#'
#' @param iter Integer; iteration number (unused but required by
#'   `sfClusterApplyLB`).
#' @param data Numeric matrix (samples x MAGs).
#' @param k Integer; kNN parameter for MI estimation.
#' @param use_mrnet Logical; if `TRUE`, apply MRNET; if `FALSE`, extract
#'   all pairwise MI edges.
#'
#' @return A named list with `MRNET` (data.frame) and `dCor`
#'   (data.frame).
#'
#' @keywords internal
.bootstrap_single <- function(iter, data, k, use_mrnet) {
  resampled <- data[sample(nrow(data), replace = TRUE), ]

  mi_matrix <- compute_mi_matrix(resampled, k = k)

  if (use_mrnet) {
    mi_edges <- compute_mrnet_network(mi_matrix)
  } else {
    mi_edges <- .extract_mi_edges(mi_matrix)
  }

  dcor_matrix <- compute_dcor_matrix(resampled)
  gc()

  # Extract dCor edges from upper triangle
  n <- nrow(dcor_matrix)
  source_vec <- character()
  target_vec <- character()
  weight_vec <- numeric()

  for (i in seq_len(n - 1L)) {
    js <- (i + 1L):n
    source_vec <- c(source_vec, rownames(dcor_matrix)[js])
    target_vec <- c(target_vec, rep(colnames(dcor_matrix)[i], length(js)))
    weight_vec <- c(weight_vec, dcor_matrix[js, i])
  }

  dcor_edges <- data.frame(
    source = source_vec,
    target = target_vec,
    weight = weight_vec,
    stringsAsFactors = FALSE
  )
  dcor_edges <- dcor_edges[dcor_edges$weight > 0, ]

  list(MRNET = mi_edges, dCor = dcor_edges)
}

# -------------------------------------------------------------------

#' Bootstrap network inference
#'
#' Runs `n_bootstrap` bootstrap iterations in parallel using
#' \pkg{snowfall}. Each iteration resamples rows, computes MI and dCor
#' networks, and returns edge lists. Results are collected into a single
#' named list.
#'
#' @param data A numeric matrix or data.frame with observations (samples)
#'   as rows and variables (MAGs) as columns.
#' @param n_bootstrap Integer; number of bootstrap iterations.
#' @param k Integer; kNN parameter for [compute_mi_matrix()].
#' @param cpus Integer; number of parallel workers. For large datasets
#'   (>1000 MAGs), consider using more CPUs (e.g. `cpus = 14L`).
#' @param use_mrnet Logical; if `TRUE` (default), applies MRNET to the MI
#'   matrix for sparse edge selection. If `FALSE`, uses all pairwise MI
#'   edges directly.
#' @param save_path Optional directory path. When non-`NULL`, each
#'   bootstrap result is saved as an RDS file for recovery.
#'
#' @return A named list with elements `MRNET_1`, `dCor_1`, `MRNET_2`,
#'   `dCor_2`, etc.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' res <- bootstrap_networks(t(fancy_tiny_clr), n_bootstrap = 2, k = 5L, cpus = 1)
#' }
#'
#' @importFrom snowfall sfInit sfStop sfExport sfClusterApplyLB sfLibrary sfCpus
bootstrap_networks <- function(
    data,
    n_bootstrap = 100L,
    k = 5L,
    cpus = 4L,
    use_mrnet = TRUE,
    save_path = NULL
) {
  data <- as.matrix(data)
  n_bootstrap <- as.integer(n_bootstrap)
  k <- as.integer(k)

  .init_snowfall(cpus)
  on.exit(snowfall::sfStop(), add = TRUE)

  # Assign internal functions to local variables for export to workers
  .bs_data <- data
  .bs_k <- k
  .bs_use_mrnet <- use_mrnet
  .bs_compute_mi_matrix <- compute_mi_matrix
  .bs_compute_mrnet_network <- compute_mrnet_network
  .bs_extract_mi_edges <- .extract_mi_edges
  .bs_compute_dcor_matrix <- compute_dcor_matrix

  snowfall::sfExport(
    ".bs_data", ".bs_k", ".bs_use_mrnet",
    ".bs_compute_mi_matrix", ".bs_compute_mrnet_network",
    ".bs_extract_mi_edges", ".bs_compute_dcor_matrix",
    local = TRUE
  )

  worker_fn <- function(iter) {
    resampled <- .bs_data[sample(nrow(.bs_data), replace = TRUE), ]

    mi_matrix <- .bs_compute_mi_matrix(resampled, k = .bs_k)

    if (.bs_use_mrnet) {
      mi_edges <- .bs_compute_mrnet_network(mi_matrix)
    } else {
      mi_edges <- .bs_extract_mi_edges(mi_matrix)
    }

    dcor_matrix <- .bs_compute_dcor_matrix(resampled)
    gc()

    n <- nrow(dcor_matrix)
    source_vec <- character()
    target_vec <- character()
    weight_vec <- numeric()

    for (i in seq_len(n - 1L)) {
      js <- (i + 1L):n
      source_vec <- c(source_vec, rownames(dcor_matrix)[js])
      target_vec <- c(target_vec, rep(colnames(dcor_matrix)[i], length(js)))
      weight_vec <- c(weight_vec, dcor_matrix[js, i])
    }

    dcor_edges <- data.frame(
      source = source_vec,
      target = target_vec,
      weight = weight_vec,
      stringsAsFactors = FALSE
    )
    dcor_edges <- dcor_edges[dcor_edges$weight > 0, ]

    list(MRNET = mi_edges, dCor = dcor_edges)
  }

  raw <- snowfall::sfClusterApplyLB(seq_len(n_bootstrap), worker_fn)

  # Flatten into named list
  all_results <- list()
  for (i in seq_along(raw)) {
    result <- raw[[i]]

    if (!is.null(save_path)) {
      saveRDS(result, file.path(save_path,
                                paste0("bootstrap_result_", i, ".rds")))
    }

    names(result) <- c(paste0("MRNET_", i), paste0("dCor_", i))
    all_results <- c(all_results, result)
  }

  message("All bootstrap iterations processed. Total iterations: ",
          length(raw))

  all_results
}
