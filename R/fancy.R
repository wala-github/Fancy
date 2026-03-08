#' Run the full Fancy hybrid network inference pipeline
#'
#' Convenience wrapper that chains preprocessing, bootstrap network inference,
#' edge scoring, and thresholding into a single call. Uses MI (via MRNET) and
#' distance correlation with bootstrap resampling to produce a ranked,
#' thresholded edge list.
#'
#' @param data Numeric matrix or data.frame of CLR-transformed abundances with
#'   **samples as rows and MAGs as columns**. Pass `t(clr_normalize(counts))`
#'   or use the Samples x MAGs orientation directly.
#' @param n_bootstrap Integer; number of bootstrap iterations (default 100).
#' @param k Integer; kNN parameter for MI estimation (default 5).
#' @param find_k Logical; if `TRUE`, run [find_optimal_k()] to choose k
#'   automatically via the elbow method (default `FALSE`).
#' @param k_values Integer vector of k values to evaluate when
#'   `find_k = TRUE`.
#' @param cpus Integer; number of parallel workers (default 4).
#' @param use_mrnet Logical; if `TRUE` (default), apply MRNET to the MI
#'   matrix for sparse edge selection.
#' @param w1 Numeric weight for edge frequency in the hybrid score
#'   (default 0.3).
#' @param w2 Numeric weight for dCor stability in the hybrid score
#'   (default 0.7).
#' @param score_type Character; `"multiplicative"` (default) or
#'   `"additive"`. See [hybrid_score()] for details.
#' @param dcor_rescue_percentile Numeric between 0 and 1, or `NULL`
#'   (default). When set, edges with dCor stability above this
#'   percentile receive a minimum EdgeFrequency floor in
#'   multiplicative mode. See [hybrid_score()].
#' @param threshold_method Character; one of `"quantile"`, `"score"`, or
#'   `"top_n"` (default `"quantile"`).
#' @param threshold_value Numeric; meaning depends on `threshold_method`
#'   (default 0.7, i.e. top 30 percent of edges for quantile method).
#' @param save_path Optional directory path for saving bootstrap results.
#' @param verbose Logical; if `TRUE` (default), emit progress messages.
#'
#' @return A list with class `"fancy"` containing:
#'   \describe{
#'     \item{`edges`}{The final thresholded edge data.frame.}
#'     \item{`all_edges`}{All scored edges before thresholding (for
#'       re-thresholding).}
#'     \item{`all_edges_unfiltered`}{All 4950 scored pairs before the
#'       hard frequency/stability filter. Useful for inspecting edges
#'       that were dropped by the pre-filter.}
#'     \item{`k`}{The k value used.}
#'     \item{`n_bootstrap`}{Number of bootstrap iterations run.}
#'     \item{`params`}{List of all parameters for reproducibility.}
#'   }
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' result <- fancy(t(fancy_tiny_clr), n_bootstrap = 2, cpus = 1)
#' head(result$edges)
#' plot(result)
#' }
fancy <- function(
    data,
    n_bootstrap = 100L,
    k = 5L,
    find_k = FALSE,
    k_values = c(2L, 3L, 4L, 5L, 6L, 7L, 8L, 10L, 15L, 18L, 20L),
    cpus = 4L,
    use_mrnet = TRUE,
    w1 = 0.3,
    w2 = 0.7,
    score_type = "multiplicative",
    dcor_rescue_percentile = NULL,
    threshold_method = "quantile",
    threshold_value = 0.7,
    save_path = NULL,
    verbose = TRUE
) {
  data <- as.matrix(data)

  # --- Step 1: Optionally find optimal k ---
  if (find_k) {
    if (verbose) message("Finding optimal k via elbow method...")
    k_result <- find_optimal_k(data, k_values = k_values, cpus = cpus)
    if (!is.null(k_result$best_k)) {
      k <- k_result$best_k
      if (verbose) message("Selected k = ", k)
    } else {
      if (verbose) message("Could not determine elbow; using k = ", k)
    }
  }

  # --- Step 2: Bootstrap network inference ---
  if (verbose) message("Running ", n_bootstrap, " bootstrap iterations...")
  bootstrap_results <- bootstrap_networks(
    data,
    n_bootstrap = n_bootstrap,
    k = k,
    cpus = cpus,
    use_mrnet = use_mrnet,
    save_path = save_path
  )

  # --- Step 3: Compute edge frequencies ---
  if (verbose) message("Computing MRNET edge frequencies...")
  mrnet_freq <- compute_edge_frequencies(bootstrap_results, "MRNET")

  if (verbose) message("Computing dCor edge frequencies...")
  dcor_freq <- compute_edge_frequencies(bootstrap_results, "dCor")

  # --- Step 4: Compute hybrid scores ---
  if (verbose) message("Computing hybrid scores...")
  edge_metrics <- compute_hybrid_scores(mrnet_freq, dcor_freq,
                                        return_unfiltered = TRUE)

  # Extract unfiltered edges (all pairs before hard-filter)
  unfiltered_metrics <- attr(edge_metrics, "unfiltered")
  attr(edge_metrics, "unfiltered") <- NULL

  # --- Step 5: Apply hybrid scoring formula ---
  scored <- hybrid_score(edge_metrics, w1 = w1, w2 = w2,
                         score_type = score_type,
                         dcor_rescue_percentile = dcor_rescue_percentile)
  all_edges_unfiltered <- hybrid_score(unfiltered_metrics, w1 = w1, w2 = w2,
                                       score_type = score_type,
                                       dcor_rescue_percentile = dcor_rescue_percentile)

  # --- Step 6: Threshold edges ---
  if (verbose) message("Thresholding edges (method = '", threshold_method,
                       "', value = ", threshold_value, ")...")
  edges <- threshold_edges(scored, method = threshold_method,
                           value = threshold_value)

  if (verbose) message("Done. Retained ", nrow(edges), " of ", nrow(scored),
                       " edges.")

  structure(
    list(
      edges                = edges,
      all_edges            = scored,
      all_edges_unfiltered = all_edges_unfiltered,
      k                    = k,
      n_bootstrap          = n_bootstrap,
      params               = list(
        n_bootstrap      = n_bootstrap,
        k                = k,
        find_k           = find_k,
        cpus             = cpus,
        use_mrnet        = use_mrnet,
        w1               = w1,
        w2               = w2,
        score_type       = score_type,
        dcor_rescue_percentile = dcor_rescue_percentile,
        threshold_method = threshold_method,
        threshold_value  = threshold_value
      )
    ),
    class = "fancy"
  )
}
