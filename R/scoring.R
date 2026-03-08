#' Compute edge frequencies from bootstrap results
#'
#' Aggregates edge lists across bootstrap iterations for a given method.
#' Canonicalises edge direction via [.remove_directionality()] and
#' computes frequency, mean weight, and standard deviation per edge pair.
#'
#' @param bootstrap_results Named list as returned by
#'   [bootstrap_networks()].
#' @param method Character; `"MRNET"` or `"dCor"`.
#'
#' @return A data.frame with columns `source`, `target`, `freq`, `mean`,
#'   `sd`.
#'
#' @importFrom dplyr group_by summarise n arrange desc
#' @importFrom tidyr separate
#' @importFrom rlang .data
#' @importFrom stats sd
#' @keywords internal
compute_edge_frequencies <- function(bootstrap_results, method) {
  valid <- bootstrap_results[grepl(paste0("^", method, "_"),
                                   names(bootstrap_results))]
  n_bootstrap <- length(valid)

  all_edges <- do.call(rbind, valid)
  all_edges <- .remove_directionality(all_edges)

  all_edges$pair <- paste0(all_edges$source, ":::", all_edges$target)

  if (method == "MRNET") {
    pair_summary <- dplyr::group_by(all_edges, .data$pair)
    pair_summary <- dplyr::summarise(
      pair_summary,
      freq  = dplyr::n() / n_bootstrap / 2,
      mean  = mean(.data$weight),
      sd    = stats::sd(.data$weight),
      .groups = "drop"
    )
  } else {
    pair_summary <- dplyr::group_by(all_edges, .data$pair)
    pair_summary <- dplyr::summarise(
      pair_summary,
      freq  = 1,
      mean  = mean(.data$weight),
      sd    = stats::sd(.data$weight),
      .groups = "drop"
    )
  }

  pair_summary <- tidyr::separate(
    pair_summary, "pair",
    into = c("source", "target"),
    sep  = ":::"
  )

  pair_summary
}

# -------------------------------------------------------------------

#' Compute hybrid scores from MRNET/MI and dCor frequencies
#'
#' Full-joins MRNET (or raw MI) and dCor edge frequency tables, computes
#' dCor stability (mean / sd), robust-scales stability to \[0, 1\] using
#' the 20th--99th percentile range, and filters edges based on frequency
#' and stability thresholds.
#'
#' @param mrnet_freq data.frame from [compute_edge_frequencies()] with
#'   method `"MRNET"`.
#' @param dcor_freq data.frame from [compute_edge_frequencies()] with
#'   method `"dCor"`.
#' @param return_unfiltered Logical; if `TRUE`, the unfiltered edge
#'   table (before hard-filter) is attached as an `"unfiltered"`
#'   attribute on the returned data.frame (default `FALSE`).
#'
#' @return A data.frame with columns `source`, `target`,
#'   `EdgeFrequency`, `mean.dcor`, `sd.dcor`, `Stability.dcor`,
#'   `Stability.dcor.scaled`.
#'
#' @importFrom dplyr full_join mutate if_else filter
#' @importFrom stats quantile
#' @keywords internal
compute_hybrid_scores <- function(mrnet_freq, dcor_freq,
                                  return_unfiltered = FALSE) {
  hybrid <- dplyr::full_join(
    mrnet_freq, dcor_freq,
    by = c("source", "target"),
    relationship = "many-to-many"
  )

  hybrid <- dplyr::mutate(
    hybrid,
    EdgeFrequency = dplyr::if_else(is.na(.data$freq.x), 0, .data$freq.x),
    mean.dcor     = dplyr::if_else(is.na(.data$mean.y), 0, .data$mean.y),
    sd.dcor       = dplyr::if_else(is.na(.data$sd.y),   1, .data$sd.y),
    Stability.dcor = dplyr::if_else(
      is.infinite(.data$mean.dcor / .data$sd.dcor),
      .data$mean.dcor / 0.001,
      .data$mean.dcor / .data$sd.dcor
    )
  )

  # Robust scaling to [0, 1]
  q20 <- as.numeric(stats::quantile(hybrid$Stability.dcor, 0.20))
  q99 <- as.numeric(stats::quantile(hybrid$Stability.dcor, 0.99))

  hybrid <- dplyr::mutate(
    hybrid,
    Stability.dcor.scaled = pmin(1, pmax(0,
      (.data$Stability.dcor - q20) / (q99 - q20)))
  )

  # Drop intermediate join columns before filtering
  hybrid <- dplyr::mutate(
    hybrid,
    freq.x = NULL, freq.y = NULL,
    mean.x = NULL, sd.x   = NULL,
    mean.y = NULL, sd.y   = NULL
  )

  # Save unfiltered snapshot
  if (return_unfiltered) {
    unfiltered <- hybrid
  }

  # Filtering
  dcor_cut <- as.numeric(
    stats::quantile(hybrid$Stability.dcor, 0.8)
  )

  hybrid <- dplyr::filter(
    hybrid,
    .data$EdgeFrequency > 0.6 |
      (.data$EdgeFrequency >= 0.3 &
       .data$EdgeFrequency <= 0.6 &
       .data$Stability.dcor > dcor_cut)
  )

  if (return_unfiltered) {
    attr(hybrid, "unfiltered") <- unfiltered
  }

  hybrid
}

# -------------------------------------------------------------------

#' Compute hybrid edge score
#'
#' Scores edges by combining MRNET edge frequency with dCor stability.
#' Two scoring modes are available:
#' \describe{
#'   \item{`"multiplicative"` (default)}{
#'     \deqn{HybridScore = EdgeFrequency^{w1} \times
#'       Stability.dcor.scaled^{w2}}{HybridScore = EdgeFrequency^w1 *
#'       Stability.dcor.scaled^w2}
#'     Edges need support from both components; a near-zero value in
#'     either component pulls the score toward zero.}
#'   \item{`"additive"`}{
#'     \deqn{HybridScore = w1 \times EdgeFrequency + w2 \times
#'       Stability.dcor.scaled}{HybridScore = w1 * EdgeFrequency +
#'       w2 * Stability.dcor.scaled}
#'     Allows edges with strong dCor but weak MRNET frequency (or
#'     vice versa) to retain a meaningful score.}
#' }
#'
#' When `dcor_rescue_percentile` is set (multiplicative mode only),
#' edges whose `Stability.dcor.scaled` exceeds the given percentile
#' have their `EdgeFrequency` raised to a minimum floor (the 30th
#' percentile of non-zero EdgeFrequency values). This prevents
#' near-zero MRNET frequency from suppressing edges that dCor finds
#' reliably.
#'
#' @param edge_metrics A data.frame containing at least `EdgeFrequency`
#'   and `Stability.dcor.scaled` columns (as returned by
#'   [compute_hybrid_scores()]).
#' @param w1 Numeric weight for edge frequency (default 0.3).
#' @param w2 Numeric weight for dCor stability (default 0.7).
#' @param score_type Character; `"multiplicative"` (default) or
#'   `"additive"`.
#' @param dcor_rescue_percentile Numeric between 0 and 1, or `NULL`
#'   (default). When set, edges with `Stability.dcor.scaled` above
#'   this percentile receive a minimum `EdgeFrequency` floor.
#'   Only applies to multiplicative scoring.
#'
#' @return The input data.frame with an added `HybridScore` column.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' result <- fancy(t(fancy_tiny_clr), n_bootstrap = 2, cpus = 1)
#' scored <- hybrid_score(result$all_edges, w1 = 0.3, w2 = 0.7)
#' head(scored)
#' }
#'
#' @importFrom stats quantile
hybrid_score <- function(
    edge_metrics,
    w1 = 0.3,
    w2 = 0.7,
    score_type = "multiplicative",
    dcor_rescue_percentile = NULL
) {
  score_type <- match.arg(score_type, c("multiplicative", "additive"))

  ef   <- edge_metrics$EdgeFrequency
  stab <- edge_metrics$Stability.dcor.scaled

  if (score_type == "additive") {
    edge_metrics$HybridScore <- w1 * ef + w2 * stab
  } else {
    # dCor rescue: floor EdgeFrequency for high-stability edges
    if (!is.null(dcor_rescue_percentile)) {
      stab_thresh <- as.numeric(
        stats::quantile(stab, dcor_rescue_percentile, na.rm = TRUE)
      )
      nonzero_ef <- ef[ef > 0]
      ef_floor <- if (length(nonzero_ef) > 0) {
        as.numeric(stats::quantile(nonzero_ef, 0.3, na.rm = TRUE))
      } else {
        0.3
      }
      rescue <- stab >= stab_thresh & ef < ef_floor
      ef[rescue] <- ef_floor
    }

    edge_metrics$HybridScore <- (ef ^ w1) * (stab ^ w2)
  }

  edge_metrics
}

# -------------------------------------------------------------------

#' Threshold scored edges
#'
#' Filters a scored edge list using one of three methods:
#' \describe{
#'   \item{`"quantile"`}{Retain edges above the given quantile of
#'     `HybridScore`. Default `value = 0.7` retains the top 30
#'     percent.}
#'   \item{`"score"`}{Retain edges with `HybridScore >= value`.}
#'   \item{`"top_n"`}{Retain the top `value` edges by `HybridScore`.}
#' }
#'
#' @param scored_edges A data.frame with a `HybridScore` column.
#' @param method Character; one of `"quantile"`, `"score"`, `"top_n"`.
#' @param value Numeric; meaning depends on `method`.
#'
#' @return A data.frame of retained edges, sorted by descending
#'   `HybridScore`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' result <- fancy(t(fancy_tiny_clr), n_bootstrap = 2, cpus = 1)
#' top <- threshold_edges(result$all_edges, method = "top_n", value = 50)
#' nrow(top)
#' }
#'
#' @importFrom stats quantile
threshold_edges <- function(
    scored_edges,
    method = "quantile",
    value = 0.7
) {
  method <- match.arg(method, c("quantile", "score", "top_n"))

  if (method == "quantile") {
    cutoff <- as.numeric(
      stats::quantile(scored_edges$HybridScore, probs = value, na.rm = TRUE)
    )
    out <- scored_edges[scored_edges$HybridScore >= cutoff, ]
  } else if (method == "score") {
    out <- scored_edges[scored_edges$HybridScore >= value, ]
  } else {
    # top_n
    value <- as.integer(value)
    ord <- order(scored_edges$HybridScore, decreasing = TRUE)
    out <- scored_edges[ord[seq_len(min(value, nrow(scored_edges)))], ]
  }

  out <- out[order(out$HybridScore, decreasing = TRUE), ]
  rownames(out) <- NULL
  out
}
