#' Plot hybrid score distribution for a fancy result
#'
#' S3 plot method for objects returned by [fancy()]. Draws a histogram of
#' `HybridScore` across all scored edges with a vertical line at the
#' threshold cutoff. The number of retained edges is annotated.
#'
#' @param x An object of class `"fancy"` as returned by [fancy()].
#' @param ... Additional arguments passed to [graphics::hist()].
#'
#' @return Invisibly returns `x`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' result <- fancy(t(fancy_tiny_clr), n_bootstrap = 2, cpus = 1)
#' plot(result)
#' }
#'
#' @importFrom graphics hist abline text par
plot.fancy <- function(x, ...) {
  scores <- x$all_edges$HybridScore
  method <- x$params$threshold_method
  value  <- x$params$threshold_value

  # Determine the actual cutoff used
  if (method == "quantile") {
    cutoff <- as.numeric(stats::quantile(scores, probs = value, na.rm = TRUE))
  } else if (method == "score") {
    cutoff <- value
  } else {
    # top_n: cutoff is the minimum HybridScore in the retained edges
    cutoff <- min(x$edges$HybridScore, na.rm = TRUE)
  }

  old_par <- graphics::par(mar = c(5, 4, 4, 2) + 0.1)
  on.exit(graphics::par(old_par), add = TRUE)

  graphics::hist(
    scores,
    breaks = 30,
    main = "Hybrid Score Distribution",
    xlab = "HybridScore",
    col = "grey80",
    border = "white",
    ...
  )

  graphics::abline(v = cutoff, col = "red", lwd = 2, lty = 2)

  n_kept <- nrow(x$edges)
  n_total <- nrow(x$all_edges)
  label <- paste0(n_kept, "/", n_total, " edges retained")
  graphics::text(cutoff, graphics::par("usr")[4] * 0.9, labels = label,
                 pos = 4, col = "red", cex = 0.8)

  invisible(x)
}

#' Plot elbow curve for k selection
#'
#' Takes the list returned by [find_optimal_k()] and plots the elbow curve
#' (k vs average MI). When `$best_k` is present, marks the elbow point with
#' a distinct symbol and label.
#'
#' @param k_results A list as returned by [find_optimal_k()], containing
#'   `$results` (data.frame with `k` and `avg_mi` columns) and optionally
#'   `$best_k`.
#'
#' @return Invisibly returns `k_results`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' k_res <- find_optimal_k(t(fancy_tiny_clr), k_values = c(3L, 5L), cpus = 1)
#' plot_k_elbow(k_res)
#' }
#'
#' @importFrom graphics legend plot points
plot_k_elbow <- function(k_results) {
  df <- k_results$results

  graphics::plot(
    df$k, df$avg_mi,
    type = "b",
    pch = 19,
    xlab = "k (neighbours)",
    ylab = "Average MI",
    main = "Elbow Curve for k Selection"
  )

  if (!is.null(k_results$best_k)) {
    best_row <- df[df$k == k_results$best_k, ]
    graphics::points(best_row$k, best_row$avg_mi, pch = 8, cex = 2,
                     col = "red", lwd = 2)
    graphics::text(best_row$k, best_row$avg_mi,
                   labels = paste0("k = ", best_row$k),
                   pos = 3, col = "red", cex = 0.9)
  }

  invisible(k_results)
}

#' Plot the co-abundance network
#'
#' Draws a network graph from a `"fancy"` result using igraph. Nodes are
#' coloured by phylum (via [phyla_palette()]) and labelled with genus names.
#' Edge width is proportional to `HybridScore`. When `community = TRUE`,
#' modularity-based community detection is applied and clusters are shown
#' as shaded convex hulls around groups of co-abundant MAGs.
#'
#' @param fancy_result An object of class `"fancy"` as returned by
#'   [fancy()].
#' @param taxonomy A data.frame with MAG IDs as row names and at least
#'   `Phyla` and `Genus` columns.
#' @param top_n Integer or `NULL`; if set, only the top N edges by
#'   `HybridScore` are plotted (default `NULL` = all thresholded edges).
#' @param layout Character; igraph layout algorithm name (default
#'   `"fruchterman.reingold"`).
#' @param community Logical; if `TRUE`, run modularity-based community
#'   detection and draw shaded convex hulls around clusters (default
#'   `FALSE`).
#'
#' @return Invisibly returns the igraph graph object.
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' data(fancy_tiny_taxonomy)
#' result <- fancy(t(fancy_tiny_clr), n_bootstrap = 2, cpus = 1)
#' plot_network(result, fancy_tiny_taxonomy)
#' }
plot_network <- function(
    fancy_result,
    taxonomy,
    top_n = NULL,
    layout = "fruchterman.reingold",
    community = FALSE
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for plot_network(). ",
         "Install it with install.packages(\"igraph\").",
         call. = FALSE)
  }

  edge_df <- fancy_result$edges

  if (!is.null(top_n)) {
    top_n <- as.integer(top_n)
    ord <- order(edge_df$HybridScore, decreasing = TRUE)
    edge_df <- edge_df[ord[seq_len(min(top_n, nrow(edge_df)))], ]
  }

  g <- igraph::graph_from_data_frame(
    edge_df[, c("source", "target", "HybridScore")],
    directed = FALSE
  )

  # Node colours from phyla
  node_ids <- igraph::V(g)$name
  node_phyla <- ifelse(
    node_ids %in% rownames(taxonomy),
    taxonomy[node_ids[node_ids %in% rownames(taxonomy)], "Phyla"][
      match(node_ids, rownames(taxonomy))
    ],
    "Unknown"
  )
  palette <- phyla_palette(node_phyla)
  node_colours <- palette[node_phyla]

  # Node labels from Genus
  node_labels <- ifelse(
    node_ids %in% rownames(taxonomy),
    taxonomy[node_ids[node_ids %in% rownames(taxonomy)], "Genus"][
      match(node_ids, rownames(taxonomy))
    ],
    node_ids
  )

  # Edge widths scaled to [0.5, 4]
  scores <- igraph::E(g)$HybridScore
  if (max(scores) > min(scores)) {
    edge_widths <- 0.5 + 3.5 * (scores - min(scores)) /
      (max(scores) - min(scores))
  } else {
    edge_widths <- rep(2, length(scores))
  }

  # Layout
  layout_fun <- switch(
    layout,
    fruchterman.reingold = igraph::layout_with_fr,
    kamada.kawai         = igraph::layout_with_kk,
    circle               = igraph::layout_in_circle,
    igraph::layout_with_fr
  )
  # set.seed(42)
  if (identical(layout_fun, igraph::layout_with_fr)) {
    coords <- layout_fun(g, niter = 1000,
                         area = igraph::vcount(g)^2.8)
  } else {
    coords <- layout_fun(g)
  }

  if (isTRUE(community)) {
    comm <- igraph::cluster_fast_greedy(g)
    n_comm <- length(comm)

    # Push communities apart to reduce hull overlap
    mem <- igraph::membership(comm)
    for (ci in seq_len(n_comm)) {
      idx <- which(mem == ci)
      if (length(idx) > 1) {
        cx <- mean(coords[idx, 1])
        cy <- mean(coords[idx, 2])
        coords[idx, 1] <- cx + (coords[idx, 1] - cx) * 0.75
        coords[idx, 2] <- cy + (coords[idx, 2] - cy) * 0.75
      }
    }

    # Generate translucent hull colours
    hull_pal <- grDevices::hcl.colors(n_comm, palette = "Set 2", alpha = 0.15)

    # plot.communities dispatches via S3 on the communities object
    graphics::plot(
      comm, g,
      layout           = coords,
      vertex.color     = node_colours,
      vertex.label     = node_labels,
      vertex.label.cex = 0.55,
      vertex.label.color = "black",
      vertex.size      = 8,
      edge.width       = edge_widths,
      edge.color       = grDevices::adjustcolor("grey50", alpha.f = 0.6),
      col              = hull_pal,
      mark.expand      = 12,
      main             = "Fancy Co-Abundance Network"
    )
    graphics::legend("bottomleft", bty = "n", cex = 0.7,
                     legend = paste(n_comm, "communities detected"))
  } else {
    igraph::plot.igraph(
      g,
      layout      = coords,
      vertex.color = node_colours,
      vertex.label = node_labels,
      vertex.label.cex = 0.55,
      vertex.label.color = "black",
      vertex.size = 8,
      edge.width  = edge_widths,
      edge.color  = grDevices::adjustcolor("grey50", alpha.f = 0.6),
      main        = "Fancy Co-Abundance Network"
    )
  }

  # Phyla legend
  unique_phyla <- sort(unique(node_phyla))
  unique_phyla <- unique_phyla[unique_phyla != "Unknown"]
  leg_cols <- palette[unique_phyla]
  graphics::legend("topright", legend = unique_phyla, fill = leg_cols,
                   cex = 1.1, bty = "n", title = "Phyla", ncol = 1)

  invisible(g)
}
