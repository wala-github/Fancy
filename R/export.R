#' Assign colours to phyla
#'
#' Returns a named character vector mapping each unique phylum to a
#' colourblind-friendly hex colour. Colours are assigned in alphabetical
#' order of phylum names so the mapping is deterministic.
#'
#' @param phyla Character vector of phylum names (duplicates are OK).
#'
#' @return A named character vector of hex colours, one per unique phylum
#'   (sorted alphabetically).
#'
#' @export
#'
#' @examples
#' phyla_palette(c("Bacillota", "Pseudomonadota", "Bacillota"))
phyla_palette <- function(phyla) {
  unique_phyla <- sort(unique(phyla))
  n <- length(unique_phyla)

  # Curated 12-colour palette (colourblind-friendly, based on Brewer Set3/Paired)
  base_colours <- c(
    "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
    "#FF7F00", "#A65628", "#F781BF", "#999999",
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3"
  )

  if (n <= length(base_colours)) {
    cols <- base_colours[seq_len(n)]
  } else {
    cols <- grDevices::hcl.colors(n, palette = "Dynamic")
  }

  stats::setNames(cols, unique_phyla)
}

#' Export a fancy result for Cytoscape
#'
#' Writes node and edge tables as tab-separated files suitable for import
#' into Cytoscape. The node table includes taxonomy columns and a hex
#' colour derived from [phyla_palette()].
#'
#' @param fancy_result An object of class `"fancy"` as returned by
#'   [fancy()].
#' @param taxonomy A data.frame with MAG IDs as row names and at least a
#'   `Phyla` column. Additional taxonomy ranks (Domain, Class, Order,
#'   Family, Genus, Species) are included when present.
#' @param file_prefix Character; path prefix for the output files
#'   (default `"fancy_network"`). Two files are written:
#'   `{prefix}_edges.tsv` and `{prefix}_nodes.tsv`.
#' @param edges Character; which edge set to export. `"thresholded"`
#'   (default) uses `fancy_result$edges`; `"all"` uses
#'   `fancy_result$all_edges`.
#'
#' @return Invisibly, a list with elements `$nodes` and `$edges`
#'   (data.frames matching the written files).
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(fancy_tiny_clr)
#' data(fancy_tiny_taxonomy)
#' result <- fancy(t(fancy_tiny_clr), n_bootstrap = 2, cpus = 1)
#' out <- export_cytoscape(result, fancy_tiny_taxonomy,
#'                         file_prefix = file.path(tempdir(), "fancy_net"))
#' head(out$edges)
#' head(out$nodes)
#' }
export_cytoscape <- function(
    fancy_result,
    taxonomy,
    file_prefix = "fancy_network",
    edges = "thresholded"
) {
  edges <- match.arg(edges, c("thresholded", "all"))

  edge_df <- if (edges == "thresholded") {
    fancy_result$edges
  } else {
    fancy_result$all_edges
  }

  # Build edge table: source, target, HybridScore
  edge_out <- data.frame(
    source     = edge_df$source,
    target     = edge_df$target,
    HybridScore = edge_df$HybridScore,
    stringsAsFactors = FALSE
  )

  # Identify MAGs present in edges
  node_ids <- sort(unique(c(edge_out$source, edge_out$target)))

  # Standard taxonomy column order
  tax_ranks <- c("Domain", "Phyla", "Class", "Order", "Family", "Genus",
                 "Species")
  available_ranks <- intersect(tax_ranks, colnames(taxonomy))

  # Build node table
  node_out <- data.frame(id = node_ids, stringsAsFactors = FALSE)

  # Add taxonomy columns for nodes that have taxonomy info
  for (rank in available_ranks) {
    node_out[[rank]] <- ifelse(
      node_ids %in% rownames(taxonomy),
      taxonomy[node_ids[node_ids %in% rownames(taxonomy)], rank][
        match(node_ids, rownames(taxonomy))
      ],
      NA_character_
    )
  }

  # Add colour based on Phyla
  if ("Phyla" %in% available_ranks) {
    all_phyla <- taxonomy[node_ids[node_ids %in% rownames(taxonomy)], "Phyla"]
    palette <- phyla_palette(all_phyla)
    node_out$Color <- palette[node_out$Phyla]
  }

  # Write files
  utils::write.table(edge_out, file = paste0(file_prefix, "_edges.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(node_out, file = paste0(file_prefix, "_nodes.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)

  invisible(list(nodes = node_out, edges = edge_out))
}
