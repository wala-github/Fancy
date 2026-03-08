#' Filter MAGs by coverage and prevalence
#'
#' Applies two filtering steps from the preprocessing pipeline:
#' 1. **Coverage masking:** Binarises the coverage table at `min_coverage`
#'    and multiplies element-wise with the count table, zeroing out counts
#'    for MAG-sample pairs with insufficient genome coverage.
#' 2. **Prevalence filter:** Retains only MAGs present (non-zero after
#'    masking) in at least `min_samples` samples.
#'
#' @note The original `02_Preprocessing.Rmd` referenced `count_table_cov`
#'   on line 111, but the variable was actually `count_table_cov_stats`
#'   (line 105).
#'   This function eliminates that confusion.
#'
#' @param count_table A data.frame of read counts with MAGs as rows and
#'   samples as columns.
#' @param coverage_table A data.frame of covered-fraction values with the
#'   same dimensions and names as `count_table`.
#' @param min_coverage Numeric; minimum covered fraction to retain a
#'   count.
#' @param min_samples Integer; minimum number of samples in which a MAG
#'   must be present after coverage masking.
#'
#' @return A data.frame of filtered counts (MAGs as rows, samples as
#'   columns).
#' @export
#' @examples
#' \donttest{
#' data(fancy_tiny_counts)
#' data(fancy_tiny_coverage)
#' filt <- filter_mags(fancy_tiny_counts, fancy_tiny_coverage)
#' dim(filt)
#' }
#'
#' @importFrom dplyr mutate across
filter_mags <- function(
    count_table,
    coverage_table,
    min_coverage = 0.3,
    min_samples = 50L
) {
  # Step 1: coverage masking — binarise, then multiply with counts
  mask <- as.data.frame(
    lapply(coverage_table, function(col) ifelse(col > min_coverage, 1, 0))
  )
  masked_counts <- as.data.frame(
    mapply(`*`, mask, count_table, SIMPLIFY = FALSE)
  )
  rownames(masked_counts) <- rownames(count_table)

  # Step 2: prevalence filter — keep MAGs present in >= min_samples samples
  row_presence <- rowSums(mask)
  keep <- row_presence >= min_samples
  masked_counts[keep, , drop = FALSE]
}

# -------------------------------------------------------------------

#' CLR-normalise a count table
#'
#' Applies centred log-ratio (CLR) transformation via
#' [compositions::clr()], adding a pseudocount first to handle zeros.
#' Equivalent to the original `clr(as.matrix(selected_count + 1))`.
#'
#' @note In the original `02_Preprocessing.Rmd` the CLR output was saved
#'   as `normalised_count`, but every downstream script (05_, 06_, 07_)
#'   references `clr_transform_df`.  The saved RData was renamed at some
#'   point; this function's output corresponds to both names.
#'
#' @param count_table A data.frame or matrix of counts (MAGs as rows,
#'   samples as columns).
#' @param pseudocount Numeric value added to all counts before
#'   transformation. Default `1`.
#'
#' @return A data.frame of CLR-transformed values preserving the original
#'   row and column names.
#' @export
#' @examples
#' \donttest{
#' data(fancy_tiny_counts)
#' data(fancy_tiny_coverage)
#' filt <- filter_mags(fancy_tiny_counts, fancy_tiny_coverage)
#' clr_df <- clr_normalize(filt)
#' dim(clr_df)
#' }
#'
#' @importFrom compositions clr
clr_normalize <- function(count_table, pseudocount = 1) {
  mat <- as.matrix(count_table) + pseudocount
  clr_mat <- compositions::clr(mat)
  out <- as.data.frame(unclass(clr_mat))
  rownames(out) <- rownames(count_table)
  out
}
