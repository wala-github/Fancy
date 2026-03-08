#' Canonicalise edge direction
#'
#' Swaps `source` and `target` so the lexicographically smaller node is
#' always in the `source` column.
#'
#' @param df A data.frame with at least `source` and `target` columns
#'   (character).
#'
#' @return The same data.frame with `source <= target` for every row.
#'
#' @keywords internal
.remove_directionality <- function(df) {
  df$source <- as.character(df$source)
  df$target <- as.character(df$target)

  idx <- df$source > df$target

  tmp <- df$source[idx]
  df$source[idx] <- df$target[idx]
  df$target[idx] <- tmp

  df
}
