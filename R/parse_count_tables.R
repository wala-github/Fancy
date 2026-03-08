# ---- Sample-name cleaning ------------------------------------------------

#' Clean sample names from TSV column headers
#'
#' Removes the `X` prefix that [utils::read.csv()] adds to column names
#' starting with a digit, replaces `.` with `-` to restore the original
#' sample IDs, and strips format-specific suffixes such as
#' `.QUALITY_PASSED_R1.fastq.*`.
#'
#' @param x A character vector of sample names (typically `colnames(df)`).
#' @param remove_prefix Logical; remove leading `X` added by R? Default
#'   `TRUE`.
#' @param dot_to_dash Logical; replace `.` with `-`? Default `TRUE`.
#' @param suffix_patterns Character vector of regex patterns to strip from
#'   the end of each name.
#'
#' @return A character vector of cleaned names.
#' @export
#' @examples
#' clean_sample_names(c("X12.34.AB.QUALITY_PASSED_R1.fastq.Read.Count"))
clean_sample_names <- function(
    x,
    remove_prefix = TRUE,
    dot_to_dash = TRUE,
    suffix_patterns = c(
      "\\.QUALITY_PASSED_R1\\.fastq\\.Read\\.Count$",
      "\\.QUALITY_PASSED_R1\\.fastq\\.Covered\\.Fraction$",
      "\\.QUALITY_PASSED_R1\\.fastq\\.Relative\\.Abundance.*$"
    )
) {
  for (pat in suffix_patterns) {
    x <- gsub(pat, "", x)
  }
  if (remove_prefix) {
    x <- gsub("^X", "", x)
  }
  if (dot_to_dash) {
    x <- gsub("\\.", "-", x)
  }
  x
}

# ---- Internal: read one column from a set of per-sample TSVs -------------

#' Read a single column from each TSV and merge by MAG ID
#'
#' @param file_paths Character vector of file paths.
#' @param column Integer; which data column to select (1-based, after the
#'   row-names column).
#' @param id_col Character; name used for the MAG-ID join column.
#'
#' @return A data.frame with MAG IDs in the first column and one sample
#'   column per file.
#'
#' @importFrom purrr map reduce
#' @importFrom dplyr full_join
#' @importFrom tibble rownames_to_column
#' @importFrom utils read.csv
#' @keywords internal
.read_column_from_tsvs <- function(file_paths, column, id_col = "mags_ids") {
  purrr::map(file_paths, function(fp) {
    file_data <- utils::read.csv(fp, header = TRUE, sep = "\t",
                                 row.names = 1)
    file_data <- file_data[, column, drop = FALSE]
    tibble::rownames_to_column(file_data, var = id_col)
  }) |>
    purrr::reduce(dplyr::full_join, by = id_col)
}

# ---- Parse count tables --------------------------------------------------

#' Parse per-sample TSV files into combined MAG-by-sample matrices
#'
#' Reads a directory of per-sample TSV files (as produced by read-mapping
#' pipelines) and merges them into MAG-by-sample matrices for read counts,
#' relative abundance, and covered fraction.
#'
#' Each TSV is expected to have MAG identifiers as row names and at least
#' three data columns (count, relative abundance, covered fraction).  The
#' first row in the raw files is often a duplicated header and is dropped by
#' default (`drop_first_row = TRUE`).
#'
#' @param path Directory containing the TSV files.
#' @param pattern Glob pattern passed to [base::list.files()].
#' @param columns Named integer vector mapping result names to column
#'   indices.  Defaults to
#'   `c(count = 1L, relative_abundance = 2L, covered_fraction = 3L)`.
#'   Set to `NULL` or a subset to read only the matrices you need.
#' @param clean_names Logical; apply [clean_sample_names()] to column
#'   names?
#' @param drop_first_row Logical; drop the first data row (header
#'   duplication artifact in raw TSVs)?
#'
#' @return A named list of data.frames, one per entry in `columns`.  Each
#'   data.frame has MAG IDs as row names and samples as columns.
#' @export
#' @examples
#' \donttest{
#' tables <- parse_count_tables("path/to/tsv_dir")
#' count   <- tables$count
#' cov_frac <- tables$covered_fraction
#' }
#'
#' @importFrom dplyr select
parse_count_tables <- function(
    path,
    pattern = "*.tsv",
    columns = c(
      count = 1L,
      relative_abundance = 2L,
      covered_fraction = 3L
    ),
    clean_names = TRUE,
    drop_first_row = TRUE
) {
  file_paths <- list.files(path, pattern = pattern, full.names = TRUE)
  if (length(file_paths) == 0L) {
    stop("No files matching '", pattern, "' found in ", path)
  }

  result <- lapply(stats::setNames(columns, names(columns)), function(col) {
    combined <- .read_column_from_tsvs(file_paths, column = col)

    if (drop_first_row && nrow(combined) > 1L) {
      combined <- combined[-1L, ]
    }

    # Set MAG IDs as rownames
    id_col <- colnames(combined)[1L]
    rownames(combined) <- combined[[id_col]]
    combined[[id_col]] <- NULL

    if (clean_names) {
      colnames(combined) <- clean_sample_names(colnames(combined))
    }

    combined
  })

  result
}
