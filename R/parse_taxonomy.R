#' Parse a GTDB-format taxonomy file
#'
#' Reads a GTDB taxonomy TSV and cleans every rank by stripping rank
#' prefixes (`d__`, `p__`, etc.) and removing sub-designation suffixes on
#' Phyla and Genus (e.g. `Bacillota_C` becomes `Bacillota`).  Species are
#' simplified by removing the genus portion before the space and stripping
#' trailing digits after `"sp"`.
#'
#' When `fill_missing = TRUE` (the default):
#' * Family `"NA"` is replaced with the Order value.
#' * Genus `"NA"` or `""` is replaced with the Family value.
#' * Species `""` is replaced with `"sp"`.
#'
#' @param file Path to a GTDB taxonomy TSV file.
#' @param strain_col Character; name of the column containing MAG / strain
#'   identifiers (used as row names). Default `"Strain"`.
#' @param fill_missing Logical; fill missing taxonomy levels as described
#'   above? Default `TRUE`.
#'
#' @return A data.frame with `strain_col` values as row names and columns:
#'   Domain, Phyla, Class, Order, Family, Genus, Species.
#' @export
#' @examples
#' \donttest{
#' tax <- parse_gtdb_taxonomy("taxonomy.MAGS.2.tsv")
#' head(tax)
#' }
#'
#' @importFrom readr read_tsv
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr mutate if_else
#' @importFrom rlang .data
#' @importFrom stringr str_remove str_remove_all
parse_gtdb_taxonomy <- function(
    file,
    strain_col = "Strain",
    fill_missing = TRUE
) {
  tax <- readr::read_tsv(file, show_col_types = FALSE) |>
    tibble::column_to_rownames(strain_col)

  # Strip rank prefixes and sub-designation suffixes
  tax <- tax |>
    dplyr::mutate(
      Domain  = stringr::str_remove_all(.data$Domain,  ".*d__"),
      Phyla   = stringr::str_remove_all(.data$Phyla,   ".*p__"),
      Phyla   = stringr::str_remove_all(.data$Phyla,   "_..*"),
      Class   = stringr::str_remove_all(.data$Class,   ".*c__"),
      Order   = stringr::str_remove_all(.data$Order,   ".*o__"),
      Family  = stringr::str_remove_all(.data$Family,  ".*f__"),
      Genus   = stringr::str_remove_all(.data$Genus,   ".*g__"),
      Genus   = stringr::str_remove_all(.data$Genus,   "_..*"),
      Species = stringr::str_remove_all(.data$Species, ".*s__")
    )

  # Simplify Species: remove genus portion before space, strip trailing

  # digits after "sp"
  tax <- tax |>
    dplyr::mutate(
      Species = stringr::str_remove(.data$Species, ".*[[:space:]]"),
      Species = stringr::str_remove(.data$Species, "(?<=sp)\\d.*")
    )

  # Fill missing values

  if (fill_missing) {
    tax <- tax |>
      dplyr::mutate(
        Family  = dplyr::if_else(.data$Family == "NA",
                                 .data$Order, .data$Family),
        Genus   = dplyr::if_else(.data$Genus == "NA",
                                 .data$Family, .data$Genus),
        Genus   = dplyr::if_else(.data$Genus == "",
                                 .data$Family, .data$Genus),
        Species = dplyr::if_else(.data$Species == "",
                                 "sp", .data$Species)
      )
  }

  as.data.frame(tax)
}
