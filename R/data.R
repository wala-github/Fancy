# ---- Tiny dataset (100 MAGs x 321 samples) --------------------------------

#' Tiny example CLR-transformed abundance table
#'
#' Pre-computed centred log-ratio (CLR) transformed abundance table for
#' 100 MAGs across 321 rumen samples. The CLR was computed on the full
#' 2178-MAG compositional reference before subsetting, preserving the
#' proper geometric mean. This is the primary input for network inference
#' with [fancy()].
#'
#' MAGs 1--4 are known nonlinear interaction pairs:
#' MAGs 1--2 (Bulleidia vs AC2028) show threshold competitive exclusion,
#' and MAGs 3--4 (RUG023 vs Cryptobacteroides) show an L-shaped
#' nonlinear relationship.
#'
#' @format A data.frame with 100 rows (MAGs) and 321 columns (samples).
#'   Values are CLR-transformed (continuous, centred around zero).
#'
#' @source Subset from a published rumen metagenome study via
#'   `data-raw/make_datasets.R` with `set.seed(42)`.
#'
#' @examples
#' data(fancy_tiny_clr)
#' dim(fancy_tiny_clr)
#'
#' @docType data
"fancy_tiny_clr"

#' Tiny example count table
#'
#' Raw read-count table for the same 100 MAGs and 321 samples as
#' [fancy_tiny_clr]. Provided for coverage filtering and exploratory
#' analysis; network inference should use the pre-CLR matrix.
#'
#' @format A data.frame with 100 rows (MAGs) and 321 columns (samples).
#'   Values are non-negative integers.
#'
#' @source Subset from a published rumen metagenome study via
#'   `data-raw/make_datasets.R` with `set.seed(42)`.
#'
#' @examples
#' data(fancy_tiny_counts)
#' dim(fancy_tiny_counts)
#'
#' @docType data
"fancy_tiny_counts"

#' Tiny example covered-fraction table
#'
#' Genome covered-fraction table matching [fancy_tiny_clr].
#' Values represent the fraction of each MAG's reference genome
#' covered by mapped reads in each sample.
#'
#' @format A data.frame with 100 rows (MAGs) and 321 columns (samples).
#'   Values are numeric fractions in \[0, 1\].
#'
#' @source Subset from a published rumen metagenome study via
#'   `data-raw/make_datasets.R` with `set.seed(42)`.
#'
#' @examples
#' data(fancy_tiny_coverage)
#' summary(unlist(fancy_tiny_coverage))
#'
#' @docType data
"fancy_tiny_coverage"

#' Tiny example GTDB taxonomy table
#'
#' GTDB taxonomy for the 100 MAGs in [fancy_tiny_clr].
#' MAGs 1--2 are a nonlinear threshold-exclusion pair
#' (Bulleidia vs AC2028, both Bacillota), and MAGs 3--4 are a
#' nonlinear L-shaped pair (RUG023 / Spirochaetota vs
#' Cryptobacteroides / Bacteroidota). The 100 MAGs span 9 phyla.
#'
#' @format A data.frame with 100 rows (MAGs as row names) and 7
#'   columns: Domain, Phyla, Class, Order, Family, Genus, Species.
#'
#' @source Subset from a published rumen metagenome study via
#'   `data-raw/make_datasets.R` with `set.seed(42)`.
#'
#' @examples
#' data(fancy_tiny_taxonomy)
#' table(fancy_tiny_taxonomy$Phyla)
#'
#' @docType data
"fancy_tiny_taxonomy"

#' Tiny example sample metadata
#'
#' Per-sample metadata for the 321 samples in [fancy_tiny_clr].
#'
#' @format A data.frame with 321 rows and 3 columns:
#' \describe{
#'   \item{Sample}{Sample identifier (anonymised).}
#'   \item{Breed}{Cattle breed (Aberdeen_Angus_X, Charolais_X,
#'     Limousin_X, or Luing).}
#'   \item{CH4}{Methane emissions in grams per day (continuous,
#'     range approximately 8--36).}
#' }
#'
#' @source Subset from a published rumen metagenome study via
#'   `data-raw/make_datasets.R` with `set.seed(42)`.
#'
#' @examples
#' data(fancy_tiny_metadata)
#' hist(fancy_tiny_metadata$CH4)
#'
#' @docType data
"fancy_tiny_metadata"
