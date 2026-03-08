## ---- Build 100 MAGs x 321 samples real subset ----
## Expands the 100-sample vignette dataset to all 321 available samples.
## CLR values come from 02_preprocessed.RData (normalised_count), which was
## CLR-transformed on the full 2178-MAG reference.
##
## Source files (all in D:/Fancy/):
##   02_preprocessed.RData          -> normalised_count (CLR, 2178x346),
##                                     selected_count  (raw, 2178x346),
##                                     metadata        (321 rows)
##   00._Raw-selected-321samples2178vars.RData -> metadata (321 rows)
##   Covered_fraction.csv           -> coverage (14262 x 346)
##   data-raw/real_subset_100x100.rds -> MAG IDs, taxonomy, known_mags

cat("Loading existing 100x100 subset for MAG IDs and taxonomy...\n")
dat100 <- readRDS("data-raw/real_subset_100x100.rds")
mag_ids <- dat100$mags          # 100 original MAG IDs
taxonomy <- dat100$taxonomy     # 100 x 7
known_mags <- dat100$known_mags # 4 known nonlinear MAG IDs
all_taxonomy <- dat100$all_taxonomy  # full taxonomy table

cat("Loading 02_preprocessed.RData...\n")
load("02_preprocessed.RData")
# Provides: normalised_count (CLR, 2178 x 346), selected_count (raw, 2178 x 346), metadata (321 x 6)

## Identify the 321 samples present in both metadata and normalised_count
samples_321 <- intersect(metadata$Sample, colnames(normalised_count))
stopifnot(length(samples_321) == 321L)

cat("Subsetting CLR matrix: 100 MAGs x 321 samples...\n")
clr_321 <- normalised_count[mag_ids, samples_321]
stopifnot(identical(dim(clr_321), c(100L, 321L)))

cat("Subsetting raw count matrix: 100 MAGs x 321 samples...\n")
counts_321 <- selected_count[mag_ids, samples_321]
stopifnot(identical(dim(counts_321), c(100L, 321L)))

cat("Loading coverage data...\n")
cov_raw <- read.csv("Covered_fraction.csv", row.names = 1, check.names = FALSE)
## Coverage CSV column names match normalised_count column names
cov_samples <- intersect(samples_321, colnames(cov_raw))
cov_mags    <- intersect(mag_ids, rownames(cov_raw))
cat("  Coverage available for", length(cov_mags), "of 100 MAGs and",
    length(cov_samples), "of 321 samples\n")

## Build coverage matrix; fill missing MAGs/samples with NA
coverage_321 <- matrix(NA_real_, nrow = 100, ncol = 321,
                       dimnames = list(mag_ids, samples_321))
coverage_321[cov_mags, cov_samples] <- as.matrix(cov_raw[cov_mags, cov_samples])

cat("Building metadata for 321 samples...\n")
meta_321 <- metadata[match(samples_321, metadata$Sample), ]
## Strip data_origin column (privacy: contains "Rainer")
meta_321$data_origin <- NULL
rownames(meta_321) <- NULL

## Keep only Sample, Breed, CH4 (sufficient for demo)
meta_321 <- meta_321[, c("Sample", "Breed", "CH4")]

cat("Assembling final list...\n")
result <- list(
  clr        = clr_321,
  counts     = counts_321,
  coverage   = coverage_321,
  taxonomy   = taxonomy,
  metadata   = meta_321,
  mags       = mag_ids,
  samples    = samples_321,
  known_mags = known_mags,
  all_taxonomy = all_taxonomy
)

out_path <- "data-raw/real_subset_100x321.rds"
saveRDS(result, out_path)
cat("Saved:", out_path, "\n")
cat("CLR dims:", paste(dim(result$clr), collapse = " x "), "\n")
cat("Metadata dims:", paste(dim(result$metadata), collapse = " x "), "\n")
cat("Done.\n")
