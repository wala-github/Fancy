# Fancy

**Hybrid Network Inference for Microbiome Co-Abundance Analysis**

Fancy (Frequency And Nonlinear Correlation hYbrid) infers microbial co-abundance
networks by combining k-nearest neighbour mutual information (via MRNET) and 
distance correlation into a single hybrid edge score. Bootstrap resampling 
provides frequency-based confidence, and a weighted scoring function ranks 
edges for downstream module analysis. Designed for metagenome-assembled
genome (MAG) count tables.

## Installation

```r
# From GitHub
devtools::install_github("wala-github/Fancy")

# From Bioconductor (once accepted)
BiocManager::install("Fancy")
```

## Package structure

```
Fancy/
├── R/                 Core source code (13 files)
│   ├── fancy.R            Main pipeline wrapper
│   ├── bootstrap.R        Parallel bootstrap network inference (snowfall)
│   ├── mi_network.R       Mutual information matrices (kNN estimator)
│   ├── dcor_network.R     Distance correlation matrices (energy::dcor)
│   ├── scoring.R          Edge frequency and hybrid score computation
│   ├── preprocess.R       MAG filtering by coverage and prevalence
│   ├── parse_count_tables.R  Sample name cleaning
│   ├── parse_taxonomy.R   GTDB taxonomy parsing
│   ├── plot.R             S3 plot method for fancy objects
│   ├── export.R           Cytoscape export with phyla colour palette
│   ├── utils.R            Internal helpers
│   ├── data.R             Dataset documentation
│   └── Fancy-package.R    Package-level documentation
│
├── man/               roxygen2-generated documentation (31 .Rd files)
│
├── data/              Bundled example data
│   └── fancy_tiny.rda     Contains 5 datasets (100 MAGs x 321 samples):
│                          fancy_tiny_clr, fancy_tiny_counts,
│                          fancy_tiny_coverage, fancy_tiny_taxonomy,
│                          fancy_tiny_metadata
│
├── inst/extdata/      Reserved for external example data files
│
├── tests/             testthat test suite
│   ├── testthat.R         Test harness
│   └── testthat/          Test scripts
│
├── vignettes/         Package vignettes
│   └── Fancy_workflow.Rmd HTML vignette (BiocStyle) walking through
│                          the full hybrid network inference pipeline
│
├── DESCRIPTION        Package metadata and dependencies
└── NAMESPACE          Exports and imports (14 exported functions)
```

## Quick start

```r
library(Fancy)

# Run the full pipeline on a CLR-normalised count table
# (~15 minutes on the example data with 100 bootstraps and 4 CPUs;
#  runtime scales with dataset size and number of bootstraps)
result <- fancy(t(fancy_tiny_clr), n_bootstrap = 100, cpus = 4)

# Inspect the result
plot(result)

# Export for Cytoscape
data(fancy_tiny_taxonomy)
export_cytoscape(result, fancy_tiny_taxonomy, file_prefix = "my_network")
```

## License

GPL (>= 3)
