# Fancy

**Frequency and Nonlinear Correlation Hybrid Network Inference**

Fancy (Frequency And Nonlinear Correlation hYbrid) infers association
networks by combining k-nearest neighbour mutual information (via MRNET) and
distance correlation into a single hybrid edge score. Bootstrap resampling
provides frequency-based confidence, and a weighted scoring function ranks
edges for downstream module analysis. While designed for microbial
co-abundance networks from metagenome-assembled genome (MAG) count tables,
Fancy works with any numeric feature matrix — including gene co-expression,
metabolite correlations, and cross-omics (e.g. gene-metabolite) associations.

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

## Why Fancy?

Standard correlation methods (e.g. Pearson, Spearman) miss nonlinear
microbial interactions. Fancy's hybrid score captures patterns that
linear metrics cannot, including threshold effects, context-dependent
relationships, and feedback loops between MAGs.

| Threshold competitive exclusion | L-shaped nonlinear relationship | Context-dependent interaction |
|:---:|:---:|:---:|
| ![](man/figures/scatter-pair1-1.png) | ![](man/figures/scatter-pair2-1.png) | ![](man/figures/scatter-pair3-1.png) |
| Bulleidia vs AC2028: Pearson r = -0.47 but the relationship is nonlinear with a clear threshold. Fancy hybrid score = 0.40. | RUG023 vs Cryptobacteroides: Pearson r = -0.31 underestimates a strong L-shaped dependency. Fancy hybrid score = 0.25. | CAG-791 vs Methanobrevibacter: Pearson r = -0.02 (no linear signal), yet Fancy hybrid score = 0.32 reveals a hidden association. |

Points are coloured by methane emission level (CH4 g/day: blue = low,
yellow = mid, red = high), highlighting environment-dependent structure
within these interactions.

## Quick start

```r
library(Fancy)

# Run the full pipeline on a CLR-normalised count table
# (~15 minutes on the example data with 100 bootstraps and 4 CPUs;
#  runtime scales with dataset size and number of bootstraps)
result <- fancy(t(fancy_tiny_clr), n_bootstrap = 100, cpus = 4)

# Inspect the result
# Show the hybrid score distribution
plot(result)

# Export for Cytoscape
data(fancy_tiny_taxonomy)
export_cytoscape(result, fancy_tiny_taxonomy, file_prefix = "my_network")
```

## License

GPL (>= 3)

## Citation

If you use this software, please cite:

Fancy v0.99.0  
DOI: https://doi.org/10.5281/zenodo.20003547
