# PhyloIgSeq

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**PhyloIgSeq** provides tools for microbiota diversity analysis and immunoglobulin (Ig) coating quantification from 16S sequencing data. It wraps and extends [phyloseq](https://joey711.github.io/phyloseq/) workflows with ordination, compositional analysis, and a suite of Ig coating scores — including a novel sliding Z-score — designed for use in both interactive pipelines and web applications.

## Installation

PhyloIgSeq is currently available from GitHub only. Install it with:

```r
# install.packages("remotes")
remotes::install_github("vlad-antipin/phyloigseq")
```

One dependency ([speedyseq](https://github.com/mikemc/speedyseq)) is also GitHub-only and will be installed automatically.

Bioconductor dependencies (phyloseq, microbiome, ComplexHeatmap) must be present. If they are not already installed:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"))
```

## Overview

A typical analysis proceeds in the following stages:

1. **Import** — load a `phyloseq` object produced by a DADA2 pipeline.
2. **Filter** — apply sample- and taxa-level filters.
3. **Alpha diversity** — compute richness, Shannon entropy, and related metrics.
4. **Beta diversity** — run constrained and unconstrained ordination (PCA, MDS, CA, RDA, CCA) with permutation tests.
5. **Compositional analysis** — generate barplots, heatmaps, and phylogenetic trees.
6. **Ig coating** — build a `PhyloIgSeq` object and compute coating scores (Palm, Kau, probability-index variants, sliding Z-score).
7. **Export** — write metadata, taxonomy, and abundance tables for downstream regression modelling.

## Usage

```r
library(PhyloIgSeq)

# --- Alpha diversity ---------------------------------------------------------
alpha_plot <- plot_alpha_diversity(
  physeq    = my_phyloseq,
  measures  = c("Observed", "Shannon"),
  group_var = "condition"
)

# --- Beta diversity ----------------------------------------------------------
beta_res <- compute_beta_diversity(
  physeq   = my_phyloseq,
  method   = "PCoA",
  distance = "bray"
)

# --- Ig coating --------------------------------------------------------------
igseq <- new("PhyloIgSeq",
  ig_coating              = coating_df,
  positive_fraction_name  = "IgA+",
  first_negative_fraction_name = "IgA-"
)

scores <- compute_ig_score(
  method = "kau",
  pos    = igseq@ig_coating$pos_counts,
  neg    = igseq@ig_coating$neg_counts
)
```

## Key features

| Feature | Details |
|---|---|
| Alpha diversity | Richness, Shannon, Simpson, Faith's PD via phyloseq/microbiome |
| Beta diversity | PCA, PCoA/MDS, CA, RDA, CCA; PERMANOVA and homogeneity tests |
| Visualisation | ggplot2-based plots with faceting, grouping, and interactive plotly output |
| Compositional analysis | Stacked barplots, annotated heatmaps (ComplexHeatmap), phylogenetic trees |
| Ig coating scores | Palm, Kau, probability index/ratio, purity-corrected variants |
| Sliding Z-score | Novel per-taxon coating score with ellipse coordinate output |
| Export | Excel-ready tables (openxlsx) and HTML widgets for downstream use |

## Dependencies

PhyloIgSeq requires R ≥ 3.5 and depends on packages from CRAN, Bioconductor, and GitHub:

- **Bioconductor**: phyloseq, microbiome, ComplexHeatmap
- **CRAN**: vegan, ggplot2, ggpubr, ggrepel, plotly, dplyr, tidyr, rstatix, MASS, umap, Rtsne, and others (see `DESCRIPTION`)
- **GitHub**: [speedyseq](https://github.com/mikemc/speedyseq)

## License

MIT © Vladislav Antipin


