# PhyloIgSeq

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

<!-- badges: end -->

PhyloIgSeq extends [phyloseq](https://joey711.github.io/phyloseq/) with tools for microbiota diversity analysis and IgSeq (immunoglobulin coating) quantification from 16S sequencing data.

## Why PhyloIgSeq

- **IgSeq scores, including a novel sliding Z-score.** `getPhyloIgSeq()` builds an Ig-coating profile from Ig+/Ig− (and optional pre-sort) fractions and scores it with classical Palm and Kau indices, probability-index/ratio variants, purity-corrected variants, and a sliding Z-score developed for this package (`compute_slide_z()`) that estimates a local null distribution along the abundance axis instead of assuming one global null. When a second negative fraction is available, the null can instead be modeled empirically from Ig−.1 vs Ig−.2, capturing technical variability directly rather than assuming a theoretical null.
- **A sparse-matrix engine for beta diversity.** `as_sparse_phyloseq()` keeps the OTU table as a `dgCMatrix`; `sparse_distance()` and `sparse_unifrac()` compute Bray-Curtis, Jaccard, (weighted) UniFrac and other distances directly on it, without ever densifying the table. This is a drop-in, numerically equivalent replacement for `phyloseq::distance()`/`UniFrac()` that scales to much larger ASV tables.
- **Approximate UniFrac without a real phylogeny.** `get_taxonomy_tree()` builds a `phy_tree` from the ranks in `tax_table` alone, so UniFrac-family distances remain available when no sequencing-based tree exists. This is a prototype/heuristic construction and not a substitute for a real phylogeny, treat resulting distances as approximate.
- **A companion app for reproducible analysis.** [phyloigseq-app](https://www.funkycells.com/main/index.php/lab-tools/phyloigseq) exposes the same pipeline through a point-and-click Shiny interface and exports every filter, transform, and scoring parameter used in a session to a dedicated sheet alongside the results, so an analysis can be reproduced exactly.

## Installation

PhyloIgSeq is currently available from GitHub only:

```r
# install.packages("remotes")
remotes::install_github("vlad-antipin/phyloigseq")
```

This pulls in [speedyseq](https://github.com/mikemc/speedyseq) automatically. Bioconductor dependencies (`phyloseq`, `microbiome`, `ComplexHeatmap`) must be installed separately if not already present:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("phyloseq", "microbiome", "ComplexHeatmap"))
```

## Quick start

```r
library(PhyloIgSeq)
library(phyloseq)
data(ps_16s_refinement) # bundled example: 515 taxa x 37 samples, with a phy_tree

# --- Alpha diversity ---------------------------------------------------------
alpha <- get_alpha_diversity(ps_16s_refinement, measure = "Shannon")
plot_alpha_diversity(alpha, x = "Protocol")

# --- Beta diversity (Bray-Curtis PCoA) ---------------------------------------
beta <- get_beta_diversity(ps_16s_refinement, method = "PCoA", dist = "bray")
ggplot_beta_diversity(beta)

# --- Same distances via the sparse engine, without densifying the table -----
ps_sparse <- as_sparse_phyloseq(ps_16s_refinement)
sparse_distance(ps_sparse, method = "bray")

# --- UniFrac from an approximate, taxonomy-derived tree ---------------------
phy_tree(ps_16s_refinement) <- get_taxonomy_tree(ps_16s_refinement)
sparse_unifrac(as_sparse_phyloseq(ps_16s_refinement), method = "wunifrac")
```

IgSeq scoring works on Ig+/Ig− fraction data, identified by columns in `sample_data`:

```r
igseq <- getPhyloIgSeq(
  physeq,
  sample_id_name = "sample_id",   # column identifying each biological sample
  fraction_id_name = "fraction",  # column indicating the sort fraction
  positive_fraction_name = "IgA+",
  first_negative_fraction_name = "IgA-"
)
plot_ig_score(igseq, score_name = "slide_z")
```

The lower-level scoring function operates directly on abundance vectors:

```r
compute_ig_score(method = "kau", pos = c(120, 45, 0, 300), neg = c(80, 60, 10, 250))
```

## Documentation

- `vignette("introduction", package = "PhyloIgSeq")`
- `?getPhyloIgSeq`, `?compute_slide_z`, `?sparse_distance`, `?sparse_unifrac`, `?get_taxonomy_tree`

## License

MIT © Vladislav Antipin, Martin Larsen
