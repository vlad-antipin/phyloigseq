# PhyloIgSeq pipeline

This package allows to analyze microbiota diversity in following steps:

- Import of `phyloseq` object (generated as the result of DADA2 pipeline 
from 16S sequencing)
- Filters on samples and taxa
- Alpha-diversity analysis (richness, Shannon metrics...)
- Beta-diversity analysis with (un)constrained ordination models (PCA, 
MDS, CA, RDA, CCA ...)
- Analysis of taxonomic composition with barplots, heatmaps and 
phylogenetic trees 
- Ig coating analysis with negative dispersion plots and computing 
various Ig scores one of which is a novel sliding z-score
- Export the results along with metadata, taxonomy and taxa abundance for 
straightforward application of regression models in Feature Selector App
