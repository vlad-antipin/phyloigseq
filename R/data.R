#' 16S sequencing approach comparison dataset
#'
#' A \code{phyloseq} object used to benchmark and compare different 16S
#' sequencing approaches. Contains 515 taxa across 37 samples with 9 sample
#' metadata variables, 7 taxonomic ranks, and a phylogenetic tree.
#'
#' @format A \code{\link[phyloseq]{phyloseq}} object with:
#' \describe{
#'   \item{otu_table}{515 taxa x 37 samples}
#'   \item{sample_data}{37 samples x 9 variables}
#'   \item{tax_table}{515 taxa x 7 taxonomic ranks (Kingdom to Species)}
#'   \item{phy_tree}{515 tips, 513 internal nodes}
#' }
"ps_16s_refinement"
