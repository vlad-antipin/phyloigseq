#' 16S sequencing approach comparison dataset
#'
#' A \code{phyloseq} object used to benchmark and compare different 16S
#' sequencing approaches. Contains 515 taxa across 37 samples with 10 sample
#' metadata variables, 7 taxonomic ranks, and a phylogenetic tree.
#'
#' @format A \code{\link[phyloseq]{phyloseq}} object with:
#' \describe{
#'   \item{otu_table}{515 taxa x 37 samples}
#'   \item{sample_data}{37 samples x 10 variables, including
#'     \code{LogConcentration}, the base-10 logarithm of \code{Concentration}
#'     (e.g. \code{8} for \code{"10^8"}), provided as a continuous sample
#'     variable}
#'   \item{tax_table}{515 taxa x 7 taxonomic ranks (Kingdom to Species)}
#'   \item{phy_tree}{515 tips, 513 internal nodes}
#' }
#' @examples
#' data(ps_16s_refinement)
#' ps_16s_refinement
#' phyloseq::sample_variables(ps_16s_refinement)
"ps_16s_refinement"

#' Toy IgSeq dataset
#'
#' A small, anonymized \code{phyloseq} object derived from a mouse IgA-Seq
#' experiment, intended as a toy dataset for examples and tests. Sample names and
#' \code{sample_id} values have been anonymized (\code{experiment_*} /
#' \code{sample_*}). \code{sex} and \code{age} are synthetic, randomly-generated
#' sample-level variables not present in the original data (constant across all
#' rows sharing the same \code{sample_id}). \code{batch} and \code{operator} are
#' synthetic, randomly-generated experiment-level variables, so they may vary
#' across rows even when \code{sample_id} is the same.
#'
#' @format A \code{\link[phyloseq]{phyloseq}} object with:
#' \describe{
#'   \item{otu_table}{4703 taxa x 30 samples}
#'   \item{sample_data}{30 samples x 6 variables (sample_id,
#'     sorting_fraction, sex, age, batch, operator)}
#'   \item{tax_table}{4703 taxa x 8 taxonomic ranks}
#' }
#' @examples
#' data(ps_igseq)
#' ps_igseq
#' phyloseq::sample_data(ps_igseq)
"ps_igseq"
