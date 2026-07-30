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
#' across rows even when \code{sample_id} is the same. The \code{"whole"}
#' (pre-sort) fraction row and the \code{ig_pheno} column (Ig+ frequency
#' phenotype, constant per \code{sample_id} like \code{sex}/\code{age}) are
#' likewise synthetic, added so the toy dataset can also exercise the
#' positive-fraction-only workflow (\code{positive_fraction_name} +
#' \code{presorting_fraction_name} + \code{ig_freq_name}, no negative
#' fraction) — their values carry no biological meaning.
#'
#' @format A \code{\link[phyloseq]{phyloseq}} object with:
#' \describe{
#'   \item{otu_table}{4703 taxa x 38 samples}
#'   \item{sample_data}{38 samples x 7 variables (sample_id,
#'     sorting_fraction, sex, age, batch, operator, ig_pheno)}
#'   \item{tax_table}{4703 taxa x 8 taxonomic ranks}
#' }
#' @examples
#' data(ps_igseq)
#' ps_igseq
#' phyloseq::sample_data(ps_igseq)
"ps_igseq"

#' Toy phage-display panning dataset
#'
#' A small, entirely synthetic \code{phyloseq} object simulating a phage-display
#' panning experiment, intended as a toy dataset for trying out
#' \code{\link{getPhyloIgSeq}}'s support for multiple positive fractions (see
#' the "Multiple positive fractions" section of
#' \code{vignette("igseq-scoring", package = "PhyloIgSeq")}). Each of 8
#' panning campaigns (\code{sample_id}) was sequenced at 4 stages
#' (\code{phage_fraction}): \code{"input"} (the naive/unselected input phage
#' library -- the shared negative/pre-sort fraction) and three sequential
#' output rounds (\code{"round1"}/\code{"round2"}/\code{"round3"}), each
#' progressively and increasingly enriched for a handful of campaign-specific
#' "responder" clones. Taxa represent individual phage clone sequences rather
#' than taxonomic units, so \code{tax_table} carries no real taxonomy — just a
#' single \code{sequence_id} column (identical to the taxon name).
#'
#' @format A \code{\link[phyloseq]{phyloseq}} object with:
#' \describe{
#'   \item{otu_table}{400 clones x 32 samples (8 campaigns x 4 fractions)}
#'   \item{sample_data}{32 samples x 10 variables: \code{sample_id},
#'     \code{phage_fraction} (\code{input}/\code{round1}/\code{round2}/\code{round3}),
#'     \code{target_antigen}, \code{library_type}, \code{host_species}
#'     (constant per \code{sample_id}, like \code{sex}/\code{age} in
#'     \code{ps_igseq}), \code{pct_target_specific} (a synthetic
#'     independent-assay phenotype, also constant per \code{sample_id},
#'     usable as \code{ig_freq_name} together with
#'     \code{presorting_fraction_name = "input"} to exercise the
#'     positive-fraction-only workflow), and \code{operator},
#'     \code{sequencing_batch}, \code{output_titer_log10}, \code{wash_count}
#'     (vary per fraction row within a campaign, like \code{batch}/\code{operator}
#'     in \code{ps_igseq})}
#'   \item{tax_table}{400 clones x 1 column (\code{sequence_id}); no real
#'     taxonomic ranks}
#' }
#' @examples
#' data(ps_phage_display)
#' ps_phage_display
#' phyloseq::sample_data(ps_phage_display)
#'
#' # Multiple positive fractions sharing one negative ("input") fraction:
#' pis <- getPhyloIgSeq(
#'   physeq = ps_phage_display,
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "phage_fraction",
#'   positive_fraction_name = c("round1", "round2", "round3"),
#'   first_negative_fraction_name = "input",
#'   scores = c("palm", "kau", "prob_index")
#' )
#' pis
"ps_phage_display"
