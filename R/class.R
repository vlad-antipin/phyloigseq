setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
#' PhyloIgSeq class
#'
#' An S4 class to represent the results of an IgSeq experiment, including Ig coating scores,
#' fractions, and taxonomic or sample-level metadata.
#'
#' @details
#' \code{ig_coating} always carries the identifier columns \code{taxon_id}/\code{sample_id},
#' followed by the actual Ig score columns (named in \code{score_names}), followed by
#' whichever per-fraction abundance and scoring-diagnostic columns (e.g. the fraction names
#' themselves, \code{zeros_imputed}, \code{ellipse_level}, \code{obs_change}/\code{obs_abundance}/
#' \code{null_change}/\code{null_abundance}) happened to be produced upstream. Only the columns
#' listed in \code{score_names} are Ig-coating scores; use \code{\link{get_ig_score}} to pull one
#' out without having to know which of the remaining columns are metadata.
#'
#' @slot ig_coating A data.frame containing per-taxon/sample Ig scores plus supporting metadata
#'   (see Details)
#' @slot score_names Character. Names of the \code{ig_coating} columns that are actual Ig scores
#'   (as opposed to fraction/diagnostic metadata columns)
#' @slot positive_fraction_name Character. Name of the positive Ig-coated fraction
#' @slot first_negative_fraction_name Character. Name of the main negative fraction (e.g., 90%)
#' @slot second_negative_fraction_name Character or NULL. Name of the secondary negative fraction (e.g., 10%)
#' @slot presorting_fraction_name Character or NULL. Name of the pre-sorting (whole community) fraction
#' @slot ig_freq_name Character or NULL. Name of the column containing total Ig+ frequency per sample
#' @slot ellipse_coords A data.frame or NULL. Stores coordinates for sliding Z-score ellipses
#' @slot sample_data A data.frame or NULL. Optional metadata for each sample
#' @slot tax_table A data.frame or NULL. Taxonomic information
#' @slot total_reads A data.frame or NULL. Total read counts per sample and fraction before rarefaction
#' @slot imputed_taxa List or NULL. Taxa that had zeros imputed, stored per sample
#'
#' @exportClass PhyloIgSeq
setClass(
  Class = "PhyloIgSeq",
  slots = list(
    ig_coating = "data.frame",
    score_names = "character",
    positive_fraction_name = "character",
    first_negative_fraction_name = "character", # 9/10 of the whole negative fraction for IgSeq
    second_negative_fraction_name = "characterOrNULL", # 1/10 -//-
    presorting_fraction_name = "characterOrNULL", # before sorting
    ig_freq_name = "characterOrNULL",
    ellipse_coords = "data.frameOrNULL",
    sample_data = "data.frameOrNULL",
    tax_table = "data.frameOrNULL",
    total_reads = "data.frameOrNULL",
    imputed_taxa = "listOrNULL"
  ),
  prototype = list(score_names = character(0))
)

#' @rdname PhyloIgSeq-class
#' @param object A PhyloIgSeq object.
setMethod("show", "PhyloIgSeq", function(object) {
  n_samples <- length(unique(object@ig_coating$sample_id))
  n_taxa <- length(unique(object@ig_coating$taxon_id))
  scores_label <- if (length(object@score_names) > 0) {
    paste(object@score_names, collapse = ", ")
  } else {
    "(none computed)"
  }

  fraction_names <- c(
    positive = object@positive_fraction_name,
    neg1 = object@first_negative_fraction_name,
    neg2 = object@second_negative_fraction_name,
    presort = object@presorting_fraction_name
  )

  cat("PhyloIgSeq-class Ig-coating scoring result\n")
  cat(sprintf(
    "ig_coating       Ig scores:      [ %s ] across %d sample(s), %d taxa\n",
    scores_label,
    n_samples,
    n_taxa
  ))
  if (length(fraction_names) > 0) {
    cat(
      "Fractions        ",
      paste(sprintf("%s=\"%s\"", names(fraction_names), fraction_names), collapse = ", "),
      "\n"
    )
  }
  if (!is.null(object@sample_data)) {
    cat(sprintf(
      "sample_data()    Sample Data:    [ %d samples by %d sample variables ]\n",
      nrow(object@sample_data),
      ncol(object@sample_data)
    ))
  }
  if (!is.null(object@tax_table)) {
    cat(sprintf(
      "tax_table()      Taxonomy Table: [ %d taxa by %d taxonomic ranks ]\n",
      nrow(object@tax_table),
      ncol(object@tax_table)
    ))
  }
  invisible(NULL)
})

#' Collapse a List of PhyloIgSeq objects
#'
#' Row-binds the \code{ig_coating}, \code{ellipse_coords}, \code{sample_data} and
#' \code{imputed_taxa} of a list of per-sample \code{\link{PhyloIgSeq-class}} objects (as produced
#' internally by \code{\link{getPhyloIgSeq}}) into a single object.
#'
#' @param phyloigseq_list A list of \code{PhyloIgSeq} objects. Elements that are not a
#'   \code{PhyloIgSeq} object are skipped.
#'
#' @return A single \code{PhyloIgSeq} object. Its \code{score_names} is the union of
#'   \code{score_names} across all input objects; \code{tax_table} and the fraction-name slots are
#'   not combined here (left at their defaults) and are set by the caller.
#'
#' @examples
#' pis_1 <- new(
#'   "PhyloIgSeq",
#'   ig_coating = data.frame(taxon_id = 1:2, sample_id = "s1", slide_z = c(0.5, -0.2)),
#'   score_names = "slide_z",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg"
#' )
#' pis_2 <- new(
#'   "PhyloIgSeq",
#'   ig_coating = data.frame(taxon_id = 1:2, sample_id = "s2", slide_z = c(1.1, 0.3)),
#'   score_names = "slide_z",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg"
#' )
#' collapsePhyloIgSeq(list(pis_1, pis_2))
#'
#' @export
collapsePhyloIgSeq <- function(phyloigseq_list) {
  ig_coating <- data.frame()
  sample_data <- data.frame()
  ellipse_coords <- data.frame()
  imputed_taxa <- list()
  score_names <- character(0)
  for (phyloigseq_obj in phyloigseq_list) {
    if (is(phyloigseq_obj, "PhyloIgSeq")) {
      # bind_rows() matches columns by name, fills in NA for missing columns
      ig_coating <- bind_rows(ig_coating, phyloigseq_obj@ig_coating)
      ellipse_coords <- bind_rows(ellipse_coords, phyloigseq_obj@ellipse_coords)
      sample_data <- bind_rows(sample_data, phyloigseq_obj@sample_data)
      imputed_taxa <- c(imputed_taxa, phyloigseq_obj@imputed_taxa)
      score_names <- union(score_names, phyloigseq_obj@score_names)
    }
  }

  return(new(
    Class = "PhyloIgSeq",
    ig_coating = ig_coating,
    score_names = score_names,
    ellipse_coords = ellipse_coords,
    sample_data = sample_data,
    tax_table = NULL,
    imputed_taxa = imputed_taxa
  ))
}

#' Compute Ig Scores from a Phyloseq Object
#'
#' This function computes various Ig scores based on sample fraction data.
#'
#' @param physeq A `phyloseq` object containing raw count data.
#' @param taxrank Character or NULL. Taxonomic rank to agglomerate to before computing scores.
#' @param sample_id_name Name of the column identifying unique samples.
#' @param sample_ids Optional. A character vector of sample IDs to subset.
#' @param fraction_id_name Name of the column indicating fraction (e.g., pos, neg).
#' @param rarefy_by_sample Logical. Rarefy read counts across fractions within each sample.
#' @param transform_by_sample Transformation method (e.g., "identity", "log").
#' @param positive_fraction_name Name of the positive fraction.
#' @param first_negative_fraction_name Name of the first negative fraction.
#' @param second_negative_fraction_name Optional. Name of a second negative fraction.
#' @param presorting_fraction_name Optional. Name of the presorting fraction.
#' @param ig_freq_name Optional. Column name for Ig frequency, if precomputed.
#' @param zero_treatment How to handle zeros ("no_zero", "pseudocount", etc.).
#' @param window_size Integer. Window size for smoothing.
#' @param empirical_null_distribution Logical. Whether to estimate null distribution.
#' @param confidence_levels Optional. Confidence levels for scoring.
#' @param scores Vector of score names to compute.
#' @param taxon_id_source How to derive the `taxon_id` used throughout `ig_coating`/`tax_table`:
#'   `"sequential"` (default) renumbers taxa as fresh sequential integers, recoverable via
#'   `tax_table$taxon_name`; `"original"` uses `physeq`'s own (possibly `taxrank`-agglomerated,
#'   and not necessarily ASV-level if `physeq` was already agglomerated upstream) taxa name
#'   directly as `taxon_id`, matching the identifiers shown by [group_sorted_samples()] when
#'   called directly (e.g. for a single-sample preview).
#'
#' @return A \code{\link{PhyloIgSeq-class}} object. Its \code{ig_coating} slot holds one row per
#'   taxon/sample with the requested \code{scores} as columns (also recorded in \code{score_names})
#'   plus supporting fraction/diagnostic columns; see \code{\link{get_ig_score}} to retrieve a
#'   single score without dealing with the rest of \code{ig_coating}.
#'
#' @examples
#' data(ps_igseq)
#' pis <- getPhyloIgSeq(
#'   physeq = ps_igseq,
#'   sample_ids = c("sample_1", "sample_2", "sample_3"),
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2",
#'   scores = c("slide_z", "palm", "kau")
#' )
#' pis
#' get_ig_score(pis, score_name = "slide_z", sample_ids = "sample_1")
#'
#' @export
getPhyloIgSeq <- function(
  physeq, # containing raw counts
  taxrank = NULL,
  sample_id_name, # name identifying each sample_id (individual for
  # which sorting was performed)
  sample_ids = NULL, # if NULL, all samples taken
  fraction_id_name, # name of the fraction column
  rarefy_by_sample = TRUE, # inside each sample, rarefy abundances for each
  # fraction (so that fractions of the same sample have the same total sum of reads)
  transform_by_sample = "identity",
  positive_fraction_name = "pos",
  first_negative_fraction_name = "neg",
  second_negative_fraction_name = NULL,
  presorting_fraction_name = NULL, # before sorting
  ig_freq_name = NULL,
  zero_treatment = "no_zero",
  window_size = 50,
  empirical_null_distribution = TRUE,
  confidence_levels = NULL,
  scores = IG_SCORES,
  # TODO: purity corrected scores
  taxon_id_source = c("sequential", "original")
) {
  taxon_id_source <- match.arg(taxon_id_source)
  if (is.null(second_negative_fraction_name) & empirical_null_distribution) {
    warning(
      "No second negative fraction furnished, cannot model empirical null ( Ig-.1 vs Ig-.2) distribution...\n"
    )
    empirical_null_distribution <- FALSE
  }
  if (!is.null(taxrank)) {
    physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank])
  }

  original_taxa_names <- taxa_names(physeq)
  if (taxon_id_source == "sequential") {
    # To make matching and computation easier, don't carry long strings as taxa names
    taxon_ids <- seq_along(original_taxa_names)
    names(original_taxa_names) <- taxon_ids # keep the mapping
    taxa_names(physeq) <- taxon_ids
  } else {
    # "original": keep physeq's own (possibly taxrank-agglomerated) taxa_names
    # as taxon_id directly instead of renumbering — original_taxa_names becomes
    # an identity mapping so downstream lookups (e.g. `taxon_name = original_taxa_names[taxon_ids_to_keep]`
    # below) still work unchanged.
    names(original_taxa_names) <- original_taxa_names
  }

  all_fraction_names <- c(
    positive_fraction_name,
    first_negative_fraction_name,
    second_negative_fraction_name,
    presorting_fraction_name
  )

  grouped_data <-
    group_sorted_samples(
      physeq = physeq,
      taxrank = NULL, # already done above
      sample_id_name = sample_id_name,
      sample_ids = sample_ids,
      fraction_id_name = fraction_id_name,
      fraction_ids = all_fraction_names,
      rarefy_by_sample = rarefy_by_sample,
      # NOTE: only pos, neg1 and neg2 are rarefied
      fractions_to_rarefy = c(
        positive_fraction_name,
        first_negative_fraction_name,
        second_negative_fraction_name
      ),
      transform_by_sample = transform_by_sample
    )

  if (is.null(sample_ids)) {
    sample_ids <- names(grouped_data)
  } else {
    sample_ids <- sample_ids[sample_ids %in% names(grouped_data)]
  }

  if (length(sample_ids) == 0) {
    warning("No suitable samples left, returning NULL\n")
    return(NULL)
  }

  metadata <- phyloseq::sample_data(physeq) %>% as("data.frame")
  # to avoid potential problems if there is already a column called "sample_id"
  # which is not a name of smaple_id column indicated by a user
  names(metadata)[
    names(metadata) != sample_id_name & names(metadata) == "sample_id"
  ] <-
    paste0(
      "original___",
      names(metadata)[
        names(metadata) != sample_id_name & names(metadata) == "sample_id"
      ]
    )

  names(metadata)[names(metadata) == sample_id_name] <- "sample_id"
  # put sample_id on the first place
  metadata <- metadata[, c("sample_id", setdiff(names(metadata), "sample_id"))]

  phyloigseq_list <- list()

  total_reads <- data.frame()

  for (sample_id in sample_ids) {
    # Keep track of total presorting read counts for each sample (used to compute relative abundances later)
    if (!is.null(presorting_fraction_name)) {
      total_reads <- rbind(
        total_reads,
        data.frame(
          sample_id = sample_id,
          total_reads = if (
            all(is.na(grouped_data[[sample_id]][[presorting_fraction_name]]))
          ) {
            NA
          } else {
            sum(
              grouped_data[[sample_id]][[presorting_fraction_name]],
              na.rm = TRUE
            )
          }
        )
      )
    }
    # Retrieve sample  metadata
    sam_metadata_df <-
      metadata[
        metadata[["sample_id"]] == sample_id &
          metadata[[fraction_id_name]] %in% all_fraction_names, ,
        drop = FALSE
      ]
    sam_metadata_row <- data.frame(matrix(NA, nrow = 1, ncol = ncol(metadata)))
    names(sam_metadata_row) <- names(metadata)

    for (var_name in names(metadata)) {
      unique_values <- unique(sam_metadata_df[[var_name]])

      if (length(unique_values) == 1) {
        sam_metadata_row[[var_name]] <- unique_values
      }
    }
    sam_metadata_row <- sam_metadata_row[
      ,
      names(sam_metadata_row) != fraction_id_name
    ]

    # Handle zeros
    present_fraction_names <- all_fraction_names[
      all_fraction_names %in% colnames(grouped_data[[sample_id]])
    ]

    zero_imputation_result <-
      impute_zeros(
        data = grouped_data[[sample_id]],
        # Don't impute zeros in other fractions!
        fraction_names = intersect(
          present_fraction_names,
          c(
            positive_fraction_name,
            first_negative_fraction_name,
            second_negative_fraction_name
          )
        ),
        method = zero_treatment
      )

    ig_coating <- zero_imputation_result$data %>%
      select(all_of(unique(c("taxon_id", "sample_id", present_fraction_names))))

    ig_coating$zeros_imputed <- ig_coating$taxon_id %in%
      zero_imputation_result$imputed_taxa

    if (nrow(ig_coating) == 0) {
      warning(paste0(
        sample_id,
        " excluded: no taxa left after zero treatment\n"
      ))
      next
    }

    # Compute scores
    if ("slide_z" %in% scores) {
      slide_z_result <-
        get_slide_z(
          sorted_sample_df = ig_coating,
          positive_fraction_name = positive_fraction_name,
          first_negative_fraction_name = first_negative_fraction_name,
          second_negative_fraction_name = second_negative_fraction_name,
          window_size = window_size,
          empirical_null_distribution = empirical_null_distribution,
          confidence_levels = confidence_levels,
          imputed_taxa = zero_imputation_result$imputed_taxa
        )

      ig_coating$slide_z <- slide_z_result$slide_z
      ig_coating$ellipse_level <- slide_z_result$ellipse_level
      if (prod(dim(slide_z_result$ma_coords)) != 0) {
        ig_coating <- cbind(
          ig_coating,
          slide_z_result$ma_coords[, c(
            "obs_change",
            "obs_abundance",
            if (empirical_null_distribution) {
              c("null_change", "null_abundance")
            }
          )]
        )
      }
      ellipse_coords <- slide_z_result$ellipse_coords
    } else {
      ellipse_coords <- data.frame()
    }

    # Other Ig scores:
    for (score in scores[scores != "slide_z"]) {
      ig_coating[[score]] <-
        compute_ig_score(
          method = score,
          pos = ig_coating[[positive_fraction_name]],
          neg = ig_coating[[first_negative_fraction_name]],
          pre = if (!is.null(presorting_fraction_name)) {
            ig_coating[[presorting_fraction_name]]
          },
          ig_freq = if (!is.null(ig_freq_name)) {
            sam_metadata_row[[ig_freq_name]]
          }
          # TODO: purity corrected scores:
          # pos_purity = pos_purity, # P(real Ig+ | Ig+ fraction)
          # neg_impurity = neg_impurity, # P(real Ig+ | Ig- fraction)
          # pos_fraction = pos_fraction, # P(Ig+ fraction)
          # neg_fraction = neg_fraction
        )
    }

    # Ig scores are the columns users look for first; keep them right after the
    # taxon_id/sample_id identifiers, ahead of fraction/diagnostic columns.
    score_names_present <- intersect(scores, names(ig_coating))
    ig_coating <- ig_coating %>%
      relocate(all_of(score_names_present), .after = "sample_id")

    imputed_taxa <- list()
    imputed_taxa[[sample_id]] <- zero_imputation_result$imputed_taxa

    phyloigseq_list[[sample_id]] <-
      new(
        Class = "PhyloIgSeq",
        ig_coating = ig_coating,
        score_names = score_names_present,
        positive_fraction_name = positive_fraction_name,
        first_negative_fraction_name = first_negative_fraction_name,
        second_negative_fraction_name = second_negative_fraction_name,
        ellipse_coords = ellipse_coords,
        sample_data = sam_metadata_row,
        tax_table = NULL,
        imputed_taxa = imputed_taxa
      )
  }

  phyloigseq_obj <- collapsePhyloIgSeq(phyloigseq_list)

  if (prod(dim(total_reads)) != 0) {
    names(total_reads)[2] <- presorting_fraction_name
    phyloigseq_obj@total_reads <- total_reads[
      total_reads$sample_id %in% phyloigseq_obj@ig_coating$sample_id,
    ]
  } else {
    phyloigseq_obj@total_reads <- NULL
  }

  phyloigseq_obj@positive_fraction_name <- positive_fraction_name
  phyloigseq_obj@first_negative_fraction_name <- first_negative_fraction_name
  phyloigseq_obj@second_negative_fraction_name <- second_negative_fraction_name
  phyloigseq_obj@presorting_fraction_name <- presorting_fraction_name
  phyloigseq_obj@ig_freq_name <- ig_freq_name

  taxon_ids_to_keep <- unique(phyloigseq_obj@ig_coating$taxon_id)
  tax_table <- as.matrix(phyloseq::tax_table(physeq)@.Data)[
    taxon_ids_to_keep,
  ] %>%
    as.data.frame()

  # In the hierarchical order
  tax_table <- cbind(
    tax_table,
    data.frame(
      taxon_name = original_taxa_names[taxon_ids_to_keep],
      taxon_id = taxon_ids_to_keep
    )
  )
  rownames(tax_table) <- NULL
  # TODO: keep only taxa that are left
  phyloigseq_obj@tax_table <- tax_table

  return(phyloigseq_obj)
}

#' Get an Ig Score from a PhyloIgSeq Object
#'
#' Retrieves a single Ig score from \code{ig_coating}, optionally restricted to a subset of taxa
#' and/or samples, without having to know which of \code{ig_coating}'s other columns are metadata.
#'
#' @param phyloigseq_obj A \code{\link{PhyloIgSeq-class}} object.
#' @param score_name Character. Name of the score to retrieve; must be one of
#'   \code{phyloigseq_obj@score_names}.
#' @param taxa_ids Optional. A vector of taxon IDs to restrict to; \code{NULL} (the default) keeps
#'   all taxa.
#' @param sample_ids Optional. A vector of sample IDs to restrict to; \code{NULL} (the default)
#'   keeps all samples.
#'
#' @return A data frame with columns \code{taxon_id}, \code{sample_id} and \code{score_name},
#'   filtered to \code{taxa_ids}/\code{sample_ids} if given.
#'
#' @examples
#' data(ps_igseq)
#' pis <- getPhyloIgSeq(
#'   physeq = ps_igseq,
#'   sample_ids = c("sample_1", "sample_2", "sample_3"),
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2",
#'   scores = c("slide_z", "palm", "kau")
#' )
#' get_ig_score(pis, score_name = "palm", sample_ids = c("sample_1", "sample_2"))
#'
#' @export
get_ig_score <- function(
  phyloigseq_obj,
  score_name,
  taxa_ids = NULL,
  sample_ids = NULL
) {
  if (!is(phyloigseq_obj, "PhyloIgSeq")) {
    stop("`phyloigseq_obj` must be a PhyloIgSeq object")
  }
  if (!score_name %in% phyloigseq_obj@score_names) {
    stop(
      "`score_name` must be one of: ",
      paste(phyloigseq_obj@score_names, collapse = ", ")
    )
  }

  result <- phyloigseq_obj@ig_coating[, c("taxon_id", "sample_id", score_name)]
  if (!is.null(taxa_ids)) {
    result <- result[result$taxon_id %in% taxa_ids, ]
  }
  if (!is.null(sample_ids)) {
    result <- result[result$sample_id %in% sample_ids, ]
  }
  rownames(result) <- NULL
  result
}

# Mimics seq_table() function with these differences:
# - naming of variables and columns
# - preprocessing (rarefaction + transformation) integrated in the function
# - only those taxa that have zero counts for all fractions are excluded
# - no restriction of names or number of fractions (but they have to be unique for each sample)
