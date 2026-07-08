setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("vectorOrNULL", c("vector", "NULL"))
setClassUnion("listOrNULL", c("list", "NULL"))
#' PhyloIgSeq class
#'
#' An S4 class to represent the results of an Ig-Seq experiment, including Ig coating scores,
#' fractions, and taxonomic or sample-level metadata.
#'
#' @slot ig_coating A data.frame containing calculated Ig scores per taxon/sample
#' @slot positive_fraction_name Character. Name of the positive Ig-coated fraction
#' @slot first_negative_fraction_name Character. Name of the main negative fraction (e.g., 90%)
#' @slot second_negative_fraction_name Character or NULL. Name of the secondary negative fraction (e.g., 10%)
#' @slot presorting_fraction_name Character or NULL. Name of the pre-sorting (whole community) fraction
#' @slot ig_freq_name Character or NULL. Name of the column containing total Ig+ frequency per sample
#' @slot ellipse_coords A data.frame or NULL. Stores coordinates for sliding Z-score ellipses
#' @slot sample_data A data.frame or NULL. Optional metadata for each sample
#' @slot tax_table A data.frame or NULL. Taxonomic information
#' @slot phyloseq_sample_ids Vector or NULL. Correspondence between phyloseq sample IDs and sample IDs used in \code{ig_coating}
#' @slot total_reads A data.frame or NULL. Total read counts per sample and fraction before rarefaction
#' @slot imputed_taxa List or NULL. Taxa that had zeros imputed, stored per sample
#'
#' @exportClass PhyloIgSeq
setClass(
  Class = "PhyloIgSeq",
  slots = list(
    ig_coating = "data.frame",
    positive_fraction_name = "character",
    first_negative_fraction_name = "character", # 9/10 of the whole negative fraction for IgSeq
    second_negative_fraction_name = "characterOrNULL", # 1/10 -//-
    presorting_fraction_name = "characterOrNULL", # before sorting
    ig_freq_name = "characterOrNULL",
    ellipse_coords = "data.frameOrNULL",
    sample_data = "data.frameOrNULL",
    tax_table = "data.frameOrNULL",
    phyloseq_sample_ids = "vectorOrNULL", # corresp. btw phyloseq_sam_id and sample_id
    #  you'll need this for exports from PhyloIgSeq!
    total_reads = "data.frameOrNULL",
    imputed_taxa = "listOrNULL"
  )
)

#' Collapse a List of PhyloIgSeq objects
#' @export
collapsePhyloIgSeq <- function(phyloigseq_list) {
  ig_coating <- data.frame()
  sample_data <- data.frame()
  ellipse_coords <- data.frame()
  imputed_taxa <- list()
  for (phyloigseq_obj in phyloigseq_list) {
    if (class(phyloigseq_obj) == "PhyloIgSeq") {
      # bind_rows() matches columns by name, fills in NA for missing columns
      ig_coating <- bind_rows(ig_coating, phyloigseq_obj@ig_coating)
      ellipse_coords <- bind_rows(ellipse_coords, phyloigseq_obj@ellipse_coords)
      sample_data <- bind_rows(sample_data, phyloigseq_obj@sample_data)
      imputed_taxa <- c(imputed_taxa, phyloigseq_obj@imputed_taxa)
    }
  }

  return(new(
    Class = "PhyloIgSeq",
    ig_coating = ig_coating,
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
#'
#' @return A data frame or list with computed scores per sample.
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
  scores = IG_SCORES
  # TODO: purity corrected scores
) {
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
  # To make matching and computation easier, don't carry long strings as taxa names
  taxon_ids <- seq_along(original_taxa_names)
  names(original_taxa_names) <- taxon_ids # keep the mapping
  taxa_names(physeq) <- taxon_ids

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
          metadata[[fraction_id_name]] %in% all_fraction_names,
        ,
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
    sam_metadata_row <- sam_metadata_row[,
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

    imputed_taxa <- list()
    imputed_taxa[[sample_id]] <- zero_imputation_result$imputed_taxa

    phyloigseq_list[[sample_id]] <-
      new(
        Class = "PhyloIgSeq",
        ig_coating = ig_coating,
        positive_fraction_name = positive_fraction_name,
        first_negative_fraction_name = first_negative_fraction_name,
        second_negative_fraction_name = first_negative_fraction_name,
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

  # a mapping (named vector) from "sample id" from otu_table (e.g. sam1_presorting)
  # to sample_id of Igseq (e.g. sam1)
  phyloigseq_obj@phyloseq_sample_ids <- as.matrix(sample_data(physeq))[,
    sample_id_name
  ]

  return(phyloigseq_obj)
}

# Mimics seq_table() function with these differences:
# - naming of variables and columns
# - preprocessing (rarefaction + transformation) integrated in the function
# - only those taxa that have zero counts for all fractions are excluded
# - no restriction of names or number of fractions (but they have to be unique for each sample)
