# Warn that a sample is being excluded from grouping, and why.
.warn_excluded <- function(sample_id, reason, detail = NULL) {
  msg <- paste0(sample_id, " excluded: ", reason)
  if (!is.null(detail)) {
    msg <- paste0(msg, paste0(detail, collapse = ", "))
  }
  warning(msg, "\n")
}

#' Group Ig Fractions by Sample
#'
#' Reshapes a [phyloseq::phyloseq-class] object with samples split across
#' IgSeq sort fractions (e.g. Ig+/Ig-/pre-sort) into one abundance data frame
#' per biological sample, with each fraction as a column. Optionally
#' agglomerates taxa at a given rank first, and rarefies/transforms
#' abundances within each sample so its fractions become comparable.
#'
#' A sample is excluded (with a `warning()`) if: a fraction id is duplicated
#' within it, it has fewer than two fractions, at least one fraction has zero
#' total reads, or its post-processing count table contains `NA`s.
#'
#' @param physeq A [phyloseq::phyloseq-class] object with raw counts.
#' @param taxrank Taxonomic rank to agglomerate to via [tax_glom()] before
#'   grouping, or `NULL` (default) to keep taxa as-is.
#' @param sample_id_name Name of the `sample_data` column identifying the
#'   biological sample (individual) each sort fraction belongs to.
#' @param sample_ids Which values of `sample_id_name` to include, or `NULL`
#'   (default) to include all.
#' @param fraction_id_name Name of the `sample_data` column identifying the
#'   sort fraction (e.g. `"pos"`/`"neg1"`/`"neg2"`).
#' @param fraction_ids Which values of `fraction_id_name` to include, or
#'   `NULL` (default) to include all.
#' @param rarefy_by_sample If `TRUE` (default), rarefy fraction abundances
#'   within each sample so all its fractions share the same total read count.
#' @param fractions_to_rarefy Which fractions to rarefy when
#'   `rarefy_by_sample = TRUE`, or `NULL` (default) to rarefy all of them.
#' @param transform_by_sample Abundance transformation applied within each
#'   sample after rarefaction: `"identity"` (default, no transformation) or
#'   `"compositional"`.
#'
#' @return A named list of data frames, one per included `sample_id`. Each
#'   data frame has one row per taxon and columns `sample_id`, `taxon_id`,
#'   and one column per retained fraction (its abundance).
#'
#' @examples
#' data(ps_igseq)
#' grouped <- group_sorted_samples(
#'   physeq = ps_igseq,
#'   sample_id_name = "sample_id",
#'   sample_ids = c("sample_1", "sample_2"),
#'   fraction_id_name = "sorting_fraction",
#'   fraction_ids = c("pos", "neg1", "neg2")
#' )
#' names(grouped)
#' head(grouped[["sample_1"]])
#'
#' @export
group_sorted_samples <- function(
  physeq, # containing raw counts
  taxrank = NULL,
  sample_id_name, # name identifying each sample_id (individual for
  # which sorting was performed)
  sample_ids = NULL, # if NULL, all samples are taken
  fraction_id_name, # name of the fraction column
  fraction_ids = NULL, # if NULL, all fractions are taken
  rarefy_by_sample = TRUE, # inside each sample, rarefy abundances for each
  # fraction (so that fractions of the same sample have the same total sum of reads)
  fractions_to_rarefy = NULL, # if NULL, all are rarefied
  transform_by_sample = "identity" # transformation to apply to abundances
  # separately inside each sample
) {
  if (!transform_by_sample %in% c("identity", "compositional")) {
    transform_by_sample <- "identity"
  }
  .check_phyloseq(physeq)

  if (nsamples(physeq) == 0) {
    stop("No samples are present")
  }

  physeq <- .prune_by_fraction(physeq, fraction_id_name, fraction_ids)

  if (!is.null(taxrank)) {
    physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank])
  }

  # By default, keep abundances as-is (no transformation)
  if (
    is.null(transform_by_sample) ||
      !transform_by_sample %in% c("identity", "compositional")
  ) {
    warning(
      "You should use counts or relative abundance, setting transformation to 'identity'...\n"
    )
    transform_by_sample <- "identity"
  }

  full_sample_data <- sample_data(physeq) %>% as("data.frame")

  # Taxa are assumed to be by rows for this analysis
  # TODO: optimize for sparse otu table
  if (!taxa_are_rows(physeq)) {
    abundance_table <- t(otu_table(physeq)) %>% as("matrix")
  } else {
    abundance_table <- otu_table(physeq) %>% as("matrix")
  }

  if (is.null(sample_ids)) {
    sample_ids <- unique(full_sample_data[, sample_id_name])
  }

  # Inside each sample (composed of multiple fractions), perform preprocessing
  # and data grouping by fraction
  sample_list <- list() # list of grouped data by sample

  for (sample_id in sample_ids) {
    #sample_list_tmp = list()
    # row_ids - ID of a unique observation ("sample" from the point of view of phyloseq) from sample_data()
    row_ids <- rownames(full_sample_data[
      full_sample_data[, sample_id_name] == sample_id,
    ])
    fraction_ids <- full_sample_data[row_ids, fraction_id_name]

    # Each fraction should be unique and therefore there must be one to one correspondence
    # unique row_ids <-> unique fraction_ids inside each sample
    if (sum(duplicated(fraction_ids)) != 0) {
      .warn_excluded(
        sample_id,
        "duplicated fraction(s): ",
        unique(fraction_ids[duplicated(fraction_ids)])
      )
      next
    }

    # Each sample has to have at least two fractions - otherwise there's nothing
    # to compare downstream
    if (length(fraction_ids) <= 1) {
      .warn_excluded(sample_id, "only one or no fraction: ", fraction_ids)
      next
    }

    # If at least one fraction doesn't have any reads for any taxon, the whole sample is excluded
    countsums_by_fraction <- colSums(abundance_table[, row_ids]) #physeq %>% prune_samples(row_ids,.) %>% sample_sums()
    if (0 %in% countsums_by_fraction) {
      .warn_excluded(
        sample_id,
        "no reads for at least one fraction: ",
        fraction_ids[countsums_by_fraction == 0]
      )
      next
    }

    taxa_counts_by_fraction <- abundance_table[, row_ids]

    # Use actual names of fractions (and not row_ids from sample_names() of phyloseq)
    colnames(taxa_counts_by_fraction) <- fraction_ids

    # Exclude taxa that have 0 reads in ALL fractions
    taxa_counts_by_fraction <- taxa_counts_by_fraction[
      rowSums(taxa_counts_by_fraction) > 0,
    ]

    # Rarefaction makes the total sum of counts equal among all fractions
    # (separately for each sample)
    # select fractions -> (rarefy) -> transform abundances -> get otu_table
    if (rarefy_by_sample) {
      fractions_to_rarefy_by_sample <-
        if (is.null(fractions_to_rarefy)) {
          fraction_ids
        } else {
          intersect(fraction_ids, fractions_to_rarefy)
        }

      taxa_counts_by_fraction[, fractions_to_rarefy_by_sample] <-
        rarefy_abundances(
          abundance_table = taxa_counts_by_fraction[,
            fractions_to_rarefy_by_sample,
            drop = FALSE
          ],
          trim_taxa = FALSE,
          taxa_are_rows = TRUE,
          silent_warnings = TRUE
        )
    }

    if (transform_by_sample != "identity") {
      taxa_counts_by_fraction <- transform_abundances(
        abundance_table = taxa_counts_by_fraction,
        transform = transform_by_sample,
        taxa_are_rows = TRUE
      )
    }

    taxa_counts_by_fraction <- as.data.frame(taxa_counts_by_fraction)

    if (anyNA(taxa_counts_by_fraction)) {
      .warn_excluded(sample_id, "NA(s) in its count table")
      next
    }

    grouped_data <- cbind(
      data.frame(sample_id = rep(sample_id, nrow(taxa_counts_by_fraction))),
      data.frame(taxon_id = rownames(taxa_counts_by_fraction)),
      taxa_counts_by_fraction
    )
    rownames(grouped_data) <- NULL
    grouped_data <- grouped_data[
      order(apply(grouped_data[, fraction_ids], 1, sum), decreasing = TRUE),
    ]
    sample_list[[sample_id]] <- grouped_data
  }

  return(
    sample_list # a list (by sample) of data frames - taxa x (sample info and taxonomy)
    # each fraction is in a separate column of each data frame
  )
}

#' Handle Zero Abundance
#'
#' Resolves zero counts in a taxa (rows) by fraction (columns) abundance data
#' frame before Ig-score computation, since several scores are undefined when
#' a fraction abundance is exactly zero. Taxa that are zero in every retained
#' fraction are always dropped first; `method` then controls how any
#' *remaining* per-fraction zeros are handled.
#'
#' @param data A data frame of taxa (rows) by fractions and metadata
#'   (columns), as produced by [group_sorted_samples()].
#' @param fraction_names Names of the columns in `data` holding fraction
#'   abundances to impute; other columns (e.g. `taxon_id`) are left as-is.
#' @param method One of `"no_zero"` (drop any taxon with a zero in any
#'   fraction), `"pseudo_count"` (add half the minimum nonzero count across
#'   all fractions to every count), `"random_pseudo_count"` (add a count
#'   drawn uniformly between a small fraction of, and, the minimum nonzero
#'   count, to every count), `"bayesian_inference"` (Bayesian multiplicative
#'   zero replacement via [zCompositions::cmultRepl()]; only zeros are
#'   modified), or `"keep_zeros"` (leave zeros as-is).
#'
#' @return A list with `data` (the input data frame, restricted to retained
#'   taxa, with the `fraction_names` columns updated per `method`) and
#'   `imputed_taxa` (the `taxon_id`s that had at least one zero fraction
#'   before imputation, or `NULL` if `data` has no `taxon_id` column, or none
#'   were imputed).
#'
#' @examples
#' data(ps_igseq)
#' grouped <- group_sorted_samples(
#'   physeq = ps_igseq,
#'   sample_id_name = "sample_id",
#'   sample_ids = c("sample_1", "sample_2"),
#'   fraction_id_name = "sorting_fraction",
#'   fraction_ids = c("pos", "neg1", "neg2")
#' )
#' result <- impute_zeros(
#'   data = grouped[["sample_1"]],
#'   fraction_names = c("pos", "neg1", "neg2"),
#'   method = "pseudo_count"
#' )
#' head(result$data)
#' head(result$imputed_taxa)
#'
#' @export
impute_zeros <- function(
  data, # dataframe taxa (rows) x fractions (cols)
  fraction_names,
  method
) {
  fraction_names <- fraction_names[fraction_names %in% colnames(data)]
  fractions <- data[, fraction_names, drop = FALSE]
  # Don't forget to keep track of taxa as rownames and apply the exclusion
  # of taxa from sorted_sample_df, using rownames
  rownames(fractions) <- 1:nrow(data)
  # Exclude taxa that are absent in every fraction in any case
  fractions <- fractions[
    !apply(fractions, 1, function(row) {
      all(row == 0)
    }),
    ,
    drop = FALSE
  ]

  rows_with_zeros <- rownames(fractions)[apply(fractions, 1, function(row) {
    any(row == 0)
  })]

  switch(
    method,
    no_zero = {
      # Exclude all taxa having zero abundance in any fraction
      fractions <- fractions[
        !apply(fractions, 1, function(row) {
          any(row == 0)
        }),
        ,
        drop = FALSE
      ]
      rows_with_zeros <- NULL
    },
    pseudo_count = {
      # Add a fixed pseudo count (a half of the minimum count observed across all fractions)
      # to ALL counts
      fractions <- fractions + min(fractions[fractions != 0]) / 2
    },
    random_pseudo_count = {
      # Add a uniformely random psedocount between nearly zero and a minimum count
      # observed across all fractions to ALL counts
      min_count <- min(fractions[fractions != 0])
      fractions <- fractions +
        matrix(
          runif(
            n = prod(dim(fractions)),
            min = min_count / 1000,
            max = min_count
          ),
          nrow = dim(fractions)[1],
          ncol = dim(fractions)[2]
        )
    },
    bayesian_inference = {
      # Use imputation of zeros with Bayesian models, ONLY ZEROS are modified
      fractions <- zCompositions::cmultRepl(
        X = fractions,
        method = "BL",
        output = "p-counts"
      )
    },
    keep_zeros = {
      # keep zeros...
    },
    stop("Wrong 'method' argument")
  )

  imputed_taxa <- NULL
  if ("taxon_id" %in% colnames(data)) {
    imputed_taxa <- data$taxon_id[1:nrow(data) %in% rows_with_zeros]
  }

  # Rownames of fractions contain original indices of row in input dataframe
  # so we can update the taxonomy for the original data (from sorted_sample_df)
  data <- data[1:nrow(data) %in% rownames(fractions), ]

  # Update only fractions' counts
  data[, fraction_names] <- fractions

  if (length(imputed_taxa) == 0) {
    imputed_taxa <- NULL
  }

  return(list(data = data, imputed_taxa = imputed_taxa))
}


#' Validate a Single Ig+ Frequency Value
#'
#' Internal. Turns one raw `sample_data` cell into a probability in `[0, 1]`, or `NA_real_`
#' (with a `warning()`) when it can't be one: absent/duplicated, non-numeric, or out of range
#' after unit conversion. Rejecting non-numerics before any arithmetic or comparison keeps a
#' character/factor column from silently reaching `compute_ig_score()`. Treating a bad value as
#' missing for that sample only -- rather than aborting -- matches `getPhyloIgSeq()`'s general
#' "skip what can't be computed, don't abort the batch" behaviour.
#'
#' @param value The raw cell, of any type.
#' @param units `"frequency"` (already a probability) or `"percent"` (divided by 100).
#' @param sample_id,fraction_label,column_name Used to identify the offender in warnings;
#'   `fraction_label` is the fraction the value was read from (`NULL` in the `"wide"` layout,
#'   where the value isn't tied to one).
#'
#' @return A single numeric in `[0, 1]`, or `NA_real_`.
#'
#' @noRd
.resolve_ig_freq_value <- function(
  value,
  units,
  sample_id,
  fraction_label = NULL,
  column_name = "ig_freq"
) {
  where <- paste0(
    "sample ",
    sample_id,
    if (!is.null(fraction_label)) paste0(", fraction ", fraction_label),
    " (column `",
    column_name,
    "`)"
  )

  # length 0: the column doesn't exist, or the sample has no row for this
  # fraction. length > 1: duplicated fraction rows. Either way there is no
  # single value to use -- and guarding here keeps `if (is.na(value))` below
  # from erroring on the `logical(0)` a missing column would produce.
  if (length(value) != 1) {
    if (length(value) > 1) {
      warning(
        "ig_freq for ",
        where,
        " resolves to ",
        length(value),
        " values instead of one; treating as NA.\n"
      )
    }
    return(NA_real_)
  }

  if (is.na(value)) {
    return(NA_real_)
  }

  if (!is.numeric(value)) {
    warning(
      "ig_freq for ",
      where,
      " is not numeric (got class '",
      class(value)[1],
      "', value '",
      value,
      "'); treating as NA.\n"
    )
    return(NA_real_)
  }

  if (identical(units, "percent")) {
    value <- value / 100
  }

  if (value < 0 || value > 1) {
    warning(
      "ig_freq for ",
      where,
      " is ",
      signif(value, 4),
      " after unit conversion (`ig_freq_units = \"",
      units,
      "\"`), outside the expected [0, 1] probability range; treating as NA.\n"
    )
    return(NA_real_)
  }

  as.numeric(value)
}


#' Resolve the Per-Fraction Ig+ Frequencies of One Sample
#'
#' Reads `getPhyloIgSeq()`'s `ig_freq_name` column out of one sample's `sample_data`
#' rows and returns the Ig+ frequency of each sort fraction, all validated the same way
#' [getPhyloIgSeq()] validates them internally.
#'
#' Exported because callers that bypass [getPhyloIgSeq()] still need the same numbers on
#' the same terms: the companion Shiny app's MA-plot preview goes straight from
#' [group_sorted_samples()] to [get_ma_coordinates()], and the purity-corrected change
#' axis needs `pos_ig_freq`/`neg1_ig_freq` resolved with the same unit conversion and
#' `[0, 1]` range checking rather than re-implemented.
#'
#' The two layouts differ in what the column is taken to mean:
#' \describe{
#'   \item{`"wide"`}{One Ig+ frequency per biological sample, repeated on (or recorded on only
#'     one of) its fraction rows -- the historical layout. Only the pre-sort frequency P(Ig+) is
#'     recoverable; the in-fraction frequencies are `NA`, since a single overall P(Ig+) says
#'     nothing about how pure either sorted fraction turned out to be.}
#'   \item{`"long"`}{The frequency is fraction-specific: each fraction's row records the Ig+
#'     frequency measured *in that fraction*. This yields all four values -- the pre-sort P(Ig+)
#'     plus the "positive fraction purity" and "negative fraction impurity" needed by the
#'     `purity_corrected_*` scores.}
#' }
#'
#' @param sam_metadata_df One sample's `sample_data` rows (one per fraction), still carrying the
#'   `fraction_id_name` column.
#' @param fraction_id_name Name of the fraction column in `sam_metadata_df`.
#' @param ig_freq_name Name of the Ig+ frequency column, or `NULL` (all four come back `NA`).
#' @param ig_freq_units,ig_freq_layout See [getPhyloIgSeq()].
#' @param positive_fraction_name The *single* positive fraction currently being scored, so that
#'   `pos_ig_freq` follows the positive fraction of the synthetic sample it belongs to.
#' @param first_negative_fraction_name,second_negative_fraction_name,presorting_fraction_name
#'   The remaining fraction names, any of which may be `NULL`.
#' @param sample_id Used to identify the sample in warnings.
#'
#' @return A one-row data frame with numeric columns `presort_ig_freq`, `pos_ig_freq`,
#'   `neg1_ig_freq` and `neg2_ig_freq`.
#'
#' @examples
#' data(ps_igseq)
#' metadata <- as(phyloseq::sample_data(ps_igseq), "data.frame")
#' resolve_ig_freqs(
#'   sam_metadata_df = metadata[metadata$sample_id == "sample_1", ],
#'   fraction_id_name = "sorting_fraction",
#'   ig_freq_name = "ig_pheno",
#'   ig_freq_units = "frequency",
#'   ig_freq_layout = "wide",
#'   positive_fraction_name = "pos",
#'   first_negative_fraction_name = "neg1"
#' )
#'
#' @export
resolve_ig_freqs <- function(
  sam_metadata_df,
  fraction_id_name,
  ig_freq_name,
  ig_freq_units,
  ig_freq_layout,
  positive_fraction_name,
  first_negative_fraction_name = NULL,
  second_negative_fraction_name = NULL,
  presorting_fraction_name = NULL,
  sample_id = NA_character_
) {
  empty <- data.frame(
    presort_ig_freq = NA_real_,
    pos_ig_freq = NA_real_,
    neg1_ig_freq = NA_real_,
    neg2_ig_freq = NA_real_
  )

  if (is.null(ig_freq_name) || !ig_freq_name %in% names(sam_metadata_df)) {
    if (!is.null(ig_freq_name)) {
      warning(
        "`ig_freq_name` = '",
        ig_freq_name,
        "' is not a column of sample_data; treating Ig+ frequencies as NA.\n"
      )
    }
    return(empty)
  }

  raw <- sam_metadata_df[[ig_freq_name]]

  if (identical(ig_freq_layout, "wide")) {
    # Historical behaviour: the value is a property of the sample, so it must be
    # the same on every fraction row (the same "identical across fractions, else
    # NA" rule the rest of the metadata is collapsed by).
    unique_values <- unique(raw)
    if (length(unique_values) > 1) {
      warning(
        "ig_freq column `",
        ig_freq_name,
        "` varies across the fractions of sample ",
        sample_id,
        " but `ig_freq_layout = \"wide\"` expects one value per sample, so it ",
        "cannot be used for this sample. If the column records the Ig+ ",
        "frequency measured separately in each fraction, use ",
        "`ig_freq_layout = \"long\"`.\n"
      )
      return(empty)
    }
    empty$presort_ig_freq <- .resolve_ig_freq_value(
      value = unique_values,
      units = ig_freq_units,
      sample_id = sample_id,
      column_name = ig_freq_name
    )
    return(empty)
  }

  # "long": one measurement per fraction, picked out by the fraction's own name.
  for_fraction <- function(fraction) {
    if (is.null(fraction)) {
      return(NA_real_)
    }
    .resolve_ig_freq_value(
      value = raw[sam_metadata_df[[fraction_id_name]] == fraction],
      units = ig_freq_units,
      sample_id = sample_id,
      fraction_label = fraction,
      column_name = ig_freq_name
    )
  }

  data.frame(
    presort_ig_freq = for_fraction(presorting_fraction_name),
    pos_ig_freq = for_fraction(positive_fraction_name),
    neg1_ig_freq = for_fraction(first_negative_fraction_name),
    neg2_ig_freq = for_fraction(second_negative_fraction_name)
  )
}
