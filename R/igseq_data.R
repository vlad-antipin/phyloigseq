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
#'   sort fraction (e.g. `"Pos"`/`"Neg1"`/`"Neg2"`).
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
#'   fraction_ids = c("Pos", "Neg1", "Neg2")
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
#'   fraction_ids = c("Pos", "Neg1", "Neg2")
#' )
#' result <- impute_zeros(
#'   data = grouped[["sample_1"]],
#'   fraction_names = c("Pos", "Neg1", "Neg2"),
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
