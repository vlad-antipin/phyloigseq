#' Group Ig Fractions by Sample
#' @export
group_sorted_samples = function(
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
    transform_by_sample = "identity"
  }
  if (class(physeq) != "phyloseq") {
    stop("Need a phyloseq object")
  }

  if (nsamples(physeq) == 0) {
    stop("No samples are present")
  }

  if (!is.null(fraction_id_name) & !is.null(fraction_ids)) {
    physeq = prune_samples(
      sample_data(physeq)[[fraction_id_name]] %in% fraction_ids,
      physeq
    )
  }

  if (!is.null(taxrank)) {
    physeq = tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) = make.unique(tax_table(physeq)[, taxrank])
  }

  # By default, keep abundances as-is (no transformation)
  if (
    is.null(transform_by_sample) ||
      !transform_by_sample %in% c("identity", "compositional")
  ) {
    warning(
      "You should use counts or relative abundance, setting transformation to 'identity'...\n"
    )
    transform_by_sample = "identity"
  }

  full_sample_data = sample_data(physeq) %>% as("data.frame")

  # Taxa are assumed to be by rows for this analysis
  if (!taxa_are_rows(physeq)) {
    abundance_table = t(otu_table(physeq)@.Data)
  } else {
    abundance_table = otu_table(physeq)@.Data
  }

  if (is.null(sample_ids)) {
    sample_ids = unique(full_sample_data[, sample_id_name])
  }

  # Inside each sample (composed of multiple fractions), perform preprocessing
  # and data grouping by fraction
  sample_list = list() # list of grouped data by sample

  for (sample_id in sample_ids) {
    #sample_list_tmp = list()
    # row_ids - ID of a unique observation ("sample" from the point of view of phyloseq) from sample_data()
    row_ids = rownames(full_sample_data[
      full_sample_data[, sample_id_name] == sample_id,
    ])
    fraction_ids = full_sample_data[row_ids, fraction_id_name]

    # Each fraction should be unique and therefore there must be one to one correspondence
    # unique row_ids <-> unique fraction_ids inside each sample
    if (sum(duplicated(fraction_ids)) != 0) {
      warning(paste0(
        sample_id,
        " excluded: duplicated fraction(s): ",
        paste0(unique(fraction_ids[duplicated(fraction_ids)]), collapse = ", "),
        "\n"
      ))
      next
    }

    # Each sample has to have at least two fractions - otherwise there's nothing
    # to compare downstream
    if (length(fraction_ids) <= 1) {
      warning(paste0(
        sample_id,
        " excluded: only one or no fraction: ",
        paste0(fraction_ids, collapse = ", "),
        "\n"
      ))
      next
    }

    # If at least one fraction doesn't have any reads for any taxon, the whole sample is excluded
    countsums_by_fraction = colSums(abundance_table[, row_ids]) #physeq %>% prune_samples(row_ids,.) %>% sample_sums()
    if (0 %in% countsums_by_fraction) {
      warning(paste0(
        sample_id,
        " excluded: no reads for at least one fraction: ",
        paste0(fraction_ids[countsums_by_fraction == 0], collapse = ", "),
        "\n"
      ))
      next
    }

    taxa_counts_by_fraction = abundance_table[, row_ids]

    # Use actual names of fractions (and not row_ids from sample_names() of phyloseq)
    colnames(taxa_counts_by_fraction) = fraction_ids

    # Exclude taxa that have 0 reads in ALL fractions
    taxa_counts_by_fraction = taxa_counts_by_fraction[
      rowSums(taxa_counts_by_fraction) > 0,
    ]

    # Rarefaction makes the total sum of counts equal among all fractions
    # (separately for each sample)
    # select fractions -> (rarefy) -> transform abundances -> get otu_table
    if (rarefy_by_sample) {
      fractions_to_rarefy_by_sample =
        if (is.null(fractions_to_rarefy)) {
          fraction_ids
        } else {
          intersect(fraction_ids, fractions_to_rarefy)
        }

      taxa_counts_by_fraction[, fractions_to_rarefy_by_sample] =
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
      taxa_counts_by_fraction = transform_abundances(
        abundance_table = taxa_counts_by_fraction,
        transform = transform_by_sample,
        taxa_are_rows = TRUE
      )
    }

    taxa_counts_by_fraction = as.data.frame(taxa_counts_by_fraction)

    if (anyNA(taxa_counts_by_fraction)) {
      warning(paste0(sample_id, " excluded: NA(s) in its count table\n"))
      next
    }

    grouped_data = cbind(
      data.frame(sample_id = rep(sample_id, nrow(taxa_counts_by_fraction))),
      data.frame(taxon_id = rownames(taxa_counts_by_fraction)),
      taxa_counts_by_fraction
    )
    rownames(grouped_data) = NULL
    grouped_data = grouped_data[
      order(apply(grouped_data[, fraction_ids], 1, sum), decreasing = TRUE),
    ]
    sample_list[[sample_id]] = grouped_data
  }

  return(
    sample_list # a list (by sample) of data frames - taxa x (sample info and taxonomy)
    # each fraction is in a separate column of each data frame
  )
}

#' Handle Zero Abundance
#' @export
impute_zeros = function(
  data, # dataframe taxa (rows) x fractions (cols)
  fraction_names,
  method
) {
  fraction_names = fraction_names[fraction_names %in% colnames(data)]
  fractions = data[, fraction_names, drop = FALSE]
  # Don't forget to keep track of taxa as rownames and apply the exclusion
  # of taxa from sorted_sample_df, using rownames
  rownames(fractions) = 1:nrow(data)
  # Exclude taxa that are absent in every fraction in any case
  fractions = fractions[
    !apply(fractions, 1, function(row) {
      all(row == 0)
    }),
    ,
    drop = FALSE
  ]

  rows_with_zeros = rownames(fractions)[apply(fractions, 1, function(row) {
    any(row == 0)
  })]

  if (method == "no_zero") {
    # Exclude all taxa having zero abundance in any fraction
    fractions = fractions[
      !apply(fractions, 1, function(row) {
        any(row == 0)
      }),
      ,
      drop = FALSE
    ]
    rows_with_zeros = NULL
  } else if (method == "pseudo_count") {
    # Add a fixed pseudo count (a half of the minimum count observed across all fractions)
    # to ALL counts
    fractions = fractions + min(fractions[fractions != 0]) / 2
  } else if (method == "random_pseudo_count") {
    # Add a uniformely random psedocount between nearly zero and a minimum count
    # observed across all fractions to ALL counts
    min_count = min(fractions[fractions != 0])
    fractions = fractions +
      matrix(
        runif(
          n = prod(dim(fractions)),
          min = min_count / 1000,
          max = min_count
        ),
        nrow = dim(fractions)[1],
        ncol = dim(fractions)[2]
      )
  } else if (method == "bayesian_inference") {
    # Use imputation of zeros with Bayesian models, ONLY ZEROS are modified
    # TODO: check the theory
    fractions = zCompositions::cmultRepl(
      X = fractions,
      method = "BL",
      output = "p-counts"
    )
  } else if (method == "keep_zeros") {
    # keep zeros...
  } else {
    stop("Wrong 'method' argument")
  }

  if ("taxon_id" %in% colnames(data)) {
    imputed_taxa = data$taxon_id[1:nrow(data) %in% rows_with_zeros]
  }

  # Rownames of fractions contain original indices of row in input dataframe
  # so we can update the taxonomy for the original data (from sorted_sample_df)
  data = data[1:nrow(data) %in% rownames(fractions), ]

  # Update only fractions' counts
  data[, fraction_names] = fractions

  if (length(imputed_taxa) == 0) {
    imputed_taxa = NULL
  }

  return(list(data = data, imputed_taxa = imputed_taxa))
}

# Implements the idea of MA plot (from microarray analysis) to abundance analysis
