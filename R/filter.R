#' Get Valid Taxranks from Phyloseq Object
#' @export
valid_taxranks <- function(
  seq_obj # phyloseq or PhyloIgSeq object
) {
  if (class(seq_obj) == "phyloseq") {
    # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
    if (taxa_are_rows(seq_obj)) {
      seq_obj <- t(seq_obj)
    }

    tax.df <- as.data.frame(tax_table(seq_obj))
  } else if (class(seq_obj) == "PhyloIgSeq") {
    tax.df <- seq_obj@tax_table
  } else {
    return(NULL)
  }

  # Get taxranks having at least one unique non-na value
  # Count number of unique non-NA and non-empty values per column
  unique.counts <- sapply(tax.df, function(col) {
    length(unique(col[!is.na(col) & col != ""]))
  })

  # Keep only those columns with >1 unique valid entry
  valid.ranks <- names(unique.counts[unique.counts > 1])

  return(valid.ranks)
}

#' Plot Read Counts by Samples/Taxa from Phyloseq Object
#' @export
plotReads <- function(physeq, min.sample.sum = NA, min.taxa.sum = NA) {
  if (is.null(min.sample.sum)) {
    min.sample.sum <- NA
  }
  if (is.null(min.taxa.sum)) {
    min.taxa.sum <- NA
  }
  if (class(physeq) != "phyloseq") {
    return(NULL)
  }

  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  # Combine sample and taxa read counts into one data frame
  combined_df <- rbind(
    data.frame(
      TotalReads = sample_sums(physeq),
      Level = "Sample",
      Threshold = min.sample.sum
    ),
    data.frame(
      TotalReads = taxa_sums(physeq),
      Level = "Taxon",
      Threshold = min.taxa.sum
    )
  )

  # Facet plot by level (sample or taxon)
  plot <-
    ggplot(combined_df, aes(x = TotalReads)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black")
  if (all(!is.na(combined_df$Threshold))) {
    plot <- plot +
      geom_vline(
        data = combined_df,
        mapping = aes(xintercept = Threshold),
        color = "darkred",
        linetype = "dashed",
        size = 1
      )
  }
  plot <- plot +
    facet_wrap(. ~ Level, scales = "free") +
    xlab("Total Reads") +
    ylab("Count") +
    ggtitle("Distribution of Total Reads by Sample and Taxon") +
    theme_minimal() +
    ggplot2::theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.title = element_text(face = "bold", hjust = 0.5)
    )

  return(plot)
}

#' Filter Phyloseq Object by Sums of Reads by Samples and by Taxa
#' @export
filterReads <- function(
  physeq,
  min.sample.sum = 100,
  min.taxa.sum = 2,
  taxa.first = TRUE
) {
  if (class(physeq) != "phyloseq") {
    return(NULL)
  }
  sparse_input <- is(otu_table(physeq), "sparse_otu_table")
  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  if (any(is.null(c(min.sample.sum, min.taxa.sum)))) {
    return(physeq)
  }
  # /!\ Attention: global assignment
  if (taxa.first) {
    keep_taxa <<- taxa_sums(physeq) >= min.taxa.sum
    if (!any(keep_taxa)) {
      warning("All taxa filtered out\n")
      return(NULL)
    }
    physeq <- subset_taxa(
      physeq,
      keep_taxa
    )

    keep_samples <<- sample_sums(physeq) >= min.sample.sum
    if (!any(keep_samples)) {
      warning("All samples filtered out\n")
      return(NULL)
    }
    physeq <- subset_samples(
      physeq,
      keep_samples
    )
  } else {
    keep_samples <<- sample_sums(physeq) >= min.sample.sum
    if (!any(keep_samples)) {
      warning("All samples filtered out\n")
      return(NULL)
    }
    physeq <- subset_samples(
      physeq,
      keep_samples
    )

    keep_taxa <<- taxa_sums(physeq) >= min.taxa.sum
    if (!any(keep_taxa)) {
      warning("All taxa filtered out\n")
      return(NULL)
    }
    physeq <- subset_taxa(
      physeq,
      keep_taxa
    )
  }
  if (sparse_input) {
    physeq <- as_sparse_phyloseq(physeq)
  }
  return(physeq)
}


# TAKEN AS-IS FROM FEATURE SELECTOR
#' Generate a vector of strings, each element of which is a logical
#' expression corresponding to a filter passed to parameters as a list
#' @export
getFilterExpression <- function(
  data,
  filter.criteria # list of filter criteria (see example)
) {
  if (!is.null(filter.criteria) & !is.null(data)) {
    filter.expression <- sapply(1:length(filter.criteria), function(i) {
      # extract the criterion (one filter)
      crit <- filter.criteria[[i]]
      # If this criterion is numeric
      if (is.numeric(data[[crit$var]])) {
        # concatenate the values, separated by a comma
        # Otherwise, if not numeric
        values.str <- paste(crit$value, collapse = ",")
        # set operator to "%in_interval%" (customized one)
        operator <- "%in_interval%"
      } else {
        # concatenate the values, separated by a comma and surrounded by quote marks
        # since those elements must be considered as strings (note: FALSE == "FALSE")
        values.str <- paste(paste0("'", crit$value, "'"), collapse = ",")
        # set operator to %in%
        operator <- "%in%"
      }

      # Generate the complete expression, surround the variable name by backticks
      # to account for possible special characters in the name (e.g. space)
      if (i == 1) {
        # don't include the logic for the first filter (it doesn't make sense)
        paste(
          crit$include,
          paste0("`", crit$var, "`"),
          operator,
          "c(",
          values.str,
          ")"
        )
      } else {
        paste(
          crit$logic,
          crit$include,
          paste0("`", crit$var, "`"),
          operator,
          "c(",
          values.str,
          ")"
        )
      }
    })

    # Concatenate all filter expressions to get the full logical expression
    filter.expression <- paste(filter.expression, collapse = " ")
    return((filter.expression))
  } else {
    return(NULL)
  }
}


#' Filter Phyloseq Object's Sample Data based on a Filter Expression
#' @export
filterSampleData <- function(
  physeq, # phyloseq object to filter
  filter.criteria = NULL # list of filter criteria (see example)
) {
  if (class(physeq) != "phyloseq") {
    return(NULL)
  }
  sparse_input <- is(otu_table(physeq), "sparse_otu_table")
  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  data <- data.frame(
    sample_data(physeq),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  filter.expression <<- getFilterExpression(
    data = data,
    filter.criteria = filter.criteria
  )
  # If filter criteria are provided
  if (!is.null(filter.expression)) {
    # Otherwise, if no filter criteria is given, return the original data
    print(filter.expression) # /!\ temporary, to control the filters
    # Parse the resulting string to get the real logical expression and
    # filter the data
    filtered.physeq <- physeq %>%
      subset_samples(eval(parse(text = filter.expression)))
  } else {
    filtered.physeq <- physeq
  }
  if (sparse_input) {
    filtered.physeq <- as_sparse_phyloseq(filtered.physeq)
  }
  return(filtered.physeq)
}


#' Filter Phyloseq Object's Taxonomic Table based on a Filter Expression
#' @export
filterTaxTable <- function(
  physeq, # phyloseq object to filter
  filter.criteria = NULL # list of filter criteria (see example)
) {
  if (class(physeq) != "phyloseq") {
    return(NULL)
  }
  sparse_input <- is(otu_table(physeq), "sparse_otu_table")
  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  data <- data.frame(
    tax_table(physeq),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # /!\ Attention: global assignment
  filter.expression <<- getFilterExpression(
    data = data,
    filter.criteria = filter.criteria
  )
  # If filter criteria are provided
  if (!is.null(filter.expression)) {
    # Otherwise, if no filter criteria is given, return the original data
    print(filter.expression) # /!\ temporary, to control the filters
    # Parse the resulting string to get the real logical expression and
    # filter the data
    filtered.physeq <- physeq %>%
      subset_taxa(eval(parse(text = filter.expression)))
  } else {
    filtered.physeq <- physeq
  }
  if (sparse_input) {
    filtered.physeq <- as_sparse_phyloseq(filtered.physeq)
  }
  return(filtered.physeq)
}


#' Plot Phylogenetic Tree from Phyloseq Object
#' @export
plot_phylo_tree <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  circular = FALSE,
  label = NULL,
  label.levels = NULL,
  remove.na.from.plot = FALSE,
  ...
) {
  if (!is.null(physeq)) {
    if (class(physeq) != "phyloseq") {
      stop("Need a phyloseq object")
    }

    if (is.null(access(physeq, "phy_tree"))) {
      stop("Phyloseq object has to contain a tree")
    }

    if (!is.null(fraction_id_name) & !is.null(fraction_ids)) {
      physeq <- prune_samples(
        sample_data(physeq)[[fraction_id_name]] %in% fraction_ids,
        physeq
      )
    }

    if (!is.null(taxrank)) {
      physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
      taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank])
    }
  }

  if (
    !is.null(label) &&
      !is.null(label.levels) &&
      !is.numeric(sample_data(physeq)[[label]])
  ) {
    keep <- keep_levels(sample_data(physeq)[[label]], label.levels)
    samples.wo.na.global <<- keep
    physeq <- prune_samples(keep, physeq)
    sample_data(physeq)[[label]] <- factorize_levels(
      sample_data(physeq)[[label]],
      label.levels
    )
  } else if (
    remove.na.from.plot &&
      !is.null(label) &&
      !is.null(sample_data(physeq)[[label]])
  ) {
    samples.wo.na.global <<- !is.na(sample_data(physeq)[[label]])
    physeq <- physeq %>% subset_samples(samples.wo.na.global)
  }

  tree <- physeq %>% plot_tree(color = label, ...)

  if (circular) {
    tree <- tree + coord_polar(theta = "y")
  }

  return(tree)
}
