#' Filter Taxa and Samples by Minimum Read Sum
#'
#' Removes taxa with fewer than `min_taxa_sum` total reads and samples with
#' fewer than `min_sample_sum` total reads. Filtering order matters, since
#' each side's sums are recomputed on whatever subset the other side's filter
#' has already produced (see `taxa_first`).
#'
#' @param physeq A `phyloseq` object.
#' @param min_sample_sum Numeric threshold; samples with total reads below
#'   this value are removed. `NULL` skips read-count filtering entirely
#'   (`physeq` is returned unfiltered, only reoriented/re-sparsified).
#' @param min_taxa_sum Numeric threshold; taxa with total reads below this
#'   value are removed. `NULL` skips read-count filtering entirely, same as
#'   `min_sample_sum = NULL`.
#' @param taxa_first Logical. If `TRUE` (default), taxa are filtered before
#'   samples, so sample sums are computed on the taxon-filtered subset; if
#'   `FALSE`, samples are filtered first and taxon sums are computed on the
#'   sample-filtered subset.
#' @param drop_zero_taxa Logical (default `TRUE`). Whichever side is filtered
#'   last can zero out taxa on the other side (e.g. a taxon confined to a
#'   sample that gets dropped for low read count). When `TRUE`, a final pass
#'   removes any taxon left with zero total reads, regardless of
#'   `min_taxa_sum` or `taxa_first`. Applied even when `min_sample_sum` or
#'   `min_taxa_sum` is `NULL`, since it addresses taxa left empty by upstream
#'   sample filtering, not the threshold-based filtering this function does.
#'
#' @return The filtered `phyloseq` object. `NULL` if `physeq` is not a
#'   `phyloseq` object, or if filtering would remove every taxon or every
#'   sample (a `warning()` is issued in that case naming which side emptied
#'   out).
#' @export
#'
#' @examples
#' data(ps_16s_refinement)
#' filtered <- filter_reads(
#'   ps_16s_refinement,
#'   min_sample_sum = 100,
#'   min_taxa_sum = 5
#' )
#' phyloseq::ntaxa(filtered)
#' phyloseq::nsamples(filtered)
filter_reads <- function(
  physeq,
  min_sample_sum = 100,
  min_taxa_sum = 2,
  taxa_first = TRUE,
  drop_zero_taxa = TRUE
) {
  if (!is(physeq, "phyloseq")) {
    return(NULL)
  }
  sparse_input <- is(otu_table(physeq), "sparse_otu_table")
  # Force taxa to be columns of otu table - like in internal phyloseq function phyloseq:::veganifyOTU()
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }

  if (is.null(min_sample_sum) || is.null(min_taxa_sum)) {
    if (isTRUE(drop_zero_taxa)) {
      physeq <- .drop_zero_taxa(physeq)
      if (is.null(physeq)) {
        return(NULL)
      }
    }
    if (sparse_input) {
      physeq <- as_sparse_phyloseq(physeq)
    }
    return(physeq)
  }

  if (taxa_first) {
    physeq <- .filter_taxa_by_sum(physeq, min_taxa_sum)
    if (is.null(physeq)) {
      return(NULL)
    }
    physeq <- .filter_samples_by_sum(physeq, min_sample_sum)
  } else {
    physeq <- .filter_samples_by_sum(physeq, min_sample_sum)
    if (is.null(physeq)) {
      return(NULL)
    }
    physeq <- .filter_taxa_by_sum(physeq, min_taxa_sum)
  }
  if (is.null(physeq)) {
    return(NULL)
  }

  if (isTRUE(drop_zero_taxa)) {
    physeq <- .drop_zero_taxa(physeq)
    if (is.null(physeq)) {
      return(NULL)
    }
  }

  if (sparse_input) {
    physeq <- as_sparse_phyloseq(physeq)
  }
  return(physeq)
}

# Keeps taxa with taxa_sums(physeq) >= min_taxa_sum, warning + returning NULL
# if that empties the table. Uses prune_taxa() (accepts a logical vector
# directly) rather than subset_taxa()'s NSE, which cannot see a variable
# local to the caller's frame.
.filter_taxa_by_sum <- function(physeq, min_taxa_sum) {
  keep_taxa <- taxa_sums(physeq) >= min_taxa_sum
  if (!any(keep_taxa)) {
    warning("All taxa filtered out\n")
    return(NULL)
  }
  prune_taxa(keep_taxa, physeq)
}

# Sample-side counterpart to .filter_taxa_by_sum().
.filter_samples_by_sum <- function(physeq, min_sample_sum) {
  keep_samples <- sample_sums(physeq) >= min_sample_sum
  if (!any(keep_samples)) {
    warning("All samples filtered out\n")
    return(NULL)
  }
  prune_samples(keep_samples, physeq)
}

# Drops taxa with exactly zero reads across all (currently present) samples,
# e.g. ghost taxa left behind once the samples that carried them are removed.
# Unlike .filter_taxa_by_sum(), this always targets strict zeros regardless
# of min_taxa_sum, so taxa aren't left dangling just because the threshold
# happens to be 0 or unset.
.drop_zero_taxa <- function(physeq) {
  keep_taxa <- taxa_sums(physeq) > 0
  if (!any(keep_taxa)) {
    warning("All taxa filtered out\n")
    return(NULL)
  }
  if (all(keep_taxa)) {
    return(physeq)
  }
  prune_taxa(keep_taxa, physeq)
}

#' Plot Distributions of Total Read Counts by Sample and by Taxon
#'
#' Histograms of per-sample and per-taxon total read counts, faceted side by
#' side, with an optional vertical reference line marking a filtering
#' threshold on each facet.
#'
#' @param physeq A `phyloseq` object.
#' @param min_sample_sum Numeric threshold to draw as a vertical reference
#'   line on the sample-level histogram. `NA` (default) or `NULL` omits it.
#' @param min_taxa_sum Numeric threshold to draw as a vertical reference line
#'   on the taxon-level histogram. `NA` (default) or `NULL` omits it.
#' @param log_scale Logical; if `TRUE`, plot the x axis (total reads) on a
#'   log10 scale. Default `FALSE`.
#'
#' @return A `ggplot` object faceted by level (`"Sample"`/`"Taxon"`). `NULL`
#'   if `physeq` is not a `phyloseq` object.
#' @export
#'
#' @examples
#' data(ps_16s_refinement)
#' plot_reads(ps_16s_refinement, min_sample_sum = 100, min_taxa_sum = 5)
plot_reads <- function(
  physeq,
  min_sample_sum = NA,
  min_taxa_sum = NA,
  log_scale = FALSE
) {
  if (is.null(min_sample_sum)) {
    min_sample_sum <- NA
  }
  if (is.null(min_taxa_sum)) {
    min_taxa_sum <- NA
  }
  if (!is(physeq, "phyloseq")) {
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
      Threshold = min_sample_sum
    ),
    data.frame(
      TotalReads = taxa_sums(physeq),
      Level = "Taxon",
      Threshold = min_taxa_sum
    )
  )

  # Facet plot by level (sample or taxon)
  plot <-
    ggplot(combined_df, aes(x = TotalReads)) +
    geom_histogram(bins = 30, fill = "#4C6A92", color = "black")
  if (all(!is.na(combined_df$Threshold))) {
    plot <- plot +
      geom_vline(
        data = combined_df,
        mapping = aes(xintercept = Threshold),
        color = "darkred",
        linetype = "dashed",
        linewidth = 1
      )
  }
  if (isTRUE(log_scale)) {
    plot <- plot + scale_x_log10()
  }
  plot <- plot +
    facet_wrap(. ~ Level, scales = "free") +
    xlab("Total Reads") +
    ylab("Count") +
    ggtitle("Distribution of Total Reads by Sample and Taxon") +
    .plot_title_theme()

  return(plot)
}
