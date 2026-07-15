#' PhyloIgSeq
#'
#' Package containing wrappers mainly around \pkg{phyloseq} and \pkg{vegan} packages
#' for downstream analysis of 16S amplicon sequencing results and allowing to further
#' analyze the Ig-Seq data, computing and plotting various Ig scores. The functions
#' from this package are used in PhyloIgSeq Shiny Web application on
#'  \url{https://www.immulab.fr/cms/index.php/team/tools/lab-tools}.
#'
#' @import openxlsx
#' @import DT
#' @import ade4
#' @import vegan
#' @import Rtsne
#' @import umap
#' @import VIM
#' @import dplyr
#' @import tidyr
#' @import ggpubr
#' @import tools
#' @import htmlwidgets
#' @import rstatix
#' @import scales
#' @import phyloseq
#' @import ggplot2
#' @import methods
#' @importFrom microbiome transform
#' @import zCompositions
#' @importFrom magrittr %>%
#' @importFrom gtools mixedsort
#' @import rlang
#' @import ggExtra
#' @importFrom car dataEllipse
#' @import sp
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation anno_text
#' @importFrom grid gpar unit
#' @import ragg
#' @import gganimate
#' @import gifski
#' @importFrom utils globalVariables
"_PACKAGE"

#' Optimal `facet_wrap()` Column Count
#'
#' Internal. Picks the number of columns for [ggplot2::facet_wrap()] that
#' minimizes empty cells, preferring wider (more columns) layouts over
#' taller ones when several column counts tie. Avoids degenerate
#' single-row/single-column results for `n > 3`.
#'
#' @param n Integer. Number of facets/panels to lay out.
#'
#' @return Integer. Recommended `ncol` for `facet_wrap()`.
#'
#' @noRd
smart_facet_ncol <- function(n) {
  if (n <= 3) {
    return(n)
  }
  divisors <- which(n %% seq_len(n) == 0)
  perfect_wide <- divisors[divisors >= sqrt(n) & divisors <= ceiling(n / 2)]
  if (length(perfect_wide) > 0) {
    return(min(perfect_wide))
  }
  nc_range <- seq(ceiling(sqrt(n)), ceiling(n / 2))
  empties <- ceiling(n / nc_range) * nc_range - n
  max(nc_range[empties == min(empties)])
}


#' Ig-Coating Score Names
#'
#' Character vector of the Ig-coating score names computable via \code{\link{compute_ig_score}}
#' and requestable through \code{\link{getPhyloIgSeq}}'s \code{scores} argument.
#'
#' @examples
#' IG_SCORES
#'
#' @export
IG_SCORES <- c("slide_z", "palm", "kau", "prob_index", "prob_ratio")
# scores to come:
# "purity_corrected_prob_index", "purity_corrected_prob_ratio")

# Prevent R CMD check NOTE about undefined global variable
utils::globalVariables("IG_SCORES")

#' Reorient a Phyloseq OTU Table to Taxa-as-Columns
#'
#' Internal. Similar to `phyloseq`'s (unexported) `veganifyOTU()`. Transposes `physeq`
#' so taxa are columns and samples are rows if not already, otherwise returns it
#' unchanged.
#'
#' @noRd
reverseASV <- function(physeq) {
  if (taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  return(physeq)
}

#' `ggplot2::geom_jitter()` with Vertical Jitter Off by Default
#'
#' Internal. Shadows [ggplot2::geom_jitter()] within this package's namespace
#' (unqualified `geom_jitter()` calls elsewhere in `R/*.R` resolve here, not to
#' `ggplot2`) to default `height = 0`, since y-axis jitter is rarely wanted for
#' this package's plots and `ggplot2::geom_jitter()`'s own default jitters both
#' axes. Falls through to `ggplot2::geom_jitter(...)` untouched if `position` is
#' explicitly supplied, since `height` conflicts with an explicit `position`.
#'
#' @noRd
geom_jitter <- function(..., height = 0) {
  if (hasArg(position)) {
    ggplot2::geom_jitter(...)
  } else {
    ggplot2::geom_jitter(..., height = height)
  }
}

#' Is `x` Within `interval`?
#'
#' Internal. `TRUE`/`FALSE` for whether `x` falls within the closed range
#' spanned by `interval` (`interval` need not be sorted; it is normalized via
#' [range()]). E.g. `3.5 %in_interval% c(3, 4)` is `TRUE`. Referenced by name
#' as a string from `filter.R`'s `eval(parse(...))`-based filter DSL, so it
#' must stay resolvable as `PhyloIgSeq:::\`%in_interval%\`` even though it is
#' not exported.
#'
#' @noRd
`%in_interval%` <- function(x, interval) {
  interval <- range(interval, na.rm = TRUE)
  x >= interval[1] & x <= interval[2]
}

#' Transform an Abundance Table
#'
#' Internal. Applies a compositional or Hellinger transform to `abundance_table`,
#' column-wise (i.e. per sample).
#'
#' @param abundance_table A numeric matrix of abundances.
#' @param transform Character. `"compositional"` (default; relative abundance,
#'   each sample's values divided by its total) or `"hellinger"` (square root
#'   of the compositional transform).
#' @param taxa_are_rows Logical. Whether `abundance_table` has taxa as rows
#'   (`TRUE`, default) or as columns (`FALSE`). The transform is always
#'   applied column-wise internally; `abundance_table` is transposed before
#'   and after when `FALSE`.
#'
#' @return A numeric matrix the same shape and orientation as `abundance_table`.
#'
#' @noRd
transform_abundances <- function(
  abundance_table,
  transform = c("compositional", "hellinger"),
  taxa_are_rows = TRUE
) {
  transform <- match.arg(transform)

  # assumes taxa are rows
  if (!taxa_are_rows) {
    abundance_table <- t(abundance_table)
  }

  abundance_table <- switch(
    transform,
    compositional = sweep(
      abundance_table,
      2,
      colSums(abundance_table),
      "/"
    ),
    hellinger = sqrt(sweep(
      abundance_table,
      2,
      colSums(abundance_table),
      "/"
    ))
  )

  if (!taxa_are_rows) {
    return(t(abundance_table))
  } else {
    return(abundance_table)
  }
}

#' Plot Rarefaction Curves
#'
#' Computes and plots rarefaction curves for all samples in a \code{phyloseq}
#' object using \code{\link[vegan]{rarecurve}}.
#'
#' @param ps A \code{\link[phyloseq]{phyloseq}} object with an OTU table and
#'   sample data. Samples with zero total counts are automatically excluded.
#' @param step Integer. Sampling-depth increment passed to
#'   \code{\link[vegan]{rarecurve}}. Smaller values yield smoother curves at
#'   the cost of computation time. Default \code{100}.
#' @param show_legend Logical. Whether to display the sample legend.
#'   Default \code{TRUE}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object with one rarefaction curve
#'   per sample, coloured by sample name.
#'
#' @seealso \code{\link[vegan]{rarecurve}}, \code{\link{plot_seq_depth}}
#'
#' @examples
#' data(ps_16s_refinement)
#' plot_rarefaction(ps_16s_refinement, step = 200)
#'
#' @export
plot_rarefaction <- function(ps, step = 100, show_legend = TRUE) {
  # Extract OTU table
  otu <- as(otu_table(ps), "matrix")

  # Ensure samples are rows
  if (taxa_are_rows(ps)) {
    otu <- t(otu)
  }

  # Remove samples with zero counts
  otu <- otu[rowSums(otu) > 0, ]

  # Compute rarefaction curves. col/lty are passed explicitly (unused by the
  # tidy = TRUE code path) to stop rarecurve() from calling par("col")/
  # par("lty") internally, which otherwise opens a graphics device as a side
  # effect even though nothing is actually plotted here.
  rare_df <- vegan::rarecurve(
    otu,
    step = step,
    xlab = "Reads",
    ylab = "Richness",
    tidy = TRUE,
    label = FALSE,
    col = 1,
    lty = 1
  )
  rare_df <- rare_df %>%
    rename(Reads = Sample, Richness = Species, Sample = Site)

  # Plot
  p <- ggplot(
    rare_df,
    aes(x = Reads, y = Richness, group = Sample, color = Sample)
  ) +
    geom_line(alpha = 0.8) +
    labs(
      x = "Sequencing depth",
      y = "Observed richness",
      title = "Rarefaction Curves"
    )

  p <- p +
    theme_minimal() +
    ggplot2::theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.title = if (show_legend) {
        element_text(face = "bold", hjust = 0.5)
      },
      legend.position = if (show_legend) {
        "right"
      } else {
        "none"
      }
    )

  return(p)
}

#' Plot Sequencing Depth
#'
#' Visualises per-sample sequencing depth from a \code{\link[phyloseq]{phyloseq}}
#' object as either a bar chart (one bar per sample) or a grouped boxplot with
#' jittered points.
#'
#' @param ps A \code{\link[phyloseq]{phyloseq}} object with an OTU table and
#'   sample data.
#' @param type Character. Plot type: \code{"bar"} (default) for a bar chart
#'   ordered by sample name, or \code{"box"} for a boxplot grouped by a
#'   metadata variable.
#' @param x_var Character. Name of a sample-data column to use as the x-axis
#'   grouping variable. Required when \code{type = "box"}, ignored otherwise.
#' @param facet_var Character or \code{NULL}. Name of a sample-data column used
#'   to facet the plot with \code{\link[ggplot2]{facet_wrap}}. Only applied
#'   when \code{type = "box"}. Default \code{NULL} (no faceting).
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso \code{\link{plot_rarefaction}}
#'
#' @examples
#' data(ps_16s_refinement)
#' plot_seq_depth(ps_16s_refinement)
#' plot_seq_depth(ps_16s_refinement, type = "box", x_var = "Protocol")
#'
#' @export
plot_seq_depth <- function(
  ps,
  type = c("bar", "box"),
  x_var = NULL,
  facet_var = NULL
) {
  type <- match.arg(type)

  # Compute sequencing depth
  depth_df <- data.frame(
    Sample = sample_names(ps),
    Depth = sample_sums(ps)
  )

  # Add metadata
  meta <- as.data.frame(sample_data(ps))
  meta$Sample <- rownames(meta)

  depth_df <- left_join(depth_df, meta, by = "Sample")

  # -----------------------
  # BARPLOT
  # -----------------------
  if (type == "bar") {
    p <- ggplot(depth_df, aes(x = Sample, y = Depth)) +
      geom_bar(stat = "identity") +
      theme_bw() +
      labs(x = "Sample", y = "Sequencing depth") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  } else {
    # -----------------------
    # BOXPLOT + JITTER
    # -----------------------
    if (is.null(x_var)) {
      stop("For boxplot, please provide x_var (independent variable).")
    }

    p <- ggplot(depth_df, aes(x = .data[[x_var]], y = Depth)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.2, alpha = 0.7) +
      theme_bw() +
      labs(x = x_var, y = "Sequencing depth")

    if (!is.null(facet_var)) {
      p <- p + facet_wrap(as.formula(paste("~", facet_var)))
    }
  }

  return(p)
}

#' Rarefy Abundances to a Common Depth by Multinomial Resampling
#'
#' Subsamples each sample's counts down to a common total (\code{common_count_sum})
#' by drawing from a multinomial distribution over taxa, so that per-sample totals
#' become directly comparable. Samples below the target depth are optionally
#' dropped rather than up-sampled.
#'
#' @param abundance_table A numeric matrix of counts, or a
#'   \code{\link[phyloseq]{otu_table}}, oriented per \code{taxa_are_rows}.
#' @param taxa_are_rows Logical. Whether \code{abundance_table} has taxa as rows
#'   and samples as columns (\code{TRUE}, default) or the reverse (\code{FALSE}).
#' @param common_count_sum Integer. Target depth every retained sample is
#'   resampled to. Defaults to the smallest nonzero sample total in
#'   \code{abundance_table}.
#' @param trim_taxa Logical. Drop taxa left with zero counts across all samples
#'   after rarefaction. Default \code{TRUE}.
#' @param trim_samples Logical. Drop samples whose total is below
#'   \code{common_count_sum} before rarefying (such samples cannot be rarefied
#'   up to that depth). Default \code{TRUE}.
#' @param silent_warnings Logical. Suppress the warning reporting how many
#'   samples/taxa were trimmed. Default \code{FALSE}.
#'
#' @return A numeric matrix, oriented as \code{abundance_table}, with sample
#'   totals rarefied to \code{common_count_sum} (and samples/taxa dropped per
#'   \code{trim_samples}/\code{trim_taxa}).
#'
#' @examples
#' data(ps_16s_refinement)
#' otu <- as(phyloseq::otu_table(ps_16s_refinement), "matrix")
#' rarefied <- rarefy_abundances(
#'   otu,
#'   taxa_are_rows = phyloseq::taxa_are_rows(ps_16s_refinement),
#'   silent_warnings = TRUE
#' )
#' rowSums(rarefied)
#'
#' @export
rarefy_abundances <-
  function(
    abundance_table,
    taxa_are_rows = TRUE,
    common_count_sum = NULL,
    trim_taxa = TRUE,
    trim_samples = TRUE,
    silent_warnings = FALSE
  ) {
    if (taxa_are_rows) {
      sample_margin <- 2
      taxa_margin <- 1
      sample_sums_fn <- colSums
      taxa_sums_fn <- rowSums
    } else {
      sample_margin <- 1
      taxa_margin <- 2
      sample_sums_fn <- rowSums
      taxa_sums_fn <- colSums
    }

    if (is.null(common_count_sum)) {
      common_count_sum <- min(setdiff(sample_sums_fn(abundance_table), 0))
    }

    dims_orig <- dim(abundance_table)

    if (trim_samples) {
      if (taxa_are_rows) {
        abundance_table <- abundance_table[,
          sample_sums_fn(abundance_table) >= common_count_sum,
          drop = FALSE
        ]
      } else {
        abundance_table <- abundance_table[
          sample_sums_fn(abundance_table) >= common_count_sum, ,
          drop = FALSE
        ]
      }
    }

    # For each sample, draw counts from a multinomial distribution over taxa
    abundance_table <-
      apply(abundance_table, sample_margin, function(sample_counts) {
        if (sum(sample_counts) == 0) {
          return(rep(0, length(sample_counts)))
        }
        rmultinom(
          1,
          size = common_count_sum,
          prob = sample_counts / sum(sample_counts)
        )[, 1]
      })

    if (sample_margin == 1) {
      # apply() to rows puts them in columns, so transpose the matrix back
      abundance_table <- t(abundance_table)
    }

    if (trim_taxa) {
      if (taxa_are_rows) {
        abundance_table <- abundance_table[
          taxa_sums_fn(abundance_table) > 0, ,
          drop = FALSE
        ]
      } else {
        abundance_table <- abundance_table[,
          taxa_sums_fn(abundance_table) > 0,
          drop = FALSE
        ]
      }
    }

    dims_rare <- dim(abundance_table)

    if (!all(dims_rare == dims_orig) & !silent_warnings) {
      warning(paste0(
        "Trimmed ",
        dims_orig[sample_margin] - dims_rare[sample_margin],
        " samples (with count sum < common_count_sum) and ",
        dims_orig[taxa_margin] - dims_rare[taxa_margin],
        " taxa (zero counts after rarefaction)",
        "\n"
      ))
    }

    return(abundance_table)
  }

#' Check Whether Abundances Look Like Counts
#'
#' Heuristically tests whether `x` holds non-negative integer (count-like) values, as
#' opposed to already-transformed (relative abundance, log, etc.) data. For large
#' inputs, only a subsample of values is checked by default (see `consider_small_part`).
#'
#' @param x A numeric matrix, `phyloseq::otu_table`, or [sparse_otu_table-class].
#' @param allow_na Logical. If `TRUE` (default), `NA` values are ignored. If `FALSE`,
#'   any `NA` makes the function return `FALSE`.
#' @param consider_small_part Logical. If `TRUE` (default), only a subsample of `x` is
#'   checked (the first 10000 stored values for a [sparse_otu_table-class], or the
#'   top-left 100x100 corner otherwise) rather than every value, for speed on large
#'   tables. Set `FALSE` to check every value.
#'
#' @return A single logical: `TRUE` if the checked values are all non-negative, finite
#'   integers (`NA`s permitted/ignored per `allow_na`), `FALSE` otherwise (including
#'   when `x` is not numeric).
#'
#' @examples
#' is_count_like(matrix(1:6, nrow = 2))
#' is_count_like(matrix(c(0.1, 0.2, 0.3, 0.4), nrow = 2))
#' is_count_like(matrix(c(1, NA, 3, 4), nrow = 2), allow_na = FALSE)
#'
#' @export
is_count_like <- function(x, allow_na = TRUE, consider_small_part = TRUE) {
  if (!is.numeric(x)) {
    return(FALSE)
  }

  if (is(x, "sparse_otu_table")) {
    # Only check stored (non-zero) values; zero is trivially count-like.
    vals <- x@sparse_data@x
    if (consider_small_part) vals <- head(vals, 10000L)
  } else {
    # Covers both otu_table and any subclass that is not sparse_otu_table.
    if (is(x, "otu_table")) {
      x <- as(x, "matrix")
    }
    if (consider_small_part) {
      x <- x[seq_len(min(100L, nrow(x))), seq_len(min(100L, ncol(x)))]
    }
    vals <- c(x)
  }

  if (allow_na) {
    vals <- vals[!is.na(vals)]
  } else if (anyNA(vals)) {
    return(FALSE)
  }

  all(vals >= 0 & is.finite(vals) & vals == floor(vals))
}

#' Make Taxonomy Names Unique At Every Taxonomic Rank
#'
#' Disambiguates taxa whose name at a given rank is a duplicate of another taxon's name
#' at that rank *without* sharing the same higher-rank lineage ("bad duplication", see
#' Details) by appending a `make.unique()`-style suffix (`.1`, `.2`, ...). Duplicate
#' names that do share the same full higher-rank lineage ("good duplication") are left
#' untouched. `NA` entries are replaced with the literal string `"NA"` first, so they
#' participate in the same duplication logic as any other value.
#'
#' @details
#' "Good duplication": low-rank names repeat because the whole lineage above them is
#' identical (e.g. two ASVs both classified as genus `Bacteroides`, species `caccae`) —
#' left unchanged, since collapsing them would lose real biological information (e.g.
#' via [tax_glom()]).
#'
#' "Bad duplication": the same low-rank name is reached via two different higher-rank
#' lineages (e.g. species `caccae` under both genus `Bacteroides` and genus
#' `Anaerostipes`) — these represent two different taxa that happen to share a species
#' name, so the second occurrence is renamed (e.g. `caccae.1`) to keep them distinguishable.
#'
#' @param taxa_table A taxonomy matrix/data frame (e.g. a `phyloseq::tax_table`), ranks
#'   as columns ordered from highest (leftmost) to lowest (rightmost), taxa as rows.
#'
#' @return `taxa_table` with `NA`s replaced by `"NA"` and "bad duplication" names at
#'   every rank column made unique.
#'
#' @examples
#' taxa_table <- data.frame(
#'   Genus = c("Bacteroides", "Bacteroides", "Anaerostipes"),
#'   Species = c("caccae", "caccae", "caccae")
#' )
#' make_unique_taxa_table(taxa_table)
#'
#' @export
make_unique_taxa_table <- function(taxa_table) {
  taxa_table[is.na(taxa_table)] <- "NA"

  # "good duplication" example:
  # low level duplications correspond to whole taxonomy duplications
  #  ASV   | ... |    Genus     | Species
  # __________________________________
  #  ASV1  | ... | Bacteroides  | caccae
  #  ASV2  | ... | Bacteroides  | caccae

  # "bad duplication" example:
  #  ASV   | ... |    Genus     | Species
  # ___________________________________
  #  ASV1  | ... | Bacteroides  | caccae
  #  ASV2  | ... | Anaerostipes | caccae

  # Will be changed to :
  #  ASV   | ... |    Genus     | Species
  # ___________________________________
  #  ASV1  | ... | Bacteroides  | caccae
  #  ASV2  | ... | Anaerostipes | caccae.1

  # Go down taxonomy hierarchy (assuming it goes from left to right in the table)
  # the topmost level of taxonomy has only "good" duplications
  for (i in 2:ncol(taxa_table)) {
    # taxonomy if defined by the current level + all higher levels
    # (here, by just one higher level since it's ensured to be unique)
    taxonomies <- apply(taxa_table[, (i - 1):i], 1, function(row) {
      paste(row, collapse = "#")
    })
    # get rid of "good duplications"
    unique_taxonomies <- unique(taxonomies)

    # identify "bad duplications" and make them unique
    duplicated_names <- sapply(
      strsplit(unique_taxonomies, "#"),
      function(taxons) {
        tail(taxons, 1)
      }
    )
    unique_names <- make.unique(duplicated_names)

    # skip if no "bad duplications" - most often the case
    if (all(unique_names == duplicated_names)) {
      next
    }

    # replace "bad duplications" (vectorized: a loop doing
    # `taxonomies == unique_taxonomies[j]` per name is O(n * length(unique_names))
    # and can take a minute-plus on large tables)
    taxonomies <- unique_names[match(taxonomies, unique_taxonomies)]

    taxa_table[, i] <- taxonomies
  }

  return(taxa_table)
}

#' Impute Missing Values With A Column's Central Tendency
#'
#' Fills `NA`s column-wise: numeric columns with their mean/median/mode (per
#' `central_tendency`), factor columns with their mode. Character columns are left
#' untouched.
#'
#' @param df A data frame (or matrix, coerced via `as.data.frame()`).
#' @param central_tendency One of `"median"` (default), `"mean"`, or `"mode"`; only
#'   applies to numeric columns (factor columns always use the mode).
#'
#' @return `df` (as a data frame) with `NA`s in numeric/factor columns imputed.
#'
#' @examples
#' df <- data.frame(a = c(1, NA, 3, 4), b = c(10, 20, NA, 40))
#' impute_with_central_tendency(df, central_tendency = "mean")
#'
#' @export
impute_with_central_tendency <- function(
  df,
  central_tendency = c("median", "mean", "mode")
) {
  central_tendency <- match.arg(central_tendency)
  df <- as.data.frame(df)
  get_mode <- function(v) {
    v <- v[!is.na(v)]
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }

  for (col in colnames(df)) {
    if (is.numeric(df[[col]])) {
      df[[col]][is.na(df[[col]])] <- switch(
        central_tendency,
        "mean" = mean(df[[col]], na.rm = TRUE),
        "median" = median(df[[col]], na.rm = TRUE),
        "mode" = get_mode(df[[col]])
      )
      # Impute categorical columns, character columns are not affected /!\
    } else if (is.factor(df[[col]])) {
      df[[col]][is.na(df[[col]])] <- get_mode(df[[col]])
    }
  }
  return(df)
}


#' Impute Data With KNN, Central Tendency, Or Zero
#'
#' Imputes `NA`s in `data_tmp` (excluding any `exceptions` columns, which are passed
#' through untouched) using one of three methods.
#'
#' @param data_tmp A data frame (or matrix) with `NA`s to impute.
#' @param exceptions Character vector of `data_tmp` column names to exclude from
#'   imputation and pass through unchanged. Defaults to `NULL` (no exceptions).
#' @param method One of `"KNN"` (default; [VIM::kNN()]), `"Central Tendency"`
#'   ([impute_with_central_tendency()]), or `"Replace NA with 0"`.
#' @param central_tendency Passed to [impute_with_central_tendency()] when
#'   `method = "Central Tendency"`: one of `"mean"`, `"median"`, or `"mode"`.
#'   Ignored otherwise.
#' @param nb_neighbors Number of neighbors passed to [VIM::kNN()]'s `k` when
#'   `method = "KNN"`. A common rule-of-thumb range is 3-10. Ignored otherwise.
#' @param add_imputation_indicators Logical, passed to [VIM::kNN()]'s `imp_var` when
#'   `method = "KNN"`: whether to add a `TRUE`/`FALSE` indicator column per imputed
#'   variable. Ignored otherwise.
#'
#' @return `data_tmp` with `NA`s in its non-`exceptions` columns imputed (`exceptions`
#'   columns first, in their original order, followed by the imputed columns).
#'
#' @examples
#' df <- data.frame(id = 1:4, a = c(1, NA, 3, 4), b = c(10, 20, NA, 40))
#' dataImpute(df, exceptions = "id", method = "Replace NA with 0")
#'
#' @export
dataImpute <- function(
  data_tmp,
  exceptions = NULL,
  method = c("KNN", "Central Tendency", "Replace NA with 0"),
  central_tendency = "median",
  nb_neighbors = 5, # /!\ Find optimal usually it's between 3 and 10
  add_imputation_indicators = FALSE
) {
  method <- match.arg(method)
  data_to_impute <- data_tmp[, !colnames(data_tmp) %in% exceptions, drop = FALSE]

  if (method == "KNN") {
    all_na_cols <- colnames(data_to_impute)[
      apply(data_to_impute, 2, function(col) all(is.na(col)))
    ]
    if (length(all_na_cols) > 0) {
      stop(
        "dataImpute: method = \"KNN\" cannot impute column(s) with no observed ",
        "values: ", paste(all_na_cols, collapse = ", "),
        ". Drop them first or choose a different `method`.",
        call. = FALSE
      )
    }
  }

  data_imputed <- switch(
    method,
    "KNN" = VIM::kNN(
      data_to_impute,
      k = nb_neighbors,
      imp_var = add_imputation_indicators
    ),
    "Central Tendency" = impute_with_central_tendency(
      data_to_impute,
      central_tendency = central_tendency
    ),
    "Replace NA with 0" = {
      data_to_impute[is.na(data_to_impute)] <- 0
      data_to_impute
    }
  )

  cbind(
    data_tmp[, colnames(data_tmp) %in% exceptions, drop = FALSE],
    data_imputed
  )
}

#' Plot Phylogenetic Tree from a Phyloseq Object
#'
#' Draws a [ggtree::ggtree()] phylogenetic tree from a `phyloseq` object's `phy_tree`
#' slot, optionally agglomerating taxa to a given rank first, restricting to a subset
#' of sort fractions, labeling tips, and coloring tips by a taxonomy column.
#'
#' @param physeq A `phyloseq` object containing a `phy_tree`.
#' @param taxrank `NULL` (default, no agglomeration) or a single taxonomic rank name
#'   present in `tax_table(physeq)` (e.g. `"Genus"`). If supplied, `physeq` is
#'   agglomerated to that rank via [tax_glom()] before plotting, and taxa are renamed
#'   to their rank value (de-duplicated with [make.unique()]).
#' @param fraction_id_name `NULL` (default) or the name of a `sample_data(physeq)`
#'   column identifying sort fractions (e.g. `"sorting_fraction"`). Used together with
#'   `fraction_ids` to restrict the tree to a subset of samples before plotting.
#' @param fraction_ids `NULL` (default) or a character vector of values of
#'   `fraction_id_name` to keep (via [phyloseq::prune_samples()]). Ignored unless
#'   `fraction_id_name` is also supplied.
#' @param layout Tree layout passed to [ggtree::ggtree()] (e.g. `"rectangular"`,
#'   `"circular"`, `"fan"`). Default `"rectangular"`.
#' @param tip_color `NULL` (default, no tip coloring) or the name of a
#'   `tax_table(physeq)` column to color tips by (via [ggtree::geom_tippoint()] and
#'   [ggsci::scale_color_igv()]). Silently ignored if the named column doesn't exist.
#' @param label_tips Logical, default `FALSE`. If `TRUE`, tip labels are drawn (via
#'   [ggtree::geom_tiplab()]); tip names are first truncated to 25 characters (with
#'   [make.unique()] to keep them distinct) so long ASV/OTU hashes don't clutter the
#'   plot.
#' @param label_size Tip label font size, passed to [ggtree::geom_tiplab()]. Only used
#'   when `label_tips = TRUE`. Default `2.5`.
#' @param ladderize One of `"left"`, `"right"`, or `NULL` (leave the tree's existing
#'   node order untouched). Default `"left"`. Passed through to [ggtree::ggtree()]'s
#'   own `ladderize`/`right` arguments.
#' @param ... Additional arguments passed on to [ggtree::ggtree()].
#'
#' @return A `ggtree`/`ggplot` object.
#'
#' @details If `phy_tree(physeq)` is unrooted, it is midpoint-rooted (via
#'   [phytools::midpoint_root()], with a `warning()`) before plotting, mirroring
#'   [sparse_unifrac()]'s handling of unrooted input. This matters for layouts that
#'   visually imply a root-to-tip direction (`"rectangular"`, `"circular"`, `"fan"`,
#'   ...); it is not relevant when using an explicitly unrooted `layout` such as
#'   `"equal_angle"`/`"daylight"`.
#'
#' @examples
#' data(ps_16s_refinement)
#' ps_phylum <- tax_glom(ps_16s_refinement, taxrank = "Phylum")
#' plot_phylo_tree(ps_phylum, layout = "circular", tip_color = "Phylum")
#' plot_phylo_tree(ps_phylum, label_tips = TRUE)
#'
#' # ps_16s_refinement's own tree is unrooted, so this midpoint-roots it first
#' # (with a warning)
#' plot_phylo_tree(ps_16s_refinement)
#'
#' @export
plot_phylo_tree <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  layout = "rectangular",
  tip_color = NULL,
  label_tips = FALSE,
  label_size = 2.5,
  ladderize = "left",
  ...
) {
  if (is.null(physeq) || !is(physeq, "phyloseq")) {
    stop("Need a phyloseq object")
  }
  if (is.null(access(physeq, "phy_tree"))) {
    stop("Phyloseq object has to contain a tree")
  }
  if (!is.null(ladderize) && !ladderize %in% c("left", "right")) {
    stop('`ladderize` must be NULL, "left", or "right"')
  }

  if (!is.null(fraction_id_name) && !is.null(fraction_ids)) {
    physeq <- prune_samples(
      sample_data(physeq)[[fraction_id_name]] %in% fraction_ids,
      physeq
    )
  }

  if (!is.null(taxrank)) {
    physeq <- tax_glom(physeq = physeq, taxrank = taxrank)
    taxa_names(physeq) <- make.unique(tax_table(physeq)[, taxrank])
  }

  # Truncate tip labels so long ASV hashes don't clutter the plot;
  # make.unique() prevents duplicate names after truncation
  if (isTRUE(label_tips)) {
    taxa_names(physeq) <- make.unique(substr(taxa_names(physeq), 1, 25))
  }

  tree_obj <- phy_tree(physeq)
  if (!ape::is.rooted(tree_obj)) {
    tree_obj <- phytools::midpoint_root(tree_obj)
    warning("Tree is unrooted, midpoint was set as root")
  }

  p <- ggtree::ggtree(
    tree_obj,
    layout = layout,
    ladderize = !is.null(ladderize),
    right = identical(ladderize, "right"),
    ...
  )

  if (isTRUE(label_tips)) {
    p <- p + ggtree::geom_tiplab(size = label_size, align = FALSE)
  }

  if (!is.null(tip_color) && tip_color %in% colnames(tax_table(physeq))) {
    tax_df <- as.data.frame(tax_table(physeq))
    tax_df$.join_label <- taxa_names(physeq)
    # %<+% can fail when p$data has a duplicate "label" column in some ggtree
    # versions; merge directly into p$data to avoid the conflict
    p$data <- dplyr::left_join(p$data, tax_df, by = c("label" = ".join_label"))
    p <- p +
      ggtree::geom_tippoint(
        ggplot2::aes(color = .data[[tip_color]]),
        size = 2
      ) +
      ggsci::scale_color_igv(na.value = "grey50", name = tip_color)
  }

  p
}
