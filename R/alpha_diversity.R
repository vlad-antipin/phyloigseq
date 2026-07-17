# Naming convention used across this file (and beta_diversity.R):
#   get_foo()  computes data for an analysis and returns a structured list
#   plot_foo() renders a ggplot from get_foo()'s output
#   foo()      thin convenience wrapper chaining get_foo() + plot_foo()

#' Sparse-aware Shannon alpha diversity
#'
#' Computes Shannon entropy per sample directly from the sparse matrix slot,
#' without materialising the full dense OTU table.
#'
#' @param physeq A \code{phyloseq} object whose OTU table is a
#'   \code{\link{sparse_otu_table-class}}.
#' @return A named numeric vector of Shannon entropy values (nats), one per sample.
#'
#' @examples
#' data(ps_16s_refinement)
#' ps_sparse <- as_sparse_phyloseq(ps_16s_refinement)
#' sparse_shannon(ps_sparse)
#'
#' @export
sparse_shannon <- function(physeq) {
  ot <- phyloseq::otu_table(physeq)
  sp <- ot@sparse_data
  if (phyloseq::taxa_are_rows(ot)) {
    sp <- Matrix::t(sp)
  }
  # sp: samples x taxa (rows = samples)
  N <- Matrix::rowSums(sp)
  sp_norm <- Matrix::Diagonal(x = 1 / N) %*% sp
  log_nz <- sp_norm
  log_nz@x <- log(sp_norm@x)
  h <- -Matrix::rowSums(sp_norm * log_nz)
  setNames(as.numeric(h), rownames(sp))
}

#' Get Alpha Diversity from a Phyloseq Object
#'
#' Computes one or more alpha-diversity measures per sample. When
#' \code{from_igseq = TRUE}, computes Ig-score significance-bucket richness
#' via [get_igseq_richness()] instead. Identity (raw counts) is a safe
#' transform for every measure; Chao1, ACE and Fisher additionally require
#' count data and are not meaningful on already-transformed (e.g.
#' compositional) abundances.
#'
#' @param physeq A \code{phyloseq} object.
#' @param from_igseq Logical. If \code{TRUE}, compute Ig-score
#'   significance-bucket richness via [get_igseq_richness()] instead of a
#'   standard diversity index; \code{proportions}/\code{low_lim}/
#'   \code{high_lim} apply, \code{transform_abundances}/\code{taxrank}/
#'   \code{measure}/\code{fraction_id_name}/\code{fraction_ids} are ignored.
#'   Default \code{FALSE}.
#' @param transform_abundances \code{NULL} or a transform name passed to
#'   [microbiome::transform()] (e.g. \code{"compositional"});
#'   \code{"identity"} (default) applies no transform.
#' @param proportions Logical, only used when \code{from_igseq = TRUE}:
#'   report bucket proportions instead of raw counts. Default \code{FALSE}.
#' @param low_lim,high_lim Numeric Ig-score cutoffs defining the "down"/"up"
#'   significance buckets when \code{from_igseq = TRUE}: taxa with a score
#'   below \code{low_lim}/above \code{high_lim} count as down-/up-regulated,
#'   the rest as not significant ("ns"). Default \code{-1.96}/\code{1.96}
#'   (two-tailed 95% Z threshold).
#' @param taxrank \code{NULL} (default, no agglomeration) or a taxonomic rank
#'   to agglomerate to via [tax_glom()] before computing diversity.
#' @param fraction_id_name,fraction_ids Restrict to a subset of samples
#'   before computing diversity: \code{fraction_id_name} names a
#'   \code{sample_data(physeq)} column, \code{fraction_ids} the values of it
#'   to keep. Both default \code{NULL} (no restriction).
#' @param measure Character vector of diversity measure(s) passed to
#'   [phyloseq::estimate_richness()]'s \code{measures} (e.g.
#'   \code{"Shannon"}, \code{c("Shannon", "Simpson")}). \code{"Shannon"}
#'   alone is computed via the faster [sparse_shannon()] when \code{physeq}'s
#'   OTU table is a \code{\link{sparse_otu_table-class}}. Default
#'   \code{"Shannon"}.
#'
#' @return A list with:
#'   \describe{
#'     \item{diversity}{A data frame with a \code{sample_id} column and one
#'       column per requested \code{measure}, one row per sample.}
#'     \item{measure}{The \code{measure} argument, unchanged: names which
#'       \code{diversity} column(s) hold plottable diversity values.}
#'     \item{sample_data}{\code{sample_data(physeq)} as a data frame with an
#'       added \code{sample_id} column, one row per sample.}
#'     \item{depth}{A data frame with \code{sample_id} and \code{depth}
#'       (\code{\link[phyloseq]{sample_sums}}) columns, one row per sample.}
#'   }
#'   When \code{from_igseq = TRUE}, this is instead exactly
#'   [get_igseq_richness()]'s return value.
#'
#' @seealso [plot_alpha_diversity()], [get_igseq_richness()]
#'
#' @examples
#' data(ps_16s_refinement)
#' alpha_div <- get_alpha_diversity(ps_16s_refinement, measure = "Shannon")
#' head(alpha_div$diversity)
#' alpha_div$measure
#'
#' @export
get_alpha_diversity <- function(
  physeq,
  from_igseq = FALSE,
  transform_abundances = "identity",
  proportions = FALSE,
  low_lim = -1.96,
  high_lim = 1.96,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  measure = "Shannon"
) {
  .check_phyloseq(physeq)
  physeq <- reverseASV(physeq)

  if (from_igseq) {
    return(get_igseq_richness(
      ps_ig_score = physeq,
      proportions = proportions,
      low_lim = low_lim,
      high_lim = high_lim
    ))
  }

  physeq <- .filter_glom_transform_physeq(
    physeq,
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    taxrank = taxrank,
    transform_abundances = transform_abundances
  )

  if (
    is(otu_table(physeq), "sparse_otu_table") && identical(measure, "Shannon")
  ) {
    shannon_vals <- sparse_shannon(physeq)
    diversity <- data.frame(
      Shannon = shannon_vals,
      row.names = names(shannon_vals)
    )
  } else {
    diversity <- estimate_richness(physeq, measures = measure)
    rownames(diversity) <- sample_names(physeq)
  }
  diversity$sample_id <- rownames(diversity)

  sample_data_df <- as(sample_data(physeq), "data.frame")
  sample_data_df$sample_id <- rownames(sample_data_df)

  list(
    diversity = diversity,
    measure = measure,
    sample_data = sample_data_df,
    depth = data.frame(
      sample_id = sample_names(physeq),
      depth = sample_sums(physeq)
    )
  )
}

#' Plot Alpha Diversity
#'
#' Renders a [get_alpha_diversity()] (or [get_igseq_richness()]) result as a
#' boxplot/violin/scatter plot, with optional grouping, faceting, and hover
#' text (for interactive use via \code{plotly::ggplotly(tooltip = "text")}).
#'
#' @param alpha_div A list as returned by [get_alpha_diversity()] or
#'   [get_igseq_richness()].
#' @param measure Which \code{alpha_div$diversity} column to plot on the
#'   y-axis. Defaults to \code{alpha_div$measure[1]}.
#' @param hover_variables Character vector of \code{alpha_div$sample_data}
#'   column names to include in the hover text (\code{"text"} aesthetic).
#'   Default \code{NULL} (no hover variables besides \code{depth}, if
#'   \code{check_depth = TRUE}).
#' @param x \code{NULL} (default; single unlabeled x position) or the name of
#'   a merged sample-data column to use as the x-axis.
#' @param x_levels \code{NULL} (default, all levels in sorted order) or a
#'   character vector giving the levels of \code{x} to keep, in plotting
#'   order.
#' @param group \code{NULL} (default) or a merged sample-data column name to
#'   color points/boxes/violins by.
#' @param group_levels Levels of \code{group} to keep, in plotting order.
#'   Default \code{NULL} (all levels, sorted).
#' @param facet_mode \code{"wrap"} (default, uses \code{facet}/\code{facet_levels})
#'   or \code{"grid"} (uses \code{facet_row}/\code{facet_col}).
#' @param facet,facet_levels Column to facet by (\code{facet_mode = "wrap"})
#'   and the levels of it to keep, in order. Both default \code{NULL}.
#' @param facet_row,facet_row_levels,facet_col,facet_col_levels Row/column
#'   facet variables and their kept levels (\code{facet_mode = "grid"}). All
#'   default \code{NULL}.
#' @param facet_labeller Passed to [ggplot2::facet_wrap()]/
#'   [ggplot2::facet_grid()]'s \code{labeller}: \code{"label_value"} (default)
#'   or \code{"label_both"}.
#' @param shape \code{NULL} (default) or a merged sample-data column name to
#'   map to point shape.
#' @param shape_levels Levels of \code{shape} to keep, in plotting order.
#'   Default \code{NULL}.
#' @param size \code{NULL} (default) or a merged sample-data column name to
#'   map to point size. Ignored if \code{check_depth = TRUE} (\code{depth} is
#'   used instead).
#' @param point_size Base point radius. Default \code{1.5}.
#' @param remove_na_from_plot Logical. Drop samples with an \code{NA} in any
#'   of \code{x}/\code{group}/\code{facet}/\code{facet_row}/\code{facet_col}/
#'   \code{shape}/\code{size} before plotting. Default \code{FALSE}.
#' @param plot_type \code{NULL} (default: \code{"scatter"} if \code{x} is
#'   numeric, otherwise \code{"boxplot"}), or one of \code{"boxplot"},
#'   \code{"violin"}, \code{"scatter"}.
#' @param stat Logical. Add a group-comparison p-value layer
#'   ([ggpubr::stat_compare_means()] for \code{"boxplot"}/\code{"violin"}, or
#'   [ggpubr::stat_cor()] for \code{"scatter"}). Default \code{FALSE}.
#' @param check_depth Logical. Map point size to sequencing depth
#'   (\code{alpha_div$depth}) and add it to the hover text. Not compatible
#'   with \code{plotly}. Default \code{FALSE}.
#' @param alpha Point/jitter transparency. Default \code{1}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso [get_alpha_diversity()]
#'
#' @examples
#' data(ps_16s_refinement)
#' alpha_div <- get_alpha_diversity(ps_16s_refinement, measure = "Shannon")
#' plot_alpha_diversity(alpha_div, x = "Protocol", group = "Protocol")
#'
#' @export
plot_alpha_diversity <- function(
  alpha_div,
  measure = NULL,
  hover_variables = NULL,
  x = NULL,
  x_levels = NULL,
  group = NULL,
  group_levels = NULL,
  facet_mode = "wrap", # "grid", "wrap"
  facet = NULL,
  facet_levels = NULL,
  facet_row = NULL,
  facet_row_levels = NULL,
  facet_col = NULL,
  facet_col_levels = NULL,
  facet_labeller = "label_value", #"label_both" or "label_value"
  shape = NULL,
  shape_levels = NULL,
  size = NULL,
  point_size = 1.5,
  remove_na_from_plot = FALSE,
  plot_type = NULL, # automatic if NULL
  stat = FALSE,
  check_depth = FALSE, # NOTE: not compatible with plotly
  alpha = 1
) {
  measure <- measure %||% alpha_div$measure[1]

  full_sample_data <- merge(
    alpha_div$diversity,
    alpha_div$sample_data,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )
  if (!is.null(alpha_div$depth)) {
    full_sample_data <- merge(
      full_sample_data,
      alpha_div$depth,
      by = "sample_id",
      all.x = TRUE,
      sort = FALSE
    )
  }

  if (is.null(facet_mode)) {
    facet_mode <- "wrap"
  }
  # Set proper data types

  full_sample_data <- apply_levels_if_valid(full_sample_data, group, group_levels)

  if (is_valid_factor(full_sample_data, shape)) {
    full_sample_data <- apply_levels(full_sample_data, shape, shape_levels)
  } else {
    shape <- NULL
  }

  full_sample_data <- apply_levels_if_valid(full_sample_data, x, x_levels)

  is_valid_facet <- is_valid_factor(full_sample_data, facet) &&
    facet_mode == "wrap"

  if (is_valid_facet) {
    full_sample_data <- apply_levels(full_sample_data, facet, facet_levels)
  }

  is_valid_facet_row <- is_valid_factor(full_sample_data, facet_row) &&
    facet_mode == "grid"

  if (is_valid_facet_row) {
    full_sample_data <- apply_levels(
      full_sample_data,
      facet_row,
      facet_row_levels
    )
  }

  is_valid_facet_col <- is_valid_factor(full_sample_data, facet_col) &&
    facet_mode == "grid"

  if (is_valid_facet_col) {
    full_sample_data <- apply_levels(
      full_sample_data,
      facet_col,
      facet_col_levels
    )
  }

  # Check whether depth is present
  check_depth <- check_depth && !is.null(full_sample_data$depth)

  # Remove all NA's from plot data (labels, facets or shape) by removing samples
  # having NA for at least one of graphical parameters (x-axis variable and group variable)
  # remove these samples from sample data
  if (remove_na_from_plot) {
    full_sample_data <- remove_nas(
      full_sample_data,
      c(x, group, facet, facet_row, facet_col, shape, size)
    )
  }

  # Prepare hover information about samples based on sample_data
  hover_variables <- c(
    hover_variables,
    if (check_depth) {
      "depth"
    }
  )
  full_sample_data$hover.text <-
    get_hover_text(full_sample_data, hover_variables)

  # Handle various plot types

  is_continuous_x <- !is.null(x) && is.numeric(full_sample_data[[x]])

  if (is.null(plot_type)) {
    if (is_continuous_x) {
      plot_type <- "scatter"
    } else {
      plot_type <- "boxplot"
    }
  }

  if (is.null(x)) {
    mapping <- aes(x = 0, y = .data[[measure]])
  } else {
    mapping <- aes(x = .data[[x]], y = .data[[measure]])
  }

  point_mapping <- aes()

  if (!is.null(group)) {
    point_mapping <- modifyList(
      point_mapping,
      aes(color = .data[[group]])
    )
  }

  if (length(hover_variables) > 0) {
    point_mapping <- modifyList(
      point_mapping,
      aes(text = .data[["hover.text"]])
    )
  }

  if (!is.null(shape)) {
    point_mapping <- modifyList(
      point_mapping,
      aes(shape = .data[[shape]])
    )
  }

  has_size_aes <- FALSE
  if (check_depth) {
    point_mapping <- modifyList(
      point_mapping,
      aes(size = .data[["depth"]])
    )
    has_size_aes <- TRUE
  } else if (!is.null(size)) {
    point_mapping <- modifyList(
      point_mapping,
      aes(size = .data[[size]])
    )
    has_size_aes <- TRUE
  }

  if (plot_type == "boxplot") {
    plot_layer_fn <- ggplot2::geom_boxplot
  }

  if (plot_type == "violin") {
    plot_layer_fn <- ggplot2::geom_violin
  }

  if (plot_type == "scatter") {
    plot_layer_fn <- if (has_size_aes) {
      function() {
        list(
          ggplot2::geom_point(point_mapping),
          ggplot2::scale_size(range = c(point_size * 0.5, point_size * 3))
        )
      }
    } else {
      function() {
        ggplot2::geom_point(point_mapping, size = point_size)
      }
    }
  }

  if (plot_type %in% c("boxplot", "violin")) {
    plt <- ggplot(full_sample_data, mapping) +
      plot_layer_fn(outlier.shape = NA)
  } else {
    plt <- ggplot(full_sample_data, mapping) +
      plot_layer_fn()
  }

  if (plot_type != "scatter") {
    if (has_size_aes) {
      plt <- plt +
        ggplot2::geom_jitter(point_mapping, alpha = alpha) +
        ggplot2::scale_size(range = c(point_size * 0.5, point_size * 3))
    } else {
      plt <- plt +
        ggplot2::geom_jitter(point_mapping, alpha = alpha, size = point_size)
    }
  } else {
    smooth_mapping <- aes()
    if (!is.null(group)) {
      smooth_mapping <- modifyList(smooth_mapping, aes(color = .data[[group]]))
    }
    plt <- plt +
      stat_smooth(smooth_mapping, method = "lm", alpha = 0.1)
  }

  if (stat) {
    if (plot_type != "scatter") {
      plt <- plt +
        stat_compare_means()
    } else {
      plt <- plt +
        ggpubr::stat_cor()
    }
  }

  plt <- plt + theme_minimal()

  if (is_valid_facet_row && is_valid_facet_col) {
    plt <- plt +
      facet_grid(
        rows = vars(!!sym(facet_row)),
        cols = vars(!!sym(facet_col)),
        labeller = facet_labeller
      )
  } else if (is_valid_facet_row) {
    plt <- plt +
      facet_grid(
        rows = vars(!!sym(facet_row)),
        labeller = facet_labeller
      )
  } else if (is_valid_facet_col) {
    plt <- plt +
      facet_grid(
        cols = vars(!!sym(facet_col)),
        labeller = facet_labeller
      )
  } else if (is_valid_facet) {
    plt <- plt +
      facet_wrap(
        ~ .data[[facet]],
        ncol = smart_facet_ncol(
          nlevels(factor(full_sample_data[[facet]]))
        ),
        labeller = facet_labeller
      )
  }

  plt <- plt +
    labs(title = "Alpha Diversity") +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      legend.title = element_text(face = "bold", hjust = 0.5)
    )
  if (!is.null(group) && !is.numeric(full_sample_data[[group]])) {
    plt <- plt +
      ggsci::scale_fill_npg()
  }

  if (is.null(x)) {
    plt <- plt +
      theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  return(plt)
}

#' Compute and Plot Alpha Diversity
#'
#' Thin convenience wrapper chaining [get_alpha_diversity()] and
#' [plot_alpha_diversity()] in one call.
#'
#' @inheritParams get_alpha_diversity
#' @inheritParams plot_alpha_diversity
#' @param physeq A \code{phyloseq} object.
#' @param x \code{NULL} (default; single unlabeled x position) or the name of
#'   a sample-data column to use as the x-axis.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object.
#'
#' @seealso [get_alpha_diversity()], [plot_alpha_diversity()]
#'
#' @examples
#' data(ps_16s_refinement)
#' alpha_diversity(ps_16s_refinement, x = "Protocol", group = "Protocol")
#'
#' @export
alpha_diversity <- function(
  physeq,
  taxrank = NULL,
  fraction_id_name = NULL,
  fraction_ids = NULL,
  measure = "Shannon",
  x = NULL,
  group = NULL,
  plot_type = "boxplot",
  hover_variables = NULL,
  stat = FALSE,
  check_depth = FALSE,
  point_size = 1.5,
  alpha = 1
) {
  alpha_div <- get_alpha_diversity(
    physeq = physeq,
    taxrank = taxrank,
    fraction_id_name = fraction_id_name,
    fraction_ids = fraction_ids,
    measure = measure
  )
  plot_alpha_diversity(
    alpha_div,
    measure = measure,
    hover_variables = hover_variables,
    x = x,
    group = group,
    plot_type = plot_type,
    stat = stat,
    check_depth = check_depth,
    point_size = point_size,
    alpha = alpha
  )
}

#' Get Richness Based on Ig-Score Significance
#'
#' For each sample, buckets taxa by their Ig score into "down" (below
#' \code{low_lim}), "up" (above \code{high_lim}), or "ns" (not significant,
#' in between), and counts (or reports the proportion of) taxa in each
#' bucket.
#'
#' @param ps_ig_score A \code{phyloseq} object whose OTU table holds
#'   per-taxon Ig scores (e.g. slide-Z scores). Re-oriented to taxa-as-columns
#'   internally if needed.
#' @param proportions Logical. Report bucket proportions (summing to 1 per
#'   sample) instead of raw counts. Default \code{FALSE}.
#' @param low_lim,high_lim Numeric Ig-score cutoffs. Default
#'   \code{-1.96}/\code{1.96} (two-tailed 95% Z threshold).
#'
#' @return A list with:
#'   \describe{
#'     \item{diversity}{A data frame, one row per sample per significance
#'       bucket, with \code{sample_id}, \code{significance}
#'       (\code{"down"}/\code{"ns"}/\code{"up"}), and \code{richness}
#'       columns.}
#'     \item{measure}{\code{"richness"}.}
#'     \item{sample_data}{\code{sample_data(ps_ig_score)} as a data frame
#'       with an added \code{sample_id} column, one row per sample.}
#'     \item{depth}{\code{NULL} (sequencing depth is not meaningful for
#'       Ig-score data).}
#'   }
#'
#' @seealso [get_alpha_diversity()], [plot_igseq_richness()]
#'
#' @examples
#' ig_scores <- matrix(
#'   rnorm(40, sd = 1.5),
#'   nrow = 4,
#'   dimnames = list(paste0("ASV", 1:4), paste0("S", 1:10))
#' )
#' ps_ig_score <- phyloseq::phyloseq(
#'   phyloseq::otu_table(ig_scores, taxa_are_rows = TRUE),
#'   phyloseq::sample_data(data.frame(
#'     Group = rep(c("A", "B"), 5),
#'     row.names = colnames(ig_scores)
#'   ))
#' )
#' igseq_richness <- get_igseq_richness(ps_ig_score)
#' head(igseq_richness$diversity)
#'
#' @export
get_igseq_richness <- function(
  ps_ig_score,
  proportions = FALSE,
  low_lim = -1.96,
  high_lim = 1.96
) {
  if (taxa_are_rows(ps_ig_score)) {
    ps_ig_score <- t(ps_ig_score)
  }

  counts <- matrix(
    NA_real_,
    nrow = nsamples(ps_ig_score),
    ncol = 3,
    dimnames = list(sample_names(ps_ig_score), c("down", "ns", "up"))
  )
  for (sample_id in sample_names(ps_ig_score)) {
    scores <- otu_table(ps_ig_score)[sample_id, ]
    counts[sample_id, "down"] <- sum(scores < low_lim, na.rm = TRUE)
    counts[sample_id, "up"] <- sum(scores > high_lim, na.rm = TRUE)
    counts[sample_id, "ns"] <- sum(
      scores >= low_lim & scores <= high_lim,
      na.rm = TRUE
    )
  }

  if (proportions) {
    counts <- counts / rowSums(counts, na.rm = TRUE)
  }

  diversity <- as.data.frame(counts)
  diversity$sample_id <- rownames(diversity)
  diversity <- as.data.frame(pivot_longer(
    diversity,
    cols = c("down", "ns", "up"),
    names_to = "significance",
    values_to = "richness"
  ))

  sample_data_df <- as(sample_data(ps_ig_score), "data.frame")
  sample_data_df$sample_id <- rownames(sample_data_df)

  list(
    diversity = diversity,
    measure = "richness",
    sample_data = sample_data_df,
    depth = NULL
  )
}

#' Plot Richness Based on Ig-Score Significance
#'
#' Jitter+violin plot of [get_igseq_richness()]'s bucketed richness, grouped
#' by a sample-data variable and faceted by significance bucket.
#'
#' @param igseq_richness A list as returned by [get_igseq_richness()].
#' @param group Sample-data column name to use as the x-axis grouping
#'   variable.
#' @param color Sample-data column name to color points/violins by.
#' @param exclude_ns Logical. Drop the "ns" (not significant) bucket before
#'   plotting. Default \code{FALSE}.
#'
#' @return A \code{\link[ggplot2]{ggplot}} object, faceted by significance
#'   bucket.
#'
#' @seealso [get_igseq_richness()]
#'
#' @examples
#' ig_scores <- matrix(
#'   rnorm(40, sd = 1.5),
#'   nrow = 4,
#'   dimnames = list(paste0("ASV", 1:4), paste0("S", 1:10))
#' )
#' ps_ig_score <- phyloseq::phyloseq(
#'   phyloseq::otu_table(ig_scores, taxa_are_rows = TRUE),
#'   phyloseq::sample_data(data.frame(
#'     Group = rep(c("A", "B"), 5),
#'     row.names = colnames(ig_scores)
#'   ))
#' )
#' igseq_richness <- get_igseq_richness(ps_ig_score)
#' plot_igseq_richness(igseq_richness, group = "Group", color = "Group")
#'
#' @export
plot_igseq_richness <- function(igseq_richness, group, color, exclude_ns = FALSE) {
  full_data <- merge(
    igseq_richness$diversity,
    igseq_richness$sample_data,
    by = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  if (exclude_ns) {
    full_data <- full_data[full_data$significance != "ns", ]
  }

  full_data[[group]] <- factor(
    full_data[[group]],
    levels = gtools::mixedsort(unique(full_data[[group]]))
  )

  ggplot(
    full_data,
    aes(
      x = .data[[group]],
      y = .data[[igseq_richness$measure]],
      color = .data[[color]]
    )
  ) +
    geom_jitter() +
    geom_violin(alpha = 0.1) +
    facet_grid(~significance) +
    theme_minimal()
}
