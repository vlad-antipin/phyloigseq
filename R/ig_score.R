#' Compute Ig-Coating Scores From Sort-Fraction Abundances
#'
#' Computes a per-taxon immunoglobulin(Ig)-coating score from Ig+/Ig- (and optionally
#' pre-sort) fraction abundances, using one of several closed-form indices. Unlike
#' [get_slide_z()]'s sliding Z-score, these are classical, non-windowed scores computed
#' directly from the pos/neg(/pre) abundance vectors. See `IG_SCORES` for the set of
#' method names used elsewhere in the package's pipeline.
#'
#' `"purity_corrected_prob_index"`/`"purity_corrected_prob_ratio"` are experimental:
#' their derivation (see the inline comments in the source) is not yet fully verified and
#' they are not part of `IG_SCORES`/the main [getPhyloIgSeq()] pipeline. They remain
#' available here for direct use while under development.
#'
#' @param method One of `"palm"`, `"kau"`, `"prob_index"`, `"prob_ratio"`,
#'   `"purity_corrected_prob_index"`, `"purity_corrected_prob_ratio"`. See Details.
#' @param pos Numeric vector of Ig+ fraction counts or relative abundances, one value per
#'   taxon.
#' @param neg Numeric vector of Ig- fraction counts or relative abundances, same length
#'   as `pos`. Required by `"palm"`, `"kau"`, `"prob_ratio"`, and the `purity_corrected_*`
#'   methods.
#' @param pre Numeric vector of pre-sort counts or relative abundances, same length as
#'   `pos`. Required by `"prob_index"` and the `purity_corrected_*` methods.
#' @param ig_freq Single probability in `[0, 1]`, the total frequency of real Ig+
#'   bacteria (P(Ig+)). Required by `"prob_index"`/`"prob_ratio"`.
#' @param pos_purity Single probability in `[0, 1]`, P(real Ig+ | Ig+ fraction). Required
#'   by the `purity_corrected_*` methods.
#' @param neg_impurity Single probability in `[0, 1]`, P(real Ig+ | Ig- fraction).
#'   Required by the `purity_corrected_*` methods.
#' @param pos_fraction Single probability in `[0, 1]`, P(Ig+ fraction). Required by the
#'   `purity_corrected_*` methods.
#' @param neg_fraction Single probability in `[0, 1]`, P(Ig- fraction). Required by the
#'   `purity_corrected_*` methods.
#'
#' @return A numeric vector of the same length as `pos`, one score per taxon. `NaN`/
#'   infinite values (e.g. from a zero-abundance taxon in a denominator) are converted to
#'   `NA`.
#'
#' @examples
#' pos <- c(50, 30, 20, 5)
#' neg <- c(5, 10, 40, 45)
#' pre <- c(20, 20, 20, 40)
#'
#' compute_ig_score(method = "palm", pos = pos, neg = neg)
#' compute_ig_score(method = "kau", pos = pos, neg = neg)
#' compute_ig_score(method = "prob_index", pos = pos, pre = pre, ig_freq = 0.3)
#' compute_ig_score(method = "prob_ratio", pos = pos, neg = neg, ig_freq = 0.3)
#'
#' @export
compute_ig_score <- function(
  method = c(
    "palm",
    "kau",
    "prob_index",
    "prob_ratio",
    "purity_corrected_prob_index",
    "purity_corrected_prob_ratio"
  ),
  # vectors of the same length
  # corresponding to abundance of each taxon in:
  pos, # Ig+ fraction counts or rel. abundance
  neg = NULL, # Ig- fraction -//-
  pre = NULL, # before sorting -//-
  # Probabilities in [0,1]:
  ig_freq = NULL, # total frequency of Ig+ bacteria
  pos_purity = NULL, # P(real Ig+ | Ig+ fraction)
  neg_impurity = NULL, # P(real Ig+ | Ig- fraction)
  pos_fraction = NULL, # P(Ig+ fraction)
  neg_fraction = NULL # P(Ig- fraction)
) {
  method <- match.arg(method)

  # Transform counts to relative abundances
  if (is.numeric(pos)) {
    # P(taxon | Ig+ fraction)
    pos_abund <- pos / sum(pos)
  } else {
    pos_abund <- NA
  }

  if (is.numeric(neg)) {
    # P(taxon | Ig- fraction)
    neg_abund <- neg / sum(neg)
  } else {
    neg_abund <- NA
  }
  if (is.numeric(pre)) {
    # P(taxon) - presorting
    pre_abund <- pre / sum(pre)
  } else {
    pre_abund <- NA
  }

  if (!is.numeric(ig_freq)) {
    ig_freq <- NA
  }

  # TODO: Sequences tend to be sequenced at equal molar ratios, so even without rarefaction the counts are not comparable !
  # so user should provide P(Ig+ fraction) themselves based on phenotyping!
  # if( sum(pos) == sum(neg) ){
  #   warning("It looks like fraction counts have been rarefied to the same abundance, it would not allow to correctly eastimate P(Ig+ fraction)\n")
  # }
  # if(is.null(pos_fraction)){pos_fraction = pos / (pos + neg)} # P(Ig+ fraction) # or pos / pre?
  # if(is.null(neg_fraction)){neg_fraction = neg / (pos + neg)} # P(Ig- fraction)
  if (is.null(pos_fraction)) {
    pos_fraction <- NA
  }
  if (is.null(neg_fraction)) {
    neg_fraction <- NA
  }

  if (is.null(pos_purity)) {
    pos_purity <- NA
  }
  if (is.null(neg_impurity)) {
    neg_impurity <- NA
  }

  score <- switch(
    method,
    palm = pos_abund / neg_abund,
    kau = {
      # minus to negate the fact that log10(pos_abund*neg_abund) is negative
      -log2(pos_abund / neg_abund) / log10(pos_abund * neg_abund)
    },
    prob_index = pos_abund * ig_freq / pre_abund, # P(Ig+ | taxon) = P(taxon | Ig+) * P(Ig+) / P(taxon)
    prob_ratio = log2(pos_abund * ig_freq / (neg_abund * (1 - ig_freq))),
    purity_corrected_prob_index = {
      # TODO: verify the proof!!!
      # Here, we discriminate `real Ig+` and `Ig+ fraction` events, assuming that
      # the purity of a fraction is independent of a taxon abundance in this fraction
      # i.e. independence of `real Ig+` and `taxon` events conditional on knowing
      # the fraction (Ig+/-):
      #
      # P(taxon & real Ig+ | fraction) = P(taxon | fraction) * P(real Ig+ | fraction)
      # flu -> temp -> test => P(test & flu | temp) = P(flu | temp) + P(test | temp)
      # which corresponds to this Bayesian network:
      # taxon -> fraction -> real Ig+ (knowing fraction breaks dependence between `taxon` and `real Ig+`)
      # BUT: theoretically, the correct network should be:
      # taxon -> real Ig+ -> fraction
      # But this will require conditioning on `real Ig+` and knowing P(taxon | real Ig+)
      # which is technically unavailable :(
      #
      # So this assumption is not totally justified theoretically, but is probably
      # still better in case of low fraction purity (?) - to verify!
      #
      # P(real Ig+ | taxon) =
      #   (   P(taxon | Ig+ fraction) * P(real Ig+ | Ig+ fraction) * P(Ig+ fraction)
      #     + P(taxon | Ig- fraction) * P(real Ig+ | Ig- fraction) * P(Ig- fraction)
      #    ) / P(taxon)
      (pos_abund *
        pos_purity *
        pos_fraction +
        neg_abund * neg_impurity * neg_fraction) /
        pre_abund
    },
    purity_corrected_prob_ratio = {
      # TODO: verify!
      prob <- pos_abund *
        pos_purity *
        pos_fraction +
        neg_abund * neg_impurity * neg_fraction
      log2(prob / (1 - prob))
    }
  )

  score[is.nan(score) | is.infinite(score)] <- NA

  return(score)
}


#' Compute a Jitter Offset for Off-Axis Points
#'
#' Shared by [plot_slide_z()] and [plot_ma()]: points that can't be placed at
#' their true x-position (e.g. an imputed/undefined log-abundance) are instead
#' jittered in a band just past the low end of the plotted x-range.
#'
#' @param values Numeric vector the x-range is computed from (the abundance
#'   values of the points that *do* get plotted at their true position).
#'
#' @return A list with `width` (the jitter band's width, `range(values) / 6`)
#'   and `x` (the band's center, three widths below `min(values)`).
#' @noRd
.jitter_offset <- function(values) {
  jitter_width <- diff(range(values)) / 6
  list(width = jitter_width, x = min(values) - jitter_width * 3)
}

#' Truncate Long Values for a ggplot Hover Tooltip
#'
#' @param x Vector coerced to character; values longer than `max_chars` are
#'   cut short with a trailing `"..."` (taxon names/IDs can be full ASV/OTU
#'   sequences, which would otherwise bloat the plot object and the tooltip).
#' @param max_chars Single integer, the maximum number of characters to keep.
#'
#' @return A character vector the same length as `x`.
#' @noRd
.truncate_for_tooltip <- function(x, max_chars) {
  x <- as.character(x)
  ifelse(
    is.na(x) | nchar(x) <= max_chars,
    x,
    paste0(substr(x, 1, max_chars), "...")
  )
}

#' Build Per-Taxon Hover Tooltips for `plot_slide_z()`
#'
#' @param ig_df `phyloigseq_obj@ig_coating` (or a subset), one row per
#'   taxon/sample, with a `slide_z` column.
#' @param tax_table `phyloigseq_obj@tax_table`, or `NULL`/empty if not
#'   available.
#' @param max_chars Passed to `.truncate_for_tooltip()`.
#'
#' @return A character vector the same length as `nrow(ig_df)`.
#' @noRd
.slide_z_tooltip <- function(ig_df, tax_table, max_chars) {
  if (length(tax_table) == 0) {
    return(paste0("slide_z: ", round(ig_df$slide_z, digits = 3)))
  }

  tax_cols <- colnames(tax_table)
  matched_tax <- tax_table[
    match(ig_df$taxon_id, tax_table$taxon_id),
    tax_cols,
    drop = FALSE
  ]

  hover_lines <- lapply(tax_cols, function(col) {
    paste0(
      col,
      ": ",
      .truncate_for_tooltip(matched_tax[[col]], max_chars),
      "<br>"
    )
  })

  paste0(
    do.call(paste0, hover_lines),
    "slide_z: ",
    round(ig_df$slide_z, digits = 3)
  )
}

#' Flag Imputed (`sample_id`, `taxon_id`) Pairs
#'
#' @param imputed_taxa `phyloigseq_obj@imputed_taxa`: a list of `taxon_id`
#'   vectors, named by `sample_id`.
#'
#' @return A character vector of `"sample_id taxon_id"` keys, one per
#'   imputed pair, suitable for `%in%` against `paste(sample_id, taxon_id)`.
#' @noRd
.imputed_taxa_lookup <- function(imputed_taxa) {
  if (length(imputed_taxa) == 0) {
    return(character(0))
  }
  unlist(lapply(names(imputed_taxa), function(sample_id) {
    taxa <- imputed_taxa[[sample_id]]
    if (length(taxa) == 0) {
      return(character(0))
    }
    paste(sample_id, taxa)
  }))
}

#' Plot Sliding Z-Scores
#'
#' Draws the sliding-Z MA-plot for the sample(s) in a `PhyloIgSeq` object built with `"slide_z"`
#' among its `scores`: log-abundance on the x-axis, log-ratio on the y-axis, colored/sized by
#' whether `|slide_z|` exceeds `z_alpha2`. Optionally overlays the null distribution (empirical or
#' theoretical, see [get_slide_z()]) and the confidence ellipses from
#' `phyloigseq_obj@ellipse_coords`. Taxa with an imputed zero (`phyloigseq_obj@imputed_taxa`) are
#' drawn jittered just past the plot's x-range instead of at their (undefined) true abundance, as
#' in [plot_ma()].
#'
#' @param phyloigseq_obj A [PhyloIgSeq-class] object with `"slide_z"` in `score_names` (i.e. built
#'   by [getPhyloIgSeq()] with `"slide_z"` among `scores`).
#' @param sample_ids Optional character vector of `sample_id`s to restrict/order the plot to;
#'   `NULL` (default) plots every sample in `ig_coating`, faceted by `sample_id` if there is more
#'   than one.
#' @param empirical_null_distribution Logical. If `TRUE` (default), overlay the empirical null
#'   distribution (`null_abundance`/`null_change`, from a second negative fraction). Falls back to
#'   `FALSE` with a `warning()` if `ig_coating` doesn't carry those columns.
#' @param z_alpha2 Single number, the `|slide_z|` significance threshold used to color/size points
#'   `"signif"` vs `"ns"` and to label the legend/axis. Default `1.96` (two-sided 95% threshold).
#' @param signif_colors Character vector of length 2, passed to `scale_color_manual()`. Matched to
#'   the alphabetically-sorted significance categories (`"ns"` first, then `"signif"`). Defaults to
#'   two colors from `ggsci::pal_npg()`.
#' @param ellipses Logical. If `TRUE` (default), overlay the confidence ellipse boundaries from
#'   `phyloigseq_obj@ellipse_coords`. Falls back to `FALSE` with a `warning()` if none are present.
#' @param tooltip_max_chars Single integer, the maximum number of characters of each `tax_table`
#'   value (e.g. `taxon_name`, which may be a full ASV/OTU sequence) shown in a point's hover
#'   tooltip before truncating with `"..."`. Default `40`.
#'
#' @return A `ggplot` object (one panel per `sample_id` if more than one is plotted).
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
#'   scores = c("slide_z", "palm", "kau"),
#'   confidence_levels = c(0.95, 0.99)
#' )
#' plot_slide_z(pis)
#' plot_slide_z(pis, sample_ids = "sample_1", ellipses = FALSE)
#'
#' @export
plot_slide_z <- function(
  phyloigseq_obj,
  sample_ids = NULL, # if NULL, all samples are plotted
  empirical_null_distribution = TRUE,
  z_alpha2 = 1.96,
  signif_colors = c(ggsci::pal_npg()(2)[2], ggsci::pal_npg()(2)[1]),
  ellipses = TRUE,
  tooltip_max_chars = 40
) {
  ig_df <- phyloigseq_obj@ig_coating
  ellipse_df <- phyloigseq_obj@ellipse_coords
  if (
    empirical_null_distribution &
      !all(c("null_abundance", "null_change") %in% colnames(ig_df))
  ) {
    warning(
      "No MA coordinates to plot empirical null distribution are furnished, using observed pos-neg distribution as null distribution...\n"
    )
    empirical_null_distribution <- FALSE
  }

  if (!empirical_null_distribution) {
    ig_df$null_abundance <- ig_df$obs_abundance
    ig_df$null_change <- ig_df$obs_change
  }

  if (is.null(ellipse_df) || prod(dim(ellipse_df)) == 0) {
    warning("No ellipse coordinates furnished\n")
    ellipses <- FALSE
  }

  if (!is.null(sample_ids)) {
    ig_df <- ig_df[ig_df$sample_id %in% sample_ids, ]
    ig_df$sample_id <- factor(ig_df$sample_id, levels = sample_ids)
    ellipse_df <- ellipse_df[ellipse_df$sample_id %in% sample_ids, ]
    ellipse_df$sample_id <- factor(ellipse_df$sample_id, levels = sample_ids)
  }

  ig_df$tooltip <- .slide_z_tooltip(
    ig_df,
    phyloigseq_obj@tax_table,
    tooltip_max_chars
  )

  is_imputed <- paste(ig_df$sample_id, ig_df$taxon_id) %in%
    .imputed_taxa_lookup(phyloigseq_obj@imputed_taxa)

  ig_df_imputed <- ig_df[is_imputed, , drop = FALSE]
  ig_df <- ig_df[!is_imputed, , drop = FALSE]

  stat_imputed <- ifelse(
    ig_df_imputed$slide_z >= z_alpha2 | ig_df_imputed$slide_z <= -z_alpha2,
    "signif",
    "ns"
  )
  stat <- ifelse(
    ig_df$slide_z >= z_alpha2 | ig_df$slide_z <= -z_alpha2,
    "signif",
    "ns"
  )

  plt <- ggplot(ig_df)

  jitter <- .jitter_offset(c(ig_df$null_abundance, ig_df$obs_abundance))
  jitter_width <- jitter$width
  jitter_x <- jitter$x
  plt <- plt +
    geom_jitter(
      data = ig_df_imputed,
      aes(x = jitter_x, y = null_change, text = tooltip),
      color = "darkgray",
      alpha = 0.6,
      position = position_jitter(width = jitter_width, height = 0)
    ) +
    geom_jitter(
      data = ig_df_imputed,
      aes(
        x = jitter_x,
        y = obs_change,
        size = stat_imputed,
        color = stat_imputed,
        # shape = stat_imputed,
        text = tooltip
      ),
      alpha = 0.8,
      position = position_jitter(width = jitter_width, height = 0)
    )

  if (empirical_null_distribution) {
    plt <- plt +
      geom_point(
        aes(x = null_abundance, y = null_change, text = tooltip),
        color = "darkgrey",
        alpha = 0.6
      )
  }

  plt <- plt +
    geom_point(
      aes(obs_abundance, obs_change, size = stat, color = stat, text = tooltip),
      alpha = 0.8
    ) +
    # ggsci::scale_color_npg()+
    scale_color_manual(values = signif_colors) +
    scale_size_discrete(range = c(1.5, 3))

  if (nrow(ellipse_df) > 0 & ellipses) {
    plt <- plt +
      geom_path(
        data = ellipse_df,
        aes(x = x, y = y, group = ellipse_level),
        color = "darkgrey",
        linetype = 2
      )
  }

  sample_id_levels <- unique(ig_df$sample_id)
  if (length(sample_id_levels) > 1) {
    plt <- plt +
      facet_wrap(. ~ sample_id)
  }

  fraction1 <- if (!is.null(phyloigseq_obj@positive_fraction_name)) {
    phyloigseq_obj@positive_fraction_name
  } else {
    "fraction_{1}"
  }
  fraction2 <- if (!is.null(phyloigseq_obj@first_negative_fraction_name)) {
    phyloigseq_obj@first_negative_fraction_name
  } else {
    "fraction_{2}"
  }
  plt <- plt +
    labs(
      x = latex2exp::TeX(paste0(
        "Log-Abundance: $\\log_{10}\\left(\\",
        fraction1,
        "\\cdot\\",
        fraction2,
        "\\right)$"
      )),
      y = latex2exp::TeX(paste0(
        "Log-Ratio: $\\log_{2}\\left(\\frac{\\",
        fraction1,
        "}{\\",
        fraction2,
        "}\\right)$"
      )),
      color = paste0("|sliding Z| >", z_alpha2),
      size = paste0("|sliding Z| >", z_alpha2),
      title = paste0(
        "Sliding Z Score",
        if (length(sample_id_levels) == 1) {
          paste0(" of ", sample_id_levels)
        }
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.title = element_text(face = "bold", hjust = 0.5)
    )
  return(plt)
}

#' Score-Specific Plot Boundaries For `plot_ig_score()`
#'
#' Internal. Returns the significance thresholds (`left_lim`/`right_lim`), color-scale midpoint,
#' and legal value range (`left_boundary`/`right_boundary`) used by `plot_ig_score()` to color
#' points/boxes and draw reference lines, one set per `score_name`. Errors for any `score_name`
#' not in `IG_SCORES` (i.e. not one of `"slide_z"`/`"kau"`/`"prob_ratio"`/`"palm"`/`"prob_index"`)
#' rather than leaving these unbound.
#'
#' @noRd
.ig_score_boundary <- function(score_name, z_alpha2) {
  switch(
    score_name,
    slide_z = list(
      left_lim = -z_alpha2,
      right_lim = z_alpha2,
      midpoint = 0,
      left_boundary = -Inf,
      right_boundary = Inf
    ),
    kau = ,
    prob_ratio = list(
      left_lim = 0,
      right_lim = 0,
      midpoint = 0,
      left_boundary = -Inf,
      right_boundary = Inf
    ),
    palm = list(
      left_lim = 1,
      right_lim = 1,
      midpoint = 1,
      left_boundary = 0,
      right_boundary = Inf
    ),
    prob_index = list(
      left_lim = 0.5,
      right_lim = 0.5,
      midpoint = 0.5,
      left_boundary = 0,
      right_boundary = 1
    ),
    stop(
      "`plot_ig_score()` has no known plotting boundary for score_name = '",
      score_name,
      "'. Supported names: 'slide_z', 'kau', 'prob_ratio', 'palm', 'prob_index'."
    )
  )
}

#' Central-Tendency Dispatch Shared by `agglomPhyloIgSeq()` and `.ig_score_agglomerate()`
#'
#' Internal. Collapses a numeric vector `x` to a single value via `method`: `"mean"`,
#' `"median"`, or `"weight_by_abund"` (a `weights`-weighted mean; `weights` is required in
#' that case). `NA`s are always dropped from `x` (and `weights`, for `"weight_by_abund"`).
#'
#' @noRd
.central_tendency <- function(x, method, weights = NULL) {
  switch(
    method,
    mean = mean(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    weight_by_abund = {
      if (is.null(weights)) {
        stop("Need abundance fraction to compute weighted score")
      }
      weighted.mean(x, weights, na.rm = TRUE)
    },
    stop("wrong agglomeration method")
  )
}

#' Two-Level Central-Tendency Agglomeration For `plot_ig_score()`
#'
#' Internal. Collapses `plot_data` (one row per taxon/sample) to one row per
#' (`taxrank_score`, `group_score`) combination via `score_agglom_fn`, either simultaneously
#' (`first_score_agglom_for_each = "both"`) or in two stages: first within each sample (or each
#' taxon), then across the remaining dimension. `"sample"`/`"taxon"`/`"both"` only give different
#' results when `score_agglom_fn` is non-linear across the intermediate step (e.g. `"median"`).
#' Central tendency itself is computed via the shared `.central_tendency()` (`score_agglom_fn`
#' is always `"mean"`/`"median"` here, never `"weight_by_abund"`).
#'
#' @noRd
.ig_score_agglomerate <- function(
  plot_data,
  score_name,
  score_agglom_fn,
  taxrank_score,
  taxrank_facet,
  group_score,
  group_facet,
  first_score_agglom_for_each
) {
  keep_cols <- unique(c(
    taxrank_score,
    taxrank_facet,
    group_score,
    group_facet,
    "agglom_score"
  ))

  if (first_score_agglom_for_each == "both") {
    plot_data %>%
      # Compute central tendency (mean, median) in a sample group and taxrank
      # - for all values in intersections formed by these groupings
      group_by(.data[[group_score]], .data[[taxrank_score]]) %>%
      mutate(
        agglom_score = .central_tendency(.data[[score_name]], score_agglom_fn)
      ) %>%
      select(all_of(keep_cols)) %>%
      distinct() %>%
      ungroup()
  } else {
    first_id <- if (first_score_agglom_for_each == "sample") {
      "sample_id"
    } else {
      "taxon_id"
    }
    first_agglom <- if (first_score_agglom_for_each == "sample") {
      taxrank_score
    } else {
      group_score
    }

    plot_data %>%
      # First, compute central tendency separately for each individual sample (or taxon)
      # grouping by a taxrank (or sample group)
      group_by(.data[[first_id]], .data[[first_agglom]]) %>%
      mutate(
        agglom_score = .central_tendency(.data[[score_name]], score_agglom_fn)
      ) %>%
      # ungroup and get rid of duplications, otherwise the central tendency will be false!
      ungroup() %>%
      select(all_of(keep_cols)) %>%
      distinct() %>%
      # Then, compute central tendency for each sample group (or taxrank), based on central
      # tendencies for each sample per taxrank (or for each taxon per sample group)
      group_by(.data[[group_score]], .data[[taxrank_score]]) %>%
      mutate(
        agglom_score = .central_tendency(agglom_score, score_agglom_fn)
      ) %>%
      ungroup()
  }
}

#' Pairwise `taxrank_score` Comparisons With Enough Data For `plot_ig_score()`
#'
#' Internal. Builds the `comparisons` list passed to [ggpubr::stat_compare_means()]: every pair
#' of `taxrank_score` levels that have at least 2 (already-agglomerated) data points in *every*
#' facet panel (`taxrank_facet`/`group_facet`) they appear in, not just overall — a level with
#' plenty of data summed across facets but only 1 point in one specific panel is excluded, since
#' a comparison can't be drawn from a single point in that panel. Returns `NULL` when fewer than
#' 2 levels qualify (nothing to compare).
#'
#' @noRd
.ig_score_valid_comparisons <- function(
  plot_data,
  taxrank_score,
  taxrank_facet,
  group_facet
) {
  facet_vars <- c(taxrank_facet, group_facet)
  count_cols <- unique(c(taxrank_score, facet_vars))

  group_counts <- plot_data %>%
    count(across(all_of(count_cols)), name = "n")

  if (length(facet_vars) > 0) {
    # Worst case (smallest) count across only the facet panels where this taxrank_score level
    # actually appears - a panel it's absent from isn't a "too few points" problem.
    group_counts <- group_counts %>%
      group_by(across(all_of(taxrank_score))) %>%
      summarise(n = min(n), .groups = "drop")
  }

  valid_groups <- group_counts[[taxrank_score]][group_counts$n >= 2]
  if (length(valid_groups) > 1) {
    combn(valid_groups, 2, simplify = FALSE)
  } else {
    NULL
  }
}

#' Shared `facet_grid()` For `plot_ig_score()`
#'
#' Internal. `taxrank_facet` and `group_facet` swap between rows/columns when `transpose = TRUE`,
#' matching the non-faceted axis swap done by `coord_flip()` elsewhere in `plot_ig_score()`.
#'
#' @noRd
.ig_score_facet_grid <- function(taxrank_facet, group_facet, transpose, scales) {
  row_var <- if (transpose) group_facet else taxrank_facet
  col_var <- if (transpose) taxrank_facet else group_facet
  facet_grid(
    rows = if (!is.null(row_var)) vars(.data[[row_var]]),
    cols = if (!is.null(col_var)) vars(.data[[col_var]]),
    scales = scales,
    space = "free"
  )
}

#' Shared `theme()` Blocks For `plot_ig_score()`
#'
#' Internal. `.ig_score_base_theme()` is used by both the bubbleplot and boxplot/violin render
#' paths; `.ig_score_transpose_theme()` is layered on top of it when `transpose = TRUE` (after
#' `coord_flip()`), targeting the axis text position that applies post-flip.
#'
#' @noRd
.ig_score_base_theme <- function() {
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(face = "bold", hjust = 0.5),
    legend.direction = "horizontal",
    axis.text.y.left = element_text(angle = 0, hjust = 1),
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
}

#' @noRd
.ig_score_transpose_theme <- function() {
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    legend.title = element_text(face = "bold", hjust = 0.5),
    legend.direction = "horizontal",
    axis.text.x.bottom = element_text(angle = 45, hjust = 1),
    strip.text.y.right = element_text(angle = 0, hjust = 0),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
}

#' Bubbleplot Render Path For `plot_ig_score()`
#'
#' Internal. One tile per (`group_score`, `taxrank_score`) combination, sized by
#' `abs(agglom_score)` and colored on a diverging scale centered at `boundary$midpoint`.
#'
#' @noRd
.ig_score_bubbleplot <- function(
  plot_data,
  group_score,
  taxrank_score,
  taxrank_facet,
  group_facet,
  transpose,
  signif_colors,
  boundary,
  score_agglom_fn,
  score_name
) {
  ggplot(
    plot_data,
    aes(
      x = .data[[group_score]],
      y = .data[[taxrank_score]],
      size = abs(agglom_score),
      fill = agglom_score
    )
  ) +
    geom_point(pch = 21) +
    scale_fill_gradient2(
      high = signif_colors[1],
      low = signif_colors[2],
      midpoint = boundary$midpoint,
      limits = c(
        max(
          boundary$left_boundary,
          boundary$midpoint - max(abs(plot_data$agglom_score - boundary$midpoint))
        ),
        min(
          boundary$right_boundary,
          boundary$midpoint + max(abs(plot_data$agglom_score - boundary$midpoint))
        )
      ),
      guide = guide_colourbar(title.position = "top", title.hjust = 0.5)
    ) +
    guides(size = "none") +
    .ig_score_facet_grid(taxrank_facet, group_facet, transpose, scales = "free") +
    theme_minimal() +
    labs(x = NULL, y = NULL, fill = paste(score_agglom_fn, score_name)) +
    .ig_score_base_theme()
}

#' Boxplot/Violin Render Path For `plot_ig_score()`
#'
#' Internal. One box/violin per `taxrank_score` level, jittered points colored by whether
#' `agglom_score` crosses `boundary$left_lim`/`right_lim`, with optional pairwise significance
#' brackets (`add_stats`).
#'
#' @noRd
.ig_score_boxplot <- function(
  plot_data,
  plot_type,
  taxrank_score,
  taxrank_facet,
  group_facet,
  transpose,
  signif_colors,
  boundary,
  add_stats,
  valid_comparisons,
  score_agglom_fn,
  score_name
) {
  plot_data$point_color <- ifelse(
    plot_data$agglom_score > boundary$right_lim,
    "high",
    ifelse(plot_data$agglom_score < boundary$left_lim, "low", "ns")
  )

  plt <- ggplot(data = plot_data, aes(x = agglom_score, y = .data[[taxrank_score]]))
  plt <- plt +
    if (plot_type == "violin") geom_violin() else geom_boxplot(outliers = FALSE)

  if (add_stats) {
    plt <- plt +
      stat_compare_means(
        method = "wilcox.test",
        comparisons = valid_comparisons,
        label = "p.signif",
        hide.ns = TRUE,
        p.adjust.method = "BH"
      )
  }

  plt <- plt +
    geom_jitter(
      alpha = 0.5,
      aes(
        size = abs(agglom_score - mean(boundary$left_lim, boundary$right_lim)),
        color = point_color
      )
    ) +
    scale_color_manual(
      values = c("high" = signif_colors[1], "low" = signif_colors[2], "ns" = "darkgrey")
    ) +
    guides(size = "none", color = "none") +
    scale_size_continuous(range = c(1, 2))

  if (mean(boundary$left_lim, boundary$right_lim) == 0) {
    plt <- plt +
      scale_x_continuous(
        limits = c(
          -max(abs(plot_data$agglom_score)),
          max(abs(plot_data$agglom_score))
        )
      )
  }

  plt +
    geom_vline(xintercept = unique(c(boundary$left_lim, boundary$right_lim)), linetype = 2) +
    .ig_score_facet_grid(
      taxrank_facet,
      group_facet,
      transpose,
      scales = if (transpose) "free_x" else "free_y"
    ) +
    theme_minimal() +
    labs(x = paste(score_agglom_fn, score_name), y = NULL) +
    .ig_score_base_theme()
}

#' Plot an Ig Score by Taxon and Sample Group
#'
#' Renders an `ig_coating` score (`score_name`) agglomerated by taxon (`taxrank_score`) and
#' sample group (`group_score`), either as a bubble plot (one tile per taxon x group, sized/
#' colored by the agglomerated score) or as a boxplot/violin plot (one box/violin per taxon,
#' jittered points colored by whether they cross `score_name`'s significance boundary).
#'
#' @param phyloigseq_obj A `PhyloIgSeq` object, e.g. from [getPhyloIgSeq()].
#' @param plot_type One of `"boxplot"`, `"violin"`, `"bubbleplot"`.
#' @param score_name Name of an `ig_coating` column to plot (one of `IG_SCORES`, e.g.
#'   `"slide_z"`). See Details for the boundary/significance threshold used per score.
#' @param taxrank_score Column of `tax_table` (or `"taxon_id"`) to agglomerate the score by and
#'   show on the taxon axis.
#' @param taxrank_facet Optional column of `tax_table` to facet by, in addition to
#'   `taxrank_score`.
#' @param group_score Column of `sample_data` (or `"sample_id"`) to agglomerate the score by and
#'   show on the sample-group axis.
#' @param group_facet Optional column of `sample_data` to facet by, in addition to `group_score`.
#' @param score_agglom_fn One of `"mean"`, `"median"`: central-tendency function used to
#'   agglomerate `score_name` within a group.
#' @param first_score_agglom_for_each One of `"sample"`, `"taxon"`, `"both"`: whether to
#'   agglomerate first within each sample (`group_score` level) and then across taxa, first
#'   within each taxon (`taxrank_score` level) and then across samples, or simultaneously across
#'   both dimensions. Only matters when `score_agglom_fn = "median"` (non-linear); with `"mean"`
#'   all three give the same result.
#' @param z_alpha2 Two-tailed significance threshold; only used when `score_name == "slide_z"`.
#' @param exclude_na Logical. Drop rows with `NA` `score_name` values before plotting.
#' @param transpose Logical. Flip the taxon/group axes (`coord_flip()`).
#' @param signif_colors Length-2 color vector, `c(high, low)`, for points/boundary crossings.
#' @param add_stats Logical. Add pairwise Wilcoxon significance brackets (via
#'   [ggpubr::stat_compare_means()]) between `taxrank_score` levels with enough data in every
#'   facet panel they appear in. Only applies to `plot_type = "boxplot"`/`"violin"`; silently
#'   disabled (no brackets, no error) if fewer than 2 levels qualify.
#'
#' @details
#' `score_name`'s plotting boundary (significance thresholds, color-scale midpoint, and legal
#' value range) is looked up internally for the 5 names in `IG_SCORES` (`"slide_z"`, `"kau"`,
#' `"prob_ratio"`, `"palm"`, `"prob_index"`) and errors for any other `score_name` — this covers
#' every score currently produced by [getPhyloIgSeq()]/[compute_ig_score()], but a custom/future
#' score name needs its own boundary added internally (`.ig_score_boundary()`) before it can be
#' plotted here.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' phyloigseq_obj <- getPhyloIgSeq(
#'   ps_igseq,
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   scores = c("slide_z", "palm")
#' )
#' plot_ig_score(phyloigseq_obj, plot_type = "boxplot", score_name = "slide_z")
#' plot_ig_score(phyloigseq_obj, plot_type = "bubbleplot", score_name = "palm")
#'
#' @export
plot_ig_score <- function(
  phyloigseq_obj,
  plot_type = c("boxplot", "violin", "bubbleplot"),
  score_name = "slide_z",
  taxrank_score = "taxon_id", # taxrank level to agglomerate the Ig score
  taxrank_facet = NULL, # taxrank for faceting
  group_score = "sample_id", # sample group to agglomerate the Ig score
  group_facet = NULL, # sample group for facetting
  score_agglom_fn = c("mean", "median"),
  first_score_agglom_for_each = c("sample", "taxon", "both"),
  z_alpha2 = 1.96, # in case of z score
  exclude_na = TRUE,
  transpose = FALSE,
  signif_colors = ggsci::pal_npg()(2),
  add_stats = TRUE
) {
  plot_type <- match.arg(plot_type)
  score_agglom_fn <- match.arg(score_agglom_fn)
  first_score_agglom_for_each <- match.arg(first_score_agglom_for_each)

  if (score_name == "slide_z" && is.null(z_alpha2)) {
    z_alpha2 <- 1.96
  }

  # TODO: add a possibility for multiple faceting (e.g. with timepoints)

  # Agglomeration of Ig score should be either with 1. median 2. mean 3. weighted average by
  # abundance.
  # TODO: 3.

  # TODO: check the correctness of data agglomeration and averaging

  # Merge score, sample and taxonomic data together
  plot_data <- phyloigseq_obj@ig_coating[, c(
    "taxon_id",
    "sample_id",
    score_name
  )] %>%
    merge(
      phyloigseq_obj@sample_data[,
        unique(c("sample_id", group_score, group_facet)),
        drop = FALSE
      ],
      by = "sample_id"
    ) %>%
    merge(
      phyloigseq_obj@tax_table[,
        unique(c("taxon_id", taxrank_score, taxrank_facet)),
        drop = FALSE
      ],
      by = "taxon_id"
    )

  plot_data <- .ig_score_agglomerate(
    plot_data,
    score_name = score_name,
    score_agglom_fn = score_agglom_fn,
    taxrank_score = taxrank_score,
    taxrank_facet = taxrank_facet,
    group_score = group_score,
    group_facet = group_facet,
    first_score_agglom_for_each = first_score_agglom_for_each
  )

  if (exclude_na) {
    plot_data <- na.omit(plot_data)
  }

  valid_comparisons <- NULL
  if (add_stats) {
    valid_comparisons <- .ig_score_valid_comparisons(
      plot_data,
      taxrank_score = taxrank_score,
      taxrank_facet = taxrank_facet,
      group_facet = group_facet
    )
    add_stats <- !is.null(valid_comparisons)
  }

  boundary <- .ig_score_boundary(score_name, z_alpha2)

  plt <- if (plot_type == "bubbleplot") {
    .ig_score_bubbleplot(
      plot_data,
      group_score = group_score,
      taxrank_score = taxrank_score,
      taxrank_facet = taxrank_facet,
      group_facet = group_facet,
      transpose = transpose,
      signif_colors = signif_colors,
      boundary = boundary,
      score_agglom_fn = score_agglom_fn,
      score_name = score_name
    )
  } else {
    .ig_score_boxplot(
      plot_data,
      plot_type = plot_type,
      taxrank_score = taxrank_score,
      taxrank_facet = taxrank_facet,
      group_facet = group_facet,
      transpose = transpose,
      signif_colors = signif_colors,
      boundary = boundary,
      add_stats = add_stats,
      valid_comparisons = valid_comparisons,
      score_agglom_fn = score_agglom_fn,
      score_name = score_name
    )
  }

  if (transpose) {
    plt <- plt +
      coord_flip() +
      .ig_score_transpose_theme()
  }

  plt +
    labs(
      title = "Ig Score by Taxa and Sample Groups",
      subtitle = paste0(
        "agglomerated: ",
        if (first_score_agglom_for_each == "sample") {
          paste0(
            "first for each ",
            group_score,
            " and then for ",
            taxrank_score
          )
        } else if (first_score_agglom_for_each == "taxon") {
          paste0(
            "first for each ",
            taxrank_score,
            " and then for ",
            group_score
          )
        } else {
          paste0("simultaneously for ", group_score, " and ", taxrank_score)
        }
      )
    )
}

#' Agglomerate `ig_coating` Scores To A Taxonomic Rank, Optionally Filtering By Abundance
#'
#' Central-tendency-agglomerates every score in `phyloigseq_obj@ig_coating` (each column
#' also listed in `phyloigseq_obj@score_names`) from one row per taxon/sample to one row per
#' `taxrank` level/sample, and (optionally) drops low-abundance taxa afterward. Used, for
#' example, before feeding a `PhyloIgSeq` object's scores into [PhyloIgSeq_to_phyloseq()] at a
#' coarser taxonomic resolution than the raw ASV/OTU it was scored at.
#'
#' @param phyloigseq_obj A `PhyloIgSeq` object, e.g. from [getPhyloIgSeq()].
#' @param abundance_fraction Optional name of an `ig_coating` fraction-abundance column (e.g.
#'   the positive/negative fraction name passed to [getPhyloIgSeq()]) to agglomerate (by sum)
#'   alongside the scores and, if supplied, to filter taxa by (see `abundance_quantile`/
#'   `min_rel_abundance`) and to weight by when `agglom_method = "weight_by_abund"`. Defaults
#'   to `NULL`, meaning no abundance column is carried through and no abundance filtering is
#'   applied.
#' @param taxrank Column of `phyloigseq_obj@tax_table` to agglomerate to (e.g. `"Genus"`).
#'   Defaults to `NULL`, meaning `"taxon_id"` (no agglomeration across taxa — only samples'
#'   fraction-abundance replicates, if any, get summarized).
#' @param make_unique_taxonomy Logical. Resolve non-unique taxonomy (e.g. shared genus names
#'   across different lineages) via [make_unique_taxa_table()] before agglomerating, so taxa
#'   that only differ upstream of `taxrank` aren't incorrectly merged together.
#' @param agglom_method One of `"mean"`, `"median"`, `"weight_by_abund"` (mean weighted by
#'   `abundance_fraction`, which must then be supplied). Defaults to `NULL`, meaning
#'   `"median"` with a `warning()` — pass an explicit value to silence it.
#' @param abundance_quantile Single probability in `[0, 1]`; taxa with (agglomerated)
#'   `abundance_fraction` below this quantile (per sample) are dropped. Only used when
#'   `abundance_fraction` is supplied. Defaults to `NULL`, meaning `0` (no quantile filter).
#' @param min_rel_abundance Single proportion in `[0, 1]`; taxa with (agglomerated)
#'   `abundance_fraction` below `min_rel_abundance * total_reads` (per sample) are dropped —
#'   see Details for how `total_reads` is determined. Only used when `abundance_fraction` is
#'   supplied. Defaults to `NULL`, meaning `0` (no relative-abundance filter).
#'
#' @details
#' `min_rel_abundance` is evaluated against a per-sample `total_reads`: when
#' `phyloigseq_obj@total_reads` was populated by [getPhyloIgSeq()] (i.e. a
#' `presorting_fraction_name` was supplied there) *and* `abundance_fraction` is that same
#' pre-sort fraction, `total_reads` is the true whole-fraction total. Otherwise it falls back
#' to summing `abundance_fraction` over the taxa remaining in `ig_coating` at this point (after
#' `taxrank` agglomeration, before the `abundance_quantile`/`min_rel_abundance` filter) — an
#' approximation, not a true whole-fraction total, if taxa were already dropped upstream.
#'
#' @return `phyloigseq_obj`, with `ig_coating` (and `score_names`) replaced by the
#'   agglomerated/filtered version, `tax_table` collapsed to the ranks up to and including
#'   `taxrank`, and `total_reads` (if present) restricted to the samples remaining in
#'   `ig_coating`.
#'
#' @examples
#' phyloigseq_obj <- getPhyloIgSeq(
#'   ps_igseq,
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   scores = c("slide_z", "palm")
#' )
#' agglom <- agglomPhyloIgSeq(
#'   phyloigseq_obj,
#'   abundance_fraction = "Pos",
#'   taxrank = "Genus",
#'   agglom_method = "median",
#'   abundance_quantile = 0.1,
#'   min_rel_abundance = 0.001
#' )
#' head(agglom@ig_coating)
#'
#' @export
agglomPhyloIgSeq <- function(
  phyloigseq_obj,
  abundance_fraction = NULL,
  taxrank = NULL,
  make_unique_taxonomy = TRUE,
  agglom_method = NULL,
  abundance_quantile = NULL,
  min_rel_abundance = NULL
) {
  scores <- intersect(colnames(phyloigseq_obj@ig_coating), IG_SCORES)
  scores <- scores[
    !vapply(
      scores,
      function(score) all(is.na(phyloigseq_obj@ig_coating[[score]])),
      logical(1)
    )
  ]

  if (is.null(taxrank)) {
    taxrank <- "taxon_id"
  }

  if (is.null(abundance_quantile)) {
    abundance_quantile <- 0
  }

  if (is.null(min_rel_abundance)) {
    min_rel_abundance <- 0
  }

  if (is.null(agglom_method)) {
    agglom_method <- "median"
    warning("agglomeration method is set to median")
  }

  if (make_unique_taxonomy) {
    phyloigseq_obj@tax_table <- PhyloIgSeq::make_unique_taxa_table(
      phyloigseq_obj@tax_table
    )
  }

  ig_coating_agglom <- .agglomerate_ig_coating(
    ig_coating = phyloigseq_obj@ig_coating,
    tax_table = phyloigseq_obj@tax_table,
    scores = scores,
    taxrank = taxrank,
    abundance_fraction = abundance_fraction,
    agglom_method = agglom_method
  )

  ig_coating_agglom <- .filter_agglomerated_by_abundance(
    ig_coating_agglom = ig_coating_agglom,
    total_reads = phyloigseq_obj@total_reads,
    taxrank = taxrank,
    scores = scores,
    abundance_fraction = abundance_fraction,
    abundance_quantile = abundance_quantile,
    min_rel_abundance = min_rel_abundance
  )

  phyloigseq_obj@ig_coating <- ig_coating_agglom
  phyloigseq_obj@score_names <- scores

  names(phyloigseq_obj@ig_coating)[
    names(phyloigseq_obj@ig_coating) == taxrank
  ] <- "taxon_id"

  # Update taxonomy
  phyloigseq_obj@tax_table <- phyloigseq_obj@tax_table[
    ,
    seq_len(which(colnames(phyloigseq_obj@tax_table) == taxrank))
  ] %>%
    distinct()

  if (!"taxon_id" %in% colnames(phyloigseq_obj@tax_table)) {
    phyloigseq_obj@tax_table$taxon_id <- phyloigseq_obj@tax_table[, taxrank]
  }

  if (!is.null(phyloigseq_obj@total_reads)) {
    phyloigseq_obj@total_reads <- phyloigseq_obj@total_reads[
      phyloigseq_obj@total_reads$sample_id %in%
        phyloigseq_obj@ig_coating$sample_id,
    ]
  }

  return(phyloigseq_obj)
}

#' Build The Per-`taxrank`/Sample Agglomerated `ig_coating` For `agglomPhyloIgSeq()`
#'
#' Internal. Agglomerates every column of `scores` (and, if supplied, `abundance_fraction`) to
#' one row per (`taxrank`, `sample_id`) via `.central_tendency()`/summation. Does not yet
#' attach `total_reads` or apply the abundance filter — see `.filter_agglomerated_by_abundance()`.
#'
#' @noRd
.agglomerate_ig_coating <- function(
  ig_coating,
  tax_table,
  scores,
  taxrank,
  abundance_fraction,
  agglom_method
) {
  agglom_data <- ig_coating %>%
    select(all_of(c("sample_id", "taxon_id", scores, abundance_fraction)))

  if (taxrank != "taxon_id") {
    agglom_data <- merge(
      agglom_data,
      tax_table[, unique(c("taxon_id", taxrank)), drop = FALSE],
      by = "taxon_id"
    )
  }

  agglom_data <- agglom_data %>%
    group_by(sample_id, .data[[taxrank]])

  for (score in scores) {
    agglom_data <- agglom_data %>%
      mutate(
        !!score := .central_tendency(
          .data[[score]],
          agglom_method,
          weights = if (is.null(abundance_fraction)) {
            NULL
          } else {
            .data[[abundance_fraction]]
          }
        )
      )
  }

  if (!is.null(abundance_fraction)) {
    agglom_data <- agglom_data %>%
      mutate(
        !!abundance_fraction := if (all(is.na(.data[[abundance_fraction]]))) {
          NA
        } else {
          sum(.data[[abundance_fraction]], na.rm = TRUE)
        }
      )
  }

  agglom_data %>%
    ungroup() %>%
    select(all_of(c(taxrank, "sample_id", scores, abundance_fraction))) %>%
    distinct()
}

#' Attach `total_reads` and Apply The Abundance Filter For `agglomPhyloIgSeq()`
#'
#' Internal. No-op when `abundance_fraction` is `NULL`. Otherwise attaches a per-sample
#' `total_reads` (from `total_reads` when it covers `abundance_fraction`, else a fallback sum
#' — see `agglomPhyloIgSeq()`'s `@details`) and drops taxa below `abundance_quantile`/
#' `min_rel_abundance`.
#'
#' @noRd
.filter_agglomerated_by_abundance <- function(
  ig_coating_agglom,
  total_reads,
  taxrank,
  scores,
  abundance_fraction,
  abundance_quantile,
  min_rel_abundance
) {
  if (is.null(abundance_fraction)) {
    return(ig_coating_agglom)
  }

  has_matching_total_reads <-
    !is.null(total_reads) && abundance_fraction %in% names(total_reads)

  if (has_matching_total_reads) {
    total_reads_df <- total_reads
    names(total_reads_df)[2] <- "total_reads"
    ig_coating_agglom <- merge(ig_coating_agglom, total_reads_df, by = "sample_id") %>%
      group_by(sample_id)
  } else {
    ig_coating_agglom <- ig_coating_agglom %>%
      group_by(sample_id) %>%
      mutate(total_reads = sum(.data[[abundance_fraction]], na.rm = TRUE))
  }

  ig_coating_agglom %>%
    filter(
      .data[[abundance_fraction]] >=
        quantile(.data[[abundance_fraction]], abundance_quantile, na.rm = TRUE),
      .data[[abundance_fraction]] >= min_rel_abundance * total_reads
    ) %>%
    ungroup() %>%
    select(all_of(c(taxrank, "sample_id", scores, abundance_fraction)))
}

#' Convert Ig Scores to Wide (Sample x Taxon) Format
#'
#' Pivots a long-format `ig_coating` (or agglomerated `ig_coating`) data frame — one row
#' per sample/taxon/score — into one wide sample-by-taxon matrix per score, suitable for
#' export or for feeding into ordination/heatmap functions that expect a taxa-by-samples
#' or samples-by-taxa table.
#'
#' @param ig_coating_agglom A data frame with (at least) `sample_id`, `taxon_id`, and one
#'   column per score, such as the `ig_coating` slot of a [PhyloIgSeq-class] object (see
#'   [getPhyloIgSeq()]/[agglomPhyloIgSeq()]).
#' @param scores Character vector of score column names to pivot, one wide data frame
#'   produced per score. Defaults to `NULL`, meaning every column of `ig_coating_agglom`
#'   also present in `IG_SCORES`.
#' @param shared_by Minimum fraction (in `[0, 1]`) of samples that must have a non-`NA`
#'   value for a taxon for that taxon's column to be kept in the output; taxa below this
#'   threshold are dropped. Defaults to `NULL`, meaning `0` (keep every taxon).
#'
#' @return A named list, one element per entry of `scores`, each a wide data frame with
#'   one row per `sample_id` and one column per `taxon_id` (plus `sample_id` itself)
#'   holding that score's values.
#'
#' @examples
#' ig_coating_agglom <- data.frame(
#'   sample_id = rep(paste0("sample_", 1:3), each = 2),
#'   taxon_id = rep(c("taxon_1", "taxon_2"), times = 3),
#'   slide_z = rnorm(6),
#'   palm = runif(6)
#' )
#' to_wider_ig_score(ig_coating_agglom, scores = c("slide_z", "palm"))
#'
#' @export
to_wider_ig_score <- function(
  ig_coating_agglom,
  scores = NULL,
  shared_by = NULL # what fraction of samples has this taxon
) {
  if (is.null(shared_by)) {
    shared_by <- 0
  }

  if (is.null(scores)) {
    scores <- intersect(colnames(ig_coating_agglom), IG_SCORES)
  }

  score_list <- list()
  for (score in scores) {
    score_list[[score]] <- ig_coating_agglom %>%
      select(all_of(c("sample_id", "taxon_id", score))) %>%
      pivot_wider(
        names_from = taxon_id,
        values_from = !!score
      ) %>%
      ungroup() %>%
      select(where(~ mean(!is.na(.)) >= shared_by))
  }

  return(score_list)
}

#' Extract `sample_data`/`tax_table` from a `PhyloIgSeq` object for `phyloseq` export
#'
#' Shared by both [PhyloIgSeq_to_phyloseq()] code paths: drops the `sample_id` column
#' from `sample_data` (moved to rownames instead, as `phyloseq` expects) and the
#' `taxon_id` column from `tax_table` (likewise moved to rownames), after dropping any
#' `tax_table` row with a missing `taxon_id`.
#'
#' @param phyloigseq_obj A `PhyloIgSeq` object.
#' @return A list with elements `sample_data` and `tax_table`, both matrices/data
#'   frames ready to pass to `phyloseq::sample_data()`/`phyloseq::tax_table()`.
#' @noRd
.phyloigseq_export_metadata <- function(phyloigseq_obj) {
  sample_data_ig_score <- phyloigseq_obj@sample_data[,
    !colnames(phyloigseq_obj@sample_data) %in% "sample_id",
    drop = FALSE
  ]
  rownames(sample_data_ig_score) <- phyloigseq_obj@sample_data$sample_id

  tax_table_ig_score <- as.matrix(phyloigseq_obj@tax_table)
  tax_table_ig_score <- tax_table_ig_score[
    !is.na(tax_table_ig_score[, "taxon_id"]), ,
    drop = FALSE
  ]
  rownames(tax_table_ig_score) <- tax_table_ig_score[, "taxon_id"]
  tax_table_ig_score <- tax_table_ig_score[,
    colnames(tax_table_ig_score) != "taxon_id",
    drop = FALSE
  ]

  list(sample_data = sample_data_ig_score, tax_table = tax_table_ig_score)
}

#' SVD/`incomplete_otu_table` path for [PhyloIgSeq_to_phyloseq()]
#'
#' Builds an [incomplete_otu_table-class] directly from the observed
#' (sample, taxon, score) triples, without ever materialising a dense samples-by-taxa
#' matrix, then fits a low-rank `softImpute::softImpute()` approximation used to impute
#' on demand.
#'
#' @inheritParams PhyloIgSeq_to_phyloseq
#' @return A `phyloseq` object whose `otu_table` is an [incomplete_otu_table-class].
#' @noRd
.PhyloIgSeq_to_phyloseq_svd <- function(
  phyloigseq_obj,
  score_name,
  shared_by,
  svd_rank,
  svd_lambda
) {
  eff_shared_by <- if (is.null(shared_by)) 0 else shared_by

  ig_obs <- stats::na.omit(
    phyloigseq_obj@ig_coating[, c("sample_id", "taxon_id", score_name)]
  )

  samp_lvls <- phyloigseq_obj@sample_data$sample_id
  taxa_lvls <- unique(ig_obs$taxon_id)

  if (eff_shared_by > 0) {
    n_samp <- length(samp_lvls)
    obs_frac <- tapply(ig_obs$sample_id, ig_obs$taxon_id, function(s) {
      length(unique(s)) / n_samp
    })
    taxa_lvls <- names(obs_frac[obs_frac >= eff_shared_by])
    ig_obs <- ig_obs[ig_obs$taxon_id %in% taxa_lvls, , drop = FALSE]
  }

  i_idx <- match(ig_obs$sample_id, samp_lvls)
  j_idx <- match(ig_obs$taxon_id, taxa_lvls)

  X_inc <- softImpute::Incomplete(i_idx, j_idx, ig_obs[[score_name]])
  dimnames(X_inc) <- list(samp_lvls, taxa_lvls)

  # Dense NA matrix for softImpute (passing Incomplete directly is slower)
  n_s <- length(samp_lvls)
  n_t <- length(taxa_lvls)
  X_dense <- matrix(
    NA_real_,
    nrow = n_s,
    ncol = n_t,
    dimnames = list(samp_lvls, taxa_lvls)
  )
  if (length(X_inc@x) > 0L) {
    col_idx <- rep(seq_len(n_t), diff(X_inc@p))
    row_idx <- X_inc@i + 1L
    X_dense[cbind(row_idx, col_idx)] <- X_inc@x
  }

  col_means <- colMeans(X_dense, na.rm = TRUE)
  col_means[is.nan(col_means)] <- 0
  X_centered <- sweep(X_dense, 2, col_means, "-")

  fit <- softImpute::softImpute(
    X_centered,
    rank.max = svd_rank,
    lambda = svd_lambda,
    type = "svd"
  )

  ot_ig <- incomplete_otu_table(
    X_inc = X_inc,
    svd_fit = list(u = fit$u, d = fit$d, v = fit$v),
    col_means = col_means
  )

  metadata <- .phyloigseq_export_metadata(phyloigseq_obj)

  phyloseq(
    ot_ig,
    phyloseq::sample_data(metadata$sample_data),
    phyloseq::tax_table(metadata$tax_table)
  )
}

#' Dense-matrix path for [PhyloIgSeq_to_phyloseq()]
#'
#' Pivots `score_name` to a dense samples-by-taxa matrix via [to_wider_ig_score()] and
#' (unless `imputation_method` is `NULL`) imputes its `NA`s with [dataImpute()].
#'
#' @inheritParams PhyloIgSeq_to_phyloseq
#' @return A `phyloseq` object with a standard dense `otu_table`.
#' @noRd
.PhyloIgSeq_to_phyloseq_dense <- function(
  phyloigseq_obj,
  score_name,
  shared_by,
  imputation_method,
  central_tendency,
  nb_neighbors
) {
  igseq_df <- PhyloIgSeq::to_wider_ig_score(
    ig_coating_agglom = phyloigseq_obj@ig_coating,
    scores = score_name,
    shared_by = shared_by
  )[[score_name]]

  # "OTU table" - Ig scores instead of abundances
  otu_table_ig_score <- as.matrix(igseq_df[,
    !colnames(igseq_df) %in% c("sample_id", "NA")
  ])
  rownames(otu_table_ig_score) <- igseq_df$sample_id
  if (!is.null(imputation_method)) {
    otu_table_ig_score <- PhyloIgSeq::dataImpute(
      otu_table_ig_score,
      method = imputation_method,
      central.tendency = central_tendency,
      nb.neighbors = nb_neighbors
    )
  }

  metadata <- .phyloigseq_export_metadata(phyloigseq_obj)

  phyloseq(
    otu_table(otu_table_ig_score, taxa_are_rows = FALSE),
    sample_data(metadata$sample_data),
    tax_table(metadata$tax_table)
  )
}

#' Convert Ig Scores From A PhyloIgSeq Object To A phyloseq Object
#'
#' Reshapes one score column of `phyloigseq_obj@ig_coating` (long: one row per
#' sample/taxon) into a samples-by-taxa "OTU table" holding score values instead of
#' abundances, then wraps it together with `phyloigseq_obj`'s `sample_data`/`tax_table`
#' into a standard `phyloseq` object — so downstream `phyloseq`-ecosystem analyses
#' (ordination, heatmaps, etc.) can be run directly on Ig-coating scores. Used, for
#' example, to feed a `PhyloIgSeq` object's scores into [get_alpha_diversity()]/
#' [get_beta_diversity()].
#'
#' Two independent code paths exist because `imputation_method = "SVD"` uses a
#' fundamentally different representation of the score matrix than the others:
#' \itemize{
#'   \item `"SVD"` never densifies. Missing (sample, taxon) score entries are kept as
#'     structurally missing (not zero) in an [incomplete_otu_table-class] backed by
#'     `softImpute::Incomplete()`, and `softImpute::softImpute()` fits a low-rank
#'     approximation used to impute on demand. Suitable when most (sample, taxon) pairs
#'     are unobserved, e.g. a raw ASV/OTU-level score matrix.
#'   \item Any other value (`NULL`, `"KNN"`, `"Central Tendency"`, `"Replace NA with 0"`,
#'     the values [dataImpute()] understands) pivots to a dense samples-by-taxa matrix
#'     via [to_wider_ig_score()] first, then (unless `NULL`) imputes the resulting `NA`s
#'     with [dataImpute()].
#' }
#'
#' @param phyloigseq_obj A `PhyloIgSeq` object, e.g. from [getPhyloIgSeq()]/
#'   [agglomPhyloIgSeq()].
#' @param score_name Name of one column of `phyloigseq_obj@ig_coating` (e.g.
#'   `"slide_z"`) to convert into the returned object's `otu_table`.
#' @param shared_by Minimum fraction (in `[0, 1]`) of samples that must have a non-`NA`
#'   `score_name` value for a taxon to be kept; taxa below this threshold are dropped.
#'   Defaults to `NULL`, meaning `0` (keep every taxon with at least one observed
#'   value).
#' @param imputation_method One of `NULL` (no imputation — the returned `otu_table`
#'   keeps its missing entries), `"SVD"` (see Details), `"KNN"`, `"Central Tendency"`,
#'   or `"Replace NA with 0"` (the latter three implemented by [dataImpute()]).
#'   Defaults to `NULL`.
#' @param central_tendency Passed through to [dataImpute()]'s `central.tendency` when
#'   `imputation_method = "Central Tendency"`: one of `"mean"`, `"median"`, `"mode"`.
#'   Ignored otherwise.
#' @param nb_neighbors Passed through to [dataImpute()]'s `nb.neighbors` when
#'   `imputation_method = "KNN"`: number of neighbors used by `VIM::kNN()`. Ignored
#'   otherwise. Defaults to `5`, a common rule-of-thumb starting point for k-NN
#'   imputation — tune for your data.
#' @param svd_rank `rank.max` passed to `softImpute::softImpute()`: maximum rank of the
#'   low-rank approximation. Only used when `imputation_method = "SVD"`. Defaults to
#'   `50L`, an arbitrary-but-generous default — lower it for small taxa counts or to
#'   speed up fitting, tune higher only if cross-validation supports it.
#' @param svd_lambda `lambda` (nuclear-norm regularization strength) passed to
#'   `softImpute::softImpute()`. Only used when `imputation_method = "SVD"`. Defaults
#'   to `1`, `softImpute::softImpute()`'s own default — tune for your data.
#'
#' @return A `phyloseq` object whose `otu_table` holds `score_name` values (samples as
#'   rows): an [incomplete_otu_table-class] when `imputation_method = "SVD"`, otherwise
#'   a standard dense `phyloseq::otu_table`.
#'
#' @examples
#' data(ps_igseq)
#' pis <- getPhyloIgSeq(
#'   physeq = ps_igseq,
#'   sample_id_name = "sample_id",
#'   fraction_id_name = "sorting_fraction",
#'   positive_fraction_name = "Pos",
#'   first_negative_fraction_name = "Neg1",
#'   second_negative_fraction_name = "Neg2",
#'   scores = c("slide_z", "palm")
#' )
#' ps_svd <- PhyloIgSeq_to_phyloseq(
#'   pis,
#'   score_name = "slide_z",
#'   imputation_method = "SVD",
#'   svd_rank = 5L
#' )
#' ps_dense <- PhyloIgSeq_to_phyloseq(
#'   pis,
#'   score_name = "palm",
#'   imputation_method = "Replace NA with 0"
#' )
#'
#' @export
PhyloIgSeq_to_phyloseq <- function(
  phyloigseq_obj,
  score_name,
  shared_by = NULL,
  imputation_method = NULL,
  central_tendency = NULL,
  nb_neighbors = 5,
  svd_rank = 50L,
  svd_lambda = 1
) {
  if (identical(imputation_method, "SVD")) {
    .PhyloIgSeq_to_phyloseq_svd(
      phyloigseq_obj = phyloigseq_obj,
      score_name = score_name,
      shared_by = shared_by,
      svd_rank = svd_rank,
      svd_lambda = svd_lambda
    )
  } else {
    .PhyloIgSeq_to_phyloseq_dense(
      phyloigseq_obj = phyloigseq_obj,
      score_name = score_name,
      shared_by = shared_by,
      imputation_method = imputation_method,
      central_tendency = central_tendency,
      nb_neighbors = nb_neighbors
    )
  }
}
