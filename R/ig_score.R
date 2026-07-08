#' Compute Other Ig Scores
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

  if (method == "palm") {
    score <- pos_abund / neg_abund
  } else if (method == "kau") {
    # minus to negate the fact that log10(pos_abund*neg_abund) is negative
    score <- -log2(pos_abund / neg_abund) / log10(pos_abund * neg_abund)
  } else if (method == "prob_index") {
    score <- pos_abund * ig_freq / pre_abund # P(Ig+ | taxon) = P(taxon | Ig+) * P(Ig+) / P(taxon)
  } else if (method == "prob_ratio") {
    score <- log2(pos_abund * ig_freq / (neg_abund * (1 - ig_freq)))
  } else if (method == "purity_corrected_prob_index") {
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
    score <- (pos_abund *
      pos_purity *
      pos_fraction +
      neg_abund * neg_impurity * neg_fraction) /
      pre_abund
  } else if (method == "purity_corrected_prob_ratio") {
    # TODO: verify!
    prob <- pos_abund *
      pos_purity *
      pos_fraction +
      neg_abund * neg_impurity * neg_fraction
    score <- log2(prob / (1 - prob))
  } else {
    stop("Wrong score type")
  }

  score[is.nan(score) | is.infinite(score)] <- NA

  return(score)
}


#' Plot Slide Z Score
#' @export
plot_slide_z <- function(
  phyloigseq_obj,
  sample_ids = NULL, # if NULL, all samples are plotted
  empirical_null_distribution = TRUE,
  z_alpha2 = 1.96,
  signif_colors = c(ggsci::pal_npg()(2)[2], ggsci::pal_npg()(2)[1]),
  ellipses = TRUE
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

  if (length(phyloigseq_obj@tax_table) > 0) {
    hover.text.all <- rep(NA, nrow(ig_df))
    for (i in 1:nrow(ig_df)) {
      # it makes sure that taxa names match
      taxon_id <- ig_df$taxon_id[i]
      hover.text <- ""
      # TODO: for now original taxon name is not kept
      values <- phyloigseq_obj@tax_table[
        phyloigseq_obj@tax_table$taxon_id == taxon_id,
        !colnames(phyloigseq_obj@tax_table) %in% c("taxon_id", "taxon_name"),
        drop = FALSE
      ]
      for (variable in colnames(values)) {
        value <- values[[variable]]
        hover.text <- paste0(hover.text, variable, ": ", value, "<br>")
      }
      hover.text.all[i] <- paste0(
        hover.text,
        paste("slide Z: ", round(ig_df$slide_z[i], digits = 3))
      )
    }
    ig_df$tooltip <- hover.text.all
  } else {
    ig_df <- mutate(
      ig_df,
      tooltip = paste("<br>slide_z: ", round(slide_z, digits = 3))
    )
  }

  ig_df_non_imputed <- data.frame()
  ig_df_imputed <- data.frame()

  for (sample_id in unique(ig_df$sample_id)) {
    ig_df_non_imputed <- rbind(
      ig_df_non_imputed,
      ig_df[
        ig_df$sample_id %in%
          sample_id &
          !ig_df$taxon_id %in% phyloigseq_obj@imputed_taxa[[sample_id]],
        ,
        drop = FALSE
      ]
    )
    ig_df_imputed <- rbind(
      ig_df_imputed,
      ig_df[
        ig_df$sample_id %in%
          sample_id &
          ig_df$taxon_id %in% phyloigseq_obj@imputed_taxa[[sample_id]],
        ,
        drop = FALSE
      ]
    )
  }
  stat_imputed <- ifelse(
    ig_df_imputed$slide_z >= z_alpha2 | ig_df_imputed$slide_z <= -z_alpha2,
    "signif",
    "ns"
  )

  ig_df <- ig_df_non_imputed
  stat <- ifelse(
    ig_df$slide_z >= z_alpha2 | ig_df$slide_z <= -z_alpha2,
    "signif",
    "ns"
  )

  plt <- ggplot(ig_df)

  jitter_width <- diff(range(c(ig_df$null_abundance, ig_df$obs_abundance))) / 6
  jitter_x <- min(c(ig_df$null_abundance, ig_df$obs_abundance)) -
    jitter_width * 3
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
        #shape = stat_imputed,
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
    #ggsci::scale_color_npg()+
    scale_color_manual(values = signif_colors) +
    scale_size_discrete(range = c(1.5, 3))

  if (length(ellipse_df) != 0 & ellipses) {
    plt <- plt +
      geom_path(
        data = ellipse_df,
        aes(x = x, y = y, group = ellipse_level),
        color = "darkgrey",
        linetype = 2
      )
  }

  if (length(unique(ig_df$sample_id)) > 1) {
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
        if (length(unique(ig_df$sample_id)) == 1) {
          paste0(" of ", unique(ig_df$sample_id))
        }
      )
    ) +
    theme_minimal() +
    ggplot2::theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.title = element_text(face = "bold", hjust = 0.5)
    )
  return(plt)
}

#' Plot Any Ig score
#' @export
plot_ig_score <- function(
  phyloigseq_obj,
  plot_type = c("boxplot", "bubbleplot")[1],
  score_name = "slide_z",
  taxrank_score = "taxon_id", # taxrank level to agglomerate the Ig score
  taxrank_facet = NULL, # taxrank for faceting
  group_score = "sample_id", # sample group to agglomerate the Ig score
  group_facet = NULL, # sample group for facetting
  score_agglom_fn = c("mean", "median")[1],
  first_score_agglom_for_each = c("sample", "taxon", "both")[1],
  z_alpha2 = 1.96, # in case of z score
  exclude_na = TRUE,
  transpose = FALSE,
  signif_colors = ggsci::pal_npg()(2),
  add_stats = TRUE
) {
  if (score_name == "slide_z" & is.null(z_alpha2)) {
    z_alpha2 <- 1.96
  }

  # TODO: add a possibility for multiple faceting (e.g. with timepoints)
  # FIXME: fix the mess with faceting when a x or y is not unique in each facet

  # Agglomeration of Ig score should be either with 1. median 2. mean 3. weighted average by abundance.
  # TODO: 3.

  # TODO: check the correctness of data agglomeration and averaging
  if (!score_agglom_fn %in% c("mean", "median")) {
    stop("`score_agglom_fn` should be either 'mean' or 'median'")
  }
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

  if (first_score_agglom_for_each == "both") {
    plot_data <- plot_data %>%
      # Compute central tendency (mean, median) in a sample group and taxrank
      # - for all values in intersections formed by these groupings
      group_by(.data[[group_score]], .data[[taxrank_score]]) %>%
      mutate(
        agglom_score = get(score_agglom_fn)(.data[[score_name]], na.rm = TRUE)
      ) %>%
      select(all_of(unique(c(
        taxrank_score,
        taxrank_facet,
        group_score,
        group_facet,
        "agglom_score"
      )))) %>%
      distinct() %>%
      ungroup()
  } else if (first_score_agglom_for_each %in% c("sample", "taxon")) {
    if (first_score_agglom_for_each == "sample") {
      first_id <- "sample_id"
      first_agglom <- taxrank_score
    } else if (first_score_agglom_for_each == "taxon") {
      first_id <- "taxon_id"
      first_agglom <- group_score
    }

    plot_data <- plot_data %>%
      # First, compute central tendency separately for each individual sample (or taxon)
      # grouping by a taxrank (or sample group)
      group_by(.data[[first_id]], .data[[first_agglom]]) %>%
      mutate(
        agglom_score = get(score_agglom_fn)(.data[[score_name]], na.rm = TRUE)
      ) %>%
      # ungroup and get rid of duplications, otherwise the central tendency will be false!
      ungroup() %>%
      select(all_of(unique(c(
        taxrank_score,
        taxrank_facet,
        group_score,
        group_facet,
        "agglom_score"
      )))) %>%
      distinct() %>%
      # Then, compute central tendency for each sample group (or taxrank), based on central tendencies for
      # each sample per taxrank ( or for each taxon per sample group)
      group_by(.data[[group_score]], .data[[taxrank_score]]) %>%
      mutate(
        agglom_score = get(score_agglom_fn)(agglom_score, na.rm = TRUE)
      ) %>%
      ungroup()
  } else {
    stop(
      "`first_score_agglom_by` argument should be 'sample', 'taxon' or 'both'"
    )
  }

  if (exclude_na) {
    plot_data <- na.omit(plot_data)
  }

  if (add_stats) {
    # FIXME: didn't account for group_score
    group_counts <- plot_data %>%
      count(.data[[taxrank_score]]) %>%
      filter(n >= 2)
    valid_groups <- group_counts[[taxrank_score]]
    if (length(valid_groups) > 1) {
      valid_comparisons <- combn(valid_groups, 2, simplify = FALSE)
    } else {
      add_stats <- FALSE
    }
  }

  if (score_name == "slide_z") {
    left_lim <- -z_alpha2
    right_lim <- z_alpha2
    midpoint <- 0
    left_boundary <- -Inf
    right_boundary <- +Inf
  } else if (score_name %in% c("kau", "prob_ratio")) {
    left_lim <- 0
    right_lim <- 0
    midpoint <- 0
    left_boundary <- -Inf
    right_boundary <- +Inf
  } else if (score_name == "palm") {
    left_lim <- 1
    right_lim <- 1
    midpoint <- 1
    left_boundary <- 0
    right_boundary <- +Inf
  } else if (score_name == "prob_index") {
    left_lim <- 0.5
    right_lim <- 0.5
    midpoint <- 0.5
    left_boundary <- 0
    right_boundary <- 1
  }

  if (plot_type == "bubbleplot") {
    plt <-
      ggplot(
        plot_data,
        aes(
          x = .data[[group_score]],
          y = .data[[taxrank_score]],
          size = abs(agglom_score),
          fill = agglom_score
        )
      ) +
      geom_point(pch = 21)

    plt <- plt +
      scale_fill_gradient2(
        high = signif_colors[1],
        low = signif_colors[2],
        midpoint = midpoint,
        limits = c(
          max(
            left_boundary,
            midpoint - max(abs(plot_data$agglom_score - midpoint))
          ),
          min(
            right_boundary,
            midpoint + max(abs(plot_data$agglom_score - midpoint))
          )
        ),
        guide = guide_colourbar(title.position = "top", title.hjust = 0.5)
      )

    plt <- plt +
      guides(size = "none")

    if (!transpose) {
      plt <- plt +
        facet_grid(
          rows = if (!is.null(taxrank_facet)) {
            vars(.data[[taxrank_facet]])
          },
          cols = if (!is.null(group_facet)) {
            vars(.data[[group_facet]])
          },
          scales = "free",
          space = "free"
        )
    } else {
      plt <- plt +
        facet_grid(
          cols = if (!is.null(taxrank_facet)) {
            vars(.data[[taxrank_facet]])
          },
          rows = if (!is.null(group_facet)) {
            vars(.data[[group_facet]])
          },
          scales = "free",
          space = "free"
        )
    }

    plt <- plt +
      theme_minimal() +
      labs(x = NULL, y = NULL, fill = paste(score_agglom_fn, score_name)) +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.direction = "horizontal",
        axis.text.y.left = element_text(angle = 0, hjust = 1),
        strip.text.y.right = element_text(angle = 0, hjust = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
  } else if (plot_type %in% c("boxplot", "violin")) {
    plot_data$point_color <- ifelse(
      plot_data$agglom_score > right_lim,
      "high",
      ifelse(plot_data$agglom_score < left_lim, "low", "ns")
    )

    plt <-
      ggplot(
        data = plot_data,
        aes(x = agglom_score, y = .data[[taxrank_score]])
      )
    if (plot_type == "violin") {
      plt <- plt + geom_violin()
    } else {
      plt <- plt + geom_boxplot(outliers = FALSE)
    }

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
          size = abs(agglom_score - mean(left_lim, right_lim)),
          color = point_color
        )
      ) +
      scale_color_manual(
        values = c(
          "high" = signif_colors[1],
          "low" = signif_colors[2],
          "ns" = "darkgrey"
        )
      ) +
      guides(size = "none", color = "none") +
      scale_size_continuous(range = c(1, 2))

    if (mean(left_lim, right_lim) == 0) {
      plt <- plt +
        scale_x_continuous(
          limits = c(
            -max(abs(plot_data$agglom_score)),
            max(abs(plot_data$agglom_score))
          )
        )
    }

    plt <- plt +
      geom_vline(xintercept = unique(c(left_lim, right_lim)), linetype = 2)

    if (!transpose) {
      plt <- plt +
        facet_grid(
          rows = if (!is.null(taxrank_facet)) {
            vars(.data[[taxrank_facet]])
          },
          cols = if (!is.null(group_facet)) {
            vars(.data[[group_facet]])
          },
          scales = "free_y",
          space = "free"
        )
    } else {
      plt <- plt +
        facet_grid(
          cols = if (!is.null(taxrank_facet)) {
            vars(.data[[taxrank_facet]])
          },
          rows = if (!is.null(group_facet)) {
            vars(.data[[group_facet]])
          },
          scales = "free_x",
          space = "free"
        )
    }

    plt <- plt +
      theme_minimal() +
      labs(x = paste(score_agglom_fn, score_name), y = NULL) +
      theme(
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 10, hjust = 0.5),
        legend.title = element_text(face = "bold", hjust = 0.5),
        legend.direction = "horizontal",
        axis.text.y.left = element_text(angle = 0, hjust = 1),
        strip.text.y.right = element_text(angle = 0, hjust = 0),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      )
  } else {
    stop("Plot type should be 'boxplot', 'violin' or 'bubbleplot'")
  }

  if (transpose) {
    plt <- plt +
      coord_flip() +
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
  plt <- plt +
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
  return(plt)
}

#' Agglomerate Ig Scores by Taxrank
#' @export
agglomPhyloIgSeq <- function(
  phyloigseq_obj,
  abundance_fraction = NULL,
  taxrank = NULL,
  make_unique_taxonomy = TRUE,
  agglom_method = NULL, # mean, median or weighted average by abundance
  # TODO: add total reads by sample (in the whole fraction)
  abundance_quantile = NULL,
  min_rel_abundance = NULL
) {
  scores <- intersect(colnames(phyloigseq_obj@ig_coating), IG_SCORES)

  for (score in scores) {
    if (all(is.na(phyloigseq_obj@ig_coating[[score]]))) {
      scores <- scores[scores != score]
    }
  }

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

  agglom_fn <- function(score, abundances, agglom_method) {
    if (agglom_method == "mean") {
      return(mean(score, na.rm = TRUE))
    } else if (agglom_method == "median") {
      return(median(score, na.rm = TRUE))
    } else if (agglom_method == "weight_by_abund") {
      if (is.null(abundances)) {
        stop("Need abundance fraction to compute weighted score")
      }
      return(weighted.mean(score, abundances, na.rm = TRUE))
    } else {
      stop("wrong agglomeration method")
    }
  }

  if (make_unique_taxonomy) {
    phyloigseq_obj@tax_table <- PhyloIgSeq::make_unique_taxa_table(
      phyloigseq_obj@tax_table
    )
  }

  # Agglomerate scores
  ig_coating_agglom <- phyloigseq_obj@ig_coating %>%
    select(all_of(c("sample_id", "taxon_id", scores, abundance_fraction))) %>%
    {
      if (taxrank == "taxon_id") {
        .
      } else {
        merge(
          .,
          phyloigseq_obj@tax_table[,
            unique(c("taxon_id", taxrank)),
            drop = FALSE
          ],
          by = "taxon_id"
        )
      }
    } %>%
    group_by(sample_id, .data[[taxrank]])

  for (score in scores) {
    ig_coating_agglom <- ig_coating_agglom %>%
      mutate(
        !!score := agglom_fn(
          .data[[score]],
          if (is.null(abundance_fraction)) {
            NULL
          } else {
            .data[[abundance_fraction]]
          },
          agglom_method
        )
      )
  }

  ig_coating_agglom <- ig_coating_agglom %>%
    {
      if (!is.null(abundance_fraction)) {
        mutate(
          .,
          !!abundance_fraction := if (all(is.na(.data[[abundance_fraction]]))) {
            NA
          } else {
            sum(.data[[abundance_fraction]], na.rm = TRUE)
          }
        )
      } else {
        .
      }
    } %>%
    ungroup() %>%
    select(all_of(c(taxrank, "sample_id", abundance_fraction, scores))) %>%
    distinct() %>%
    {
      if (!is.null(abundance_fraction)) {
        if (
          !is.null(phyloigseq_obj@total_reads) &&
            abundance_fraction %in% names(phyloigseq_obj@total_reads)
        ) {
          total_reads_df <- phyloigseq_obj@total_reads
          names(total_reads_df)[2] <- "total_reads"
          merge(., total_reads_df, by = "sample_id") %>%
            #mutate(., total_reads = phyloigseq_obj@total_reads[[abundance_fraction]][phyloigseq_obj@total_reads[["sample_id"]] == sample_id] )
            group_by(., sample_id)
        } else {
          group_by(., sample_id) %>%
            mutate(
              .,
              total_reads = sum(.data[[abundance_fraction]], na.rm = TRUE)
            )
        }
      } else {
        .
      }
    } %>%
    {
      if (!is.null(abundance_fraction)) {
        filter(
          .,
          .data[[abundance_fraction]] >=
            quantile(
              .data[[abundance_fraction]],
              abundance_quantile,
              na.rm = TRUE
            ),
          .data[[abundance_fraction]] >= min_rel_abundance * total_reads
        )
      } else {
        .
      }
    } %>%
    ungroup() %>%
    select(all_of(c(taxrank, "sample_id", abundance_fraction, scores)))

  phyloigseq_obj@ig_coating <- ig_coating_agglom

  names(phyloigseq_obj@ig_coating)[
    names(phyloigseq_obj@ig_coating) == taxrank
  ] <- "taxon_id"

  # Update taxonomy
  phyloigseq_obj@tax_table <- phyloigseq_obj@tax_table[,
    1:which(colnames(phyloigseq_obj@tax_table) == taxrank)
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

#' Convert Ig Scores to Longer Format
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

#' Convert Ig scores from PhyloIgSeq object to phyloseq
#' @export
PhyloIgSeq_to_phyloseq <-
  function(
    phyloigseq_obj,
    score_name,
    shared_by = NULL, # what fraction of samples has this taxon
    imputation_method = NULL,
    central_tendency = NULL, # "mean", "median" or "mode"
    nb_neighbors = 5,
    svd_rank = 50L, # rank.max for softImpute (SVD path only)
    svd_lambda = 1 # regularization lambda for softImpute (SVD path only)
  ) {
    # --- SVD path: build incomplete_otu_table without materialising a dense matrix ---
    if (identical(imputation_method, "SVD")) {
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

      # Sample data
      sample_data_ig_score <- phyloigseq_obj@sample_data[,
        !colnames(phyloigseq_obj@sample_data) %in% "sample_id",
        drop = FALSE
      ]
      rownames(sample_data_ig_score) <- phyloigseq_obj@sample_data$sample_id

      # Taxonomy table
      tax_table_ig_score <- phyloigseq_obj@tax_table %>% as.matrix()
      tax_table_ig_score <- tax_table_ig_score[
        !is.na(tax_table_ig_score[, "taxon_id"]),
        ,
        drop = FALSE
      ]
      rownames(tax_table_ig_score) <- tax_table_ig_score[, "taxon_id"]
      tax_table_ig_score <- tax_table_ig_score[,
        colnames(tax_table_ig_score) != "taxon_id",
        drop = FALSE
      ]

      return(phyloseq(
        ot_ig,
        phyloseq::sample_data(sample_data_ig_score),
        phyloseq::tax_table(tax_table_ig_score)
      ))
    }

    # --- Legacy path (NULL / "KNN" / "Central Tendency" / "Replace NA with 0") ---
    igseq_df <-
      PhyloIgSeq::to_wider_ig_score(
        ig_coating_agglom = phyloigseq_obj@ig_coating,
        scores = score_name,
        shared_by = shared_by
      )[[score_name]]

    # "OTU table" - Ig scores instead of abundances
    otu_table_ig_score <- igseq_df
    otu_table_ig_score <- as.matrix(otu_table_ig_score[,
      !colnames(otu_table_ig_score) %in% c("sample_id", "NA")
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

    # Sample data
    sample_data_ig_score <- phyloigseq_obj@sample_data[,
      !colnames(phyloigseq_obj@sample_data) %in% c("sample_id"),
      drop = FALSE
    ]
    rownames(sample_data_ig_score) <- phyloigseq_obj@sample_data$sample_id

    # Taxonomy table
    tax_table_ig_score <- phyloigseq_obj@tax_table %>% as.matrix()
    tax_table_ig_score <- tax_table_ig_score[
      !is.na(tax_table_ig_score[, "taxon_id"]),
    ]
    rownames(tax_table_ig_score) <- tax_table_ig_score[, "taxon_id"]
    tax_table_ig_score <- tax_table_ig_score[,
      colnames(tax_table_ig_score) != "taxon_id"
    ]

    # Put all in one phyloseq object
    return(phyloseq(
      otu_table(otu_table_ig_score, taxa_are_rows = FALSE),
      sample_data(sample_data_ig_score),
      tax_table(tax_table_ig_score)
    ))
  }
